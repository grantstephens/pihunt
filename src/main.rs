use clap::{Parser, Subcommand, ValueEnum};
use pihunt::basis::{Columns, Extra, Shape};
use pihunt::log::{Kind, Record};
use pihunt::report::{describe, formula};
use pihunt::runner::{Ctx, load_done, run_chain};
use pihunt::shard::Shard;
use pihunt::{config, log, plan, report, verify};
use rayon::prelude::*;
use rug::Integer;
use std::path::PathBuf;
use std::process::ExitCode;
use std::sync::atomic::{AtomicUsize, Ordering};

#[derive(Parser)]
#[command(version, about = "Hunt for BBP-type formulas for pi with PSLQ")]
struct Cli {
    #[command(subcommand)]
    cmd: Cmd,
}

/// Which n-th-digit algorithm to use. Default is `Thm2` (Gourdon Theorem 2, the chunked
/// remainder-tree algorithm in `docs/nthdigit-theorem2.md`): about 25x faster than Theorem 1
/// at n = 10^6 and ~60x at 10^7, with peak memory in the tens of MiB (69 MiB at 10^7, see
/// `docs/nthdigit.md`), MPFR-verified through 10^7. `Thm1` (Theorem 1, `O(log² n)` memory,
/// `O(n²)`-ish time) stays available as the strict-memory option: a flat ~5 MiB at any
/// position, at the cost of hours instead of minutes at 10^7.
#[derive(Clone, Copy, Debug, ValueEnum, PartialEq, Eq)]
enum Method {
    Thm1,
    Thm2,
}

#[derive(Subcommand)]
enum Cmd {
    /// Run every job in a batch file, appending results to its output log.
    Run {
        batch: PathBuf,
        /// Run only shard k of N (1-based), e.g. "2/4". Splits the output log too.
        #[arg(long)]
        shard: Option<Shard>,
    },
    /// Dry run: show how big a batch is without running it.
    Plan {
        batch: PathBuf,
        /// Report the job count for shard k of N (1-based), e.g. "2/4".
        #[arg(long)]
        shard: Option<Shard>,
    },
    /// Re-verify every hit in a results log at 2x precision.
    Verify { results: PathBuf },
    /// Markdown exclusion report over one or more results logs, to stdout.
    Report {
        #[arg(required = true)]
        results: Vec<PathBuf>,
    },
    /// Print `count` decimal digits of pi starting at position `pos` (1-based; position 1
    /// is the '1' in 3.14159...).
    Digit {
        pos: u64,
        #[arg(long, default_value_t = 10)]
        count: usize,
        /// Which algorithm to use; see [`Method`].
        #[arg(long, value_enum, default_value_t = Method::Thm2)]
        method: Method,
        /// Theorem 2 only: ART chunk modulus size in bits (memory budget). Default scales as
        /// `~4*sqrt(pos)` decimal digits (`pihunt::nthdigit2::default_mem_bits`), the "headline"
        /// `m ∝ √n` case from `docs/nthdigit-theorem2.md` §6.1.
        #[arg(long)]
        mem: Option<u64>,
    },
    /// Stream decimal digits of pi forever, in independently-computed blocks (memory never
    /// grows), printing each block to stdout as soon as it's ready.
    Stream {
        #[arg(long, default_value_t = 1)]
        from: u64,
        #[arg(long, default_value_t = 10)]
        block: usize,
        /// Stop after this many blocks. Hidden: for tests only.
        #[arg(long, hide = true)]
        blocks: Option<u64>,
        /// Which algorithm to use; see [`Method`].
        #[arg(long, value_enum, default_value_t = Method::Thm2)]
        method: Method,
        /// Theorem 2 only: ART chunk modulus size in bits. See `digit --mem`.
        #[arg(long)]
        mem: Option<u64>,
    },
}

fn main() -> ExitCode {
    let result = match Cli::parse().cmd {
        Cmd::Run { batch, shard } => run(batch, shard),
        Cmd::Plan { batch, shard } => show_plan(batch, shard),
        Cmd::Verify { results } => reverify(results),
        Cmd::Report { results } => write_report(results),
        Cmd::Digit {
            pos,
            count,
            method,
            mem,
        } => digit_cmd(pos, count, method, mem),
        Cmd::Stream {
            from,
            block,
            blocks,
            method,
            mem,
        } => stream_cmd(from, block, blocks, method, mem),
    };
    result.unwrap_or_else(|e| {
        eprintln!("error: {e}");
        ExitCode::FAILURE
    })
}

fn run(path: PathBuf, shard: Option<Shard>) -> Result<ExitCode, String> {
    let batch = config::load(&path)?;
    let ctx = Ctx {
        batch: &batch.name,
        finder: batch.defaults.finder.finder(),
        max_columns: batch.defaults.max_columns,
        escalate: batch.defaults.escalate,
    };
    let mut jobs = plan::plan(&batch);
    if let Some(shard) = shard {
        jobs.retain(|j| shard.owns(j, ctx.finder.name()));
    }
    let output = shard.map_or_else(|| batch.output.clone(), |s| s.output_path(&batch.output));
    let threads = batch
        .threads
        .unwrap_or_else(|| std::thread::available_parallelism().map_or(1, |n| n.get()));
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .map_err(|e| e.to_string())?;
    let done = load_done(&output)?;
    let (tx, writer) = log::spawn_writer(&output).map_err(|e| e.to_string())?;
    let (finished, ran, resumed) = (
        AtomicUsize::new(0),
        AtomicUsize::new(0),
        AtomicUsize::new(0),
    );
    let total = jobs.len();
    eprintln!(
        "{}: {total} jobs, {} finder, on {threads} threads → {} ({} attempts already logged)",
        batch.name,
        ctx.finder.name(),
        output.display(),
        done.len()
    );

    let records: Vec<Record> = pool.install(|| {
        jobs.par_iter()
            .map_with(tx, |tx, j| {
                let chain = run_chain(j, &ctx, &done, &mut |rec| {
                    ran.fetch_add(1, Ordering::Relaxed);
                    tx.send(rec.clone()).expect("writer thread alive");
                });
                resumed.fetch_add(chain.resumed, Ordering::Relaxed);
                let rec = chain.last;
                let n = finished.fetch_add(1, Ordering::Relaxed) + 1;
                let esc = if rec.escalated_from.is_some() {
                    ", escalated"
                } else {
                    ""
                };
                eprintln!(
                    "[{n}/{total}] {} → {:?} ({} ms, {} digits{esc})",
                    describe(&rec),
                    rec.outcome,
                    rec.elapsed_ms,
                    rec.params.precision_digits
                );
                rec
            })
            .collect()
    });
    writer
        .join()
        .expect("writer thread panicked")
        .map_err(|e| e.to_string())?;
    println!(
        "ran {} attempts, resumed {} attempts from the log",
        ran.into_inner(),
        resumed.into_inner()
    );
    summarise(&records);
    Ok(ExitCode::SUCCESS)
}

fn show_plan(path: PathBuf, shard: Option<Shard>) -> Result<ExitCode, String> {
    let batch = config::load(&path)?;
    let jobs = plan::plan(&batch);
    let finder = batch.defaults.finder.finder();
    let max_cols = batch.defaults.max_columns;
    let (skipped, run): (Vec<_>, Vec<_>) = jobs.iter().partition(|j| j.shape.columns() > max_cols);
    println!(
        "{}: {} jobs ({} to run, {} skipped over max_columns {max_cols}), {} finder, escalate up to {}x digits",
        batch.name,
        jobs.len(),
        run.len(),
        skipped.len(),
        finder.name(),
        1u32 << batch.defaults.escalate
    );
    if let (Some(n), Some(lo), Some(hi)) = (
        run.iter().map(|j| j.shape.columns()).max(),
        run.iter().map(|j| j.digits).min(),
        run.iter().map(|j| j.digits).max(),
    ) {
        println!("largest n = {n}, precision {lo}..={hi} digits");
    }
    if let Some(shard) = shard {
        let owned = jobs.iter().filter(|j| shard.owns(j, finder.name())).count();
        println!(
            "shard {}/{}: {owned} jobs → {}",
            shard.k,
            shard.n,
            shard.output_path(&batch.output).display()
        );
    }
    Ok(ExitCode::SUCCESS)
}

fn write_report(paths: Vec<PathBuf>) -> Result<ExitCode, String> {
    let mut records = Vec::new();
    for path in &paths {
        records.extend(log::read(path)?);
    }
    print!("{}", report::render(&records));
    Ok(ExitCode::SUCCESS)
}

fn reverify(path: PathBuf) -> Result<ExitCode, String> {
    let records = log::read(&path)?;
    let mut failed = 0;
    let hits: Vec<&Record> = records.iter().filter(|r| r.outcome == Kind::Hit).collect();
    for rec in &hits {
        let p = &rec.params;
        let extras = p
            .extras
            .iter()
            .map(|e| e.parse::<Extra>())
            .collect::<Result<Vec<_>, _>>()?;
        let shape = Shape::new(p.base, p.period, p.degrees[0], p.degrees[1], extras);
        let relation = rec.relation.as_ref().ok_or("hit without relation")?;
        let relation = relation
            .iter()
            .map(|c| c.parse::<Integer>().map_err(|e| e.to_string()))
            .collect::<Result<Vec<_>, _>>()?;
        let cols = Columns::build(&shape, p.precision_digits);
        let (ok, r) = verify::passes(&relation, &cols.hi, p.precision_digits);
        println!(
            "{} {} residual 1e{r:.0}",
            if ok { "OK  " } else { "FAIL" },
            describe(rec)
        );
        failed += !ok as usize;
    }
    println!("{} hits checked, {failed} failed", hits.len());
    Ok(if failed == 0 {
        ExitCode::SUCCESS
    } else {
        ExitCode::FAILURE
    })
}

/// Peak resident set size in KiB, read from `/proc/self/status` (`VmHWM`). `None` if the
/// file can't be read or parsed (e.g. non-Linux) — we don't want a missing metric to fail
/// the digit computation itself.
fn peak_rss_kb() -> Option<u64> {
    let status = std::fs::read_to_string("/proc/self/status").ok()?;
    status.lines().find_map(|line| {
        line.strip_prefix("VmHWM:")?
            .split_whitespace()
            .next()?
            .parse()
            .ok()
    })
}

/// Computes `count` digits at 0-based position `n` with the requested method. `mem` (bits) is
/// only meaningful for `Thm2`; `None` picks `pihunt::nthdigit2::default_mem_bits(n)`.
fn digits_with(method: Method, n: u64, count: usize, mem: Option<u64>) -> String {
    match method {
        Method::Thm1 => pihunt::nthdigit::digits(n, count),
        Method::Thm2 => {
            let mem_bits = mem.unwrap_or_else(|| pihunt::nthdigit2::default_mem_bits(n.max(1)));
            pihunt::nthdigit2::digits(n, count, mem_bits)
        }
    }
}

fn digit_cmd(pos: u64, count: usize, method: Method, mem: Option<u64>) -> Result<ExitCode, String> {
    if pos == 0 {
        return Err("position is 1-based; use pos >= 1".to_string());
    }
    let start = std::time::Instant::now();
    let s = digits_with(method, pos - 1, count, mem);
    let ms = start.elapsed().as_millis();
    match peak_rss_kb() {
        Some(kb) => eprintln!("digit {pos} (+{count}, {method:?}): {ms} ms, peak RSS {kb} KiB"),
        None => eprintln!("digit {pos} (+{count}, {method:?}): {ms} ms"),
    }
    #[cfg(feature = "nthdigit-profile")]
    {
        let (factor_ns, loop_ns) = pihunt::nthdigit::profile_totals_ns();
        eprintln!(
            "  [profile] factor_segment: {:.1} ms (summed across threads), binomial loop: {:.1} ms",
            factor_ns as f64 / 1e6,
            loop_ns as f64 / 1e6
        );
    }
    println!("{s}");
    Ok(ExitCode::SUCCESS)
}

fn stream_cmd(
    from: u64,
    block: usize,
    blocks: Option<u64>,
    method: Method,
    mem: Option<u64>,
) -> Result<ExitCode, String> {
    use std::io::Write;
    if from == 0 {
        return Err("position is 1-based; use --from >= 1".to_string());
    }
    if block == 0 {
        return Err("--block must be >= 1".to_string());
    }
    let mut pos = from;
    let mut done = 0u64;
    loop {
        let start = std::time::Instant::now();
        let s = digits_with(method, pos - 1, block, mem);
        println!("{s}");
        std::io::stdout().flush().map_err(|e| e.to_string())?;
        let ms = start.elapsed().as_millis();
        let end = pos + block as u64 - 1;
        match peak_rss_kb() {
            Some(kb) => eprintln!("[block {done}] pos {pos}..{end} ({ms} ms, peak RSS {kb} KiB)"),
            None => eprintln!("[block {done}] pos {pos}..{end} ({ms} ms)"),
        }
        pos += block as u64;
        done += 1;
        if blocks.is_some_and(|limit| done >= limit) {
            break;
        }
    }
    Ok(ExitCode::SUCCESS)
}

fn summarise(records: &[Record]) {
    use Kind::*;
    println!("\n== summary ==");
    for kind in [
        Hit,
        Excluded,
        Inconclusive,
        Skipped,
        Suspicious,
        Spurious,
        Junk,
    ] {
        let n = records.iter().filter(|r| r.outcome == kind).count();
        if n > 0 {
            println!("{kind:?}: {n}");
        }
    }
    for rec in records.iter().filter(|r| r.outcome == Hit) {
        let tag = rec.tag.as_deref().unwrap_or("?");
        let banner = if tag == "NEW" {
            "!!!!!!!!!! NEW HIT !!!!!!!!!!\n"
        } else {
            ""
        };
        println!(
            "{banner}hit [{tag}] {}: {}",
            describe(rec),
            formula(&rec.columns, rec.relation.as_deref().unwrap_or_default())
        );
    }
}
