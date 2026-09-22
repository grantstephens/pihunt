use clap::{Parser, Subcommand};
use pihunt::basis::{Columns, Extra, Shape};
use pihunt::log::{Kind, Record};
use pihunt::report::{describe, formula};
use pihunt::runner::{Ctx, load_done, run_chain};
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

#[derive(Subcommand)]
enum Cmd {
    /// Run every job in a batch file, appending results to its output log.
    Run { batch: PathBuf },
    /// Dry run: show how big a batch is without running it.
    Plan { batch: PathBuf },
    /// Re-verify every hit in a results log at 2x precision.
    Verify { results: PathBuf },
    /// Markdown exclusion report over one or more results logs, to stdout.
    Report {
        #[arg(required = true)]
        results: Vec<PathBuf>,
    },
}

fn main() -> ExitCode {
    let result = match Cli::parse().cmd {
        Cmd::Run { batch } => run(batch),
        Cmd::Plan { batch } => show_plan(batch),
        Cmd::Verify { results } => reverify(results),
        Cmd::Report { results } => write_report(results),
    };
    result.unwrap_or_else(|e| {
        eprintln!("error: {e}");
        ExitCode::FAILURE
    })
}

fn run(path: PathBuf) -> Result<ExitCode, String> {
    let batch = config::load(&path)?;
    let jobs = plan::plan(&batch);
    let threads = batch
        .threads
        .unwrap_or_else(|| std::thread::available_parallelism().map_or(1, |n| n.get()));
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .map_err(|e| e.to_string())?;
    let done = load_done(&batch.output)?;
    let (tx, writer) = log::spawn_writer(&batch.output).map_err(|e| e.to_string())?;
    let ctx = Ctx {
        batch: &batch.name,
        finder: batch.defaults.finder.finder(),
        max_columns: batch.defaults.max_columns,
        escalate: batch.defaults.escalate,
    };
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
        batch.output.display(),
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

fn show_plan(path: PathBuf) -> Result<ExitCode, String> {
    let batch = config::load(&path)?;
    let jobs = plan::plan(&batch);
    let max_cols = batch.defaults.max_columns;
    let (skipped, run): (Vec<_>, Vec<_>) = jobs.iter().partition(|j| j.shape.columns() > max_cols);
    println!(
        "{}: {} jobs ({} to run, {} skipped over max_columns {max_cols}), {} finder, escalate up to {}x digits",
        batch.name,
        jobs.len(),
        run.len(),
        skipped.len(),
        batch.defaults.finder.finder().name(),
        1u32 << batch.defaults.escalate
    );
    if let (Some(n), Some(lo), Some(hi)) = (
        run.iter().map(|j| j.shape.columns()).max(),
        run.iter().map(|j| j.digits).min(),
        run.iter().map(|j| j.digits).max(),
    ) {
        println!("largest n = {n}, precision {lo}..={hi} digits");
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
