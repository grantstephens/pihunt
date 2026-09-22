//! Timing baseline, classic vs multilevel. Run with:
//! cargo test --release --test timing -- --ignored --nocapture

use pihunt::basis::{Extra, Shape, auto_digits};
use pihunt::job::run_job;
use pihunt::plan::Job;
use pihunt::pslq::{RelationFinder, classic::ClassicPslq, multilevel::MultilevelPslq};

/// Classic gets skipped above this many columns: it takes minutes to hours there.
const CLASSIC_MAX_N: usize = 47;

#[test]
#[ignore]
fn timing_baseline() {
    let all = [
        Extra::Pi2,
        Extra::Log2,
        Extra::Log3,
        Extra::Log5,
        Extra::Catalan,
        Extra::Zeta3,
    ];
    let shapes = [
        Shape::new(10, 8, 1, 1, vec![]),
        Shape::new(10, 12, 1, 1, all.to_vec()),
        Shape::new(10, 12, 1, 2, all[..2].to_vec()),
        Shape::new(10, 16, 1, 2, all.to_vec()),
        Shape::new(10, 20, 1, 2, all.to_vec()),
        Shape::new(10, 16, 1, 3, all.to_vec()),
        Shape::new(10, 20, 1, 3, all.to_vec()),
    ];
    let finders: [&dyn RelationFinder; 2] = [&ClassicPslq, &MultilevelPslq];
    println!("| n | digits | finder | outcome | iterations | seconds |");
    println!("|---|---|---|---|---|---|");
    for shape in shapes {
        let n = shape.columns();
        let digits = auto_digits(n, 1000);
        let job = Job {
            shape,
            coeff_bound: 1000,
            digits,
            gamma: 1.16,
            max_iterations: 10_000_000,
        };
        for f in finders {
            if f.name() == "classic" && n > CLASSIC_MAX_N {
                println!("| {n} | {digits} | classic | (skipped) | | |");
                continue;
            }
            let rec = run_job(&job, "timing", f, 200);
            println!(
                "| {n} | {digits} | {} | {:?} | {} | {:.2} |",
                f.name(),
                rec.outcome,
                rec.iterations,
                rec.elapsed_ms as f64 / 1000.0
            );
        }
    }
}
