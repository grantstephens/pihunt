use pihunt::basis::{Extra, Shape};
use pihunt::config::parse;
use pihunt::plan::{Job, plan};

const GRID: &str = r#"
name   = "t"
output = "results/t.jsonl"

[defaults]
coeff_bound = 1000

[search]
mode    = "grid"
bases   = [10, 16]
periods = { from = 2, to = 4 }
degrees = [[1, 1], [1, 2]]
extras  = [[], ["log2", "log5"], ["log5", "log2"]]
"#;

fn sample(method: &str, seed: u32) -> String {
    GRID.replace("mode    = \"grid\"", "mode    = \"sample\"")
        + &format!("\n[sample]\nmethod = \"{method}\"\ncount = 20\nseed = {seed}\n")
}

#[test]
fn grid_is_cartesian_and_deduplicated() {
    // Two extras sets are the same after sorting, so 2 * 3 * 2 * 2 unique jobs.
    let jobs = plan(&parse(GRID).unwrap());
    assert_eq!(jobs.len(), 24);
    assert_eq!(jobs[0].shape, Shape::new(10, 2, 1, 1, vec![]));
    assert_eq!(jobs[0].digits, 50 + (3.0f64 * 3.0 * 1.5).ceil() as u32);
}

#[test]
fn sample_modes_are_deterministic_and_in_range() {
    for method in ["sobol", "lhs"] {
        let a = plan(&parse(&sample(method, 1)).unwrap());
        let b = plan(&parse(&sample(method, 1)).unwrap());
        assert_eq!(a, b, "{method}");
        assert!(!a.is_empty() && a.len() <= 20);
        for j in &a {
            assert!([10, 16].contains(&j.shape.base));
            assert!((2..=4).contains(&j.shape.period));
        }
    }
}

fn job(extras: Vec<Extra>) -> Job {
    Job {
        shape: Shape::new(10, 4, 1, 2, extras),
        coeff_bound: 1000,
        digits: 120,
        gamma: 1.16,
        max_iterations: 10,
    }
}

#[test]
fn job_ids_are_stable_and_precise() {
    let a = job(vec![Extra::Log2, Extra::Log5]);
    let b = job(vec![Extra::Log5, Extra::Log2]);
    assert_eq!(a.id("classic"), b.id("classic"));
    assert_eq!(a.id("classic").len(), 32);
    assert_ne!(a.id("classic"), a.id("multilevel"));
    let variants = [
        Job {
            digits: 121,
            ..a.clone()
        },
        Job {
            coeff_bound: 999,
            ..a.clone()
        },
        Job {
            gamma: 1.2,
            ..a.clone()
        },
        Job {
            max_iterations: 11,
            ..a.clone()
        },
        Job {
            shape: Shape::new(100, 4, 1, 2, vec![Extra::Log2, Extra::Log5]),
            ..a.clone()
        },
        Job {
            shape: Shape::new(10, 5, 1, 2, vec![Extra::Log2, Extra::Log5]),
            ..a.clone()
        },
        Job {
            shape: Shape::new(10, 4, 1, 1, vec![Extra::Log2, Extra::Log5]),
            ..a.clone()
        },
        Job {
            shape: Shape::new(10, 4, 1, 2, vec![Extra::Log2]),
            ..a.clone()
        },
    ];
    for v in variants {
        assert_ne!(a.id("classic"), v.id("classic"), "{v:?}");
    }
}
