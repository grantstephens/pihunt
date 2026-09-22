use pihunt::basis::Extra;
use pihunt::config::{Method, Mode, Precision, parse};

const GOOD: &str = r#"
name   = "t"
output = "results/t.jsonl"

[defaults]
coeff_bound = 1000

[search]
mode    = "grid"
bases   = [10, 16]
periods = { from = 2, to = 4 }
degrees = [[1, 1], [1, 2]]
extras  = [[], ["log2", "log5"]]
"#;

#[test]
fn parses_with_defaults() {
    let b = parse(GOOD).unwrap();
    assert_eq!(b.defaults.gamma, 1.16);
    assert_eq!(b.defaults.max_iterations, 1_000_000);
    assert_eq!(b.defaults.max_columns, 80);
    assert_eq!(b.defaults.precision_digits, Precision::Auto);
    assert_eq!(b.search.mode, Mode::Grid);
    assert_eq!(b.search.extras[1], vec![Extra::Log2, Extra::Log5]);
    assert_eq!(b.threads, None);
}

#[test]
fn parses_fixed_precision_and_sample() {
    let text = GOOD
        .replace(
            "coeff_bound = 1000",
            "coeff_bound = 1000\nprecision_digits = 300",
        )
        .replace("mode    = \"grid\"", "mode    = \"sample\"")
        + "\n[sample]\nmethod = \"lhs\"\ncount = 10\nseed = 7\n";
    let b = parse(&text).unwrap();
    assert_eq!(b.defaults.precision_digits, Precision::Fixed(300));
    assert_eq!(b.sample.unwrap().method, Method::Lhs);
}

#[test]
fn rejects_bad_configs() {
    let cases = [
        GOOD.replace("coeff_bound = 1000", "coeff_bound = 1000\ngamma = 1.1"),
        GOOD.replace("\"log5\"", "\"log7\""),
        GOOD.replace("from = 2, to = 4", "from = 5, to = 4"),
        GOOD.replace("mode    = \"grid\"", "mode    = \"sample\""),
        GOOD.replace(
            "coeff_bound = 1000",
            "coeff_bound = 1000\nprecision_digits = 30",
        ),
        GOOD.replace(
            "coeff_bound = 1000",
            "coeff_bound = 1000\nprecision_digits = \"lots\"",
        ),
        GOOD.replace("[[1, 1], [1, 2]]", "[[2, 1]]"),
        GOOD.replace("bases   = [10, 16]", "bases   = []"),
        GOOD.replace("name   = \"t\"", "name   = \"t\"\ntypo = 1"),
    ];
    for (i, c) in cases.iter().enumerate() {
        assert!(parse(c).is_err(), "case {i} should fail:\n{c}");
    }
}
