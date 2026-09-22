//! Batch file parsing and validation.

use crate::basis::Extra;
use crate::pslq::{RelationFinder, classic::ClassicPslq, multilevel::MultilevelPslq};
use serde::Deserialize;
use std::path::{Path, PathBuf};

#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Batch {
    pub name: String,
    pub output: PathBuf,
    /// Defaults to all cores.
    pub threads: Option<usize>,
    pub defaults: Defaults,
    pub search: Search,
    pub sample: Option<Sample>,
}

#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Defaults {
    pub coeff_bound: u64,
    #[serde(default)]
    pub precision_digits: Precision,
    #[serde(default = "default_gamma")]
    pub gamma: f64,
    #[serde(default = "default_max_iterations")]
    pub max_iterations: u64,
    #[serde(default = "default_max_columns")]
    pub max_columns: usize,
    #[serde(default)]
    pub finder: FinderName,
    /// Inconclusive jobs are retried at 2x, 4x, ... digits up to this many times.
    #[serde(default = "default_escalate")]
    pub escalate: u32,
}

#[derive(Debug, Clone, Copy, PartialEq, Default, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum FinderName {
    #[default]
    Classic,
    Multilevel,
}

impl FinderName {
    pub fn finder(self) -> &'static dyn RelationFinder {
        match self {
            FinderName::Classic => &ClassicPslq,
            FinderName::Multilevel => &MultilevelPslq,
        }
    }
}

fn default_gamma() -> f64 {
    1.16
}
fn default_max_iterations() -> u64 {
    1_000_000
}
fn default_max_columns() -> usize {
    80
}
fn default_escalate() -> u32 {
    2
}

/// `"auto"` or a fixed number of decimal digits.
#[derive(Debug, Clone, Copy, PartialEq, Default, Deserialize)]
#[serde(try_from = "toml::Value")]
pub enum Precision {
    #[default]
    Auto,
    Fixed(u32),
}

impl TryFrom<toml::Value> for Precision {
    type Error = String;
    fn try_from(v: toml::Value) -> Result<Self, String> {
        match v {
            toml::Value::String(s) if s == "auto" => Ok(Precision::Auto),
            toml::Value::Integer(i) if i > 0 && i <= u32::MAX as i64 => {
                Ok(Precision::Fixed(i as u32))
            }
            other => Err(format!(
                "precision_digits must be \"auto\" or a positive integer, got {other}"
            )),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Mode {
    Grid,
    Sample,
}

#[derive(Debug, Clone, Copy, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Range {
    pub from: u32,
    pub to: u32,
}

#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Search {
    pub mode: Mode,
    pub bases: Vec<u32>,
    pub periods: Range,
    /// Each entry is an inclusive degree range [s_lo, s_hi].
    pub degrees: Vec<[u32; 2]>,
    /// Each entry is one set of extras; use `[[]]` for none.
    pub extras: Vec<Vec<Extra>>,
}

#[derive(Debug, Clone, Copy, PartialEq, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Method {
    Sobol,
    Lhs,
}

#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Sample {
    pub method: Method,
    pub count: u32,
    pub seed: u32,
}

pub fn load(path: &Path) -> Result<Batch, String> {
    let text = std::fs::read_to_string(path).map_err(|e| format!("{}: {e}", path.display()))?;
    parse(&text)
}

pub fn parse(text: &str) -> Result<Batch, String> {
    let batch: Batch = toml::from_str(text).map_err(|e| e.to_string())?;
    validate(&batch)?;
    Ok(batch)
}

fn validate(b: &Batch) -> Result<(), String> {
    let d = &b.defaults;
    let s = &b.search;
    if d.gamma <= (4.0f64 / 3.0).sqrt() {
        return Err(format!(
            "gamma must be > sqrt(4/3) ≈ 1.1547, got {}",
            d.gamma
        ));
    }
    if d.escalate > 4 {
        return Err(format!(
            "escalate must be <= 4 (16x digits), got {}",
            d.escalate
        ));
    }
    if d.coeff_bound < 2 {
        return Err("coeff_bound must be >= 2".into());
    }
    if let Precision::Fixed(p) = d.precision_digits
        && p < 60
    {
        return Err(format!("precision_digits must be >= 60, got {p}"));
    }
    if s.bases.is_empty() || s.bases.iter().any(|&b| b < 2) {
        return Err("bases must be non-empty and every base >= 2".into());
    }
    if s.periods.from < 1 || s.periods.from > s.periods.to {
        return Err(format!(
            "periods must satisfy 1 <= from <= to, got {:?}",
            s.periods
        ));
    }
    if s.degrees.is_empty() || s.degrees.iter().any(|[lo, hi]| *lo < 1 || lo > hi) {
        return Err(
            "degrees must be non-empty and every [lo, hi] must satisfy 1 <= lo <= hi".into(),
        );
    }
    if s.extras.is_empty() {
        return Err("extras must be non-empty; use [[]] for no extras".into());
    }
    if b.threads == Some(0) {
        return Err("threads must be >= 1".into());
    }
    match (s.mode, &b.sample) {
        (Mode::Sample, None) => return Err("mode = \"sample\" needs a [sample] section".into()),
        (Mode::Sample, Some(smp)) if smp.count == 0 || smp.count > 65_536 => {
            return Err("sample.count must be in 1..=65536".into());
        }
        _ => {}
    }
    Ok(())
}
