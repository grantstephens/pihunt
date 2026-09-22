//! Append-only JSONL results log.

use serde::{Deserialize, Serialize};
use std::fs::OpenOptions;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;
use std::sync::mpsc::{Sender, channel};
use std::thread::JoinHandle;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum Kind {
    Hit,
    Junk,
    Suspicious,
    Spurious,
    Excluded,
    Inconclusive,
    Skipped,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Params {
    pub base: u32,
    pub period: u32,
    pub degrees: [u32; 2],
    pub extras: Vec<String>,
    pub coeff_bound: u64,
    pub precision_digits: u32,
    pub gamma: f64,
    pub max_iterations: u64,
    pub finder: String,
    pub algo_version: u32,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct DroppedRecord {
    pub column: String,
    pub relation: Vec<String>,
    pub tag: Option<String>,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Verify {
    pub passed: bool,
    pub residual_log10: f64,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Record {
    pub job_id: String,
    pub batch: String,
    pub pihunt_version: String,
    pub started: String,
    pub elapsed_ms: u64,
    pub params: Params,
    pub columns: Vec<String>,
    pub dropped: Vec<DroppedRecord>,
    pub outcome: Kind,
    /// Exclusion bound reached, formatted like "1.2e5". None for skipped jobs.
    pub bound: Option<String>,
    pub iterations: u64,
    /// Aligned with `columns`; present for hit/junk/suspicious/spurious.
    pub relation: Option<Vec<String>>,
    pub verify: Option<Verify>,
    /// "known" or "NEW" for hits.
    pub tag: Option<String>,
    /// Why a job was inconclusive or skipped.
    pub note: Option<String>,
    /// Job ID of the inconclusive attempt this one retries at higher precision.
    #[serde(default)]
    pub escalated_from: Option<String>,
}

/// Spawn the single writer thread. Drop the sender to finish; join to surface IO errors.
pub fn spawn_writer(
    path: &Path,
) -> std::io::Result<(Sender<Record>, JoinHandle<std::io::Result<()>>)> {
    if let Some(dir) = path.parent().filter(|d| !d.as_os_str().is_empty()) {
        std::fs::create_dir_all(dir)?;
    }
    let mut file = OpenOptions::new().create(true).append(true).open(path)?;
    let (tx, rx) = channel::<Record>();
    let handle = std::thread::spawn(move || {
        for rec in rx {
            let line = serde_json::to_string(&rec).expect("records always serialise");
            writeln!(file, "{line}")?;
            file.flush()?;
        }
        Ok(())
    });
    Ok((tx, handle))
}

pub fn read(path: &Path) -> Result<Vec<Record>, String> {
    let file = std::fs::File::open(path).map_err(|e| format!("{}: {e}", path.display()))?;
    BufReader::new(file)
        .lines()
        .enumerate()
        .filter(|(_, l)| !matches!(l, Ok(s) if s.trim().is_empty()))
        .map(|(i, l)| {
            let l = l.map_err(|e| e.to_string())?;
            serde_json::from_str(&l).map_err(|e| format!("{}:{}: {e}", path.display(), i + 1))
        })
        .collect()
}
