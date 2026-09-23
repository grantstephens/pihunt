//! Splits a batch's jobs across `N` machines by hashing each job's base ID.
//!
//! Shard assignment must survive edits to the batch file that reorder jobs, so it is keyed on
//! `Job::id`, not position: the first 16 hex chars of the base job's ID (escalation attempt 0,
//! i.e. the un-escalated `Job` as produced by `plan::plan`) parsed as `u64`, reduced mod `N`.
//! Using the base ID keeps a whole escalation chain in one shard, since `run_chain` escalates
//! internally after a job has already been assigned.

use crate::plan::Job;
use std::path::{Path, PathBuf};
use std::str::FromStr;

/// Shard `k` of `n`, both 1-based (`1 <= k <= n`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Shard {
    pub k: u32,
    pub n: u32,
}

impl Shard {
    /// True if `job` (the base, un-escalated job) belongs to this shard.
    pub fn owns(&self, job: &Job, finder: &str) -> bool {
        let id = job.id(finder);
        let hash = u64::from_str_radix(&id[..16], 16).expect("job id is hex");
        hash % self.n as u64 == (self.k - 1) as u64
    }

    /// `<stem>.shard-k-of-n.<ext>` next to `base`, so shards never share a log file.
    pub fn output_path(&self, base: &Path) -> PathBuf {
        let dir = base.parent().unwrap_or_else(|| Path::new(""));
        let stem = base.file_stem().unwrap_or_default();
        let mut name = stem.to_os_string();
        name.push(format!(".shard-{}-of-{}", self.k, self.n));
        if let Some(ext) = base.extension() {
            name.push(".");
            name.push(ext);
        }
        dir.join(name)
    }
}

impl FromStr for Shard {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let (k, n) = s
            .split_once('/')
            .ok_or_else(|| format!("shard must be \"k/N\", got {s:?}"))?;
        let k: u32 = k
            .parse()
            .map_err(|_| format!("shard k must be a positive integer, got {k:?}"))?;
        let n: u32 = n
            .parse()
            .map_err(|_| format!("shard N must be a positive integer, got {n:?}"))?;
        if n == 0 {
            return Err("shard N must be >= 1".to_string());
        }
        if k == 0 || k > n {
            return Err(format!("shard k must satisfy 1 <= k <= N, got {k}/{n}"));
        }
        Ok(Shard { k, n })
    }
}
