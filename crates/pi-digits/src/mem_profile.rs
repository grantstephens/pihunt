//! Ad hoc heap-accounting instrumentation for the `nthdigit2` memory-composition
//! investigation (see `src/nthdigit2.rs` module docs' "Memory" section and
//! `docs/nthdigit-theorem2.md` §7). Enable with `--features mem-profile`.
//!
//! This installs a counting global allocator (tracking live and peak bytes across the whole
//! process) and exposes `checkpoint`/`reset_peak`/`peak_bytes`/`live_bytes` for call sites to
//! log progress. When the feature is disabled, the global allocator is untouched (plain
//! system allocator) and every function here compiles to a no-op — this module is not part of
//! the shipped `nthdigit2::digits` cost model.
//!
//! A per-call `Layout` counting allocator can't distinguish "which named structure" a byte
//! belongs to, so callers that want a breakdown by structure (e.g. `nthdigit2::c_part`'s
//! `lucas_small` vs `cofactor` Vecs) additionally compute those structures' own heap footprint
//! directly (`.capacity() * size_of::<T>()`) and log that alongside a `checkpoint`. The global
//! peak still catches anything neither of us thought to size explicitly (rayon per-thread
//! working sets, `rug::Integer` bignums, glibc arena overhead, etc).

#[cfg(feature = "mem-profile")]
mod imp {
    use std::alloc::{GlobalAlloc, Layout, System};
    use std::sync::atomic::{AtomicI64, AtomicUsize, Ordering};

    static CURRENT: AtomicI64 = AtomicI64::new(0);
    static PEAK: AtomicUsize = AtomicUsize::new(0);

    fn bump(delta: i64) {
        let cur = CURRENT.fetch_add(delta, Ordering::Relaxed) + delta;
        if cur > 0 {
            PEAK.fetch_max(cur as usize, Ordering::Relaxed);
        }
    }

    pub struct CountingAlloc;

    unsafe impl GlobalAlloc for CountingAlloc {
        unsafe fn alloc(&self, layout: Layout) -> *mut u8 {
            let ptr = unsafe { System.alloc(layout) };
            if !ptr.is_null() {
                bump(layout.size() as i64);
            }
            ptr
        }
        unsafe fn dealloc(&self, ptr: *mut u8, layout: Layout) {
            unsafe { System.dealloc(ptr, layout) };
            bump(-(layout.size() as i64));
        }
        unsafe fn alloc_zeroed(&self, layout: Layout) -> *mut u8 {
            let ptr = unsafe { System.alloc_zeroed(layout) };
            if !ptr.is_null() {
                bump(layout.size() as i64);
            }
            ptr
        }
        unsafe fn realloc(&self, ptr: *mut u8, layout: Layout, new_size: usize) -> *mut u8 {
            let new_ptr = unsafe { System.realloc(ptr, layout, new_size) };
            if !new_ptr.is_null() {
                bump(new_size as i64 - layout.size() as i64);
            }
            new_ptr
        }
    }

    #[global_allocator]
    static GLOBAL: CountingAlloc = CountingAlloc;

    /// Logs a labelled live/peak-so-far snapshot to stderr.
    pub fn checkpoint(label: &str) {
        let cur = CURRENT.load(Ordering::Relaxed).max(0) as usize;
        let peak = PEAK.load(Ordering::Relaxed);
        eprintln!(
            "[mem-profile] {label}: live={:.2} MiB peak-so-far={:.2} MiB",
            cur as f64 / (1024.0 * 1024.0),
            peak as f64 / (1024.0 * 1024.0)
        );
    }

    /// Resets the running peak to the current live size (so a later `peak_bytes()` reports the
    /// peak *since this call*, e.g. to isolate one phase's contribution).
    pub fn reset_peak() {
        let cur = CURRENT.load(Ordering::Relaxed).max(0) as usize;
        PEAK.store(cur, Ordering::Relaxed);
    }

    pub fn peak_bytes() -> usize {
        PEAK.load(Ordering::Relaxed)
    }

    pub fn live_bytes() -> usize {
        CURRENT.load(Ordering::Relaxed).max(0) as usize
    }
}

#[cfg(feature = "mem-profile")]
pub use imp::{checkpoint, live_bytes, peak_bytes, reset_peak};

#[cfg(not(feature = "mem-profile"))]
pub fn checkpoint(_label: &str) {}
#[cfg(not(feature = "mem-profile"))]
pub fn reset_peak() {}
#[cfg(not(feature = "mem-profile"))]
pub fn peak_bytes() -> usize {
    0
}
#[cfg(not(feature = "mem-profile"))]
pub fn live_bytes() -> usize {
    0
}
