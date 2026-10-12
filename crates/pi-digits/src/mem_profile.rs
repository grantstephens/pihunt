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
//! working sets, glibc arena overhead, etc) — **except `rug::Integer`/`Float` bignums under the
//! `gmp` backend**, which were invisible until [`install_gmp_hook`]: this workspace links
//! *system* GMP (see the root `Cargo.toml`'s comment on `gmp-mpfr-sys`), so GMP's own `malloc`
//! calls never go through Rust's `#[global_allocator]`. [`install_gmp_hook`] routes them through
//! the same counter via `mp_set_memory_functions`. Measured gap this closed at n=1e7: process RSS
//! 19.9 MiB vs. 12.5 MiB tracked before the hook — see `nthdigit2.rs`'s "Memory" section.

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

    #[cfg(feature = "gmp")]
    mod gmp_hook {
        use super::bump;
        use gmp_mpfr_sys::gmp;
        use std::alloc::{GlobalAlloc, Layout, System, handle_alloc_error};
        use std::ffi::c_void;
        use std::sync::Once;

        // GMP doesn't report the alignment it needs; 16 matches what `malloc` guarantees on
        // every platform this project targets, which is what GMP's own default allocator uses.
        const ALIGN: usize = 16;

        fn layout(size: usize) -> Layout {
            Layout::from_size_align(size, ALIGN).expect("GMP allocation size overflows isize")
        }

        extern "C" fn alloc(size: usize) -> *mut c_void {
            let l = layout(size);
            let ptr = unsafe { System.alloc(l) };
            if ptr.is_null() {
                handle_alloc_error(l);
            }
            bump(size as i64);
            ptr as *mut c_void
        }

        unsafe extern "C" fn realloc(
            ptr: *mut c_void,
            old_size: usize,
            new_size: usize,
        ) -> *mut c_void {
            let new_ptr = unsafe { System.realloc(ptr as *mut u8, layout(old_size), new_size) };
            if new_ptr.is_null() {
                handle_alloc_error(layout(new_size));
            }
            bump(new_size as i64 - old_size as i64);
            new_ptr as *mut c_void
        }

        unsafe extern "C" fn free(ptr: *mut c_void, size: usize) {
            unsafe { System.dealloc(ptr as *mut u8, layout(size)) };
            bump(-(size as i64));
        }

        static INSTALLED: Once = Once::new();

        /// Routes GMP's own allocations through this module's counter — see the module docs'
        /// note on why `rug::Integer`/`Float` were invisible without this. Idempotent; call
        /// before the computation you want to profile starts (allocations made before the
        /// call aren't retroactively counted, but anything live at call time keeps whatever
        /// size GMP last told its *previous* allocator, which is fine for a peak-since-start
        /// measurement taken from near the start of a run).
        pub fn install() {
            INSTALLED.call_once(|| unsafe {
                gmp::set_memory_functions(Some(alloc), Some(realloc), Some(free));
            });
        }
    }

    #[cfg(feature = "gmp")]
    pub use gmp_hook::install as install_gmp_hook;
    #[cfg(not(feature = "gmp"))]
    pub fn install_gmp_hook() {}
}

#[cfg(feature = "mem-profile")]
pub use imp::{checkpoint, install_gmp_hook, live_bytes, peak_bytes, reset_peak};

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
#[cfg(not(feature = "mem-profile"))]
pub fn install_gmp_hook() {}
