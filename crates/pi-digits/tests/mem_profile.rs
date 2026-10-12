//! `install_gmp_hook` closes a real blind spot (see `src/mem_profile.rs`'s module docs): GMP's
//! own `malloc` calls, made by `rug::Integer`/`Float` under the `gmp` backend, never touch
//! Rust's `#[global_allocator]` unless hooked via `mp_set_memory_functions`.
//!
//! Needs its own process (a separate integration test file, not a `#[cfg(test)]` module in
//! `src/`): the tracked counter is process-global, and `cargo test`'s unit-test binary runs this
//! crate's many other `rug::Integer`-using tests concurrently in one process, which pollutes any
//! before/after delta this test would otherwise check. A separate file gets its own test binary
//! (a separate process), so it's the only thing touching the counter.
//!
//! Gated on `mem-profile`+`gmp` so a default `cargo test` (neither enabled) skips this file
//! instead of running it against meaningless no-op stubs.
#![cfg(all(feature = "mem-profile", feature = "gmp"))]

use pi_digits::mem_profile::{install_gmp_hook, live_bytes, reset_peak};
use rug::Integer;

#[test]
fn gmp_allocations_are_invisible_until_hooked_then_tracked() {
    const SIGNAL_BITS: u32 = 2_000_000; // ~250 KB

    let mut untracked = Integer::from(1);
    reset_peak();
    let before_hook = live_bytes();
    untracked <<= SIGNAL_BITS;
    assert_eq!(
        live_bytes(),
        before_hook,
        "GMP allocation moved the counter before the hook was installed"
    );
    drop(untracked);

    install_gmp_hook();
    reset_peak();
    let before = live_bytes();
    let mut tracked = Integer::from(1);
    tracked <<= SIGNAL_BITS;
    let after = live_bytes();
    assert!(
        after >= before + 200_000,
        "expected the ~250 KB allocation to be tracked: before={before} after={after}"
    );
    drop(tracked);
    assert!(
        live_bytes() < after,
        "freeing the bignum should reduce live bytes"
    );
}
