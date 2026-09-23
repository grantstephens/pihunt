//! Sequential fallback machinery for the crate's `#[cfg(feature = "parallel")]` rayon call
//! sites, so `nthdigit.rs` and `nthdigit2.rs` keep exactly one code path per call site instead
//! of duplicating whole functions for the non-parallel case (no threads on WASM; also usable
//! single-threaded on native). Every macro here picks between a rayon call and its plain
//! `std::iter` equivalent behind the same `parallel` feature gate `bignum.rs` uses for `gmp`
//! vs `pure`.
//!
//! The certification-sensitive accumulators these wrap (`u128::wrapping_add`, term-count
//! sums, `HashMap::extend`) are all associative and commutative, so a sequential single-pass
//! `fold` and a parallel tree of per-chunk `fold`s-then-`reduce` compute the identical result;
//! only the term *count* need match, which it does since both paths visit every item exactly
//! once.

/// `into_par_iter()` when `parallel` is on, plain `into_iter()` otherwise.
#[cfg(feature = "parallel")]
macro_rules! maybe_into_par_iter {
    ($e:expr) => {
        rayon::iter::IntoParallelIterator::into_par_iter($e)
    };
}
#[cfg(not(feature = "parallel"))]
macro_rules! maybe_into_par_iter {
    ($e:expr) => {
        ($e).into_iter()
    };
}

/// `par_iter()` when `parallel` is on, plain `iter()` otherwise.
#[cfg(feature = "parallel")]
macro_rules! maybe_par_iter {
    ($e:expr) => {
        rayon::iter::IntoParallelRefIterator::par_iter($e)
    };
}
#[cfg(not(feature = "parallel"))]
macro_rules! maybe_par_iter {
    ($e:expr) => {
        ($e).iter()
    };
}

/// `rayon::join(a, b)` when `parallel` is on, sequential `(a(), b())` otherwise.
#[cfg(feature = "parallel")]
macro_rules! maybe_join {
    ($a:expr, $b:expr) => {
        rayon::join($a, $b)
    };
}
#[cfg(not(feature = "parallel"))]
macro_rules! maybe_join {
    ($a:expr, $b:expr) => {
        (($a)(), ($b)())
    };
}

/// Closes a `.map(...)` chain: `.reduce(identity, op)` when `parallel` is on (rayon combines
/// per-thread partials with `op`), `.fold(identity(), op)` otherwise (one running total). Same
/// identity closure and combining op either way.
#[cfg(feature = "parallel")]
macro_rules! maybe_reduce {
    ($iter:expr, $identity:expr, $op:expr) => {
        $iter.reduce($identity, $op)
    };
}
#[cfg(not(feature = "parallel"))]
macro_rules! maybe_reduce {
    ($iter:expr, $identity:expr, $op:expr) => {
        $iter.fold(($identity)(), $op)
    };
}

/// Rayon's current worker-thread count (`parallel` on), or `1` (no threads otherwise).
#[cfg(feature = "parallel")]
pub(crate) fn current_threads() -> u64 {
    rayon::current_num_threads().max(1) as u64
}
#[cfg(not(feature = "parallel"))]
pub(crate) fn current_threads() -> u64 {
    1
}

pub(crate) use maybe_into_par_iter;
pub(crate) use maybe_join;
pub(crate) use maybe_par_iter;
pub(crate) use maybe_reduce;
