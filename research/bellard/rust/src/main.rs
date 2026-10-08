//! Rust port of Bellard's 1997 `pi.c` and Feb 1997 `pi1.c` (n-th decimal digit of pi, O(n^2)
//! time, tiny memory).
//!
//! usage: bellard-pi <n> [--formula pi|pi1|pi1-simd4] [--variant naive|fast] [--threads T]
//!
//! * `--formula pi` (default): Bellard's original Jan 1997 formula (`pi.c`).
//! * `--formula pi1`: the Feb 1997 variant using Gosper's series (`pi1.c`), ~4x faster.
//! * `--variant naive`: line-for-line port (u64 `%` mulmod, one extended-Euclid inverse per k).
//! * `--variant fast` (default): Montgomery u32 arithmetic, batched modular inversion
//!   (Montgomery's trick, one Euclid per block instead of per k), sieved primes, one thread per
//!   core over the primes. For `pi1`, the single prime `a == 2` can't use Montgomery arithmetic
//!   (it needs an odd modulus) and silently falls back to the naive per-k computation just for
//!   that one prime -- negligible cost, see [`prime_frac_fast_gosper`]'s docs.
//!
//! All print the same digits as the corresponding C program (position n has the same meaning:
//! digits n.. of pi's fractional part, i.e. `pi.c 1` / `pi1.c 1` both print 141592653).
use std::sync::atomic::{AtomicUsize, Ordering};

fn inv_mod(x: i64, y: i64) -> i64 {
    let (mut u, mut v, mut c, mut a) = (x, y, 1i64, 0i64);
    loop {
        let q = v / u;
        let t = c;
        c = a - q * c;
        a = t;
        let t = u;
        u = v - q * u;
        v = t;
        if u == 0 {
            break;
        }
    }
    a %= y;
    if a < 0 { a + y } else { a }
}

fn mul_mod(a: u64, b: u64, m: u64) -> u64 {
    a * b % m
}

fn pow_mod(mut a: u64, mut b: u64, m: u64) -> u64 {
    let mut r = 1 % m;
    loop {
        if b & 1 == 1 {
            r = mul_mod(r, a, m);
        }
        b >>= 1;
        if b == 0 {
            return r;
        }
        a = mul_mod(a, a, m);
    }
}

fn primes_upto(limit: usize) -> Vec<u32> {
    let mut sieve = vec![true; limit + 1];
    let mut out = Vec::new();
    for i in 2..=limit {
        if sieve[i] {
            out.push(i as u32);
            let mut j = i * i;
            while j <= limit {
                sieve[j] = false;
                j += i;
            }
        }
    }
    out
}

/// Same `vmax`/`av` as pi.c.
fn prime_power(a: u32, n2: u32) -> (u32, u32) {
    let vmax = ((n2 as f64).ln() / (a as f64).ln()) as u32;
    let mut av = 1u32;
    for _ in 0..vmax {
        av *= a;
    }
    (vmax, av)
}

/// Strip factors of `a` from the running k and (2k-1) terms; shared by both variants.
struct Stepper {
    a: u32,
    v: i32,
    kq: u32,
    kq2: u32,
}

impl Stepper {
    fn new(a: u32) -> Self {
        Stepper {
            a,
            v: 0,
            kq: 1,
            kq2: 1,
        }
    }
    #[inline(always)]
    fn step(&mut self, k: u32) -> (u32, u32) {
        let a = self.a;
        let mut t = k;
        if self.kq >= a {
            loop {
                t /= a;
                self.v -= 1;
                if t % a != 0 {
                    break;
                }
            }
            self.kq = 0;
        }
        self.kq += 1;
        let tn = t;

        let mut t = 2 * k - 1;
        if self.kq2 >= a {
            if self.kq2 == a {
                loop {
                    t /= a;
                    self.v += 1;
                    if t % a != 0 {
                        break;
                    }
                }
            }
            self.kq2 -= a;
        }
        self.kq2 += 2;
        (tn, t)
    }
}

fn prime_frac_naive(a: u32, n_pos: u32, big_n: u32) -> f64 {
    let (vmax, av) = prime_power(a, 2 * big_n);
    let (vmax, av) = (vmax as i32, av as u64);
    let (mut s, mut num, mut den) = (0u64, 1u64, 1u64);
    let mut st = Stepper::new(a);
    for k in 1..=big_n {
        let (tn, td) = st.step(k);
        num = mul_mod(num, tn as u64, av);
        den = mul_mod(den, td as u64, av);
        if st.v > 0 {
            let mut t = inv_mod(den as i64, av as i64) as u64;
            t = mul_mod(t, num, av);
            t = mul_mod(t, k as u64, av);
            for _ in st.v..vmax {
                t = mul_mod(t, a as u64, av);
            }
            s += t;
            if s >= av {
                s -= av;
            }
        }
    }
    let t = pow_mod(10, n_pos as u64 - 1, av);
    let s = mul_mod(s, t, av);
    s as f64 / av as f64
}

#[derive(Clone, Copy)]
struct Mont {
    n: u32,
    ninv: u32, // -n^-1 mod 2^32
    r2: u32,   // R^2 mod n
    r3: u32,   // R^3 mod n
}

impl Mont {
    fn new(n: u32) -> Self {
        debug_assert!(n & 1 == 1 && n < (1 << 31));
        let mut x = n;
        for _ in 0..5 {
            x = x.wrapping_mul(2u32.wrapping_sub(n.wrapping_mul(x)));
        }
        let nn = n as u64;
        let r = (1u64 << 32) % nn;
        let r2 = r * r % nn;
        let r3 = r2 * r % nn;
        Mont {
            n,
            ninv: x.wrapping_neg(),
            r2: r2 as u32,
            r3: r3 as u32,
        }
    }
    /// a * b / R mod n, for a < n, b < 2^32. The conditional subtract at the end used to be a
    /// plain `if u >= n { u - n } else { u }` -- data-dependent, ~50/50, and *not* something
    /// LLVM was already turning into a branchless sequence: switching to [`reduce_once`] measured
    /// a real ~12-21% speedup (see `research/bellard/README.md`).
    #[inline(always)]
    fn mul(&self, a: u32, b: u32) -> u32 {
        let t = a as u64 * b as u64;
        let m = (t as u32).wrapping_mul(self.ninv);
        let u = ((t + m as u64 * self.n as u64) >> 32) as u32;
        reduce_once(u, self.n)
    }
}

/// Branchless `if x >= m { x - m } else { x }`, for `x < 2*m`. Used by [`Mont::mul`]'s
/// conditional subtract: that one data-dependent branch, alone, measured a real ~12-21% win
/// over the plain `if` (LLVM wasn't already generating a branchless sequence for it) --
/// see `research/bellard/README.md`.
///
/// Deliberately *not* used for the `fast` variants' batch accumulator (`s += t; if s >= av { s
/// -= av }` in `flush!`), even though it's the same pattern: tried it, measured it, and it was a
/// small but reproducible *regression* (~1.5-2%) there instead -- LLVM apparently already
/// handles that specific, simpler call site well, and the extra arithmetic this function adds
/// just cost more than the branch it replaced. Kept as a documented negative result rather than
/// applying the trick everywhere on the assumption that it must help.
#[inline(always)]
fn reduce_once(x: u32, m: u32) -> u32 {
    let diff = x.wrapping_sub(m);
    let mask = (diff as i32 >> 31) as u32;
    diff.wrapping_add(m & mask)
}

const BLK: usize = 512;

fn prime_frac_fast(a: u32, n_pos: u32, big_n: u32) -> f64 {
    let (vmax, av) = prime_power(a, 2 * big_n);
    let vmax = vmax as i32;
    let mo = Mont::new(av);
    // pa[e] = a^e * R mod av
    let r = ((1u64 << 32) % av as u64) as u64;
    let mut pa = vec![0u32; vmax as usize + 1];
    let mut p = r;
    for e in pa.iter_mut() {
        *e = p as u32;
        p = p * a as u64 % av as u64;
    }

    let (mut num, mut den) = (1u32, 1u32);
    let mut s = 0u32; // sum of t * R^-1
    let mut st = Stepper::new(a);
    let (mut ba, mut bb, mut bk, mut be) = ([0u32; BLK], [0u32; BLK], [0u32; BLK], [0u8; BLK]);
    let mut pre = [0u32; BLK];
    let mut cnt = 0usize;

    macro_rules! flush {
        () => {
            if cnt > 0 {
                pre[0] = bb[0];
                for i in 1..cnt {
                    pre[i] = mo.mul(pre[i - 1], bb[i]);
                }
                // Montgomery form of (prod b)^-1
                let inv = inv_mod(pre[cnt - 1] as i64, av as i64) as u32;
                let mut ip = mo.mul(inv, mo.r3);
                for i in (0..cnt).rev() {
                    let ib = if i == 0 { ip } else { mo.mul(ip, pre[i - 1]) };
                    if i > 0 {
                        ip = mo.mul(ip, bb[i]);
                    }
                    let mut t = mo.mul(ba[i], ib); // num/den
                    t = mo.mul(t, bk[i]); // * k / R
                    t = mo.mul(t, pa[be[i] as usize]); // * a^e
                    s += t;
                    if s >= av {
                        s -= av;
                    }
                }
                cnt = 0;
            }
        };
    }

    for k in 1..=big_n {
        let (tn, td) = st.step(k);
        num = mo.mul(num, tn);
        den = mo.mul(den, td);
        if st.v > 0 {
            ba[cnt] = num;
            bb[cnt] = mo.mul(den, mo.r2); // Montgomery form of den
            bk[cnt] = k;
            be[cnt] = (vmax - st.v).max(0) as u8;
            cnt += 1;
            if cnt == BLK {
                flush!();
            }
        }
    }
    flush!();
    let s = mo.mul(s, mo.r2); // undo the R^-1
    let t = pow_mod(10, n_pos as u64 - 1, av as u64);
    let s = mul_mod(s as u64, t, av as u64);
    s as f64 / av as f64
}

// ---------------------------------------------------------------------------------
// pi1.c (Feb 1997, Gosper's series): pi = sum_k (25k-3) / (C(3k,k) * 2^(k-1))
// ---------------------------------------------------------------------------------

/// Inverse of `u` mod odd `v`, via the binary extended-Euclidean algorithm `pi1.c` uses for
/// every prime except 2 (the general `inv_mod` above stays correct there too, but `pi1.c` uses
/// this one specifically; ported faithfully rather than simplified, since the point of the
/// `naive` variant is to mirror the C line for line). `pi1.c`'s version enters its inner loop via
/// a `goto` that skips the first halving step when `u` is odd; `enter_at_y4` reproduces that by
/// skipping exactly one iteration of the halving branch.
fn inv_mod2(u: i64, v: i64) -> i64 {
    let mut u1 = 1i64;
    let mut u3 = u;
    let mut v1 = v;
    let mut v3 = v;
    let (mut t1, mut t3);
    let mut enter_at_y4;
    if u & 1 != 0 {
        t1 = 0;
        t3 = -v;
        enter_at_y4 = true;
    } else {
        t1 = 1;
        t3 = u;
        enter_at_y4 = false;
    }
    loop {
        loop {
            if !enter_at_y4 {
                if t1 & 1 == 0 {
                    t1 >>= 1;
                    t3 >>= 1;
                } else {
                    t1 = (t1 + v) >> 1;
                    t3 >>= 1;
                }
            }
            enter_at_y4 = false;
            if t3 & 1 != 0 {
                break;
            }
        }
        if t3 >= 0 {
            u1 = t1;
            u3 = t3;
        } else {
            v1 = v - t1;
            v3 = -t3;
        }
        t1 = u1 - v1;
        t3 = u3 - v3;
        if t1 < 0 {
            t1 += v;
        }
        if t3 == 0 {
            break;
        }
    }
    u1
}

/// `pi1.c`'s per-prime fraction: four running terms (`2k`, `2k-1`, `3*(3k-1)`, `3k-2`) instead of
/// `pi.c`'s two, plus a special `a == 2` branch (Gosper's `2^(k-1)` denominator needs `v` to start
/// at `-n_pos` instead of `0`, and `vmax` grown by `big_n - n_pos`, since the base-2 part is
/// otherwise handled the same way `pi.c` handles base 10 directly). `DIVN` in `pi1.c` becomes the
/// `divn!` macro below: advance a running counter, and once it wraps past `a`, strip every factor
/// of `a` out of the term (adjusting `v`).
fn prime_frac_naive_gosper(a: u32, n_pos: u32, big_n: u32) -> f64 {
    let is_two = a == 2;
    let mut vmax = (((3 * big_n) as f64).ln() / (a as f64).ln()) as i32;
    if is_two {
        vmax += big_n as i32 - n_pos as i32;
        if vmax <= 0 {
            return 0.0;
        }
    }
    let mut av: u64 = 1;
    for _ in 0..vmax {
        av *= a as u64;
    }

    let (mut kq1, mut kq2, mut kq3, mut kq4): (i32, i32, i32, i32) = (0, -1, -3, -2);
    let mut v: i32;
    let mut num: u64;
    let mut den: u64 = 1;
    let mut s: u64 = 0;

    if is_two {
        num = 1;
        v = -(n_pos as i32);
    } else {
        num = pow_mod(2, n_pos as u64, av);
        v = 0;
    }

    macro_rules! divn {
        ($t:expr, $vinc:expr, $kq:expr, $kqinc:expr) => {{
            $kq += $kqinc;
            if $kq >= a as i32 {
                while $kq >= a as i32 {
                    $kq -= a as i32;
                }
                if $kq == 0 {
                    loop {
                        $t /= a as u64;
                        v += $vinc;
                        if $t % a as u64 != 0 {
                            break;
                        }
                    }
                }
            }
        }};
    }

    for k in 1..=big_n as u64 {
        let mut t = 2 * k;
        divn!(t, -1, kq1, 2);
        num = mul_mod(num, t, av);

        let mut t = 2 * k - 1;
        divn!(t, -1, kq2, 2);
        num = mul_mod(num, t, av);

        let mut t = 3 * (3 * k - 1);
        divn!(t, 1, kq3, 9);
        den = mul_mod(den, t, av);

        let mut t = 3 * k - 2;
        divn!(t, 1, kq4, 3);
        if !is_two {
            t *= 2;
        } else {
            v += 1;
        }
        den = mul_mod(den, t, av);

        if v > 0 {
            let mut t = if !is_two {
                inv_mod2(den as i64, av as i64) as u64
            } else {
                inv_mod(den as i64, av as i64) as u64
            };
            t = mul_mod(t, num, av);
            for _ in v..vmax {
                t = mul_mod(t, a as u64, av);
            }
            let t1 = 25 * k - 3;
            t = mul_mod(t, t1, av);
            s += t;
            if s >= av {
                s -= av;
            }
        }
    }

    let t = pow_mod(5, n_pos as u64 - 1, av);
    let s = mul_mod(s, t, av);
    s as f64 / av as f64
}

/// Montgomery/batched-inversion/threaded version of [`prime_frac_naive_gosper`], mirroring
/// [`prime_frac_fast`]'s structure (same `Mont`, same `pa[e] = a^e * R mod av` table, same
/// prefix-product batched-inverse trick) but with four running terms instead of two, and
/// `bk[i]` holding `25k-3` (the per-term scalar Gosper's formula multiplies by) instead of
/// `pi.c`'s bare `k`.
///
/// `a == 2` falls back to [`prime_frac_naive_gosper`] instead: Montgomery multiplication needs
/// an odd modulus (it inverts the modulus mod 2^32), and `av = 2^vmax` for that branch never is.
/// There's exactly one prime equal to 2 among the O(N/log N) primes up to `3*big_n`, so falling
/// back to the naive per-k computation for just that one prime costs nothing measurable.
fn prime_frac_fast_gosper(a: u32, n_pos: u32, big_n: u32) -> f64 {
    if a == 2 {
        return prime_frac_naive_gosper(2, n_pos, big_n);
    }
    let vmax = (((3 * big_n) as f64).ln() / (a as f64).ln()) as i32;
    let mut av: u32 = 1;
    for _ in 0..vmax {
        av *= a;
    }
    let mo = Mont::new(av);
    // pa[e] = a^e * R mod av
    let r = (1u64 << 32) % av as u64;
    let mut pa = vec![0u32; vmax as usize + 1];
    let mut p = r;
    for e in pa.iter_mut() {
        *e = p as u32;
        p = p * a as u64 % av as u64;
    }

    let mut num: u32 = pow_mod(2, n_pos as u64, av as u64) as u32;
    let mut den: u32 = 1;
    let mut s: u32 = 0; // sum of t * R^-1
    let (mut kq1, mut kq2, mut kq3, mut kq4): (i32, i32, i32, i32) = (0, -1, -3, -2);
    let mut v: i32 = 0;
    let (mut ba, mut bb, mut bk, mut be) = ([0u32; BLK], [0u32; BLK], [0u32; BLK], [0u8; BLK]);
    let mut pre = [0u32; BLK];
    let mut cnt = 0usize;

    macro_rules! flush {
        () => {
            if cnt > 0 {
                pre[0] = bb[0];
                for i in 1..cnt {
                    pre[i] = mo.mul(pre[i - 1], bb[i]);
                }
                // Montgomery form of (prod b)^-1 -- inv_mod2, not inv_mod: every prime reaching
                // this function is odd (a==2 is routed away before this point).
                let inv = inv_mod2(pre[cnt - 1] as i64, av as i64) as u32;
                let mut ip = mo.mul(inv, mo.r3);
                for i in (0..cnt).rev() {
                    let ib = if i == 0 { ip } else { mo.mul(ip, pre[i - 1]) };
                    if i > 0 {
                        ip = mo.mul(ip, bb[i]);
                    }
                    let mut t = mo.mul(ba[i], ib); // num/den
                    t = mo.mul(t, bk[i]); // * (25k-3) / R
                    t = mo.mul(t, pa[be[i] as usize]); // * a^e
                    s += t;
                    if s >= av {
                        s -= av;
                    }
                }
                cnt = 0;
            }
        };
    }

    macro_rules! divn {
        ($t:expr, $vinc:expr, $kq:expr, $kqinc:expr) => {{
            $kq += $kqinc;
            if $kq >= a as i32 {
                while $kq >= a as i32 {
                    $kq -= a as i32;
                }
                if $kq == 0 {
                    loop {
                        $t /= a;
                        v += $vinc;
                        if $t % a != 0 {
                            break;
                        }
                    }
                }
            }
        }};
    }

    for k in 1..=big_n {
        let mut t = 2 * k;
        divn!(t, -1, kq1, 2);
        num = mo.mul(num, t);

        let mut t = 2 * k - 1;
        divn!(t, -1, kq2, 2);
        num = mo.mul(num, t);

        let mut t = 3 * (3 * k - 1);
        divn!(t, 1, kq3, 9);
        den = mo.mul(den, t);

        let mut t = 3 * k - 2;
        divn!(t, 1, kq4, 3);
        t *= 2; // a != 2 always here, so this always applies (pi1.c's `if (a != 2) t *= 2`)
        den = mo.mul(den, t);

        if v > 0 {
            ba[cnt] = num;
            bb[cnt] = mo.mul(den, mo.r2); // Montgomery form of den
            bk[cnt] = 25 * k - 3;
            be[cnt] = (vmax - v).max(0) as u8;
            cnt += 1;
            if cnt == BLK {
                flush!();
            }
        }
    }
    flush!();
    let s = mo.mul(s, mo.r2); // undo the R^-1
    let t = pow_mod(5, n_pos as u64 - 1, av as u64);
    let s = mul_mod(s as u64, t, av as u64);
    s as f64 / av as f64
}

// ---------------------------------------------------------------------------------
// pi1-simd4: AVX2, 4 *different primes* (each its own modulus) processed in lockstep through
// the num/den update chain -- the part of prime_frac_fast_gosper that runs on *every* k
// regardless of whether that k contributes to the sum, and so is the dominant cost (see
// research/bellard/README.md's "redesign" entry). A separate `--formula`, not a `pi1` variant,
// on purpose: this must never change what `pi`/`pi1`'s own numbers measure.
// ---------------------------------------------------------------------------------

#[cfg(target_arch = "x86_64")]
mod simd4 {
    use std::arch::x86_64::*;

    /// Four lanes, each a u32 value held in the low 32 bits of its own 64-bit slot -- the
    /// layout `_mm256_mul_epu32` wants (it reads each 64-bit lane's low 32 bits and widens the
    /// product to the full 64 bits; no interleaving needed, unlike packing 4 values into 8x u32
    /// sub-lanes would require).
    #[inline(always)]
    fn pack(x: [u32; 4]) -> __m256i {
        // `_mm256_set_epi64x` (an earlier version of this used it) lowers to a multi-instruction
        // insert sequence -- this is called on every `k`, so it was worth the one instruction
        // AVX2 actually has for "4x u32 -> 4x zero-extended u64" (vpmovzxdq), via a 128-bit load
        // of the 4 values packed contiguously plus a single widen.
        unsafe { _mm256_cvtepu32_epi64(_mm_loadu_si128(x.as_ptr() as *const __m128i)) }
    }

    #[inline(always)]
    fn unpack(x: __m256i) -> [u32; 4] {
        let mut out = [0u64; 4];
        unsafe { _mm256_storeu_si256(out.as_mut_ptr() as *mut __m256i, x) };
        [out[0] as u32, out[1] as u32, out[2] as u32, out[3] as u32]
    }

    /// Four-lane Montgomery multiply: `a[i]*b[i]/R mod n[i]`, same REDC as [`super::Mont::mul`],
    /// each lane its own modulus -- REDC has no cross-lane dependency, so batching it across
    /// different moduli is free (the hard part of SIMD-ing this algorithm is elsewhere: the
    /// *factor-stripping* has per-lane-divergent control flow, which is why only this part --
    /// the uniform, always-executes-every-k part -- is vectorised; see [`Lane4`]).
    #[target_feature(enable = "avx2")]
    unsafe fn montmul4(a: __m256i, b: __m256i, n: __m256i, ninv: __m256i) -> __m256i {
        let t = _mm256_mul_epu32(a, b); // 4x a*b, exact u64 (a,b < 2^32)
        // m = (t as u32) * ninv, mod 2^32: mullo_epi32 on the low-u32 sub-lane of each
        // 64-bit slot gives exactly this. The "odd" sub-lane product it also computes is
        // never read -- mul_epu32 below only looks at bits [31:0] of each 64-bit lane.
        let m = _mm256_mullo_epi32(t, ninv);
        let mn = _mm256_mul_epu32(m, n); // 4x m*n, exact u64
        let sum = _mm256_add_epi64(t, mn); // divisible by 2^32 per lane, by REDC construction
        let u = _mm256_srli_epi64(sum, 32); // per-lane >> 32: zero-extended u32 result
        // Branchless conditional subtract, same trick as the scalar `reduce_once` (measured
        // a real win there; see research/bellard/README.md) -- done across all 4 lanes via a
        // signed-compare sign mask (AVX2 has no 64-bit arithmetic right shift to lift the
        // sign bit directly; `cmpgt` against zero gets the same mask).
        let diff = _mm256_sub_epi64(u, n);
        let neg = _mm256_cmpgt_epi64(_mm256_setzero_si256(), diff);
        _mm256_add_epi64(diff, _mm256_and_si256(n, neg))
    }

    /// Montgomery constants for 4 lanes at once -- just [`montmul4`]'s `n`/`ninv` operands,
    /// pre-packed. Used only for the hot, uniform num/den update; the batched-inverse flush
    /// below reuses the existing scalar [`super::Mont`] per lane instead (same lane, same
    /// modulus, but that part is cheap and per-lane-divergent, so there's nothing to gain from
    /// vectorising it too -- see this module's top comment).
    struct Mont4 {
        n_v: __m256i,
        ninv_v: __m256i,
    }

    impl Mont4 {
        fn new(n: [u32; 4]) -> Self {
            let mut ninv = [0u32; 4];
            for (lane, &nn) in n.iter().enumerate() {
                debug_assert!(nn & 1 == 1 && nn < (1 << 31));
                let mut x = nn;
                for _ in 0..5 {
                    x = x.wrapping_mul(2u32.wrapping_sub(nn.wrapping_mul(x)));
                }
                ninv[lane] = x.wrapping_neg();
            }
            Mont4 {
                n_v: pack(n),
                ninv_v: pack(ninv),
            }
        }
        #[target_feature(enable = "avx2")]
        unsafe fn mul(&self, a: __m256i, b: __m256i) -> __m256i {
            unsafe { montmul4(a, b, self.n_v, self.ninv_v) }
        }
    }

    /// One lane's batched-inverse flush -- identical logic to [`super::prime_frac_fast_gosper`]'s
    /// `flush!` macro, extracted to a plain function since this part stays scalar per lane (see
    /// the module doc comment for why).
    #[allow(clippy::too_many_arguments)]
    fn flush_lane(
        mo: &super::Mont,
        pa: &[u32],
        av: u32,
        cnt: &mut usize,
        ba: &[u32; super::BLK],
        bb: &[u32; super::BLK],
        bk: &[u32; super::BLK],
        be: &[u8; super::BLK],
        s: &mut u32,
    ) {
        if *cnt == 0 {
            return;
        }
        let mut pre = [0u32; super::BLK];
        pre[0] = bb[0];
        for i in 1..*cnt {
            pre[i] = mo.mul(pre[i - 1], bb[i]);
        }
        let inv = super::inv_mod2(pre[*cnt - 1] as i64, av as i64) as u32;
        let mut ip = mo.mul(inv, mo.r3);
        for i in (0..*cnt).rev() {
            let ib = if i == 0 { ip } else { mo.mul(ip, pre[i - 1]) };
            if i > 0 {
                ip = mo.mul(ip, bb[i]);
            }
            let mut t = mo.mul(ba[i], ib);
            t = mo.mul(t, bk[i]);
            t = mo.mul(t, pa[be[i] as usize]);
            *s += t;
            if *s >= av {
                *s -= av;
            }
        }
        *cnt = 0;
    }

    /// The vectorised counterpart of [`super::prime_frac_fast_gosper`]: four *different* odd
    /// primes processed in lockstep. `primes4` must all be odd (the caller routes `a == 2`
    /// elsewhere, same as the scalar function does) -- debug-checked, not re-handled here.
    ///
    /// Only the num/den update (two terms each, every k, regardless of whether k contributes)
    /// is vectorised: that's the part with no per-lane-divergent control flow, and the part that
    /// runs unconditionally for every one of `big_n` values of k, which is why it's the
    /// dominant cost (see `research/bellard/README.md`'s "redesign" entry for the reasoning).
    /// Factor-stripping (`DIVN`) and the batched-inverse extraction stay scalar, one lane at a
    /// time, reusing the same per-prime bookkeeping shape `prime_frac_fast_gosper` already uses.
    #[target_feature(enable = "avx2")]
    unsafe fn frac4_avx2(primes4: [u32; 4], n_pos: u32, big_n: u32) -> [f64; 4] {
        let mut vmax = [0i32; 4];
        let mut av = [0u32; 4];
        for lane in 0..4 {
            let a = primes4[lane];
            debug_assert!(
                a != 2 && a % 2 == 1,
                "frac4_avx2 requires odd primes, got {a}"
            );
            vmax[lane] = (((3 * big_n) as f64).ln() / (a as f64).ln()) as i32;
            let mut avv = 1u32;
            for _ in 0..vmax[lane] {
                avv *= a;
            }
            av[lane] = avv;
        }

        let mo_scalar: [super::Mont; 4] = std::array::from_fn(|lane| super::Mont::new(av[lane]));
        let mo4 = Mont4::new(av);

        // pa[lane][e] = a[lane]^e * R[lane] mod av[lane], same construction as the scalar fast
        // variants' own `pa` table.
        let pa: [Vec<u32>; 4] = std::array::from_fn(|lane| {
            let r = (1u64 << 32) % av[lane] as u64;
            let mut v = Vec::with_capacity(vmax[lane] as usize + 1);
            let mut p = r;
            for _ in 0..=vmax[lane] {
                v.push(p as u32);
                p = p * primes4[lane] as u64 % av[lane] as u64;
            }
            v
        });

        let num: [u32; 4] =
            std::array::from_fn(|lane| super::pow_mod(2, n_pos as u64, av[lane] as u64) as u32);
        // Carried forward as vectors across the whole k loop, never round-tripped through scalar
        // form unless a lane actually needs to push a contribution (measured: packing/unpacking
        // on every k, even when nothing needed it that iteration, was enough overhead to erase
        // the gain from vectorising the multiplies in the first place -- see
        // research/bellard/README.md's account of this).
        let mut numv = pack(num);
        let mut denv = pack([1u32; 4]);

        let mut kq1 = [0i32; 4];
        let mut kq2 = [-1i32; 4];
        let mut kq3 = [-3i32; 4];
        let mut kq4 = [-2i32; 4];
        let mut v = [0i32; 4];

        let mut ba = [[0u32; super::BLK]; 4];
        let mut bb = [[0u32; super::BLK]; 4];
        let mut bk = [[0u32; super::BLK]; 4];
        let mut be = [[0u8; super::BLK]; 4];
        let mut cnt = [0usize; 4];
        let mut s = [0u32; 4];

        for k in 1..=big_n {
            let mut t_num1 = [0u32; 4];
            let mut t_num2 = [0u32; 4];
            let mut t_den1 = [0u32; 4];
            let mut t_den2 = [0u32; 4];

            for lane in 0..4 {
                let a = primes4[lane] as i32;

                let mut t = 2 * k;
                kq1[lane] += 2;
                if kq1[lane] >= a {
                    while kq1[lane] >= a {
                        kq1[lane] -= a;
                    }
                    if kq1[lane] == 0 {
                        loop {
                            t /= primes4[lane];
                            v[lane] -= 1;
                            if t % primes4[lane] != 0 {
                                break;
                            }
                        }
                    }
                }
                t_num1[lane] = t;

                let mut t = 2 * k - 1;
                kq2[lane] += 2;
                if kq2[lane] >= a {
                    while kq2[lane] >= a {
                        kq2[lane] -= a;
                    }
                    if kq2[lane] == 0 {
                        loop {
                            t /= primes4[lane];
                            v[lane] -= 1;
                            if t % primes4[lane] != 0 {
                                break;
                            }
                        }
                    }
                }
                t_num2[lane] = t;

                let mut t = 3 * (3 * k - 1);
                kq3[lane] += 9;
                if kq3[lane] >= a {
                    while kq3[lane] >= a {
                        kq3[lane] -= a;
                    }
                    if kq3[lane] == 0 {
                        loop {
                            t /= primes4[lane];
                            v[lane] += 1;
                            if t % primes4[lane] != 0 {
                                break;
                            }
                        }
                    }
                }
                t_den1[lane] = t;

                let mut t = 3 * k - 2;
                kq4[lane] += 3;
                if kq4[lane] >= a {
                    while kq4[lane] >= a {
                        kq4[lane] -= a;
                    }
                    if kq4[lane] == 0 {
                        loop {
                            t /= primes4[lane];
                            v[lane] += 1;
                            if t % primes4[lane] != 0 {
                                break;
                            }
                        }
                    }
                }
                t_den2[lane] = t * 2; // a != 2 always here, so this always applies
            }

            numv = unsafe { mo4.mul(numv, pack(t_num1)) };
            numv = unsafe { mo4.mul(numv, pack(t_num2)) };
            denv = unsafe { mo4.mul(denv, pack(t_den1)) };
            denv = unsafe { mo4.mul(denv, pack(t_den2)) };

            if v.iter().any(|&x| x > 0) {
                let num = unpack(numv);
                let den = unpack(denv);
                for lane in 0..4 {
                    if v[lane] > 0 {
                        let c = cnt[lane];
                        ba[lane][c] = num[lane];
                        bb[lane][c] = mo_scalar[lane].mul(den[lane], mo_scalar[lane].r2);
                        bk[lane][c] = 25 * k - 3;
                        be[lane][c] = (vmax[lane] - v[lane]).max(0) as u8;
                        cnt[lane] += 1;
                        if cnt[lane] == super::BLK {
                            flush_lane(
                                &mo_scalar[lane],
                                &pa[lane],
                                av[lane],
                                &mut cnt[lane],
                                &ba[lane],
                                &bb[lane],
                                &bk[lane],
                                &be[lane],
                                &mut s[lane],
                            );
                        }
                    }
                }
            }
        }

        let mut out = [0f64; 4];
        for lane in 0..4 {
            flush_lane(
                &mo_scalar[lane],
                &pa[lane],
                av[lane],
                &mut cnt[lane],
                &ba[lane],
                &bb[lane],
                &bk[lane],
                &be[lane],
                &mut s[lane],
            );
            let s_lane = mo_scalar[lane].mul(s[lane], mo_scalar[lane].r2);
            let t = super::pow_mod(5, n_pos as u64 - 1, av[lane] as u64);
            let s_lane = super::mul_mod(s_lane as u64, t, av[lane] as u64);
            out[lane] = s_lane as f64 / av[lane] as f64;
        }
        out
    }

    /// Safe entry point: checks for AVX2 at runtime (this module compiles for any x86_64 target,
    /// but not every x86_64 CPU has AVX2) before calling the unsafe, `#[target_feature]`-gated
    /// implementation. `primes4` must be four distinct odd primes; `a == 2` is the caller's
    /// responsibility to route elsewhere (see [`super::prime_frac_fast_gosper`]'s own docs for
    /// why `a == 2` can't go through Montgomery arithmetic at all).
    pub fn frac4(primes4: [u32; 4], n_pos: u32, big_n: u32) -> [f64; 4] {
        assert!(
            is_x86_feature_detected!("avx2"),
            "pi1-simd4 requires an AVX2-capable CPU"
        );
        unsafe { frac4_avx2(primes4, n_pos, big_n) }
    }

    #[cfg(test)]
    mod tests {
        use super::*;

        /// Scalar reference REDC, copied rather than reusing `super::super::Mont` directly, so
        /// this test doesn't share a bug with the code it's checking.
        fn scalar_redc(a: u32, b: u32, n: u32, ninv: u32) -> u32 {
            let t = a as u64 * b as u64;
            let m = (t as u32).wrapping_mul(ninv);
            let u = ((t + m as u64 * n as u64) >> 32) as u32;
            if u >= n { u - n } else { u }
        }

        fn ninv_of(n: u32) -> u32 {
            let mut x = n;
            for _ in 0..5 {
                x = x.wrapping_mul(2u32.wrapping_sub(n.wrapping_mul(x)));
            }
            x.wrapping_neg()
        }

        /// Deterministic xorshift64, matching the convention used elsewhere in this repo
        /// (crates/pi-digits's tests) so this needs no extra dependency.
        struct Xorshift(u64);
        impl Xorshift {
            fn next(&mut self) -> u64 {
                self.0 ^= self.0 << 13;
                self.0 ^= self.0 >> 7;
                self.0 ^= self.0 << 17;
                self.0
            }
        }

        #[test]
        fn montmul4_matches_scalar_redc_on_random_inputs() {
            if !is_x86_feature_detected!("avx2") {
                eprintln!("skipping: no avx2 on this machine");
                return;
            }
            let mut rng = Xorshift(0xC0FF_EE15_BADD_CAFE);
            for _ in 0..20_000 {
                let mut n = [0u32; 4];
                let mut a = [0u32; 4];
                let mut b = [0u32; 4];
                let mut ninv = [0u32; 4];
                for lane in 0..4 {
                    // Odd modulus < 2^31, matching Mont::new's own precondition.
                    n[lane] = ((rng.next() as u32) | 1) & 0x7FFF_FFFF;
                    if n[lane] == 0 {
                        n[lane] = 1;
                    }
                    a[lane] = rng.next() as u32 % n[lane];
                    b[lane] = rng.next() as u32 % n[lane];
                    ninv[lane] = ninv_of(n[lane]);
                }
                let got = unsafe { unpack(montmul4(pack(a), pack(b), pack(n), pack(ninv))) };
                for lane in 0..4 {
                    let want = scalar_redc(a[lane], b[lane], n[lane], ninv[lane]);
                    assert_eq!(
                        got[lane], want,
                        "lane {lane}: a={} b={} n={} ninv={}",
                        a[lane], b[lane], n[lane], ninv[lane]
                    );
                }
            }
        }

        /// Stronger than an end-to-end digit match: each lane's own contribution from
        /// [`frac4`] must agree with calling the already-validated
        /// [`super::prime_frac_fast_gosper`] on that same prime individually. An end-to-end
        /// digit match could in principle hide two compensating bugs; this can't, since it
        /// checks lane-by-lane with no cross-lane room to cancel.
        #[test]
        fn frac4_matches_prime_frac_fast_gosper_per_lane() {
            if !is_x86_feature_detected!("avx2") {
                eprintln!("skipping: no avx2 on this machine");
                return;
            }
            for &n_pos in &[10_000u32, 20_000, 50_000] {
                let big_n = ((n_pos as f64 + 20.0) * 10f64.ln() / 13.5f64.ln()) as u32;
                let primes: Vec<u32> = super::super::primes_upto(3 * big_n as usize)
                    .into_iter()
                    .filter(|&p| p != 2)
                    .take(16)
                    .collect();
                for group in primes.chunks_exact(4) {
                    let group4: [u32; 4] = group.try_into().unwrap();
                    let got = frac4(group4, n_pos, big_n);
                    for lane in 0..4 {
                        let want = super::super::prime_frac_fast_gosper(group4[lane], n_pos, big_n);
                        assert!(
                            (got[lane] - want).abs() < 1e-9,
                            "n_pos={n_pos} big_n={big_n} lane {lane} prime={}: got={} want={}",
                            group4[lane],
                            got[lane],
                            want
                        );
                    }
                }
            }
        }
    }
}

/// `--formula pi1-simd4`: same math as [`prime_frac_fast_gosper`] (Gosper's series), but the
/// odd primes are processed four at a time via [`simd4::frac4`] instead of one at a time. `a ==
/// 2` and any primes left over after grouping into fours go through the existing scalar paths
/// ([`prime_frac_naive_gosper`] for `a == 2`, [`prime_frac_fast_gosper`] for the 0-3 leftover
/// primes) -- both already validated, nothing new to get wrong there.
fn run_pi1_simd4(n: u32, threads: usize) {
    let big_n = ((n as f64 + 20.0) * 10f64.ln() / 13.5f64.ln()) as u32;
    let all_primes: Vec<u32> = primes_upto(3 * big_n as usize)
        .into_iter()
        .filter(|&p| p >= 2)
        .collect();
    let odd: Vec<u32> = all_primes.iter().copied().filter(|&p| p != 2).collect();
    let has_two = all_primes.first() == Some(&2);

    let n_full_chunks = odd.len() / 4;
    let chunks: Vec<[u32; 4]> = (0..n_full_chunks)
        .map(|c| std::array::from_fn(|i| odd[c * 4 + i]))
        .collect();
    let remainder: Vec<u32> = odd[n_full_chunks * 4..].to_vec();

    fn parallel_sum<T: Sync>(items: &[T], threads: usize, f: impl Fn(&T) -> f64 + Sync) -> f64 {
        let next = AtomicUsize::new(0);
        let partials: Vec<f64> = std::thread::scope(|sc| {
            let hs: Vec<_> = (0..threads.max(1))
                .map(|_| {
                    sc.spawn(|| {
                        let mut acc = 0f64;
                        loop {
                            let i = next.fetch_add(1, Ordering::Relaxed);
                            if i >= items.len() {
                                return acc;
                            }
                            acc = (acc + f(&items[i])) % 1.0;
                        }
                    })
                })
                .collect();
            hs.into_iter().map(|h| h.join().unwrap()).collect()
        });
        partials.into_iter().fold(0f64, |a, b| (a + b) % 1.0)
    }

    let mut sum = if has_two {
        prime_frac_naive_gosper(2, n, big_n)
    } else {
        0.0
    };
    sum =
        (sum + parallel_sum(&chunks, threads, |&group| {
            simd4::frac4(group, n, big_n)
                .into_iter()
                .fold(0f64, |a, b| (a + b) % 1.0)
        })) % 1.0;
    sum =
        (sum + parallel_sum(&remainder, threads, |&p| {
            prime_frac_fast_gosper(p, n, big_n)
        })) % 1.0;

    println!(
        "Decimal digits of pi at position {}: {:09}",
        n,
        (sum * 1e9) as i32
    );
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let n: u32 = args
        .get(1)
        .and_then(|s| s.parse().ok())
        .filter(|&n| n > 0)
        .unwrap_or_else(|| {
            eprintln!(
                "usage: bellard-pi <n> [--formula pi|pi1|pi1-simd4] [--variant naive|fast] [--threads T]"
            );
            std::process::exit(1)
        });
    let opt = |name: &str| {
        args.iter()
            .position(|a| a == name)
            .and_then(|i| args.get(i + 1))
    };
    let formula = opt("--formula").map(String::as_str).unwrap_or("pi");
    let fast = opt("--variant").map(|v| v != "naive").unwrap_or(true);
    let threads: usize = opt("--threads")
        .and_then(|t| t.parse().ok())
        .unwrap_or_else(|| std::thread::available_parallelism().map_or(1, |x| x.get()));

    if formula == "pi1-simd4" {
        if opt("--variant").is_some_and(|v| v == "naive") {
            eprintln!(
                "error: --formula pi1-simd4 has no naive variant (it's the SIMD redesign itself)"
            );
            std::process::exit(1);
        }
        run_pi1_simd4(n, threads);
        return;
    }

    type PrimeFrac = fn(u32, u32, u32) -> f64;
    let (big_n, primes, f): (u32, Vec<u32>, PrimeFrac) = match formula {
        "pi" => {
            let big_n = ((n as f64 + 20.0) * 10f64.ln() / 2f64.ln()) as u32;
            let primes = primes_upto(2 * big_n as usize)
                .into_iter()
                .filter(|&p| p >= 3)
                .collect();
            let f = if fast {
                prime_frac_fast
            } else {
                prime_frac_naive
            };
            (big_n, primes, f)
        }
        "pi1" => {
            let big_n = ((n as f64 + 20.0) * 10f64.ln() / 13.5f64.ln()) as u32;
            let primes = primes_upto(3 * big_n as usize)
                .into_iter()
                .filter(|&p| p >= 2)
                .collect();
            let f = if fast {
                prime_frac_fast_gosper
            } else {
                prime_frac_naive_gosper
            };
            (big_n, primes, f)
        }
        other => {
            eprintln!("error: unknown --formula {other:?}, expected pi, pi1, or pi1-simd4");
            std::process::exit(1);
        }
    };
    let fracs: Vec<f64> = {
        let next = AtomicUsize::new(0);
        let mut out = vec![0f64; primes.len()];
        let chunks: Vec<Vec<(usize, f64)>> = std::thread::scope(|sc| {
            let hs: Vec<_> = (0..threads.max(1))
                .map(|_| {
                    sc.spawn(|| {
                        let mut mine = Vec::new();
                        loop {
                            let i = next.fetch_add(1, Ordering::Relaxed);
                            if i >= primes.len() {
                                return mine;
                            }
                            mine.push((i, f(primes[i], n, big_n)));
                        }
                    })
                })
                .collect();
            hs.into_iter().map(|h| h.join().unwrap()).collect()
        });
        for (i, v) in chunks.into_iter().flatten() {
            out[i] = v;
        }
        out
    };
    let mut sum = 0f64;
    for x in fracs {
        sum = (sum + x) % 1.0;
    }
    println!(
        "Decimal digits of pi at position {}: {:09}",
        n,
        (sum * 1e9) as i32
    );
}
