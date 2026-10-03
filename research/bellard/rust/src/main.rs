//! Rust port of Bellard's 1997 `pi.c` (n-th decimal digit of pi, O(n^2) time, tiny memory).
//!
//! usage: bellard-pi <n> [--variant naive|fast] [--threads T]
//!
//! * `naive`: line-for-line port (u64 `%` mulmod, one extended-Euclid inverse per k).
//! * `fast` (default): Montgomery u32 arithmetic, batched modular inversion (Montgomery's trick,
//!   one Euclid per block instead of per k), sieved primes, one thread per core over the primes.
//!
//! Both print the same digits as `pi.c` (position n has the same meaning: digits n.. of pi's
//! fractional part, i.e. `pi.c 1` prints 141592653).
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
        Stepper { a, v: 0, kq: 1, kq2: 1 }
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
        Mont { n, ninv: x.wrapping_neg(), r2: r2 as u32, r3: r3 as u32 }
    }
    /// a * b / R mod n, for a < n, b < 2^32.
    #[inline(always)]
    fn mul(&self, a: u32, b: u32) -> u32 {
        let t = a as u64 * b as u64;
        let m = (t as u32).wrapping_mul(self.ninv);
        let u = ((t + m as u64 * self.n as u64) >> 32) as u32;
        if u >= self.n { u - self.n } else { u }
    }
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

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let n: u32 = args.get(1).and_then(|s| s.parse().ok()).filter(|&n| n > 0).unwrap_or_else(|| {
        eprintln!("usage: bellard-pi <n> [--variant naive|fast] [--threads T]");
        std::process::exit(1)
    });
    let opt = |name: &str| args.iter().position(|a| a == name).and_then(|i| args.get(i + 1));
    let fast = opt("--variant").map(|v| v != "naive").unwrap_or(true);
    let threads: usize = opt("--threads")
        .and_then(|t| t.parse().ok())
        .unwrap_or_else(|| std::thread::available_parallelism().map_or(1, |x| x.get()));

    let big_n = ((n as f64 + 20.0) * 10f64.ln() / 2f64.ln()) as u32;
    let primes: Vec<u32> = primes_upto(2 * big_n as usize).into_iter().filter(|&p| p >= 3).collect();
    let f = if fast { prime_frac_fast } else { prime_frac_naive };
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
    println!("Decimal digits of pi at position {}: {:09}", n, (sum * 1e9) as i32);
}
