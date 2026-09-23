//! The only big-integer surface the digit algorithms use. Two backends: `gmp` (rug) for the
//! CLI, `pure` (dashu-int) for WASM. Keep this list minimal: every method here must exist, with
//! identical semantics, in both backends, and `tests/bignum.rs` runs against both.

#[cfg(feature = "gmp")]
mod imp {
    use rug::{Integer, ops::Pow};

    #[derive(Clone, Debug, PartialEq, Eq)]
    pub struct Big(Integer);

    impl Big {
        pub fn from_u64(v: u64) -> Big {
            Big(Integer::from(v))
        }
        pub fn zero() -> Big {
            Big(Integer::new())
        }
        pub fn one() -> Big {
            Big(Integer::from(1))
        }
        pub fn mul(&self, o: &Big) -> Big {
            Big(Integer::from(&self.0 * &o.0))
        }
        pub fn mul_u64(&self, v: u64) -> Big {
            Big(Integer::from(&self.0 * v))
        }
        pub fn add(&self, o: &Big) -> Big {
            Big(Integer::from(&self.0 + &o.0))
        }
        pub fn sub(&self, o: &Big) -> Big {
            Big(Integer::from(&self.0 - &o.0))
        }
        pub fn rem(&self, m: &Big) -> Big {
            Big(Integer::from(&self.0 % &m.0))
        }
        pub fn rem_u64(&self, m: u64) -> u64 {
            Integer::from(&self.0 % m).to_u64().expect("rem < m")
        }
        pub fn div_u64_exact(&self, d: u64) -> Big {
            debug_assert!(Integer::from(&self.0 % d) == 0, "div_u64_exact: not exact");
            Big(Integer::from(&self.0 / d))
        }
        pub fn div_u64(&self, d: u64) -> Big {
            Big(Integer::from(&self.0 / d))
        }
        pub fn to_u64(&self) -> Option<u64> {
            self.0.to_u64()
        }
        pub fn is_zero(&self) -> bool {
            self.0 == 0
        }
        pub fn bits(&self) -> u64 {
            self.0.significant_bits() as u64
        }
        pub fn shl(&self, k: u32) -> Big {
            Big(Integer::from(&self.0 << k))
        }
        pub fn pow_u64(base: u64, exp: u32) -> Big {
            Big(Integer::from(base).pow(exp))
        }
        pub fn binomial(n: u64, k: u32) -> Big {
            Big(Integer::from(n).binomial(k))
        }
        pub fn to_decimal_string(&self) -> String {
            self.0.to_string()
        }
    }
}

#[cfg(feature = "pure")]
mod imp {
    use dashu_base::BitTest;
    use dashu_int::UBig;

    #[derive(Clone, Debug, PartialEq, Eq)]
    pub struct Big(UBig);

    impl Big {
        pub fn from_u64(v: u64) -> Big {
            Big(UBig::from(v))
        }
        pub fn zero() -> Big {
            Big(UBig::ZERO)
        }
        pub fn one() -> Big {
            Big(UBig::ONE)
        }
        pub fn mul(&self, o: &Big) -> Big {
            Big(&self.0 * &o.0)
        }
        pub fn mul_u64(&self, v: u64) -> Big {
            Big(&self.0 * v)
        }
        pub fn add(&self, o: &Big) -> Big {
            Big(&self.0 + &o.0)
        }
        pub fn sub(&self, o: &Big) -> Big {
            Big(&self.0 - &o.0)
        }
        pub fn rem(&self, m: &Big) -> Big {
            Big(&self.0 % &m.0)
        }
        pub fn rem_u64(&self, m: u64) -> u64 {
            &self.0 % m
        }
        pub fn div_u64_exact(&self, d: u64) -> Big {
            debug_assert!(&self.0 % d == 0, "div_u64_exact: not exact");
            Big(&self.0 / d)
        }
        pub fn div_u64(&self, d: u64) -> Big {
            Big(&self.0 / d)
        }
        pub fn to_u64(&self) -> Option<u64> {
            u64::try_from(&self.0).ok()
        }
        pub fn is_zero(&self) -> bool {
            self.0 == UBig::ZERO
        }
        pub fn bits(&self) -> u64 {
            self.0.bit_len() as u64
        }
        pub fn shl(&self, k: u32) -> Big {
            Big(&self.0 << k as usize)
        }
        pub fn pow_u64(base: u64, exp: u32) -> Big {
            Big(UBig::from(base).pow(exp as usize))
        }
        /// `C(n, k)`, via the standard incremental identity `C(n,i+1) = C(n,i)*(n-i)/(i+1)`
        /// (exact at every step, since `C(n,i)*(n-i) = C(n,i+1)*(i+1)` identically).
        pub fn binomial(n: u64, k: u32) -> Big {
            if k as u64 > n {
                return Big::zero();
            }
            let mut acc = UBig::ONE;
            for i in 0..k as u64 {
                acc = &acc * (n - i);
                acc = &acc / (i + 1);
            }
            Big(acc)
        }
        pub fn to_decimal_string(&self) -> String {
            self.0.to_string()
        }
    }
}

pub use imp::Big;
