use std::fs::File;
use std::io::Read;

use rug::{rand::RandState, Integer as Int};

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Modulus {
    pub n: Int,
}

impl Modulus {
    pub fn new(n: &Int) -> Self {
        Self { n: n.to_owned() }
    }

    pub fn add(&self, a: &Int, b: &Int) -> Int {
        Int::from(a + b).modulo(&self.n)
    }

    pub fn sub(&self, a: &Int, b: &Int) -> Int {
        Int::from(a - b).modulo(&self.n)
    }

    pub fn mul(&self, a: &Int, b: &Int) -> Int {
        Int::from(a * b).modulo(&self.n)
    }

    pub fn div(&self, a: &Int, b: &Int) -> Option<Int> {
        let i = self.inv(b)?;
        let ret = Int::from(a * &i).modulo(&self.n);
        Some(ret)
    }

    pub fn neg(&self, a: &Int) -> Int {
        Int::from(-a).modulo(&self.n)
    }

    pub fn inv(&self, a: &Int) -> Option<Int> {
        match a.clone().invert(&self.n) {
            Ok(inverse) => {
                let one = Int::from(a * &inverse).modulo(&self.n);
                assert_eq!(one, 1);
                Some(inverse)
            }
            Err(_) => None,
        }
    }

    pub fn rand(&self) -> Int {
        // `RandState::new()` uses a fixed default seed; reseed from the OS
        // RNG so successive calls produce independent values.
        let mut rng = RandState::new();
        let mut seed_bytes = [0u8; 32];
        File::open("/dev/urandom")
            .and_then(|mut f| f.read_exact(&mut seed_bytes))
            .expect("seed Modulus::rand from /dev/urandom");
        let seed = Int::from_digits(&seed_bytes, rug::integer::Order::Msf);
        rng.seed(&seed);
        self.n.clone().random_below(&mut rng)
    }

    pub fn pow(&self, val: &Int, exp: &Int) -> Option<Int> {
        val.clone().pow_mod(exp, &self.n).ok()
    }

    pub fn has_sqrt(&self, a: &Int) -> bool {
        a.legendre(&self.n) == 1
    }

    // Square root mod n: Tonelli–Shanks (chapter 2.4).
    // Returns None on quadratic non-residues; otherwise an x with x² ≡ a mod n.
    pub fn sqrt(&self, a: &Int) -> Option<Int> {
        // Reduce a mod n first; a = 0 has sqrt 0.
        let a = self.add(a, &Int::ZERO);
        if a.is_zero() {
            return Some(Int::ZERO);
        }
        if !self.has_sqrt(&a) {
            return None;
        }

        let p = self.n.clone();
        // Fast path when p ≡ 3 (mod 4): sqrt(a) = a^((p+1)/4).
        if p.get_bit(0) && p.get_bit(1) {
            let exp = Int::from(&p + 1) / 4;
            return self.pow(&a, &exp);
        }

        // Decompose p - 1 = q · 2^s with q odd.
        let mut q: Int = p.clone() - 1;
        let mut s: u32 = 0;
        while !q.get_bit(0) {
            q >>= 1u32;
            s += 1;
        }

        // Find a non-residue z to build a primitive 2^s-th root c = z^q.
        let z = loop {
            let cand = self.rand();
            if !cand.is_zero() && !self.has_sqrt(&cand) {
                break cand;
            }
        };

        // Standard Tonelli–Shanks state: M, c, t, r (Wikipedia notation).
        let mut m_state: u32 = s;
        let mut c = z.pow_mod(&q, &p).expect("c = z^q mod p");
        let mut t = a.clone().pow_mod(&q, &p).expect("t = a^q mod p");
        let r_exp = Int::from(&q + 1) / 2;
        let mut r = a.pow_mod(&r_exp, &p).expect("r = a^((q+1)/2) mod p");

        loop {
            if t == 1 {
                return Some(r);
            }
            // Smallest i ∈ [1, M) such that t^(2^i) = 1.
            let mut i: u32 = 1;
            let mut tmp =
                t.clone().pow_mod(&Int::from(2), &p).expect("t²");
            while tmp != 1 {
                i += 1;
                if i >= m_state {
                    // Should not happen for a true QR; bail rather than loop.
                    return None;
                }
                tmp = tmp.pow_mod(&Int::from(2), &p).expect("repeated sq");
            }
            // b = c^(2^(M-i-1))
            let shift = m_state - i - 1;
            let exp = Int::from(1) << shift;
            let b = c.clone().pow_mod(&exp, &p).expect("b");
            // Update: M = i, c = b², t = t·b², r = r·b
            let b2 = b.clone().pow_mod(&Int::from(2), &p).expect("b²");
            m_state = i;
            c = b2.clone();
            t = (t * &b2).modulo(&p);
            r = (r * &b).modulo(&p);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Exercise Tonelli-Shanks (the slow `p ≡ 1 mod 4` path) — the BN254 base
    /// field prime has p ≡ 1 mod 4 in the lower bits' sense for this test, so
    /// the algorithm actually runs the Tonelli-Shanks loop. We square a known
    /// integer and verify sqrt recovers it (up to sign).
    #[test]
    fn test_sqrt_tonelli_shanks_path() {
        // Pick a prime that is ≡ 1 mod 4 to force the slow path.
        // 13 ≡ 1 mod 4. has_sqrt(4) is true, sqrt(4) ∈ {2, 11}.
        let m = Modulus::new(&Int::from(13));
        let s = m.sqrt(&Int::from(4)).expect("sqrt(4) mod 13");
        assert!(s == Int::from(2) || s == Int::from(11));
        // Another QR: sqrt(9) ∈ {3, 10}.
        let s = m.sqrt(&Int::from(9)).expect("sqrt(9) mod 13");
        assert!(s == Int::from(3) || s == Int::from(10));
        // Non-residue: sqrt(2) returns None.
        assert!(m.sqrt(&Int::from(2)).is_none());
    }

    #[test]
    fn test_rand_is_nondeterministic() {
        // Successive calls must not return the same value (with overwhelming
        // probability over a large modulus). Two consecutive calls returning
        // the same value would catch a regression to the unseeded default.
        let m = Modulus::new(&Int::from_str_radix(
            "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFEFFFFFC2F",
            16,
        ).unwrap());
        let a = m.rand();
        let b = m.rand();
        assert_ne!(a, b);
    }
}
