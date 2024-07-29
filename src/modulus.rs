use std::ops::{Add, Div, Sub};

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
        let mut rng = RandState::new();
        self.n.clone().random_below(&mut rng)
    }

    pub fn pow(&self, val: &Int, exp: &Int) -> Option<Int> {
        val.clone().pow_mod(exp, &self.n).ok()
    }

    pub fn has_sqrt(&self, a: &Int) -> bool {
        a.legendre(&self.n) == 1
    }

    // Square root mod n: Tonelli and Shanks
    // (Chapter 2.4, Listings 2.6, 2.7, 2.8)
    pub fn sqrt(&self, a: &Int) -> Option<Int> {
        // Exit immediately if nonresidue
        if !self.has_sqrt(a) {
            return None;
        }

        // Last 2 bits set? Return `x = a^((q+1)/4)`
        if self.n.get_bit(0) && self.n.get_bit(1) {
            let q = self.n.clone().add(Int::ONE).div(&Int::from(4));
            let x = self.pow(a, &q)?;
            return Some(x);
        }

        let p = self.n.clone();
        // `q = p - 1`
        let mut q = p.clone().sub(Int::ONE);
        // find number of binary zeros (first index of binary one)
        let e = q.find_one(0).expect("e");
        //  break down `p - 1` into `2^e * q`
        for _ in 0..e {
            q = q.div(2);
        }

        // find a generator
        let n = loop {
            // randomly search for nonresidue
            let x = self.rand();
            if self.has_sqrt(&x) {
                break x;
            }
        };

        // initialize working components
        let mut y = n.pow_mod(&q, &p).expect("y = n^q mod p");
        let mut r = e;
        let x = a
            .clone()
            .pow_mod(&q.clone().sub(Int::ONE).div(Int::from(2)), &p)
            .expect("x = a^((q-1)/2) mod p");
        let mut b = (a.clone() * x.clone() * x.clone()).modulo(&p);
        let mut x = (a.clone() * x.clone()).modulo(&p);

        // loop on algorithm
        loop {
            // Minimum `m` such that `b^(2^m) = 1 mod p`
            let mut m = 1;
            while m < r {
                let one = b
                    .clone()
                    .pow_mod(&Int::from(1 << (m - 1)), &p)
                    .expect("b^(2^m) = 1 mod p");
                if &one == Int::ONE {
                    break;
                }
                m += 1;
            }
            if m == r {
                unreachable!("'should never happen because `a` is quadratic residue'");
            }

            let e = r - m - 1;
            let e = 1 << (e - 1);
            let t = y.clone().pow_mod(&Int::from(e), &p).expect("t");
            y = t.clone().pow_mod(&Int::from(2), &p).expect("y");
            r = m;
            x = (x * t).modulo(&p);
            b = (b * y.clone()).modulo(&p);
            if &b == Int::ONE {
                return Some(x);
            }
        }
    }
}
