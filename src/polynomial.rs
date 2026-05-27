use rug::Integer as Int;

use crate::modulus::Modulus;

#[derive(Clone, Debug, PartialEq)]
pub struct Polynomial {
    coef: Vec<Int>,
}

#[allow(clippy::needless_range_loop)]
impl Polynomial {
    pub fn zeros(len: usize) -> Self {
        Self {
            coef: vec![Int::ZERO; len],
        }
    }

    #[allow(clippy::new_without_default)]
    pub fn new() -> Self {
        Self { coef: Vec::new() }
    }

    pub fn degree(&self) -> usize {
        self.coef.len().saturating_sub(1)
    }

    pub fn is_zero(&self) -> bool {
        self.coef.is_empty() || self.coef.iter().all(|k| k.is_zero())
    }

    pub fn trim(mut self) -> Self {
        while self.coef.last().is_some_and(|k| k.is_zero()) {
            self.coef.pop();
        }
        self
    }

    pub fn set(&mut self, pow: usize, int: Int) {
        if pow > self.degree() {
            self.coef.resize(pow + 1, Int::ZERO);
        }
        self.coef[pow] = int;
    }

    pub fn get(&self, pow: usize) -> Int {
        self.coef.get(pow).cloned().unwrap_or(Int::ZERO)
    }

    pub fn neg(&self) -> Self {
        let mut ret = self.clone();
        for k in ret.coef.iter_mut() {
            *k = Int::from(-&*k);
        }
        ret
    }

    pub fn add(&self, that: &Self) -> Self {
        let mut ret = Self::new();
        let len = self.coef.len().max(that.coef.len());
        ret.coef.resize(len, Int::ZERO);
        for i in 0..len {
            let lhs = self.coef.get(i).cloned().unwrap_or(Int::ZERO);
            let rhs = that.coef.get(i).cloned().unwrap_or(Int::ZERO);
            ret.coef[i] = lhs + rhs;
        }
        ret.trim()
    }

    pub fn sub(&self, that: &Self) -> Self {
        let that = that.neg();
        self.add(&that)
    }

    // Polynomials multiplication & reduction by irreducible polynomial over prime field
    pub fn mul(
        &self,
        that: &Self,
        irreducible: &Self,
        m: &Modulus,
    ) -> Self {
        let r = self.degree() + that.degree();
        let mut ret = Self::zeros(r + 1);

        for i in 0..=r {
            for j in 0..=i {
                let x = m.mul(&self.get(j), &that.get(i - j));
                let k = m.add(&ret.get(i), &x);
                ret.set(i, k);
            }
        }

        // Listing 8.7
        if r < irreducible.degree() {
            return ret;
        }

        let mut out = Self::zeros(irreducible.degree());
        for i in 0..irreducible.degree() {
            out.set(i, ret.get(i))
        }

        let table = irreducible.mulprep(m);
        for i in irreducible.degree()..=r {
            for j in 0..irreducible.degree() {
                let t = m.mul(&ret.get(i), &table[i][j]);
                let k = m.add(&out.get(j), &t);
                out.set(j, k);
            }
        }

        out.trim()
    }

    fn mulprep(&self, m: &Modulus) -> Vec<Vec<Int>> {
        let degree = self.degree();
        let mut table = vec![vec![Int::ZERO; degree]; degree * 2];

        // Listing 8.2
        for i in 0..degree {
            table[i][i] = Int::from(1);
        }

        // Listing 8.3
        let norm = self.normal(m);
        for j in 0..degree {
            table[degree][j] = m.neg(&norm.get(j));
        }

        // Listing 8.4
        for i in (degree + 1)..(degree * 2) {
            for j in 1..degree {
                let t = m.mul(
                    &table[degree][j],
                    &table[i - 1][degree - 1],
                );
                table[i][j] = m.add(&table[i - 1][j - 1], &t);
            }
            table[i][0] =
                m.mul(&table[degree][0], &table[i - 1][degree - 1]);
        }

        table
    }

    // Modular subtraction (coefficients reduced mod p)
    fn sub_mod(&self, that: &Self, m: &Modulus) -> Self {
        let len = self.coef.len().max(that.coef.len());
        let mut ret = Self::zeros(len);
        for i in 0..len {
            ret.set(i, m.sub(&self.get(i), &that.get(i)));
        }
        ret.trim()
    }

    // Multiply all coefficients by scalar s mod p
    fn scale(&self, s: &Int, m: &Modulus) -> Self {
        let mut ret = self.clone();
        for k in ret.coef.iter_mut() {
            *k = m.mul(k, s);
        }
        ret.trim()
    }

    // Listing 10.1: Euclidean division a/b → (quotient, remainder)
    pub fn euclid(&self, b: &Self, m: &Modulus) -> (Self, Self) {
        let mut r = self.clone().trim();
        let mut q = Self::zeros(1);

        if r.is_zero() || b.degree() > r.degree() {
            return (q, r);
        }

        while !r.is_zero() && r.degree() >= b.degree() {
            let j = r.degree() - b.degree();
            let s = m
                .div(&r.get(r.degree()), &b.get(b.degree()))
                .unwrap();
            q.set(j, s.clone());
            for i in 0..=b.degree() {
                let sb = m.mul(&s, &b.get(i));
                let ri = m.sub(&r.get(i + j), &sb);
                r.set(i + j, ri);
            }
            r = r.trim();
        }

        (q.trim(), r)
    }

    // Listing 10.2 + 10.3: GCD of two polynomials
    pub fn gcd(&self, b: &Self, m: &Modulus) -> Self {
        if self.is_zero() {
            return b.clone();
        }
        if b.is_zero() {
            return self.clone();
        }

        let (mut aw, mut bw) = if self.degree() >= b.degree() {
            (self.clone(), b.clone())
        } else {
            (b.clone(), self.clone())
        };

        while !bw.is_zero() && bw.degree() > 0 {
            let (_, r) = aw.euclid(&bw, m);
            aw = bw;
            bw = r;
        }

        if bw.is_zero() {
            aw
        } else {
            bw
        }
    }

    // Listing 10.4 + 10.5: Inverse of self modulo irreducible
    pub fn inv(&self, irreducible: &Self, m: &Modulus) -> Self {
        let mut v = irreducible.clone();
        // Reduce mod p so that a leading coefficient ≡ 0 mod p is trimmed.
        let self_reduced = {
            let mut tmp = self.clone();
            for c in tmp.coef.iter_mut() {
                *c = m.add(c, &Int::ZERO);
            }
            tmp.trim()
        };
        let mut u = self_reduced.normal(m);
        let mut w = Self::zeros(1);
        let mut y = Self::zeros(1);
        y.set(
            0,
            m.inv(&self_reduced.get(self_reduced.degree())).unwrap(),
        );
        let one = {
            let mut p = Self::zeros(1);
            p.set(0, Int::from(1));
            p
        };

        loop {
            let (q, r) = v.euclid(&u, m);
            if r.is_zero() {
                break;
            }
            let rho = m.inv(&r.get(r.degree())).unwrap();
            let r = r.normal(m);
            // t = ρ(w − q·y)
            let qy = q.mul(&y, irreducible, m);
            let t = w.sub_mod(&qy, m).scale(&rho, m);
            w = y;
            y = t;
            v = u;
            u = r.clone();
            if r == one {
                break;
            }
        }
        y
    }

    // Listing 10.6: Division self/other modulo irreducible
    pub fn div(
        &self,
        other: &Self,
        irreducible: &Self,
        m: &Modulus,
    ) -> Self {
        let inv = other.inv(irreducible, m);
        self.mul(&inv, irreducible, m)
    }

    // Listing 9.1 + 9.2: Polynomial exponentiation g^k mod irreducible
    // Uses square-and-multiply (MSB to LSB)
    pub fn pow(
        &self,
        k: &Int,
        irreducible: &Self,
        m: &Modulus,
    ) -> Self {
        let bits = k.significant_bits();
        if bits == 0 {
            let mut one = Self::zeros(1);
            one.set(0, Int::from(1));
            return one;
        }

        let mut r = self.clone();
        for i in (0..bits - 1).rev() {
            let flag = k.get_bit(i);
            // Square r, then optionally multiply by self
            let tmp = r.mul(&r, irreducible, m);
            r = if flag {
                tmp.mul(self, irreducible, m)
            } else {
                tmp
            };
        }
        r
    }

    // Listing 9.4: Raise polynomial to field prime power (x^p)
    pub fn xp(&self, irreducible: &Self, m: &Modulus) -> Self {
        self.pow(&m.n, irreducible, m)
    }

    // Listing 9.5: raise a polynomial to `(p-1)/2` where `p` is the *base*
    // prime — not the GF(p^k) Euler criterion (that exponent is `(p^k-1)/2`).
    // This is the half-power used by the Cantor–Zassenhaus split step over
    // F_p: `(x + r)^((p-1)/2) mod f` is +1 or −1 at each F_p-root of `f`,
    // which lets a random `r` separate roots into two halves.
    pub fn gpow_p2(&self, irreducible: &Self, m: &Modulus) -> Self {
        let exp = Int::from(&m.n - 1) / 2;
        self.pow(&exp, irreducible, m)
    }

    // Listing 11.4-11.7: Find an irreducible trinomial x^n + x + a₀ over F_p
    // Uses Ben-Or's algorithm: test gcd(x^(p^i) - x, f) = 1 for i = 1..n/2
    pub fn find_irreducible(n: usize, m: &Modulus) -> Option<Self> {
        let mlimt = n / 2;

        let mut x_poly = Self::zeros(2);
        x_poly.set(1, Int::from(1));

        let mut j = Int::from(2);
        while j < m.n {
            // Candidate: r = x^n + x + j
            let mut r = Self::zeros(n + 1);
            r.set(n, Int::from(1));
            r.set(1, Int::from(1));
            r.set(0, j.clone());

            let mut is_irrd = true;
            let mut xp_cur = x_poly.clone();

            for _ in 0..mlimt {
                xp_cur = xp_cur.xp(&r, m);
                let xpm1 = xp_cur.sub_mod(&x_poly, m);
                let g = xpm1.gcd(&r, m);
                if !g.is_zero() && g.degree() > 0 {
                    is_irrd = false;
                    break;
                }
            }

            if is_irrd {
                return Some(r);
            }

            j += 1;
        }

        None
    }

    // Listing 8.9
    fn normal(&self, m: &Modulus) -> Self {
        let mut ret = self.clone();
        // Reduce every coefficient mod p first, so trim drops coefficients
        // that are 0 mod p but happened to be stored as un-reduced ints.
        for c in ret.coef.iter_mut() {
            *c = m.add(c, &Int::ZERO);
        }
        ret = ret.trim();

        let k = ret.get(ret.degree());
        if ret.is_zero() || k == 1 {
            return ret;
        }

        let d = m.inv(&k).unwrap();
        for k in ret.coef.iter_mut() {
            *k = m.mul(k, &d);
        }
        ret
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_polynomial_example() {
        // Listing 8.{11,12}
        const M: i32 = 7;

        let m = Modulus::new(&Int::from(M));
        let mut p = Polynomial::zeros(5);
        p.set(4, Int::from(1));
        p.set(3, Int::from(2));
        p.set(2, Int::from(1));
        p.set(1, Int::from(3));
        p.set(0, Int::from(5));

        let table = p.mulprep(&m);
        for row in &table {
            println!("{:?}", row);
        }

        assert_eq!(table[4], vec![M - 5, M - 3, M - 1, M - 2]);
        assert_eq!(table[7], vec![0, 6, 1, 5]);
    }

    // Chapter 9 tests — tiny example: p = 43, irreducible = x² + x + 3
    fn tiny_setup() -> (Modulus, Polynomial) {
        let m = Modulus::new(&Int::from(43));
        let mut irrd = Polynomial::zeros(3);
        irrd.set(2, Int::from(1));
        irrd.set(1, Int::from(1));
        irrd.set(0, Int::from(3));
        (m, irrd)
    }

    #[test]
    fn test_poly_pow_x5() {
        // Section 9.1 worked example: x^5 mod (x²+x+3) mod 43 = x + 28
        let (m, irrd) = tiny_setup();
        let mut x = Polynomial::zeros(2);
        x.set(1, Int::from(1));

        let result = x.pow(&Int::from(5), &irrd, &m);
        assert_eq!(result.get(0), Int::from(28));
        assert_eq!(result.get(1), Int::from(1));
    }

    #[test]
    fn test_poly_pow_listing_9_3() {
        // Listing 9.3 output: (11x + 3)^25 mod (x²+x+3) mod 43 = 3x + 26
        let (m, irrd) = tiny_setup();
        let mut tst = Polynomial::zeros(2);
        tst.set(1, Int::from(11));
        tst.set(0, Int::from(3));

        let result = tst.pow(&Int::from(25), &irrd, &m);
        assert_eq!(result.get(0), Int::from(26));
        assert_eq!(result.get(1), Int::from(3));
    }

    #[test]
    fn test_poly_pow_fermat() {
        // Fermat's Little Theorem: a^(p^k - 1) = 1
        // p=43, k=2 (degree of irreducible), p^k - 1 = 1848
        let (m, irrd) = tiny_setup();
        let mut tst = Polynomial::zeros(2);
        tst.set(1, Int::from(11));
        tst.set(0, Int::from(3));

        let result = tst.pow(&Int::from(1848), &irrd, &m);
        assert_eq!(result.get(0), Int::from(1));
        assert!(result.get(1).is_zero());
    }

    #[test]
    fn test_poly_xp() {
        // x^p should work via xp() shorthand
        let (m, irrd) = tiny_setup();
        let mut x = Polynomial::zeros(2);
        x.set(1, Int::from(1));

        let via_xp = x.xp(&irrd, &m);
        let via_pow = x.pow(&Int::from(43), &irrd, &m);
        assert_eq!(via_xp, via_pow);
    }

    // Chapter 10 tests

    #[test]
    fn test_euclid() {
        // Section 10.1 example: (x³+2x²+3x+1) / (5x³+x+6) mod 7
        // Expected: q = 3, r = 2x²+4
        let m = Modulus::new(&Int::from(7));
        let mut a = Polynomial::zeros(4);
        a.set(3, Int::from(1));
        a.set(2, Int::from(2));
        a.set(1, Int::from(3));
        a.set(0, Int::from(1));

        let mut b = Polynomial::zeros(4);
        b.set(3, Int::from(5));
        b.set(1, Int::from(1));
        b.set(0, Int::from(6));

        let (q, r) = a.euclid(&b, &m);
        assert_eq!(q.get(0), Int::from(3));
        assert!(q.get(1).is_zero());
        assert_eq!(r.get(2), Int::from(2));
        assert_eq!(r.get(0), Int::from(4));
        assert!(r.get(1).is_zero());
    }

    #[test]
    fn test_gcd() {
        // Section 10.1 example: gcd(x²−x+4, x²+12x+1) mod 17
        // Factors: (x+5)(x−6) and (x+5)(x+7), so gcd ∝ (x+5)
        let m = Modulus::new(&Int::from(17));
        let mut a = Polynomial::zeros(3);
        a.set(2, Int::from(1));
        a.set(1, Int::from(16)); // -1 mod 17
        a.set(0, Int::from(4));

        let mut b = Polynomial::zeros(3);
        b.set(2, Int::from(1));
        b.set(1, Int::from(12));
        b.set(0, Int::from(1));

        let g = a.gcd(&b, &m);
        // Result is c*(x+5) for some constant c — normalize to check
        let g = g.normal(&m);
        assert_eq!(g.get(1), Int::from(1));
        assert_eq!(g.get(0), Int::from(5));
    }

    #[test]
    fn test_inv() {
        // Section 10.2 example: inv(x+17) mod (x²+x+3) mod 43 = 5x+6
        let (m, irrd) = tiny_setup();
        let mut b = Polynomial::zeros(2);
        b.set(1, Int::from(1));
        b.set(0, Int::from(17));

        let result = b.inv(&irrd, &m);
        assert_eq!(result.get(1), Int::from(5));
        assert_eq!(result.get(0), Int::from(6));

        // Verify: (x+17) * (5x+6) ≡ 1 mod (x²+x+3) mod 43
        let product = b.mul(&result, &irrd, &m);
        assert_eq!(product.get(0), Int::from(1));
        assert!(product.get(1).is_zero());
    }

    #[test]
    fn test_div() {
        // a/b * b ≡ a mod irreducible
        let (m, irrd) = tiny_setup();
        let mut a = Polynomial::zeros(2);
        a.set(1, Int::from(11));
        a.set(0, Int::from(3));

        let mut b = Polynomial::zeros(2);
        b.set(1, Int::from(1));
        b.set(0, Int::from(17));

        let quotient = a.div(&b, &irrd, &m);
        let back = quotient.mul(&b, &irrd, &m);
        assert_eq!(back.get(0), a.get(0));
        assert_eq!(back.get(1), a.get(1));
    }

    // Chapter 11 tests

    #[test]
    fn test_find_irreducible_tiny() {
        // Book's tiny example: x² + x + 3 is irreducible mod 43
        // a₀=2 is reducible (has root 24), a₀=3 is the first irreducible
        let m = Modulus::new(&Int::from(43));
        let r = Polynomial::find_irreducible(2, &m).unwrap();
        assert_eq!(r.get(2), Int::from(1));
        assert_eq!(r.get(1), Int::from(1));
        assert_eq!(r.get(0), Int::from(3));
    }

    #[test]
    fn test_find_irreducible_exercise_11_1() {
        // Exercise 11.1: x^7 + x + 7 is irreducible mod 29
        let m = Modulus::new(&Int::from(29));
        let r = Polynomial::find_irreducible(7, &m).unwrap();
        assert_eq!(r.get(7), Int::from(1));
        assert_eq!(r.get(1), Int::from(1));
        assert_eq!(r.get(0), Int::from(7));
    }

    #[test]
    fn test_poly_gpow_p2() {
        // g^((p-1)/2) is the Euler criterion (quadratic residue test)
        // For p=43: (p-1)/2 = 21
        // (11x + 3)^21 should give ±1 or 0 in some sense,
        // but more importantly (g^21)^2 = g^42 = g^(p-1)
        let (m, irrd) = tiny_setup();
        let mut tst = Polynomial::zeros(2);
        tst.set(1, Int::from(11));
        tst.set(0, Int::from(3));

        let half = tst.gpow_p2(&irrd, &m);
        let full = half.mul(&half, &irrd, &m);
        let direct = tst.pow(&Int::from(42), &irrd, &m);
        assert_eq!(full, direct);
    }
}
