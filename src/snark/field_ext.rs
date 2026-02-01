//! Extension field polynomial types for pairing-based cryptography.
//!
//! Maps to C types: POLY, POLY_POINT, POLY_CURVE

use rug::Integer as Int;
use std::fmt;

use crate::modulus::Modulus;

/// Maximum degree for extension field polynomials (matches C MAXDEGREE)
pub const MAX_DEGREE: usize = 32;

/// Extension field polynomial with coefficients in a prime field.
/// Represents elements in F_p^k where k is the extension degree.
#[derive(Clone)]
pub struct Poly {
    pub deg: usize,
    pub coef: Vec<Int>,
}

impl Default for Poly {
    fn default() -> Self {
        Self::new()
    }
}

impl Poly {
    /// Create a new zero polynomial
    pub fn new() -> Self {
        Self {
            deg: 0,
            coef: vec![Int::ZERO; MAX_DEGREE],
        }
    }

    /// Create polynomial from coefficients (lowest degree first)
    pub fn from_coefs(coefs: &[Int]) -> Self {
        let mut p = Self::new();
        for (i, c) in coefs.iter().enumerate() {
            p.coef[i] = c.clone();
        }
        p.deg = if coefs.is_empty() { 0 } else { coefs.len() - 1 };
        p.normalize();
        p
    }

    /// Create polynomial with single coefficient (constant)
    pub fn constant(c: Int) -> Self {
        let mut p = Self::new();
        p.coef[0] = c;
        p.deg = 0;
        p
    }

    /// Normalize polynomial degree (find actual highest non-zero coefficient)
    pub fn normalize(&mut self) {
        while self.deg > 0 && self.coef[self.deg].is_zero() {
            self.deg = self.deg.saturating_sub(1);
        }
    }

    /// Check if polynomial is zero
    pub fn is_zero(&self) -> bool {
        self.deg == 0 && self.coef[0].is_zero()
    }

    /// Copy polynomial
    pub fn copy_from(&mut self, other: &Poly) {
        self.deg = other.deg;
        for i in 0..=other.deg {
            self.coef[i] = other.coef[i].clone();
        }
    }

    /// Polynomial addition in extension field
    pub fn add(&self, other: &Poly, m: &Modulus) -> Poly {
        let mut result = Poly::new();
        result.deg = self.deg.max(other.deg);

        for i in 0..=result.deg {
            let a = if i <= self.deg {
                &self.coef[i]
            } else {
                &Int::ZERO
            };
            let b = if i <= other.deg {
                &other.coef[i]
            } else {
                &Int::ZERO
            };
            result.coef[i] = m.add(a, b);
        }
        result.normalize();
        result
    }

    /// Polynomial subtraction in extension field
    pub fn sub(&self, other: &Poly, m: &Modulus) -> Poly {
        let mut result = Poly::new();
        result.deg = self.deg.max(other.deg);

        for i in 0..=result.deg {
            let a = if i <= self.deg {
                &self.coef[i]
            } else {
                &Int::ZERO
            };
            let b = if i <= other.deg {
                &other.coef[i]
            } else {
                &Int::ZERO
            };
            result.coef[i] = m.sub(a, b);
        }
        result.normalize();
        result
    }

    /// Polynomial multiplication with reduction by irreducible polynomial
    pub fn mul(
        &self,
        other: &Poly,
        irrd: &Poly,
        m: &Modulus,
    ) -> Poly {
        let deg_sum = self.deg + other.deg;

        // First compute full product
        let mut temp = vec![Int::ZERO; deg_sum + 1];
        for i in 0..=self.deg {
            for j in 0..=other.deg {
                let product = m.mul(&self.coef[i], &other.coef[j]);
                temp[i + j] = m.add(&temp[i + j], &product);
            }
        }

        // Reduce by irreducible polynomial if needed
        if deg_sum < irrd.deg {
            let mut result = Poly::new();
            result.deg = deg_sum;
            for i in 0..=deg_sum {
                result.coef[i] = temp[i].clone();
            }
            return result;
        }

        // Use precomputed reduction table
        let table = irrd.mulprep(m);
        let mut result = Poly::new();
        result.deg = irrd.deg - 1;

        // Copy lower coefficients
        for i in 0..irrd.deg {
            result.coef[i] = temp[i].clone();
        }

        // Reduce higher coefficients
        for i in irrd.deg..=deg_sum {
            for j in 0..irrd.deg {
                let t = m.mul(&temp[i], &table[i - irrd.deg][j]);
                result.coef[j] = m.add(&result.coef[j], &t);
            }
        }

        result.normalize();
        result
    }

    /// Precompute reduction table for polynomial multiplication
    /// Returns table[i][j] = coefficient j of x^(deg+i) mod irrd
    fn mulprep(&self, m: &Modulus) -> Vec<Vec<Int>> {
        let deg = self.deg;
        let mut table = vec![vec![Int::ZERO; deg]; deg];

        // First row: x^deg mod irrd = -c_{deg-1}x^{deg-1} - ... - c_0
        // (normalized so leading coef is 1)
        let norm = self.normal(m);
        for j in 0..deg {
            table[0][j] = m.neg(&norm.coef[j]);
        }

        // Subsequent rows: multiply by x and reduce
        for i in 1..deg {
            // x * table[i-1] = shift right and reduce
            for j in (1..deg).rev() {
                let t = m.mul(&table[0][j], &table[i - 1][deg - 1]);
                table[i][j] = m.add(&table[i - 1][j - 1], &t);
            }
            table[i][0] = m.mul(&table[0][0], &table[i - 1][deg - 1]);
        }

        table
    }

    /// Normalize polynomial (make leading coefficient 1)
    fn normal(&self, m: &Modulus) -> Poly {
        let mut result = self.clone();
        if result.is_zero() {
            return result;
        }

        let lead = &result.coef[result.deg];
        if lead == &Int::from(1) {
            return result;
        }

        if let Some(inv) = m.inv(lead) {
            for i in 0..=result.deg {
                result.coef[i] = m.mul(&result.coef[i], &inv);
            }
        }
        result
    }

    /// Polynomial division in extension field
    pub fn div(
        &self,
        other: &Poly,
        irrd: &Poly,
        m: &Modulus,
    ) -> Option<Poly> {
        // a / b = a * b^(-1)
        let inv = other.inv(irrd, m)?;
        Some(self.mul(&inv, irrd, m))
    }

    /// Polynomial inverse in extension field using extended Euclidean algorithm
    pub fn inv(&self, irrd: &Poly, m: &Modulus) -> Option<Poly> {
        if self.is_zero() {
            return None;
        }

        // Extended GCD to find inverse
        let mut old_r = irrd.clone();
        let mut r = self.clone();
        let mut old_s = Poly::new();
        let mut s = Poly::constant(Int::from(1));

        while !r.is_zero() {
            let (q, rem) = poly_divmod(&old_r, &r, m);
            old_r = r;
            r = rem;

            let new_s = old_s.sub(&q.mul(&s, irrd, m), m);
            old_s = s;
            s = new_s;
        }

        // old_r should be constant (gcd = 1 for irreducible)
        if old_r.deg != 0 {
            return None;
        }

        // Scale by inverse of the constant
        if let Some(inv) = m.inv(&old_r.coef[0]) {
            let mut result = Poly::new();
            result.deg = old_s.deg;
            for i in 0..=old_s.deg {
                result.coef[i] = m.mul(&old_s.coef[i], &inv);
            }
            result.normalize();
            Some(result)
        } else {
            None
        }
    }

    /// Polynomial exponentiation using square-and-multiply
    pub fn pow(&self, exp: &Int, irrd: &Poly, m: &Modulus) -> Poly {
        if exp.is_zero() {
            return Poly::constant(Int::from(1));
        }

        let mut result = Poly::constant(Int::from(1));
        let mut base = self.clone();

        let bits = exp.significant_bits();
        for i in 0..bits {
            if exp.get_bit(i) {
                result = result.mul(&base, irrd, m);
            }
            if i + 1 < bits {
                base = base.mul(&base, irrd, m);
            }
        }

        result
    }
}

impl PartialEq for Poly {
    /// Compare two polynomials
    fn eq(&self, other: &Poly) -> bool {
        if self.deg != other.deg {
            return false;
        }
        for i in 0..=self.deg {
            if self.coef[i] != other.coef[i] {
                return false;
            }
        }
        true
    }
}

impl fmt::Debug for Poly {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Poly[")?;
        for i in 0..=self.deg {
            if i > 0 {
                write!(f, ", ")?;
            }
            write!(f, "{}", self.coef[i])?;
        }
        write!(f, "]")
    }
}

/// Polynomial division with remainder
fn poly_divmod(a: &Poly, b: &Poly, m: &Modulus) -> (Poly, Poly) {
    if b.is_zero() {
        panic!("Division by zero polynomial");
    }

    if a.deg < b.deg {
        return (Poly::new(), a.clone());
    }

    let mut q = Poly::new();
    let mut r = a.clone();

    let b_lead_inv = m
        .inv(&b.coef[b.deg])
        .expect("Leading coef should be invertible");

    while !r.is_zero() && r.deg >= b.deg {
        let coef = m.mul(&r.coef[r.deg], &b_lead_inv);
        let deg_diff = r.deg - b.deg;

        q.coef[deg_diff] = m.add(&q.coef[deg_diff], &coef);
        if deg_diff > q.deg {
            q.deg = deg_diff;
        }

        for i in 0..=b.deg {
            let term = m.mul(&coef, &b.coef[i]);
            r.coef[i + deg_diff] =
                m.sub(&r.coef[i + deg_diff], &term);
        }
        r.normalize();
    }

    (q, r)
}

/// Point on an extension field curve
/// Maps to C type: POLY_POINT
#[derive(Clone)]
pub struct PolyPoint {
    pub x: Poly,
    pub y: Poly,
}

impl Default for PolyPoint {
    fn default() -> Self {
        Self::new()
    }
}

impl PolyPoint {
    /// Create a new point at infinity (0, 0)
    pub fn new() -> Self {
        Self {
            x: Poly::new(),
            y: Poly::new(),
        }
    }

    /// Create point from coordinates
    pub fn from_coords(x: Poly, y: Poly) -> Self {
        Self { x, y }
    }

    /// Check if point is at infinity
    pub fn is_inf(&self) -> bool {
        self.x.is_zero() && self.y.is_zero()
    }

    /// Copy from another point
    pub fn copy_from(&mut self, other: &PolyPoint) {
        self.x.copy_from(&other.x);
        self.y.copy_from(&other.y);
    }

    /// Point negation (same x, negate y)
    pub fn neg(&self, m: &Modulus) -> PolyPoint {
        let mut neg_y = Poly::new();
        neg_y.deg = self.y.deg;
        for i in 0..=self.y.deg {
            neg_y.coef[i] = m.neg(&self.y.coef[i]);
        }
        PolyPoint {
            x: self.x.clone(),
            y: neg_y,
        }
    }
}

impl PartialEq for PolyPoint {
    /// Check equality
    fn eq(&self, other: &PolyPoint) -> bool {
        self.x.eq(&other.x) && self.y.eq(&other.y)
    }
}

impl fmt::Debug for PolyPoint {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "PolyPoint {{ x: {:?}, y: {:?} }}", self.x, self.y)
    }
}

/// Curve over extension field: y^2 = x^3 + a4*x + a6
/// Maps to C type: POLY_CURVE
#[derive(Clone)]
pub struct PolyCurve {
    pub a4: Poly,
    pub a6: Poly,
}

impl Default for PolyCurve {
    fn default() -> Self {
        Self::new()
    }
}

impl PolyCurve {
    /// Create a new curve with zero coefficients
    pub fn new() -> Self {
        Self {
            a4: Poly::new(),
            a6: Poly::new(),
        }
    }

    /// Create curve from coefficients
    pub fn from_coeffs(a4: Poly, a6: Poly) -> Self {
        Self { a4, a6 }
    }

    /// Evaluate y^2 = x^3 + a4*x + a6 at point x
    pub fn eval(&self, x: &Poly, irrd: &Poly, m: &Modulus) -> Poly {
        // x^2
        let x2 = x.mul(x, irrd, m);
        // x^3
        let x3 = x2.mul(x, irrd, m);
        // a4 * x
        let a4x = self.a4.mul(x, irrd, m);
        // x^3 + a4*x
        let sum = x3.add(&a4x, m);
        // x^3 + a4*x + a6
        sum.add(&self.a6, m)
    }

    /// Point addition on extension curve
    pub fn add(
        &self,
        p: &PolyPoint,
        q: &PolyPoint,
        irrd: &Poly,
        m: &Modulus,
    ) -> PolyPoint {
        if p.is_inf() {
            return q.clone();
        }
        if q.is_inf() {
            return p.clone();
        }

        // Check if P = -Q (same x, y sums to 0)
        let y_sum = p.y.add(&q.y, m);
        if p.x.eq(&q.x) && y_sum.is_zero() {
            return PolyPoint::new(); // Point at infinity
        }

        let (num, den) = if p.x.eq(&q.x) {
            // Point doubling: lambda = (3*x^2 + a4) / (2*y)
            let x2 = p.x.mul(&p.x, irrd, m);
            let three = Poly::constant(Int::from(3));
            let two = Poly::constant(Int::from(2));
            let num = three.mul(&x2, irrd, m).add(&self.a4, m);
            let den = two.mul(&p.y, irrd, m);
            (num, den)
        } else {
            // Point addition: lambda = (y2 - y1) / (x2 - x1)
            let num = q.y.sub(&p.y, m);
            let den = q.x.sub(&p.x, m);
            (num, den)
        };

        // lambda = num / den
        let lambda = match num.div(&den, irrd, m) {
            Some(l) => l,
            None => return PolyPoint::new(), // Division failed, return infinity
        };

        // x3 = lambda^2 - x1 - x2
        let lambda2 = lambda.mul(&lambda, irrd, m);
        let x3 = lambda2.sub(&p.x, m).sub(&q.x, m);

        // y3 = lambda * (x1 - x3) - y1
        let x_diff = p.x.sub(&x3, m);
        let y3 = lambda.mul(&x_diff, irrd, m).sub(&p.y, m);

        PolyPoint { x: x3, y: y3 }
    }

    /// Scalar multiplication on extension curve
    pub fn mul(
        &self,
        p: &PolyPoint,
        k: &Int,
        irrd: &Poly,
        m: &Modulus,
    ) -> PolyPoint {
        if p.is_inf() || k.is_zero() {
            return PolyPoint::new();
        }

        let mut result = PolyPoint::new();
        let mut base = p.clone();

        let bits = k.significant_bits();
        for i in 0..bits {
            if k.get_bit(i) {
                result = self.add(&result, &base, irrd, m);
            }
            if i + 1 < bits {
                base = self.add(&base, &base, irrd, m);
            }
        }

        result
    }
}

impl fmt::Debug for PolyCurve {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "PolyCurve {{ a4: {:?}, a6: {:?} }}",
            self.a4, self.a6
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_poly_basic() {
        let m = Modulus::new(&Int::from(7));

        let a = Poly::from_coefs(&[
            Int::from(1),
            Int::from(2),
            Int::from(3),
        ]);
        let b = Poly::from_coefs(&[Int::from(4), Int::from(5)]);

        let sum = a.add(&b, &m);
        assert_eq!(sum.coef[0], Int::from(5)); // 1 + 4 = 5
        assert_eq!(sum.coef[1], Int::from(0)); // 2 + 5 = 7 = 0 mod 7
        assert_eq!(sum.coef[2], Int::from(3)); // 3 + 0 = 3
    }

    #[test]
    fn test_poly_mul_no_reduction() {
        let m = Modulus::new(&Int::from(7));
        let irrd = Poly::from_coefs(&[
            Int::from(1),
            Int::from(0),
            Int::from(0),
            Int::from(0),
            Int::from(1),
        ]); // x^4 + 1

        let a = Poly::from_coefs(&[Int::from(1), Int::from(1)]); // 1 + x
        let b = Poly::from_coefs(&[Int::from(1), Int::from(1)]); // 1 + x

        let prod = a.mul(&b, &irrd, &m);
        // (1+x)^2 = 1 + 2x + x^2
        assert_eq!(prod.coef[0], Int::from(1));
        assert_eq!(prod.coef[1], Int::from(2));
        assert_eq!(prod.coef[2], Int::from(1));
    }
}
