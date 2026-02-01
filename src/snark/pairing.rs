//! Bilinear pairing operations for SNARK verification.
//!
//! Implements Tate and Weil pairings over extension fields.
//! Ported from pairing.c

use rug::ops::Pow;
use rug::Integer as Int;

use crate::modulus::Modulus;

use super::field_ext::{Poly, PolyCurve, PolyPoint};

/// Pairing operations over extension fields.
pub struct Pairing {
    /// Field prime modulus
    pub modulus: Modulus,
    /// Irreducible polynomial for field extension
    pub irrd: Poly,
    /// Extension field cardinality (p^k)
    pub field_card: Int,
}

impl Pairing {
    /// Create a new pairing context
    pub fn new(prime: &Int, irrd: Poly) -> Self {
        // Compute p^k where k is the extension degree
        let k = irrd.deg as u32;
        let field_card = prime.clone().pow(k);

        Self {
            modulus: Modulus::new(prime),
            irrd,
            field_card,
        }
    }

    /// Determine if a point is in G1 or G2
    /// Returns:
    ///   1 for x in group 1, y in group 1
    ///   2 for x in group 2, y in group 1
    ///   3 for x in group 1, y in group 2
    ///   4 for x in group 2, y in group 2
    pub fn g1g2(&self, p: &PolyPoint) -> i32 {
        let x_deg = p.x.deg > 0;
        let y_deg = p.y.deg > 0;

        match (x_deg, y_deg) {
            (true, true) => 4,
            (true, false) => 2,
            (false, true) => 3,
            (false, false) => 1,
        }
    }

    /// Compute Miller's h function value.
    ///
    /// Computes the line function h_{P,Q}(R) for points P, Q, R on the curve.
    /// This is the core building block of Miller's algorithm.
    fn hpq(
        &self,
        p: &PolyPoint,
        q: &PolyPoint,
        r: &PolyPoint,
        curve: &PolyCurve,
    ) -> Poly {
        let m = &self.modulus;

        // If P or Q is infinity, return 1
        if p.is_inf() || q.is_inf() {
            return Poly::constant(Int::from(1));
        }

        // Check if P + Q = 0 (vertical line case)
        let y_sum = p.y.add(&q.y, m);
        if y_sum.is_zero() {
            let x_diff = p.x.sub(&q.x, m);
            if x_diff.is_zero() {
                // P == -Q, line is x = x_P, return x_R - x_P
                return r.x.sub(&p.x, m);
            }
            // Shouldn't reach here normally
            let num = p.y.sub(&q.y, m);
            let den = x_diff;
            return self.compute_h_value(&num, &den, p, q, r);
        }

        // Compute lambda using "secure form" to avoid division by zero
        // lambda = (x_P^2 + x_P*x_Q + x_Q^2 + a4) / (y_P + y_Q)
        let xp2 = p.x.mul(&p.x, &self.irrd, m);
        let xpxq = p.x.mul(&q.x, &self.irrd, m);
        let xq2 = q.x.mul(&q.x, &self.irrd, m);

        let num = xp2.add(&xpxq, m).add(&xq2, m).add(&curve.a4, m);
        let den = p.y.add(&q.y, m);

        self.compute_h_value(&num, &den, p, q, r)
    }

    /// Helper to compute h value given numerator and denominator of slope
    fn compute_h_value(
        &self,
        num: &Poly,
        den: &Poly,
        p: &PolyPoint,
        q: &PolyPoint,
        r: &PolyPoint,
    ) -> Poly {
        let m = &self.modulus;

        // lambda = num / den
        let lambda = match num.div(den, &self.irrd, m) {
            Some(l) => l,
            None => return Poly::constant(Int::from(1)),
        };

        // Compute h = (y_R - y_P - lambda*(x_R - x_P)) / (x_R - lambda^2 + x_P + x_Q)
        let yr_yp = r.y.sub(&p.y, m);
        let xr_xp = r.x.sub(&p.x, m);
        let lambda_xr_xp = lambda.mul(&xr_xp, &self.irrd, m);
        let t = yr_yp.sub(&lambda_xr_xp, m);

        let lambda2 = lambda.mul(&lambda, &self.irrd, m);
        let b = r.x.sub(&lambda2, m).add(&p.x, m).add(&q.x, m);

        match t.div(&b, &self.irrd, m) {
            Some(h) => h,
            None => Poly::constant(Int::from(1)),
        }
    }

    /// Miller's algorithm for computing f_{m,P}(R).
    ///
    /// Computes the Miller function using the double-and-add approach.
    pub fn miller(
        &self,
        p: &PolyPoint,
        r: &PolyPoint,
        m_order: &Int,
        curve: &PolyCurve,
    ) -> Poly {
        let mod_p = &self.modulus;

        let bits = m_order.significant_bits();
        if bits == 0 {
            return Poly::constant(Int::from(1));
        }

        let mut t = p.clone();
        let mut f = Poly::constant(Int::from(1));

        // Start from second-highest bit
        for i in (0..bits - 1).rev() {
            // f = f^2 * h_{T,T}(R)
            let h = self.hpq(&t, &t, r, curve);
            f = f.mul(&f, &self.irrd, mod_p);
            f = f.mul(&h, &self.irrd, mod_p);

            // T = 2T
            t = curve.add(&t, &t, &self.irrd, mod_p);

            if m_order.get_bit(i) {
                // f = f * h_{T,P}(R)
                let h = self.hpq(&t, p, r, curve);
                f = f.mul(&h, &self.irrd, mod_p);

                // T = T + P
                t = curve.add(&t, p, &self.irrd, mod_p);
            }
        }

        f
    }

    /// Compute Weil pairing e(P, Q).
    ///
    /// Returns an element of the multiplicative group that is an m-th root of unity.
    pub fn weil(
        &self,
        p: &PolyPoint,
        q: &PolyPoint,
        s: &PolyPoint,
        m_order: &Int,
        curve: &PolyCurve,
    ) -> Poly {
        let mod_p = &self.modulus;

        // Q + S
        let q_plus_s = curve.add(q, s, &self.irrd, mod_p);

        // -S
        let minus_s = s.neg(mod_p);

        // P + (-S) = P - S
        let p_minus_s = curve.add(p, &minus_s, &self.irrd, mod_p);

        // f_P(Q+S), f_P(S), f_Q(P-S), f_Q(-S)
        let t1 = self.miller(p, &q_plus_s, m_order, curve);
        let t2 = self.miller(p, s, m_order, curve);
        let t3 = self.miller(q, &p_minus_s, m_order, curve);
        let t4 = self.miller(q, &minus_s, m_order, curve);

        // w1 = t1 / t2
        let w1 = match t1.div(&t2, &self.irrd, mod_p) {
            Some(w) => w,
            None => return Poly::constant(Int::from(1)),
        };

        // w2 = t3 / t4
        let w2 = match t3.div(&t4, &self.irrd, mod_p) {
            Some(w) => w,
            None => return Poly::constant(Int::from(1)),
        };

        // w = w1 / w2
        match w1.div(&w2, &self.irrd, mod_p) {
            Some(w) => w,
            None => Poly::constant(Int::from(1)),
        }
    }

    /// Compute Tate pairing e(P, Q).
    ///
    /// The Tate pairing is more efficient than Weil for verification.
    /// Returns an element of F_p^k / (F_p^k)^m.
    pub fn tate(
        &self,
        p: &PolyPoint,
        q: &PolyPoint,
        s: &PolyPoint,
        m_order: &Int,
        curve: &PolyCurve,
    ) -> Poly {
        let mod_p = &self.modulus;

        // Q + S
        let q_plus_s = curve.add(q, s, &self.irrd, mod_p);

        // f_P(Q+S), f_P(S)
        let t1 = self.miller(p, &q_plus_s, m_order, curve);
        let t2 = self.miller(p, s, m_order, curve);

        // t = t1 / t2
        let t = match t1.div(&t2, &self.irrd, mod_p) {
            Some(ratio) => ratio,
            None => return Poly::constant(Int::from(1)),
        };

        // Final exponentiation: t^((p^k - 1) / m)
        // This projects the result to the m-torsion subgroup
        let exp = (self.field_card.clone() - Int::from(1)) / m_order;
        t.pow(&exp, &self.irrd, mod_p)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_g1g2_classification() {
        let prime = Int::from(7);
        let irrd = Poly::from_coefs(&[
            Int::from(1),
            Int::from(0),
            Int::from(1),
        ]); // x^2 + 1

        let pairing = Pairing::new(&prime, irrd);

        // Point with both coords degree 0 -> G1
        let p1 = PolyPoint::new();
        assert_eq!(pairing.g1g2(&p1), 1);

        // Point with x degree > 0
        let mut p2 = PolyPoint::new();
        p2.x.deg = 1;
        p2.x.coef[1] = Int::from(1);
        assert_eq!(pairing.g1g2(&p2), 2);
    }
}
