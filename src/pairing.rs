//! Bilinear pairings — Miller's algorithm, Weil and Tate pairings, plus
//! the trace-based extension cardinality formula.
//!
//! Port of `pairing.c` (Chapter 15).

use rug::ops::Pow;
use rug::Integer as Int;

use crate::elliptic::{Curve, Point};
use crate::modulus::Modulus;
use crate::poly_elliptic::{PolyCurve, PolyPoint};
use crate::polynomial::Polynomial;

/// Schoof / weil-conjecture style: given trace of Frobenius `t` and embedding
/// degree `k`, return the cardinality of E(F_{p^k}).
///
/// Computes v_i with v_0 = 2, v_1 = t, v_i = t*v_{i-1} - p*v_{i-2}, then
/// returns p^k + 1 − v_k.
pub fn cardinality(t: &Int, k: u32, p: &Int) -> Int {
    let mut v0 = Int::from(2);
    let mut v1 = t.clone();
    if k == 0 {
        return Int::from(p.clone().pow(k)) + 1 - v0;
    }
    if k == 1 {
        return Int::from(p.clone().pow(1u32)) + 1 - v1;
    }
    let mut vk = v1.clone();
    for _ in 2..=k {
        vk = Int::from(t * &v1) - Int::from(p * &v0);
        v0 = v1;
        v1 = vk.clone();
    }
    Int::from(p.clone().pow(k)) + 1 - vk
}

/// Miller's line/slope function: `h_{P,Q}(R)` per AEC p. 394.
///
/// If P or Q is infinity, returns 1. If P = -Q (vertical line case),
/// returns x_R - x_P. Otherwise returns the line-through-PQ / vertical-line
/// quotient evaluated at R, used to accumulate the Miller function `f_P`.
fn hpq(
    ec: &PolyCurve,
    p: &PolyPoint,
    q: &PolyPoint,
    r: &PolyPoint,
) -> Polynomial {
    if p.is_inf() || q.is_inf() {
        return one_poly();
    }
    let pm = &ec.p;
    let irrd = &ec.irrd;

    let b = p.y.add(&q.y);
    let (num, den) = if is_zero_mod(&b, pm) {
        let x_diff = poly_sub_mod(&p.x, &q.x, pm);
        if is_zero_mod(&x_diff, pm) {
            // P == -Q: line is vertical. h = x_R - x_P.
            return poly_sub_mod(&r.x, &p.x, pm);
        }
        // Different points whose y's happen to cancel modulo something —
        // fall back to the standard (y2 - y1)/(x2 - x1) form.
        let num = poly_sub_mod(&p.y, &q.y, pm);
        (num, x_diff)
    } else {
        let t1 = p.x.mul(&p.x, irrd, pm);
        let t2 = p.x.mul(&q.x, irrd, pm);
        let t3 = q.x.mul(&q.x, irrd, pm);
        let num = t1.add(&t2).add(&t3).add(&ec.a4);
        (num, b)
    };
    let lambda = num.div(&den, irrd, pm);

    // num/den at R: (yR - yP - λ(xR - xP)) / (xR + xP + xQ - λ²).
    let ry_yp = poly_sub_mod(&r.y, &p.y, pm);
    let rx_xp = poly_sub_mod(&r.x, &p.x, pm);
    let lambda_dx = rx_xp.mul(&lambda, irrd, pm);
    let top = poly_sub_mod(&ry_yp, &lambda_dx, pm);

    let lambda2 = lambda.mul(&lambda, irrd, pm);
    let mut bot = poly_sub_mod(&r.x, &lambda2, pm);
    bot = bot.add(&p.x).add(&q.x);
    top.div(&bot, irrd, pm)
}

/// Miller's algorithm: compute the function `f_{m,P}` evaluated at R, where
/// `m` is the order of P. Used as the core of both Weil and Tate pairings.
pub fn miller(
    ec: &PolyCurve,
    p: &PolyPoint,
    r: &PolyPoint,
    m: &Int,
) -> Polynomial {
    let pm = &ec.p;
    let irrd = &ec.irrd;
    let mut t = p.clone();
    let mut f = one_poly();
    if m.significant_bits() < 2 {
        return f;
    }
    let mut bit = (m.significant_bits() - 2) as i64;
    while bit >= 0 {
        let h = hpq(ec, &t, &t, r);
        f = f.mul(&f, irrd, pm);
        f = f.mul(&h, irrd, pm);
        t = ec.add(&t, &t);
        if m.get_bit(bit as u32) {
            let h = hpq(ec, &t, p, r);
            f = f.mul(&h, irrd, pm);
            t = ec.add(&t, p);
        }
        bit -= 1;
    }
    f
}

/// Weil pairing `e_m(P, Q)`. Requires an auxiliary point S of order different
/// from m (so that the four Miller evaluations stay defined). Returns an m-th
/// root of unity in GF(p^k).
pub fn weil(
    ec: &PolyCurve,
    p: &PolyPoint,
    q: &PolyPoint,
    s: &PolyPoint,
    m: &Int,
) -> Polynomial {
    let pm = &ec.p;
    let irrd = &ec.irrd;
    let qps = ec.add(q, s);
    // -S: same x, negated y.
    let neg_s = PolyPoint::new(s.x.clone(), poly_neg(&s.y, pm));
    let p_minus_s = ec.add(p, &neg_s);

    let t1 = miller(ec, p, &qps, m);
    let t2 = miller(ec, p, s, m);
    let t3 = miller(ec, q, &p_minus_s, m);
    let t4 = miller(ec, q, &neg_s, m);

    let w1 = t1.div(&t2, irrd, pm);
    let w2 = t3.div(&t4, irrd, pm);
    w1.div(&w2, irrd, pm)
}

/// Tate pairing `t_m(P, Q)`. Cheaper than Weil — one pair of Miller calls and
/// a final exponentiation to `(p^k - 1)/m`.
pub fn tate(
    ec: &PolyCurve,
    p: &PolyPoint,
    q: &PolyPoint,
    s: &PolyPoint,
    m: &Int,
) -> Polynomial {
    let pm = &ec.p;
    let irrd = &ec.irrd;
    let qps = ec.add(q, s);
    let t1 = miller(ec, p, &qps, m);
    let t2 = miller(ec, p, s, m);
    let t = t1.div(&t2, irrd, pm);
    let pk = Int::from(pm.n.clone().pow(irrd.degree() as u32));
    let exp = Int::from(&pk - 1) / m;
    t.pow(&exp, irrd, pm)
}

/// Order of `P` on a base curve given a factor list. First factor `r` such
/// that `r*P` is infinity wins. None if none match.
pub fn get_order(ec: &Curve, p: &Point, factors: &[Int]) -> Option<Int> {
    for f in factors {
        if ec.mul(p, f).is_inf() {
            return Some(f.clone());
        }
    }
    None
}

// --- helpers ---------------------------------------------------------------

fn one_poly() -> Polynomial {
    let mut p = Polynomial::zeros(1);
    p.set(0, Int::from(1));
    p
}

fn poly_neg(p: &Polynomial, m: &Modulus) -> Polynomial {
    let mut out = Polynomial::zeros(p.degree() + 1);
    for i in 0..=p.degree() {
        out.set(i, m.neg(&p.get(i)));
    }
    out.trim()
}

fn poly_sub_mod(a: &Polynomial, b: &Polynomial, m: &Modulus) -> Polynomial {
    let n = a.degree().max(b.degree()) + 1;
    let mut out = Polynomial::zeros(n);
    for i in 0..n {
        out.set(i, m.sub(&a.get(i), &b.get(i)));
    }
    out.trim()
}

fn is_zero_mod(p: &Polynomial, m: &Modulus) -> bool {
    for i in 0..=p.degree() {
        if !m.add(&p.get(i), &Int::ZERO).is_zero() {
            return false;
        }
    }
    true
}

#[cfg(test)]
mod tests {
    use super::*;

    /// y² = x³ + 23x + 42 mod 43, with embedding degree 2. The book uses
    /// trace t = −11 and reports |E(F_{p²})| = 1815 = 3·5·11².
    fn tiny_curve() -> PolyCurve {
        let m = Modulus::new(&Int::from(43));
        let irrd = Polynomial::find_irreducible(2, &m).expect("irreducible");
        let mut a4 = Polynomial::zeros(1);
        a4.set(0, Int::from(23));
        let mut a6 = Polynomial::zeros(1);
        a6.set(0, Int::from(42));
        PolyCurve::new(a4, a6, irrd, m)
    }

    fn lift(v: i64) -> Polynomial {
        let mut p = Polynomial::zeros(1);
        p.set(0, Int::from(v));
        p
    }

    /// Build a fixed independent reference point on the curve at a different
    /// order than `tor` (tor=11, here we pick order-55 or whatever falls out).
    fn pick_aux_point(ec: &PolyCurve, tor: &Int) -> PolyPoint {
        // Sweep starting x's until embed gives a point whose order is not
        // `tor` and is invertible against `tor`.
        let factors: Vec<Int> = [3, 5, 11, 15, 33, 55, 165, 1815]
            .into_iter()
            .map(Int::from)
            .collect();
        let mut start = lift(0);
        for _ in 0..200 {
            if let Some((p, _)) = ec.embed(&start, 1) {
                if let Some(ord) = ec.order(&p, &factors) {
                    if &ord != tor {
                        return p;
                    }
                }
            }
            let mut next = lift(0);
            next.set(0, start.get(0) + 1);
            // bump x.coef[0] by 1; OK for sweep purposes
            start = next;
        }
        panic!("could not find auxiliary point of order != tor");
    }

    /// Helper: find a point of exactly the given order `m` on the curve.
    fn pick_point_of_order(ec: &PolyCurve, m: &Int) -> PolyPoint {
        let factors: Vec<Int> = [3, 5, 11, 15, 33, 55, 165, 1815]
            .into_iter()
            .map(Int::from)
            .collect();
        let mut start = lift(0);
        for _ in 0..400 {
            if let Some((p, _)) = ec.embed(&start, 1) {
                if let Some(ord) = ec.order(&p, &factors) {
                    if &ord == m {
                        return p;
                    }
                }
            }
            let mut next = lift(0);
            next.set(0, start.get(0) + 1);
            start = next;
        }
        panic!("could not find point of order {}", m);
    }

    #[test]
    fn test_cardinality_tiny() {
        // p = 43, k = 2, t = -11 → |E(F_{p²})| = 1815.
        let p = Int::from(43);
        let t = Int::from(-11);
        let card = cardinality(&t, 2, &p);
        assert_eq!(card, Int::from(1815));
    }

    #[test]
    fn test_cardinality_k1_matches_trace() {
        // |E(F_p)| = p + 1 - t.
        let p = Int::from(43);
        let t = Int::from(7);
        let card = cardinality(&t, 1, &p);
        assert_eq!(card, Int::from(43) + 1 - 7);
    }

    #[test]
    fn test_weil_pairing_alternating() {
        // Alternating: e(P, P) = 1 (since (P, P) is degenerate).
        let ec = tiny_curve();
        let tor = Int::from(11);
        let p = pick_point_of_order(&ec, &tor);
        let s = pick_aux_point(&ec, &tor);
        let w = weil(&ec, &p, &p, &s, &tor);
        assert_eq!(w.degree(), 0);
        assert_eq!(w.get(0), Int::from(1));
    }

    #[test]
    fn test_weil_pairing_is_torsion_root() {
        // For any P, Q of order m, e(P, Q)^m == 1.
        let ec = tiny_curve();
        let tor = Int::from(11);
        let p = pick_point_of_order(&ec, &tor);
        let q = pick_point_of_order(&ec, &tor);
        let s = pick_aux_point(&ec, &tor);
        let w = weil(&ec, &p, &q, &s, &tor);
        let wm = w.pow(&tor, &ec.irrd, &ec.p);
        assert_eq!(wm.degree(), 0);
        assert_eq!(wm.get(0), Int::from(1));
    }

    #[test]
    fn test_weil_pairing_bilinear() {
        // e(P, T+Q) == e(P, T) * e(P, Q).
        let ec = tiny_curve();
        let tor = Int::from(11);
        let p = pick_point_of_order(&ec, &tor);
        let q = pick_point_of_order(&ec, &tor);
        let t = pick_point_of_order(&ec, &tor);
        let s = pick_aux_point(&ec, &tor);

        let w_pq = weil(&ec, &p, &q, &s, &tor);
        let w_pt = weil(&ec, &p, &t, &s, &tor);
        let tpq = ec.add(&t, &q);
        let w_p_tpq = weil(&ec, &p, &tpq, &s, &tor);
        let prod = w_pt.mul(&w_pq, &ec.irrd, &ec.p);
        assert_eq!(w_p_tpq, prod);
    }

    #[test]
    fn test_tate_pairing_nontrivial() {
        let ec = tiny_curve();
        let tor = Int::from(11);
        let p = pick_point_of_order(&ec, &tor);
        let q = pick_point_of_order(&ec, &tor);
        let s = pick_aux_point(&ec, &tor);
        let t = tate(&ec, &p, &q, &s, &tor);
        let tm = t.pow(&tor, &ec.irrd, &ec.p);
        assert_eq!(tm.degree(), 0);
        assert_eq!(tm.get(0), Int::from(1));
    }

    #[test]
    fn test_tate_pairing_bilinear() {
        let ec = tiny_curve();
        let tor = Int::from(11);
        let p = pick_point_of_order(&ec, &tor);
        let q = pick_point_of_order(&ec, &tor);
        let t = pick_point_of_order(&ec, &tor);
        let s = pick_aux_point(&ec, &tor);

        let t_pq = tate(&ec, &p, &q, &s, &tor);
        let t_pt = tate(&ec, &p, &t, &s, &tor);
        let tpq = ec.add(&t, &q);
        let t_p_tpq = tate(&ec, &p, &tpq, &s, &tor);
        let prod = t_pt.mul(&t_pq, &ec.irrd, &ec.p);
        assert_eq!(t_p_tpq, prod);
    }

    #[test]
    fn test_base_curve_get_order_bn254() {
        // The BN254 base point has order equal to the curve order; check
        // [order]·G = O via get_order.
        let ec = crate::elliptic::curves::curve_bn254();
        let order = ec.order.clone();
        let found = get_order(&ec, &ec.base, std::slice::from_ref(&order));
        assert_eq!(found, Some(order));
    }
}
