//! Elliptic curves over an extension field GF(p^k).
//!
//! Port of `poly_eliptic.c` (Chapter 13). Mirrors `elliptic.rs` but the
//! coordinates are polynomials in `Polynomial` reduced by an irreducible.

use rug::ops::Pow;
use rug::Integer as Int;

use crate::modulus::Modulus;
use crate::polynomial::Polynomial;

/// A point on `y^2 = x^3 + a4*x + a6` over GF(p^k).
///
/// Infinity is tracked by an explicit flag — same idea as `elliptic::Point` —
/// so curves whose `(0, 0)` happens to satisfy the equation are not aliased.
#[derive(Clone, Debug, PartialEq)]
pub struct PolyPoint {
    pub x: Polynomial,
    pub y: Polynomial,
    inf: bool,
}

impl PolyPoint {
    pub fn new(x: Polynomial, y: Polynomial) -> Self {
        Self { x, y, inf: false }
    }

    pub fn inf() -> Self {
        Self {
            x: Polynomial::zeros(1),
            y: Polynomial::zeros(1),
            inf: true,
        }
    }

    pub fn is_inf(&self) -> bool {
        self.inf
    }

    /// Group classification used by the pairing routines:
    /// `1` if both coords live in the base field (G1),
    /// `>1` if at least one coord has nonzero degree (G2 / extension).
    /// Matches the encoding from `pairing.c::g1g2`.
    pub fn g1g2(&self) -> u32 {
        let xb = self.x.degree() > 0;
        let yb = self.y.degree() > 0;
        match (xb, yb) {
            (true, true) => 4,
            (true, false) => 2,
            (false, true) => 3,
            (false, false) => 1,
        }
    }
}

/// Curve y^2 = x^3 + a4*x + a6 over GF(p^k) defined by `irrd`.
#[derive(Clone, Debug)]
pub struct PolyCurve {
    pub a4: Polynomial,
    pub a6: Polynomial,
    pub irrd: Polynomial,
    pub p: Modulus,
}

impl PolyCurve {
    pub fn new(
        a4: Polynomial,
        a6: Polynomial,
        irrd: Polynomial,
        p: Modulus,
    ) -> Self {
        Self { a4, a6, irrd, p }
    }

    /// f(x) = x^3 + a4*x + a6.
    pub fn fofx(&self, x: &Polynomial) -> Polynomial {
        let x2 = x.mul(x, &self.irrd, &self.p);
        let x3 = x2.mul(x, &self.irrd, &self.p);
        let a4x = self.a4.mul(x, &self.irrd, &self.p);
        x3.add(&a4x).add(&self.a6)
    }

    /// Test if a point lies on the curve. The identity is on every curve.
    pub fn fits(&self, p: &PolyPoint) -> bool {
        if p.is_inf() {
            return true;
        }
        let lhs = p.y.mul(&p.y, &self.irrd, &self.p);
        let rhs = self.fofx(&p.x);
        coef_eq_mod(&lhs, &rhs, &self.p)
    }

    /// P + Q using the unified slope formula
    /// `λ = (x1² + x1·x2 + x2² + a4) / (y1 + y2)`,
    /// which also covers doubling when y1 + y2 ≠ 0. P + (-P) → ∞.
    pub fn add(&self, p: &PolyPoint, q: &PolyPoint) -> PolyPoint {
        if p.is_inf() {
            return q.clone();
        }
        if q.is_inf() {
            return p.clone();
        }

        let irrd = &self.irrd;
        let pm = &self.p;

        // Unified slope: λ = (x1² + x1·x2 + x2² + a4) / (y1 + y2).
        // Breaks when y1 + y2 ≡ 0; need to disambiguate P + (-P) vs a chord
        // whose endpoints just happen to have y-coords that cancel.
        let y_sum = p.y.add(&q.y);
        let (num, den) = if coef_zero_mod(&y_sum, pm) {
            let x_diff = poly_sub_mod(&p.x, &q.x, pm);
            if coef_zero_mod(&x_diff, pm) {
                // Same x, opposite y → P + (-P) = O. Covers 2-torsion
                // doubling (y = 0) as well.
                return PolyPoint::inf();
            }
            // x1 ≠ x2 but y1 = -y2 — fall back to standard chord slope
            // (y1 - y2)/(x1 - x2) = (q.y - p.y)/(q.x - p.x).
            let num = poly_sub_mod(&q.y, &p.y, pm);
            let den = poly_sub_mod(&q.x, &p.x, pm);
            (num, den)
        } else {
            let t1 = p.x.mul(&p.x, irrd, pm);
            let t2 = p.x.mul(&q.x, irrd, pm);
            let t3 = q.x.mul(&q.x, irrd, pm);
            let num = t1.add(&t2).add(&t3).add(&self.a4);
            (num, y_sum)
        };
        let lambda = num.div(&den, irrd, pm);

        let lambda2 = lambda.mul(&lambda, irrd, pm);
        let x_sum = p.x.add(&q.x);
        let rx = poly_sub_mod(&lambda2, &x_sum, pm);
        let dx = poly_sub_mod(&p.x, &rx, pm);
        let lx = dx.mul(&lambda, irrd, pm);
        let ry = poly_sub_mod(&lx, &p.y, pm);
        PolyPoint::new(rx, ry)
    }

    /// Q = k*P via double-and-add (MSB-first).
    pub fn mul(&self, p: &PolyPoint, k: &Int) -> PolyPoint {
        if p.is_inf() || k.is_zero() {
            return PolyPoint::inf();
        }
        let mut r = p.clone();
        let mut bit = k.significant_bits() - 1;
        while bit > 0 {
            r = self.add(&r, &r);
            if k.get_bit(bit - 1) {
                r = self.add(&r, p);
            }
            bit -= 1;
        }
        r
    }

    /// Find the two points on the curve at the lowest x ≥ start that yields a
    /// quadratic residue f(x). Returns (P, -P) with the lexicographically
    /// smaller y as the first element. None if `limit` x candidates exhaust.
    pub fn embed(
        &self,
        start: &Polynomial,
        limit: usize,
    ) -> Option<(PolyPoint, PolyPoint)> {
        let mut x = start.clone();
        for _ in 0..limit {
            let f = self.fofx(&x);
            if is_quad_residue(&f, &self.irrd, &self.p) {
                let y = poly_sqrt(&f, &self.irrd, &self.p)?;
                let neg_y = poly_neg(&y, &self.p);
                let p1 = PolyPoint::new(x.clone(), y.clone());
                let p2 = PolyPoint::new(x, neg_y);
                // Order so that the first point has the smaller leading
                // coefficient (matches the book's "smaller y" convention).
                if cmp_leading(&p2.y, &p1.y) < 0 {
                    return Some((p2, p1));
                }
                return Some((p1, p2));
            }
            x = ff_bump(&x, &self.irrd, &self.p);
        }
        None
    }

    /// Order of P given a list of candidate orders (factor sweep).
    /// Returns the first factor `r` such that `r*P = O`. None if none match.
    /// Mirrors `poly_get_order` from `pairing.c`.
    pub fn order(
        &self,
        p: &PolyPoint,
        factors: &[Int],
    ) -> Option<Int> {
        for f in factors {
            if self.mul(p, f).is_inf() {
                return Some(f.clone());
            }
        }
        None
    }
}

// --- helpers ----------------------------------------------------------------

/// Compare two polynomials for equality after reducing every coefficient mod p.
fn coef_eq_mod(a: &Polynomial, b: &Polynomial, m: &Modulus) -> bool {
    let n = a.degree().max(b.degree()) + 1;
    for i in 0..n {
        let ai = m.add(&a.get(i), &Int::ZERO);
        let bi = m.add(&b.get(i), &Int::ZERO);
        if ai != bi {
            return false;
        }
    }
    true
}

/// True iff every coefficient of `p`, reduced mod m, is zero.
fn coef_zero_mod(p: &Polynomial, m: &Modulus) -> bool {
    for i in 0..=p.degree() {
        if !m.add(&p.get(i), &Int::ZERO).is_zero() {
            return false;
        }
    }
    true
}

/// Subtract polynomials, reducing each coefficient mod p.
fn poly_sub_mod(
    a: &Polynomial,
    b: &Polynomial,
    m: &Modulus,
) -> Polynomial {
    let n = a.degree().max(b.degree()) + 1;
    let mut out = Polynomial::zeros(n);
    for i in 0..n {
        out.set(i, m.sub(&a.get(i), &b.get(i)));
    }
    out.trim()
}

/// Negate every coefficient (modulo m).
fn poly_neg(p: &Polynomial, m: &Modulus) -> Polynomial {
    let mut out = Polynomial::zeros(p.degree() + 1);
    for i in 0..=p.degree() {
        out.set(i, m.neg(&p.get(i)));
    }
    out.trim()
}

/// Increment the polynomial as if its coefficients were "digits base p".
/// Bumps coef[0] by 1; on overflow (== 0 mod p) carries to coef[1] etc.
/// Mirrors `FF_bump` from `poly_eliptic.c`.
fn ff_bump(
    x: &Polynomial,
    irrd: &Polynomial,
    m: &Modulus,
) -> Polynomial {
    let n = irrd.degree();
    let one = Int::from(1);
    let mut out = Polynomial::zeros(n);
    for i in 0..n {
        out.set(i, m.add(&x.get(i), &Int::ZERO));
    }
    let mut i = 0;
    while i < n {
        let v = m.add(&out.get(i), &one);
        let wrapped = v.is_zero();
        out.set(i, v);
        if !wrapped {
            break;
        }
        i += 1;
    }
    out.trim()
}

/// Lexicographic compare of leading-to-trailing coefficients.
fn cmp_leading(a: &Polynomial, b: &Polynomial) -> i32 {
    let da = a.degree();
    let db = b.degree();
    if da != db {
        return if da < db { -1 } else { 1 };
    }
    let mut i = da as isize;
    while i >= 0 {
        let av = a.get(i as usize);
        let bv = b.get(i as usize);
        if av < bv {
            return -1;
        }
        if av > bv {
            return 1;
        }
        i -= 1;
    }
    0
}

/// Euler criterion in GF(p^k): `a` is a QR iff `a^((p^k - 1)/2) == 1`.
fn is_quad_residue(
    a: &Polynomial,
    irrd: &Polynomial,
    m: &Modulus,
) -> bool {
    if coef_zero_mod(a, m) {
        return true;
    }
    let pk = m.n.clone().pow(irrd.degree() as u32);
    let exp = (pk - 1) / 2;
    let r = a.pow(&exp, irrd, m);
    r.degree() == 0 && r.get(0) == 1
}

/// Brute-force sqrt for small GF(p^k): enumerate every element of the
/// extension field and return one whose square is `a`. Last-resort fallback
/// for the tiny `F_43²` demo curve, where the Tonelli-Shanks lifted to
/// polynomials occasionally hits a state it can't reduce out of (`M = 1,
/// t = -1`) — that state requires multiplying R by sqrt(-1), which the
/// algorithm doesn't have direct access to. Brute force is `O(p^k)` so this
/// is only viable for very small fields; on production-sized curves the
/// regular T-S path covers all inputs.
fn poly_sqrt_bruteforce(
    a: &Polynomial,
    irrd: &Polynomial,
    pm: &Modulus,
) -> Option<Polynomial> {
    let k = irrd.degree();
    let p = pm.n.clone();
    let total: Int = p.clone().pow(k as u32);
    let mut idx = Int::from(0);
    while idx < total {
        let mut cand = Polynomial::zeros(k);
        let mut rem = idx.clone();
        for i in 0..k {
            let d = rem.clone() % &p;
            cand.set(i, d.clone());
            rem /= &p;
        }
        let cand = cand.trim();
        let sq = cand.mul(&cand, irrd, pm);
        if coef_eq_mod(&sq, a, pm) {
            return Some(cand);
        }
        idx += 1;
    }
    None
}

/// Square root in GF(p^k). Fast path when `p^k ≡ 3 (mod 4)`; otherwise
/// Tonelli–Shanks lifted to polynomials. Mirrors the structure of
/// [`crate::modulus::Modulus::sqrt`] — same Wikipedia notation (M, c, t, R)
/// to keep the algorithm legible. Returns None on non-residues.
///
/// When the lifted Tonelli-Shanks hits a state it can't reduce out of
/// (which happens for some QRs in very small extension fields), this falls
/// back to [`poly_sqrt_bruteforce`] rather than returning None — the demo
/// binaries (`weil`, `tate`) need every sqrt to enumerate every point.
fn poly_sqrt(
    a: &Polynomial,
    irrd: &Polynomial,
    pm: &Modulus,
) -> Option<Polynomial> {
    if let Some(s) = poly_sqrt_inner(a, irrd, pm) {
        return Some(s);
    }
    // Tonelli-Shanks bailed but `a` may still be a QR (small-field edge
    // case). Try brute force only when the field is small enough to
    // enumerate cheaply (≤ 2^20 elements).
    let pk_bits =
        pm.n.significant_bits() * irrd.degree() as u32;
    if pk_bits <= 20 && is_quad_residue(a, irrd, pm) {
        return poly_sqrt_bruteforce(a, irrd, pm);
    }
    None
}

fn poly_sqrt_inner(
    a: &Polynomial,
    irrd: &Polynomial,
    pm: &Modulus,
) -> Option<Polynomial> {
    if coef_zero_mod(a, pm) {
        return Some(Polynomial::zeros(0));
    }
    if !is_quad_residue(a, irrd, pm) {
        return None;
    }
    let pk: Int = pm.n.clone().pow(irrd.degree() as u32);

    // Fast path: pk ≡ 3 mod 4 → sqrt(a) = a^((pk + 1)/4).
    if pk.get_bit(0) && pk.get_bit(1) {
        let exp = (pk + 1) / 4;
        return Some(a.pow(&exp, irrd, pm));
    }

    // Decompose pk - 1 = q · 2^s with q odd.
    let mut q: Int = pk.clone() - 1;
    let mut s: u32 = 0;
    while !q.get_bit(0) {
        q /= 2;
        s += 1;
    }

    // Find a non-residue z; c = z^q has order exactly 2^s.
    let z = loop {
        let cand = rand_poly(irrd.degree(), pm);
        if coef_zero_mod(&cand, pm) {
            continue;
        }
        if !is_quad_residue(&cand, irrd, pm) {
            break cand;
        }
    };

    let mut m_state: u32 = s;
    let mut c = z.pow(&q, irrd, pm);
    let mut t = a.pow(&q, irrd, pm);
    let r_exp = (q + 1) / 2;
    let mut r = a.pow(&r_exp, irrd, pm);

    loop {
        if poly_is_one(&t) {
            return Some(r);
        }
        // Smallest i ∈ [1, M) such that t^(2^i) = 1.
        let mut i: u32 = 1;
        let mut tmp = t.mul(&t, irrd, pm);
        while !poly_is_one(&tmp) {
            i += 1;
            if i >= m_state {
                // True QRs guarantee i < M; bail rather than loop.
                return None;
            }
            tmp = tmp.mul(&tmp, irrd, pm);
        }
        // After the inner loop exits we need `i < M` for the shift below
        // to be non-negative. The while-guard above only fires when we
        // *enter* the loop body — if t² == 1 immediately and M == 1 we
        // get here with i == M, so re-check.
        if i + 1 > m_state {
            return None;
        }
        // b = c^(2^(M - i - 1))
        let shift = m_state - i - 1;
        let exp = Int::from(1) << shift;
        let b = c.pow(&exp, irrd, pm);
        let b2 = b.mul(&b, irrd, pm);
        m_state = i;
        c = b2.clone();
        t = t.mul(&b2, irrd, pm);
        r = r.mul(&b, irrd, pm);
    }
}

#[cfg(test)]
mod sqrt_test {
    use super::*;

    fn lift(v: i64) -> Polynomial {
        let mut p = Polynomial::zeros(1);
        p.set(0, Int::from(v));
        p
    }

    #[test]
    fn test_poly_sqrt_constant_one() {
        let m = Modulus::new(&Int::from(43));
        let irrd = Polynomial::find_irreducible(2, &m).unwrap();
        let one = lift(1);
        let s = poly_sqrt(&one, &irrd, &m).expect("sqrt(1)");
        let s2 = s.mul(&s, &irrd, &m);
        assert!(poly_is_one(&s2));
    }

    #[test]
    fn test_poly_sqrt_random_qrs() {
        // For every random non-zero a, a² is a QR and sqrt(a²) must square
        // back to a² (up to sign). Stress-test Tonelli-Shanks branches.
        let m = Modulus::new(&Int::from(43));
        let irrd = Polynomial::find_irreducible(2, &m).unwrap();
        for seed in 1..40 {
            let mut a = Polynomial::zeros(2);
            a.set(0, Int::from(seed));
            a.set(1, Int::from(seed % 7));
            let a_sq = a.mul(&a, &irrd, &m);
            let s = poly_sqrt(&a_sq, &irrd, &m)
                .expect("a² is always a QR");
            let s2 = s.mul(&s, &irrd, &m);
            // s² must equal a² mod p (coefficient-wise after reduction).
            assert!(coef_eq_mod(&s2, &a_sq, &m), "seed {seed}");
        }
    }
}

fn rand_poly(deg: usize, m: &Modulus) -> Polynomial {
    let mut p = Polynomial::zeros(deg);
    for i in 0..deg {
        p.set(i, m.rand());
    }
    p.trim()
}

fn poly_is_one(p: &Polynomial) -> bool {
    p.degree() == 0 && p.get(0) == 1
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Standard tiny curve from the book: y^2 = x^3 + 23x + 42 over F_43,
    /// extended by an irreducible of degree 2.
    fn tiny_curve() -> PolyCurve {
        let m = Modulus::new(&Int::from(43));
        let irrd =
            Polynomial::find_irreducible(2, &m).expect("irreducible");

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

    #[test]
    fn test_polypoint_inf() {
        let inf = PolyPoint::inf();
        let p = PolyPoint::new(lift(0), lift(0));
        assert!(inf.is_inf());
        assert!(!p.is_inf());
        assert_ne!(inf, p);
    }

    #[test]
    fn test_g1g2_classification() {
        // Constant coordinates → G1.
        let pp = PolyPoint::new(lift(5), lift(7));
        assert_eq!(pp.g1g2(), 1);
        // x with degree > 0 → G2 marker (2).
        let mut x = Polynomial::zeros(2);
        x.set(0, Int::from(1));
        x.set(1, Int::from(1));
        let pp = PolyPoint::new(x.clone(), lift(7));
        assert_eq!(pp.g1g2(), 2);
        // Both degrees > 0 → 4.
        let pp = PolyPoint::new(x.clone(), x);
        assert_eq!(pp.g1g2(), 4);
    }

    #[test]
    fn test_tiny_curve_embed_then_fits() {
        let ec = tiny_curve();
        // Sweep starting from x = 0 and find a point.
        let start = lift(0);
        let (p1, p2) = ec.embed(&start, 100).expect("embed");
        assert!(ec.fits(&p1));
        assert!(ec.fits(&p2));
        // The two points must share x and have opposite y.
        assert_eq!(p1.x, p2.x);
        let y_sum = p1.y.add(&p2.y);
        assert!(coef_zero_mod(&y_sum, &ec.p));
    }

    #[test]
    fn test_tiny_curve_add_inverse_is_inf() {
        let ec = tiny_curve();
        let (p1, p2) = ec.embed(&lift(0), 100).expect("embed");
        assert!(ec.add(&p1, &p2).is_inf());
    }

    #[test]
    fn test_tiny_curve_inf_is_identity() {
        let ec = tiny_curve();
        let (p, _) = ec.embed(&lift(0), 100).expect("embed");
        let inf = PolyPoint::inf();
        assert_eq!(ec.add(&p, &inf), p);
        assert_eq!(ec.add(&inf, &p), p);
    }

    #[test]
    fn test_tiny_curve_double_is_add_self() {
        let ec = tiny_curve();
        let (p, _) = ec.embed(&lift(0), 100).expect("embed");
        let dbl = ec.add(&p, &p);
        let via_mul = ec.mul(&p, &Int::from(2));
        assert_eq!(dbl, via_mul);
        assert!(ec.fits(&dbl));
    }

    #[test]
    fn test_tiny_curve_scalar_associates() {
        // (a+b)P = aP + bP
        let ec = tiny_curve();
        let (p, _) = ec.embed(&lift(0), 100).expect("embed");
        let a = Int::from(7);
        let b = Int::from(13);
        let lhs = ec.mul(&p, &(a.clone() + b.clone()));
        let rhs = ec.add(&ec.mul(&p, &a), &ec.mul(&p, &b));
        assert_eq!(lhs, rhs);
    }

    #[test]
    fn test_tiny_curve_order_of_point() {
        // The tiny curve has 1815 = 3·5·11² points on E(F_{43²}); valid
        // factor orders include 3, 5, 11, 15, 33, 55, 165, 1815.
        let ec = tiny_curve();
        let factors: Vec<Int> = [3, 5, 11, 15, 33, 55, 165, 1815]
            .into_iter()
            .map(Int::from)
            .collect();
        let (p, _) = ec.embed(&lift(0), 100).expect("embed");
        let ord = ec.order(&p, &factors).expect("order");
        // Sanity: [ord] P = O.
        assert!(ec.mul(&p, &ord).is_inf());
        // Order should be in the factor list.
        assert!(factors.contains(&ord));
    }
}
