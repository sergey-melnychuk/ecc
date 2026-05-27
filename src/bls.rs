//! BLS-style pairing signatures (Chapter 18).
//!
//! Implements `keygen`, `sign`, `verify` and an aggregation primitive on a
//! configurable pairing-friendly setup. The signature verification equation
//! is the standard `e(σ, G2) == e(H(m), PK)`, with `e` the Weil pairing.
//!
//! Port focuses on the core algorithms from `signature.c`; the more elaborate
//! multi-signature schemes (membership keys, subgroup signing) layer on top of
//! the same primitives.

use rug::Integer as Int;

use crate::elliptic::{Curve, Point};
use crate::hash::{hash, hash_n};
use crate::pairing::weil;
use crate::poly_elliptic::{PolyCurve, PolyPoint};
use crate::polynomial::Polynomial;

/// Pairing-friendly setup parameters needed for BLS.
#[derive(Clone)]
pub struct BlsSystem {
    /// Base curve E over F_p.
    pub e: Curve,
    /// Extension curve E over F_{p^k}.
    pub ex: PolyCurve,
    /// Torsion order (prime r such that r | #E(F_p)).
    pub tor: Int,
    /// Cofactor for the base curve: cobse = #E(F_p) / tor.
    pub cobse: Int,
    /// G1 generator (point of order `tor` on E).
    pub g1: Point,
    /// G2 generator (point of order `tor` on E_x, NOT just a lift of g1).
    pub g2: PolyPoint,
    /// Reference point of order ≠ tor used in Weil pairing evaluation.
    pub aux: PolyPoint,
}

/// Lift a base-curve point into the extension curve as a degree-0 PolyPoint.
fn to_g2(p: &Point) -> PolyPoint {
    let mut x = Polynomial::zeros(1);
    x.set(0, p.x.clone());
    let mut y = Polynomial::zeros(1);
    y.set(0, p.y.clone());
    PolyPoint::new(x, y)
}

impl BlsSystem {
    /// Generate a key pair. `sk` is a uniform scalar in `[1, tor)`. `pk = sk*G2`.
    pub fn keygen(&self) -> (Int, PolyPoint) {
        let m = crate::modulus::Modulus::new(&self.tor);
        let mut sk = m.rand();
        if sk.is_zero() {
            sk = Int::from(1);
        }
        let pk = self.ex.mul(&self.g2, &sk);
        (sk, pk)
    }

    /// Hash a message to a point on G1.
    ///
    /// Hashes to an integer, embeds onto E by sweeping x until f(x) is a QR,
    /// then multiplies by the cofactor so the result lands in the order-`tor`
    /// subgroup. Returns None if no valid embedding is found within `limit`
    /// candidate x values — should not happen for cryptographic-size fields.
    pub fn hash_to_g1(&self, msg: &[u8]) -> Option<Point> {
        let h = hash(msg);
        let modulus = crate::modulus::Modulus::new(&self.e.modulus);
        let x = modulus.add(&h, &Int::ZERO);
        let p = self.e.find(&x, 1024)?;
        Some(self.e.mul(&p, &self.cobse))
    }

    /// BLS sign: σ = sk · H(m).
    pub fn sign(&self, sk: &Int, msg: &[u8]) -> Option<Point> {
        let h = self.hash_to_g1(msg)?;
        Some(self.e.mul(&h, sk))
    }

    /// BLS verify: returns true iff `e(σ, G2) == e(H(m), PK)`.
    pub fn verify(&self, pk: &PolyPoint, msg: &[u8], sig: &Point) -> bool {
        let Some(h) = self.hash_to_g1(msg) else { return false; };
        let h2 = to_g2(&h);
        let sig2 = to_g2(sig);
        let w_sig = weil(&self.ex, &sig2, &self.g2, &self.aux, &self.tor);
        let w_pk = weil(&self.ex, &h2, pk, &self.aux, &self.tor);
        w_sig == w_pk
    }

    /// Aggregate two signatures by summing their G1 points: σ_agg = σ1 + σ2.
    /// Then `e(σ_agg, G2) == e(H(m), PK1+PK2)` follows from bilinearity.
    pub fn aggregate_sigs(&self, sigs: &[Point]) -> Point {
        let mut acc = Point::inf();
        for s in sigs {
            acc = self.e.add(&acc, s);
        }
        acc
    }

    /// Aggregate two public keys by summing their G2 points.
    pub fn aggregate_pks(&self, pks: &[PolyPoint]) -> PolyPoint {
        let mut acc = PolyPoint::inf();
        for k in pks {
            acc = self.ex.add(&acc, k);
        }
        acc
    }

    /// Deterministic key derivation from a passphrase (utility for demos).
    /// Hashes the bytes to a torsion scalar.
    pub fn derive_sk(&self, seed: &[u8], domain: &[u8]) -> Int {
        let modulus = crate::modulus::Modulus::new(&self.tor);
        let h = hash_n(&[domain, seed]);
        modulus.add(&h, &Int::ZERO)
    }
}

/// Build the tiny demo system used in chapters 16–18: p=43, k=2, tor=11.
pub fn tiny_system() -> BlsSystem {
    let p = crate::modulus::Modulus::new(&Int::from(43));
    let irrd =
        Polynomial::find_irreducible(2, &p).expect("irreducible deg-2");

    // Base curve y² = x³ + 23x + 42 mod 43.
    let card_e = Int::from(55); // p + 1 - t with t = -11
    let tor = Int::from(11);
    let cobse = card_e.clone() / &tor;
    // Find a generator G1 on the base curve. The curve cardinality is 55, so
    // any order-11 point works. Sweep x's looking for a point whose order is
    // exactly 11.
    let base_curve = {
        // We don't have a known base point yet; fill in a placeholder first.
        let a = Int::from(23);
        let b = Int::from(42);
        let placeholder = Point::new(Int::from(0), Int::from(0));
        Curve::new(Int::from(43), card_e.clone(), placeholder, a, b)
    };

    let mut g1 = Point::inf();
    for cand_x in 0..43 {
        let Some(p) = base_curve.find(&Int::from(cand_x), 1) else {
            continue;
        };
        // Order check: 11·P = O and 5·P ≠ O (rules out order-1 and order-5).
        let p11 = base_curve.mul(&p, &Int::from(11));
        let p5 = base_curve.mul(&p, &Int::from(5));
        if p11.is_inf() && !p5.is_inf() {
            g1 = p;
            break;
        }
    }
    assert!(!g1.is_inf(), "no order-11 G1 generator found");

    let base_curve = Curve::new(
        Int::from(43),
        card_e,
        g1.clone(),
        Int::from(23),
        Int::from(42),
    );

    // Extension curve
    let mut a4 = Polynomial::zeros(1);
    a4.set(0, Int::from(23));
    let mut a6 = Polynomial::zeros(1);
    a6.set(0, Int::from(42));
    let ext_curve = PolyCurve::new(a4, a6, irrd, p);

    // Find a G2 generator: an order-11 point with degree>0 in at least one
    // coordinate (so it lives in G2, not just lifted G1).
    let g2 = find_g2(&ext_curve, &tor);
    let aux = find_aux(&ext_curve, &tor);

    BlsSystem {
        e: base_curve,
        ex: ext_curve,
        tor,
        cobse,
        g1,
        g2,
        aux,
    }
}

fn find_g2(ec: &PolyCurve, tor: &Int) -> PolyPoint {
    let factors: Vec<Int> = [3, 5, 11, 15, 33, 55, 165, 1815]
        .into_iter()
        .map(Int::from)
        .collect();
    let mut x = Polynomial::zeros(1);
    for _ in 0..200 {
        if let Some((p, _)) = ec.embed(&x, 1) {
            if p.g1g2() > 1 {
                if let Some(ord) = ec.order(&p, &factors) {
                    if &ord == tor {
                        return p;
                    }
                }
            }
            x = bump(&p.x, ec);
        } else {
            x = bump(&x, ec);
        }
    }
    panic!("no G2 generator found");
}

fn find_aux(ec: &PolyCurve, tor: &Int) -> PolyPoint {
    let factors: Vec<Int> = [3, 5, 11, 15, 33, 55, 165, 1815]
        .into_iter()
        .map(Int::from)
        .collect();
    let mut x = Polynomial::zeros(1);
    for _ in 0..200 {
        if let Some((p, _)) = ec.embed(&x, 1) {
            if let Some(ord) = ec.order(&p, &factors) {
                if &ord != tor {
                    return p;
                }
            }
            x = bump(&p.x, ec);
        } else {
            x = bump(&x, ec);
        }
    }
    panic!("no aux point found");
}

fn bump(x: &Polynomial, ec: &PolyCurve) -> Polynomial {
    let m = &ec.p;
    let n = ec.irrd.degree();
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tiny_system_invariants() {
        let sys = tiny_system();
        // G1 is order-11 on the base curve.
        let g1_11 = sys.e.mul(&sys.g1, &sys.tor);
        assert!(g1_11.is_inf());
        // G2 is order-11 on the extension curve.
        let g2_11 = sys.ex.mul(&sys.g2, &sys.tor);
        assert!(g2_11.is_inf());
        // G2 actually lives in the extension (degree > 0 somewhere).
        assert!(sys.g2.g1g2() > 1);
        // Aux has order ≠ tor.
        let aux_11 = sys.ex.mul(&sys.aux, &sys.tor);
        assert!(!aux_11.is_inf());
    }

    #[test]
    fn test_keygen_pk_in_g2_subgroup() {
        let sys = tiny_system();
        let (_sk, pk) = sys.keygen();
        // [tor] · PK = O (PK is in the order-tor G2 subgroup).
        let pk_tor = sys.ex.mul(&pk, &sys.tor);
        assert!(pk_tor.is_inf());
    }

    #[test]
    fn test_hash_to_g1_lands_in_torsion() {
        let sys = tiny_system();
        let h = sys.hash_to_g1(b"hello world").expect("hash");
        // [tor] · H = O.
        let hh = sys.e.mul(&h, &sys.tor);
        assert!(hh.is_inf());
    }

    #[test]
    fn test_bls_sign_and_verify_roundtrip() {
        let sys = tiny_system();
        let (sk, pk) = sys.keygen();
        let msg = b"attack at dawn";
        let sig = sys.sign(&sk, msg).expect("sign");
        assert!(sys.verify(&pk, msg, &sig));
    }

    #[test]
    fn test_bls_rejects_tampered_message() {
        let sys = tiny_system();
        let (sk, pk) = sys.keygen();
        let sig = sys.sign(&sk, b"original").expect("sign");
        assert!(!sys.verify(&pk, b"tampered", &sig));
    }

    #[test]
    fn test_bls_rejects_wrong_pk() {
        let sys = tiny_system();
        let (sk1, _pk1) = sys.keygen();
        let (_sk2, pk2) = sys.keygen();
        let sig = sys.sign(&sk1, b"hi").expect("sign");
        assert!(!sys.verify(&pk2, b"hi", &sig));
    }

    #[test]
    fn test_bls_aggregate_verifies_with_summed_pk() {
        // Two signers on the SAME message: e(σ1+σ2, G2) == e(H(m), PK1+PK2).
        let sys = tiny_system();
        let (sk1, pk1) = sys.keygen();
        let (sk2, pk2) = sys.keygen();
        let msg = b"same message";
        let s1 = sys.sign(&sk1, msg).expect("sign1");
        let s2 = sys.sign(&sk2, msg).expect("sign2");
        let agg_sig = sys.aggregate_sigs(&[s1, s2]);
        let agg_pk = sys.aggregate_pks(&[pk1, pk2]);
        assert!(sys.verify(&agg_pk, msg, &agg_sig));
    }
}
