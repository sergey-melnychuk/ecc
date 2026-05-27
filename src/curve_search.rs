//! Pairing-friendly curve search (Chapter 14).
//!
//! Port of `pairing_sweep_alpha.c`. Implements algorithms 6.19, 6.20 and 6.24
//! from Freeman/Scott/Teske, "A Taxonomy of Pairing-Friendly Elliptic Curves"
//! (2010). For each embedding degree `k` and a small set of "alpha" values,
//! sweep `x` and collect candidates where the cyclotomic-style polynomial `r`
//! and the field-size polynomial `q` are both prime.
//!
//! Output: list of `(k, alpha, x, r, q, t)` tuples that name a CM curve
//! ready to be built with [`crate::cm_curve::get_curve`].

use rug::ops::Pow;
use rug::Integer as Int;

/// Embedding degrees handled by the sweep, mapped to which algorithm to use.
/// `true` means equation 6.2, `false` means equation 6.20 — same split as the
/// `algt` table in the C reference.
const K_TABLE: &[(u32, bool)] = &[
    (5, true),
    (7, false),
    (11, false),
    (13, true),
    (17, true),
    (19, false),
    (23, false),
    (29, true),
    (31, false),
];

/// One candidate curve from the sweep.
#[derive(Debug, Clone, PartialEq)]
pub struct Candidate {
    pub k: u32,
    pub alpha: Int,
    pub x: Int,
    /// Torsion order — prime.
    pub r: Int,
    /// Field prime.
    pub q: Int,
    /// Trace of Frobenius.
    pub t: Int,
    /// log2(q)/log2(r) — efficiency ratio (lower is better; 1.0 is ideal).
    pub rho: f64,
}

/// Cyclotomic-style polynomial `Φ_{4k}` evaluated at `z = α·x²`.
/// Mirrors `phi4k` in the C reference.
pub fn phi4k(k: u32, alpha: &Int, x: &Int) -> Int {
    let z: Int = (Int::from(x * x)) * alpha;
    // r = 1 - z + z^2 - z^3 + ... ± z^(k-1) (alternating signs).
    let mut acc = Int::from(1) - &z;
    let mut zpow = z.clone();
    for i in 1..(k as i64 - 1) {
        zpow *= &z;
        if i & 1 == 1 {
            acc += &zpow;
        } else {
            acc -= &zpow;
        }
    }
    acc
}

/// q(z) for algorithm 6.20: `(z^(k+1) + z^k + 4·z^((k+1)/2) + (z+1)) / 4`.
pub fn qofz_20(k: u32, alpha: &Int, x: &Int) -> Int {
    let z: Int = (Int::from(x * x)) * alpha;
    let k1 = k + 1;
    let k2 = k1 / 2;
    let t1: Int = z.clone().pow(k1);
    let t2: Int = z.clone().pow(k);
    let t3: Int = Int::from(z.clone().pow(k2)) * 4;
    let t4: Int = z.clone() + 1;
    let sum: Int = t1 + t2 + t3 + t4;
    sum / 4
}

/// q(z) for algorithm 6.2: `(z^(k+2) + 2·z^(k+1) + z^k + (z-1)²) / 4`.
pub fn qofz_2(k: u32, alpha: &Int, x: &Int) -> Int {
    let z: Int = (Int::from(x * x)) * alpha;
    let k1 = k + 1;
    let k2 = k + 2;
    let t1: Int = z.clone().pow(k2);
    let t2: Int = Int::from(z.clone().pow(k1)) * 2;
    let t3: Int = z.clone().pow(k);
    let z_minus_1: Int = z.clone() - 1;
    let t4: Int = z_minus_1.clone() * z_minus_1;
    let sum: Int = t1 + t2 + t3 + t4;
    sum / 4
}

/// Trace of Frobenius for algorithm 6.20: `t = (αx²)^((k+1)/2) + 1`.
pub fn tofz_20(k: u32, alpha: &Int, x: &Int) -> Int {
    let z: Int = (Int::from(x * x)) * alpha;
    let k1 = (k + 1) / 2;
    Int::from(z.pow(k1)) + 1
}

/// Trace of Frobenius for algorithm 6.2: `t = 1 - αx²`.
pub fn tofz_2(_k: u32, alpha: &Int, x: &Int) -> Int {
    let z: Int = (Int::from(x * x)) * alpha;
    Int::from(1) - z
}

/// Miller-Rabin probable-prime test wrapping rug's `is_probably_prime`.
fn is_prob_prime(n: &Int, reps: u32) -> bool {
    matches!(
        n.is_probably_prime(reps),
        rug::integer::IsPrime::Yes | rug::integer::IsPrime::Probably
    )
}

/// Sweep parameters `(α, x)` looking for pairing-friendly curves of
/// embedding degree `k`, with `r` no larger than `lg2_r_max` bits.
///
/// `k` must be one of {5, 7, 11, 13, 17, 19, 23, 29, 31}. Returns the list of
/// candidates found, sorted by rho ascending (most efficient first).
pub fn sweep(k: u32, lg2_r_max: u32) -> Vec<Candidate> {
    let algt = match K_TABLE.iter().find(|(kk, _)| *kk == k) {
        Some((_, a)) => *a,
        None => return Vec::new(),
    };

    // Same α candidates as the C reference: {7, 11, 15, 19, 23} + 20·a for a∈0..8.
    let atab = [7i64, 11, 15, 19, 23];
    let max = 1u64 << (lg2_r_max / 2 / (k - 1));
    let mut out = Vec::new();

    for a in 0..8 {
        for &base in &atab {
            let alphabase = base + 20 * a;
            let alpha = Int::from(alphabase);
            let mut j = 1u64;
            while j < max {
                let x = Int::from(j);
                let r = phi4k(k, &alpha, &x);
                let rsz = r.significant_bits();
                if rsz > lg2_r_max {
                    break;
                }
                if is_prob_prime(&r, 25) {
                    let (q, t) = if algt {
                        (qofz_2(k, &alpha, &x), tofz_2(k, &alpha, &x))
                    } else {
                        (qofz_20(k, &alpha, &x), tofz_20(k, &alpha, &x))
                    };
                    let qsz = q.significant_bits();
                    if qsz > 0 && is_prob_prime(&q, 25) {
                        let rho = qsz as f64 / rsz as f64;
                        out.push(Candidate {
                            k,
                            alpha: alpha.clone(),
                            x,
                            r,
                            q,
                            t,
                            rho,
                        });
                    }
                }
                j += 2;
            }
        }
    }

    out.sort_by(|a, b| a.rho.partial_cmp(&b.rho).unwrap());
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_phi4k_small() {
        // Φ_{4·5}(α·x²) at α=1, x=1, k=5: z = 1, polynomial in z = (1)(1-z+z²-z³+z⁴).
        // = 1 - 1 + 1 - 1 + 1 = 1.
        let r = phi4k(5, &Int::from(1), &Int::from(1));
        assert_eq!(r, Int::from(1));
    }

    #[test]
    fn test_q_t_relation_k5_alg20() {
        // For algorithm 6.20: q + 1 - t should be divisible by r (group has
        // a subgroup of order r). Smoke check on small inputs.
        let k = 5u32;
        let alpha = Int::from(7);
        let x = Int::from(3);
        let r = phi4k(k, &alpha, &x);
        let q = qofz_20(k, &alpha, &x);
        let t = tofz_20(k, &alpha, &x);
        let card: Int = q + 1 - t;
        assert_eq!(Int::from(&card % &r), Int::ZERO);
    }

    #[test]
    fn test_q_t_relation_k5_alg2() {
        // Same relation for algorithm 6.2.
        let k = 5u32;
        let alpha = Int::from(7);
        let x = Int::from(3);
        let r = phi4k(k, &alpha, &x);
        let q = qofz_2(k, &alpha, &x);
        let t = tofz_2(k, &alpha, &x);
        let card: Int = q + 1 - t;
        assert_eq!(Int::from(&card % &r), Int::ZERO);
    }

    #[test]
    fn test_sweep_finds_candidates() {
        // With a small budget there should still be at least one (k, α, x)
        // where both r and q are prime. We check the basic invariants on the
        // first candidate found.
        let candidates = sweep(5, 40);
        assert!(!candidates.is_empty(), "sweep produced no candidates");
        let c = &candidates[0];
        assert!(is_prob_prime(&c.r, 25));
        assert!(is_prob_prime(&c.q, 25));
        let card: Int = c.q.clone() + 1 - &c.t;
        assert_eq!(Int::from(&card % &c.r), Int::ZERO);
        assert!(c.rho >= 1.0);
    }

    #[test]
    fn test_sweep_k7() {
        let candidates = sweep(7, 50);
        assert!(!candidates.is_empty(), "no k=7 candidates");
        let c = &candidates[0];
        let card: Int = c.q.clone() + 1 - &c.t;
        assert_eq!(Int::from(&card % &c.r), Int::ZERO);
    }
}
