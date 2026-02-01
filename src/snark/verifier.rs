//! SNARK Verifier - verifies proofs using pairing equations.
//!
//! Verification checks that:
//! e(A, B) = e(alpha*G, beta*H) * e(V, gamma*H) * e(C, delta*H)
//!
//! Where V = sum(statement_j * theta_j * G) for public inputs.
//!
//! Ported from snark_verify.c

use crate::elliptic::Point;

use super::crs::Crs;
use super::field_ext::{Poly, PolyPoint};
use super::pairing::Pairing;
use super::prover::ProofRecord;
use super::system::SigSystem;

/// SNARK Verifier.
pub struct Verifier {
    /// System parameters
    sys: SigSystem,
    /// Common Reference String
    crs: Crs,
    /// Pairing context
    pairing: Pairing,
}

impl Verifier {
    /// Create a new verifier with the given parameters.
    pub fn new(sys: SigSystem, crs: Crs) -> Self {
        let pairing = Pairing::new(&sys.prime, sys.irrd.clone());
        Self { sys, crs, pairing }
    }

    /// Verify a proof record.
    ///
    /// Returns true if the proof is valid, false otherwise.
    pub fn verify(&self, record: &ProofRecord) -> bool {
        let prime_mod = self.sys.modulus();

        // Compute V = sum(statement_j * theta_j * G) for public inputs
        let mut v = Point::inf();
        for i in 0..=record.l {
            let term = self
                .sys
                .e
                .mul(&self.crs.theta_g[i], &record.statement[i]);
            v = self.sys.e.add(&v, &term);
        }

        // Generate random point S for pairing (needed for Tate pairing)
        let s = self.random_g2_point();

        // Convert G1 points to G2 for pairing
        let alpha_g_g2 = self.sys.to_g2(&self.crs.alpha_g);
        let v_g2 = self.sys.to_g2(&v);
        let c_g2 = self.sys.to_g2(&record.proof.c);
        let a_g2 = self.sys.to_g2(&record.proof.a);

        // Compute e(alpha*G, beta*H)
        let e_ab = self.pairing.tate(
            &alpha_g_g2,
            &self.crs.beta_h,
            &s,
            &self.sys.tor,
            &self.sys.ex,
        );

        // Compute e(V, gamma*H)
        let e_vh = self.pairing.tate(
            &v_g2,
            &self.crs.gamma_h,
            &s,
            &self.sys.tor,
            &self.sys.ex,
        );

        // Compute e(C, delta*H)
        let e_ch = self.pairing.tate(
            &c_g2,
            &self.crs.delta_h,
            &s,
            &self.sys.tor,
            &self.sys.ex,
        );

        // Compute e(A, B)
        let e_a_b = self.pairing.tate(
            &a_g2,
            &record.proof.b,
            &s,
            &self.sys.tor,
            &self.sys.ex,
        );

        // LHS = e(alpha*G, beta*H) * e(V, gamma*H) * e(C, delta*H)
        let lhs = e_ab.mul(&e_vh, &self.sys.irrd, &prime_mod);
        let lhs = lhs.mul(&e_ch, &self.sys.irrd, &prime_mod);

        // RHS = e(A, B)
        // Verification passes if LHS == RHS
        lhs.eq(&e_a_b)
    }

    /// Generate a random G2 point for pairing computation.
    ///
    /// The point S is used to ensure the pairing is non-degenerate.
    fn random_g2_point(&self) -> PolyPoint {
        let prime_mod = self.sys.modulus();

        // For a valid pairing, we need a point that's linearly independent
        // from P and Q. We can use a random scalar multiple of G2.
        let random_scalar = prime_mod.rand();
        self.sys.ex.mul(
            &self.sys.g2,
            &random_scalar,
            &self.sys.irrd,
            &prime_mod,
        )
    }
}

/// Simplified verification result with diagnostic info.
pub struct VerificationResult {
    /// Whether the proof verified
    pub valid: bool,
    /// e(alpha*G, beta*H)
    pub e_alpha_beta: Poly,
    /// e(V, gamma*H)
    pub e_v_gamma: Poly,
    /// e(C, delta*H)
    pub e_c_delta: Poly,
    /// e(A, B)
    pub e_a_b: Poly,
    /// LHS product
    pub lhs: Poly,
}

impl Verifier {
    /// Verify with detailed diagnostic output.
    ///
    /// Returns both the verification result and intermediate values for debugging.
    pub fn verify_verbose(
        &self,
        record: &ProofRecord,
    ) -> VerificationResult {
        let prime_mod = self.sys.modulus();

        // Compute V = sum(statement_j * theta_j * G)
        let mut v = Point::inf();
        for i in 0..=record.l {
            let term = self
                .sys
                .e
                .mul(&self.crs.theta_g[i], &record.statement[i]);
            v = self.sys.e.add(&v, &term);
        }

        let s = self.random_g2_point();

        let alpha_g_g2 = self.sys.to_g2(&self.crs.alpha_g);
        let v_g2 = self.sys.to_g2(&v);
        let c_g2 = self.sys.to_g2(&record.proof.c);
        let a_g2 = self.sys.to_g2(&record.proof.a);

        let e_alpha_beta = self.pairing.tate(
            &alpha_g_g2,
            &self.crs.beta_h,
            &s,
            &self.sys.tor,
            &self.sys.ex,
        );

        let e_v_gamma = self.pairing.tate(
            &v_g2,
            &self.crs.gamma_h,
            &s,
            &self.sys.tor,
            &self.sys.ex,
        );

        let e_c_delta = self.pairing.tate(
            &c_g2,
            &self.crs.delta_h,
            &s,
            &self.sys.tor,
            &self.sys.ex,
        );

        let e_a_b = self.pairing.tate(
            &a_g2,
            &record.proof.b,
            &s,
            &self.sys.tor,
            &self.sys.ex,
        );

        let lhs =
            e_alpha_beta.mul(&e_v_gamma, &self.sys.irrd, &prime_mod);
        let lhs = lhs.mul(&e_c_delta, &self.sys.irrd, &prime_mod);

        let valid = lhs.eq(&e_a_b);

        VerificationResult {
            valid,
            e_alpha_beta,
            e_v_gamma,
            e_c_delta,
            e_a_b,
            lhs,
        }
    }
}

#[cfg(test)]
mod tests {
    // Integration tests would go here, requiring a full setup
    // with real curve parameters
}
