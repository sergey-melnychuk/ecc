//! SNARK Prover - generates proofs for statements.
//!
//! Given a QAP, CRS, and witness values, the prover computes
//! proof elements (A, B, C) that satisfy the verification equation.
//!
//! Ported from snark_proof.c

use rug::Integer as Int;
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Read, Write};

use crate::elliptic::Point;

use super::crs::Crs;
use super::field_ext::PolyPoint;
use super::lagrange::LagrangeInterpolator;
use super::qap::Qap;
use super::system::{
    read_ext_point, read_int, read_point, write_ext_point, write_int,
    write_point, SigSystem,
};

/// SNARK proof elements.
pub struct Proof {
    /// Proof element A (G1 point)
    pub a: Point,
    /// Proof element B (G2 point in extension field)
    pub b: PolyPoint,
    /// Proof element C (G1 point)
    pub c: Point,
}

impl Proof {
    /// Create a new empty proof
    pub fn new() -> Self {
        Self {
            a: Point::inf(),
            b: PolyPoint::new(),
            c: Point::inf(),
        }
    }

    /// Load proof from binary file.
    pub fn load(filename: &str) -> io::Result<ProofRecord> {
        let file = File::open(filename)?;
        let mut reader = BufReader::new(file);

        // Read dimensions
        let mut buf = [0u8; 4];
        reader.read_exact(&mut buf)?;
        let n = i32::from_ne_bytes(buf) as usize;
        reader.read_exact(&mut buf)?;
        let m = i32::from_ne_bytes(buf) as usize;
        reader.read_exact(&mut buf)?;
        let l = i32::from_ne_bytes(buf) as usize;

        // Read public inputs (statement values)
        let mut statement = vec![Int::ZERO; l + 1];
        for i in 0..=l {
            statement[i] = read_int(&mut reader)?;
        }

        // Read proof elements
        let a = read_point(&mut reader)?;
        let c = read_point(&mut reader)?;
        let b = read_ext_point(&mut reader)?;

        Ok(ProofRecord {
            n,
            m,
            l,
            statement,
            proof: Proof { a, b, c },
        })
    }

    /// Save proof to binary file.
    pub fn save(
        &self,
        filename: &str,
        n: usize,
        m: usize,
        l: usize,
        statement: &[Int],
    ) -> io::Result<()> {
        let file = File::create(filename)?;
        let mut writer = BufWriter::new(file);

        // Write dimensions
        writer.write_all(&(n as i32).to_ne_bytes())?;
        writer.write_all(&(m as i32).to_ne_bytes())?;
        writer.write_all(&(l as i32).to_ne_bytes())?;

        // Write public inputs
        for s in statement {
            write_int(&mut writer, s)?;
        }

        // Write proof elements
        write_point(&mut writer, &self.a)?;
        write_point(&mut writer, &self.c)?;
        write_ext_point(&mut writer, &self.b)?;

        Ok(())
    }
}

impl Default for Proof {
    fn default() -> Self {
        Self::new()
    }
}

/// A proof record contains the proof and public statement values.
pub struct ProofRecord {
    /// Number of gates
    pub n: usize,
    /// Number of wires
    pub m: usize,
    /// Number of public inputs
    pub l: usize,
    /// Public statement values (a_0, ..., a_l ordered by sw)
    pub statement: Vec<Int>,
    /// The proof itself
    pub proof: Proof,
}

/// SNARK Prover.
pub struct Prover {
    /// System parameters
    sys: SigSystem,
    /// QAP parameters
    qap: Qap,
    /// Common Reference String
    crs: Crs,
}

impl Prover {
    /// Create a new prover with the given parameters.
    pub fn new(sys: SigSystem, qap: Qap, crs: Crs) -> Self {
        Self { sys, qap, crs }
    }

    /// Compute wire values for the example circuit.
    ///
    /// Given a1 (medicine), a2 (dose), a3 (patient), computes all wire values.
    /// Maps to C function: wires()
    fn compute_wires(
        &self,
        a1: &Int,
        a2: &Int,
        a3: &Int,
    ) -> Vec<Int> {
        let tor_mod = self.sys.torsion_modulus();
        let m = self.qap.m;
        let mut aj = vec![Int::ZERO; m];

        aj[0] = Int::from(1); // a0 = 1 (constant)
        aj[1] = a1.clone(); // a1 = medicine
        aj[2] = a2.clone(); // a2 = dose
        aj[3] = a3.clone(); // a3 = patient

        // Random a4
        aj[4] = tor_mod.rand();

        // a5 = a1
        aj[5] = a1.clone();

        // a6 = a2 * a3
        aj[6] = tor_mod.mul(a2, a3);

        // a7 = (a1 + a2) * a1
        let tmp = tor_mod.add(a1, a2);
        aj[7] = tor_mod.mul(&tmp, a1);

        // a8 = (a3 + a4) * a6
        let tmp = tor_mod.add(a3, &aj[4]);
        aj[8] = tor_mod.mul(&tmp, &aj[6]);

        // a9 = a8 * a7
        aj[9] = tor_mod.mul(&aj[8], &aj[7]);

        aj
    }

    /// Compute cross-term coefficients for h(x).
    ///
    /// Maps to C function: crossterms()
    fn crossterms(&self, a: &[Int]) -> Vec<Int> {
        let tor_mod = self.sys.torsion_modulus();
        let mut coef = vec![Int::ZERO; 10];

        // These are the specific cross-terms for the example circuit
        // l0-l1
        let t1 = tor_mod.mul(&a[1], &a[2]);
        let t2 = tor_mod.mul(&a[0], &a[3]);
        coef[0] = tor_mod.add(&t1, &t2);

        // l0-l2
        let t0 = tor_mod.add(&a[1], &a[2]);
        let t1 = tor_mod.mul(&t0, &a[0]);
        let t2 = tor_mod.mul(&a[1], &a[5]);
        coef[1] = tor_mod.add(&t1, &t2);

        // l0-l3
        let t1 = tor_mod.add(&a[3], &a[4]);
        let t2 = tor_mod.mul(&t1, &a[0]);
        let t3 = tor_mod.mul(&a[1], &a[6]);
        coef[2] = tor_mod.add(&t2, &t3);

        // l0-l4
        let t4 = tor_mod.mul(&a[0], &a[8]);
        let t5 = tor_mod.mul(&a[1], &a[7]);
        coef[3] = tor_mod.add(&t4, &t5);

        // l1-l2
        let t2 = tor_mod.mul(&t0, &a[2]);
        let t3 = tor_mod.mul(&a[3], &a[5]);
        coef[4] = tor_mod.add(&t2, &t3);

        // l1-l3
        let t1 = tor_mod.add(&a[3], &a[4]);
        let t2 = tor_mod.mul(&t1, &a[2]);
        let t3 = tor_mod.mul(&a[3], &a[6]);
        coef[5] = tor_mod.add(&t2, &t3);

        // l1-l4
        let t4 = tor_mod.mul(&a[3], &a[7]);
        let t5 = tor_mod.mul(&a[2], &a[8]);
        coef[6] = tor_mod.add(&t4, &t5);

        // l2-l3
        let t0 = tor_mod.add(&a[1], &a[2]);
        let t1 = tor_mod.add(&a[3], &a[4]);
        let t3 = tor_mod.mul(&t0, &a[6]);
        let t4 = tor_mod.mul(&t1, &a[5]);
        coef[7] = tor_mod.add(&t3, &t4);

        // l2-l4
        let t2 = tor_mod.mul(&t0, &a[7]);
        let t3 = tor_mod.mul(&a[5], &a[8]);
        coef[8] = tor_mod.add(&t2, &t3);

        // l3-l4
        let t2 = tor_mod.mul(&t1, &a[7]);
        let t3 = tor_mod.mul(&a[6], &a[8]);
        coef[9] = tor_mod.add(&t2, &t3);

        coef
    }

    /// Generate a proof for the given statement values.
    ///
    /// # Arguments
    /// * `a1` - Medicine number (public)
    /// * `a2` - Dose (public)
    /// * `a3` - Patient number (witness)
    pub fn prove(&self, a1: &Int, a2: &Int, a3: &Int) -> ProofRecord {
        let tor_mod = self.sys.torsion_modulus();
        let prime_mod = self.sys.modulus();
        let interp = LagrangeInterpolator::new(&self.sys.tor);

        let n = self.qap.n;
        let m = self.qap.m;
        let l = self.qap.l;

        // Compute all wire values
        let aj = self.compute_wires(a1, a2, a3);

        // Compute h(x) coefficients from cross-terms
        let cross_coef = self.crossterms(&aj);

        // Flatten htable and compute h(x) coefficients using matflat
        let htable_flat: Vec<Int> =
            self.qap.htable.iter().flatten().cloned().collect();
        let hx_coef =
            interp.matflat(&htable_flat, n - 1, &cross_coef);

        // Choose random r and s for zero-knowledge
        let r = tor_mod.rand();
        let s = tor_mod.rand();

        // Compute av = sum(aj * vj) and aw = sum(aj * wj)
        // These are the combined left/right input polynomials evaluated at z
        let vj_flat: Vec<Int> =
            self.qap.vj.iter().flatten().cloned().collect();
        let wj_flat: Vec<Int> =
            self.qap.wj.iter().flatten().cloned().collect();
        let av = interp.matflat(&vj_flat, n, &aj);
        let aw = interp.matflat(&wj_flat, n, &aj);

        // Compute proof element A = alpha*G + sum(av[i] * z^i * G) + r*delta*G
        let j = n - 1; // For indexing
        let mut proof_a = self.crs.alpha_g.clone();
        for i in 0..n {
            let term = self.sys.e.mul(&self.crs.z_g[j - i], &av[i]);
            proof_a = self.sys.e.add(&proof_a, &term);
        }
        let r_delta_g = self.sys.e.mul(&self.crs.delta_g, &r);
        proof_a = self.sys.e.add(&proof_a, &r_delta_g);

        // Compute proof element B = beta*H + sum(aw[i] * z^i * H) + s*delta*H
        let mut proof_b = self.crs.beta_h.clone();
        for i in 0..n {
            let term = self.sys.ex.mul(
                &self.crs.z_h[j - i],
                &aw[i],
                &self.sys.irrd,
                &prime_mod,
            );
            proof_b = self.sys.ex.add(
                &proof_b,
                &term,
                &self.sys.irrd,
                &prime_mod,
            );
        }
        let s_delta_h = self.sys.ex.mul(
            &self.crs.delta_h,
            &s,
            &self.sys.irrd,
            &prime_mod,
        );
        proof_b = self.sys.ex.add(
            &proof_b,
            &s_delta_h,
            &self.sys.irrd,
            &prime_mod,
        );

        // Compute proof element C
        // C = r*(beta*G + sum(aw[i]*z^i*G)) + s*A + sum(witness aj * theta_j * G) + sum(hx[i] * zt[i] * G)
        let mut proof_c = self.crs.beta_g.clone();
        for i in 0..n {
            let term = self.sys.e.mul(&self.crs.z_g[j - i], &aw[i]);
            proof_c = self.sys.e.add(&proof_c, &term);
        }
        proof_c = self.sys.e.mul(&proof_c, &r);

        // + s*A
        let s_a = self.sys.e.mul(&proof_a, &s);
        proof_c = self.sys.e.add(&proof_c, &s_a);

        // + witness terms (indexes l+1 to m-1)
        for i in (l + 1)..m {
            let wire_idx = self.qap.sw[i];
            let term =
                self.sys.e.mul(&self.crs.theta_g[i], &aj[wire_idx]);
            proof_c = self.sys.e.add(&proof_c, &term);
        }

        // + h(x)*t(x)/delta terms
        for i in 0..(n - 1) {
            let term = self
                .sys
                .e
                .mul(&self.crs.zt_g[j - 1 - i], &hx_coef[i]);
            proof_c = self.sys.e.add(&proof_c, &term);
        }

        // Extract public statement values (ordered by sw permutation)
        let mut statement = vec![Int::ZERO; l + 1];
        for i in 0..=l {
            statement[i] = aj[self.qap.sw[i]].clone();
        }

        ProofRecord {
            n,
            m,
            l,
            statement,
            proof: Proof {
                a: proof_a,
                b: proof_b,
                c: proof_c,
            },
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_proof_new() {
        let proof = Proof::new();
        assert!(proof.a.is_inf());
        assert!(proof.b.is_inf());
        assert!(proof.c.is_inf());
    }
}
