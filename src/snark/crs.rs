//! Common Reference String (CRS) for SNARK setup.
//!
//! The CRS is generated from "toxic waste" random values (alpha, beta, gamma, delta, z)
//! that must be destroyed after setup to ensure soundness.
//!
//! Ported from snark_crs.c

use rug::Integer as Int;
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Read, Write};

use crate::elliptic::Point;

use super::field_ext::PolyPoint;
use super::lagrange::LagrangeInterpolator;
use super::qap::Qap;
use super::system::{
    read_ext_point, read_point, write_ext_point, write_point,
    SigSystem,
};

/// Common Reference String for SNARK proving and verification.
#[derive(Clone)]
pub struct Crs {
    /// Number of gates
    pub n: usize,
    /// Number of wires
    pub m: usize,
    /// Number of public inputs
    pub l: usize,

    // G1 elements (base curve points)
    /// alpha * G
    pub alpha_g: Point,
    /// beta * G
    pub beta_g: Point,
    /// delta * G
    pub delta_g: Point,
    /// z^i * G for i = 0..n
    pub z_g: Vec<Point>,
    /// theta_j * G for each wire j (ordered by sw permutation)
    pub theta_g: Vec<Point>,
    /// z^i * t(z) / delta * G for i = 0..n-1
    pub zt_g: Vec<Point>,

    // G2 elements (extension curve points)
    /// beta * H
    pub beta_h: PolyPoint,
    /// delta * H
    pub delta_h: PolyPoint,
    /// gamma * H
    pub gamma_h: PolyPoint,
    /// z^i * H for i = 0..n
    pub z_h: Vec<PolyPoint>,
}

impl Crs {
    /// Create an empty CRS structure
    pub fn new(n: usize, m: usize, l: usize) -> Self {
        Self {
            n,
            m,
            l,
            alpha_g: Point::inf(),
            beta_g: Point::inf(),
            delta_g: Point::inf(),
            z_g: vec![Point::inf(); n],
            theta_g: vec![Point::inf(); m],
            zt_g: vec![Point::inf(); n - 1],
            beta_h: PolyPoint::new(),
            delta_h: PolyPoint::new(),
            gamma_h: PolyPoint::new(),
            z_h: vec![PolyPoint::new(); n],
        }
    }

    /// Generate CRS from QAP and system parameters.
    ///
    /// This is the "trusted setup" phase. The random values (alpha, beta, gamma, delta, z)
    /// must be securely destroyed after this function returns.
    pub fn generate(qap: &Qap, sys: &SigSystem) -> Self {
        let tor_mod = sys.torsion_modulus();
        let prime_mod = sys.modulus();
        let interp = LagrangeInterpolator::new(&sys.tor);

        let n = qap.n;
        let m = qap.m;
        let l = qap.l;

        // Generate toxic waste (random values)
        let alpha = tor_mod.rand();
        let beta = tor_mod.rand();
        let gamma = tor_mod.rand();
        let delta = tor_mod.rand();
        let z = tor_mod.rand();

        // Compute powers of z: z^0, z^1, ..., z^{n-1}
        let mut zpow = vec![Int::from(1); n];
        for i in 1..n {
            zpow[i] = tor_mod.mul(&zpow[i - 1], &z);
        }

        // Compute theta_j = beta * v_j(z) + alpha * w_j(z) + y_j(z)^2 for each wire j
        let mut theta = vec![Int::ZERO; m];
        for i in 0..m {
            // Evaluate v_j(z)
            let vj_z = interp.lcalc(&z, &qap.vj[i]);
            // Evaluate w_j(z)
            let wj_z = interp.lcalc(&z, &qap.wj[i]);
            // Evaluate y_j(z)
            let yj_z = interp.lcalc(&z, &qap.yj[i]);

            // theta_j = beta * v_j(z)
            theta[i] = tor_mod.mul(&vj_z, &beta);
            // + alpha * w_j(z)
            let alpha_wj = tor_mod.mul(&wj_z, &alpha);
            theta[i] = tor_mod.add(&theta[i], &alpha_wj);
            // + y_j(z)^2
            let yj_sq = tor_mod.mul(&yj_z, &yj_z);
            theta[i] = tor_mod.add(&theta[i], &yj_sq);
        }

        // Divide statement thetas by gamma, witness thetas by delta
        // Statement indexes: 0 to l (via sw permutation)
        for i in 0..=l {
            let wire_idx = qap.sw[i];
            theta[wire_idx] =
                tor_mod.div(&theta[wire_idx], &gamma).unwrap();
        }
        // Witness indexes: l+1 to m-1
        for i in (l + 1)..m {
            let wire_idx = qap.sw[i];
            theta[wire_idx] =
                tor_mod.div(&theta[wire_idx], &delta).unwrap();
        }

        // Compute t(z) / delta
        let t_z = interp.tofz(&z, &qap.list);
        let t_delta = tor_mod.div(&t_z, &delta).unwrap();

        // Create CRS points on G1
        let mut crs = Crs::new(n, m, l);

        // alpha * G
        crs.alpha_g = sys.e.mul(&sys.g1, &alpha);
        // beta * G
        crs.beta_g = sys.e.mul(&sys.g1, &beta);
        // delta * G
        crs.delta_g = sys.e.mul(&sys.g1, &delta);

        // z^i * G for i = 0..n
        for i in 0..n {
            crs.z_g[i] = sys.e.mul(&sys.g1, &zpow[i]);
        }

        // theta_j * G for each wire (ordered by sw permutation)
        for i in 0..m {
            let wire_idx = qap.sw[i];
            crs.theta_g[i] = sys.e.mul(&sys.g1, &theta[wire_idx]);
        }

        // z^i * t(z) / delta * G for i = 0..n-1
        let mut zt_coef = t_delta.clone();
        crs.zt_g[0] = sys.e.mul(&sys.g1, &zt_coef);
        for i in 1..(n - 1) {
            zt_coef = tor_mod.mul(&zt_coef, &z);
            crs.zt_g[i] = sys.e.mul(&sys.g1, &zt_coef);
        }

        // Create CRS points on G2 (extension curve)
        // beta * H
        crs.beta_h =
            sys.ex.mul(&sys.g2, &beta, &sys.irrd, &prime_mod);
        // gamma * H
        crs.gamma_h =
            sys.ex.mul(&sys.g2, &gamma, &sys.irrd, &prime_mod);
        // delta * H
        crs.delta_h =
            sys.ex.mul(&sys.g2, &delta, &sys.irrd, &prime_mod);

        // z^i * H for i = 0..n
        for i in 0..n {
            crs.z_h[i] =
                sys.ex.mul(&sys.g2, &zpow[i], &sys.irrd, &prime_mod);
        }

        crs
    }

    /// Load CRS from binary file.
    pub fn load(filename: &str) -> io::Result<Self> {
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

        let mut crs = Self::new(n, m, l);

        // Read G1 points
        crs.alpha_g = read_point(&mut reader)?;
        crs.beta_g = read_point(&mut reader)?;
        crs.delta_g = read_point(&mut reader)?;

        for i in 0..n {
            crs.z_g[i] = read_point(&mut reader)?;
        }
        for i in 0..m {
            crs.theta_g[i] = read_point(&mut reader)?;
        }
        for i in 0..(n - 1) {
            crs.zt_g[i] = read_point(&mut reader)?;
        }

        // Read G2 points
        crs.beta_h = read_ext_point(&mut reader)?;
        crs.delta_h = read_ext_point(&mut reader)?;
        crs.gamma_h = read_ext_point(&mut reader)?;

        for i in 0..n {
            crs.z_h[i] = read_ext_point(&mut reader)?;
        }

        Ok(crs)
    }

    /// Save CRS to binary file.
    pub fn save(&self, filename: &str) -> io::Result<()> {
        let file = File::create(filename)?;
        let mut writer = BufWriter::new(file);

        // Write dimensions
        writer.write_all(&(self.n as i32).to_ne_bytes())?;
        writer.write_all(&(self.m as i32).to_ne_bytes())?;
        writer.write_all(&(self.l as i32).to_ne_bytes())?;

        // Write G1 points
        write_point(&mut writer, &self.alpha_g)?;
        write_point(&mut writer, &self.beta_g)?;
        write_point(&mut writer, &self.delta_g)?;

        for i in 0..self.n {
            write_point(&mut writer, &self.z_g[i])?;
        }
        for i in 0..self.m {
            write_point(&mut writer, &self.theta_g[i])?;
        }
        for i in 0..(self.n - 1) {
            write_point(&mut writer, &self.zt_g[i])?;
        }

        // Write G2 points
        write_ext_point(&mut writer, &self.beta_h)?;
        write_ext_point(&mut writer, &self.delta_h)?;
        write_ext_point(&mut writer, &self.gamma_h)?;

        for i in 0..self.n {
            write_ext_point(&mut writer, &self.z_h[i])?;
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_crs_new() {
        let crs = Crs::new(5, 10, 4);
        assert_eq!(crs.n, 5);
        assert_eq!(crs.m, 10);
        assert_eq!(crs.l, 4);
        assert_eq!(crs.z_g.len(), 5);
        assert_eq!(crs.theta_g.len(), 10);
        assert_eq!(crs.zt_g.len(), 4);
    }
}
