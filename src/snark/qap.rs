//! Quadratic Arithmetic Program (QAP) representation.
//!
//! A QAP encodes an arithmetic circuit as polynomial constraints.
//! For a circuit with n gates and m wires, we have:
//! - v_j(x), w_j(x), y_j(x) polynomials for each wire j
//! - t(x) = (x - r_0)(x - r_1)...(x - r_{n-1}) target polynomial
//! - h(x) = (v(x) * w(x) - y(x)) / t(x) quotient polynomial
//!
//! Ported from snark_qap.c

use rug::Integer as Int;
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Read, Write};

use super::lagrange::LagrangeInterpolator;
use super::system::{read_int, write_int};

/// Quadratic Arithmetic Program parameters.
pub struct Qap {
    /// Number of gates (evaluation points)
    pub n: usize,
    /// Number of wires (variables)
    pub m: usize,
    /// Number of public (statement) inputs
    pub l: usize,
    /// Wire permutation: maps wire index to statement/witness ordering
    pub sw: Vec<usize>,
    /// Gate evaluation points (primes r_0, ..., r_{n-1})
    pub list: Vec<Int>,
    /// Left input coefficients: vj[wire][gate]
    pub vj: Vec<Vec<Int>>,
    /// Right input coefficients: wj[wire][gate]
    pub wj: Vec<Vec<Int>>,
    /// Output coefficients: yj[wire][gate]
    pub yj: Vec<Vec<Int>>,
    /// Precomputed h(x) coefficient table for L_i * L_j / t(x)
    pub htable: Vec<Vec<Int>>,
}

impl Qap {
    /// Create a new empty QAP
    pub fn new(n: usize, m: usize, l: usize) -> Self {
        Self {
            n,
            m,
            l,
            sw: vec![0; m],
            list: vec![Int::ZERO; n],
            vj: vec![vec![Int::ZERO; n]; m],
            wj: vec![vec![Int::ZERO; n]; m],
            yj: vec![vec![Int::ZERO; n]; m],
            htable: Vec::new(),
        }
    }

    /// Load QAP parameters from binary file.
    ///
    /// File format matches snark_qap.c output.
    pub fn load(filename: &str) -> io::Result<Self> {
        let file = File::open(filename)?;
        let mut reader = BufReader::new(file);

        // Read dimensions
        let mut buf = [0u8; 4];
        reader.read_exact(&mut buf)?;
        let n = i32::from_ne_bytes(buf) as usize;

        reader.read_exact(&mut buf)?;
        let m = i32::from_ne_bytes(buf) as usize;

        // Read wire permutation
        let mut sw = vec![0usize; m];
        for item in sw.iter_mut() {
            reader.read_exact(&mut buf)?;
            *item = i32::from_ne_bytes(buf) as usize;
        }

        // Read l (number of public inputs)
        reader.read_exact(&mut buf)?;
        let l = i32::from_ne_bytes(buf) as usize;

        // Read gate evaluation points
        let mut list = vec![Int::ZERO; n];
        for item in list.iter_mut() {
            *item = read_int(&mut reader)?;
        }

        // Read coefficient matrices (n * m elements each)
        let k = n * m;
        let mut vj_flat = vec![Int::ZERO; k];
        let mut wj_flat = vec![Int::ZERO; k];
        let mut yj_flat = vec![Int::ZERO; k];

        for item in vj_flat.iter_mut() {
            *item = read_int(&mut reader)?;
        }
        for item in wj_flat.iter_mut() {
            *item = read_int(&mut reader)?;
        }
        for item in yj_flat.iter_mut() {
            *item = read_int(&mut reader)?;
        }

        // Convert to 2D arrays [wire][gate]
        let mut vj = vec![vec![Int::ZERO; n]; m];
        let mut wj = vec![vec![Int::ZERO; n]; m];
        let mut yj = vec![vec![Int::ZERO; n]; m];

        for wire in 0..m {
            for gate in 0..n {
                vj[wire][gate] = vj_flat[wire * n + gate].clone();
                wj[wire][gate] = wj_flat[wire * n + gate].clone();
                yj[wire][gate] = yj_flat[wire * n + gate].clone();
            }
        }

        // Read h(x) table: (n-1)^2 * n / 2 elements
        // Actually it's n*(n-1)/2 pairs, each with n-1 coefficients
        let h_pairs = n * (n - 1) / 2;
        let h_coef_len = n - 1;
        let mut htable = vec![vec![Int::ZERO; h_coef_len]; h_pairs];

        for pair in htable.iter_mut() {
            for coef in pair.iter_mut() {
                *coef = read_int(&mut reader)?;
            }
        }

        Ok(Self {
            n,
            m,
            l,
            sw,
            list,
            vj,
            wj,
            yj,
            htable,
        })
    }

    /// Save QAP parameters to binary file.
    pub fn save(&self, filename: &str) -> io::Result<()> {
        let file = File::create(filename)?;
        let mut writer = BufWriter::new(file);

        // Write dimensions
        writer.write_all(&(self.n as i32).to_ne_bytes())?;
        writer.write_all(&(self.m as i32).to_ne_bytes())?;

        // Write wire permutation
        for &s in &self.sw {
            writer.write_all(&(s as i32).to_ne_bytes())?;
        }

        // Write l
        writer.write_all(&(self.l as i32).to_ne_bytes())?;

        // Write gate evaluation points
        for item in &self.list {
            write_int(&mut writer, item)?;
        }

        // Write coefficient matrices (flattened by wire then gate)
        for wire in 0..self.m {
            for gate in 0..self.n {
                write_int(&mut writer, &self.vj[wire][gate])?;
            }
        }
        for wire in 0..self.m {
            for gate in 0..self.n {
                write_int(&mut writer, &self.wj[wire][gate])?;
            }
        }
        for wire in 0..self.m {
            for gate in 0..self.n {
                write_int(&mut writer, &self.yj[wire][gate])?;
            }
        }

        // Write h(x) table
        for pair in &self.htable {
            for coef in pair {
                write_int(&mut writer, coef)?;
            }
        }

        Ok(())
    }

    /// Get the coefficient polynomial for wire j's left input v_j(x).
    /// Returns coefficients in decreasing power order.
    pub fn get_vj_coefs(&self, wire: usize) -> &[Int] {
        &self.vj[wire]
    }

    /// Get the coefficient polynomial for wire j's right input w_j(x).
    pub fn get_wj_coefs(&self, wire: usize) -> &[Int] {
        &self.wj[wire]
    }

    /// Get the coefficient polynomial for wire j's output y_j(x).
    pub fn get_yj_coefs(&self, wire: usize) -> &[Int] {
        &self.yj[wire]
    }
}

/// Builder for constructing a QAP from circuit constraints.
pub struct QapBuilder {
    /// Number of gates
    n: usize,
    /// Number of wires
    m: usize,
    /// Number of public inputs
    l: usize,
    /// Lagrange interpolator
    interp: LagrangeInterpolator,
    /// Gate evaluation points
    list: Vec<Int>,
    /// Wire permutation
    sw: Vec<usize>,
    /// Left coefficients by wire
    vj: Vec<Vec<Int>>,
    /// Right coefficients by wire
    wj: Vec<Vec<Int>>,
    /// Output coefficients by wire
    yj: Vec<Vec<Int>>,
}

impl QapBuilder {
    /// Create a new QAP builder.
    ///
    /// # Arguments
    /// * `n` - Number of gates (constraints)
    /// * `m` - Number of wires (variables)
    /// * `l` - Number of public inputs (statement variables)
    /// * `modulus` - Prime modulus for coefficient arithmetic
    pub fn new(n: usize, m: usize, l: usize, modulus: &Int) -> Self {
        Self {
            n,
            m,
            l,
            interp: LagrangeInterpolator::new(modulus),
            list: vec![Int::ZERO; n],
            sw: (0..m).collect(), // Identity permutation by default
            vj: vec![vec![Int::ZERO; n]; m],
            wj: vec![vec![Int::ZERO; n]; m],
            yj: vec![vec![Int::ZERO; n]; m],
        }
    }

    /// Set the gate evaluation points (should be distinct).
    pub fn set_gate_points(&mut self, points: &[Int]) {
        assert_eq!(
            points.len(),
            self.n,
            "Must provide n gate points"
        );
        self.list = points.to_vec();
    }

    /// Set the wire permutation (maps internal index to statement/witness ordering).
    pub fn set_permutation(&mut self, perm: &[usize]) {
        assert_eq!(
            perm.len(),
            self.m,
            "Must provide m permutation values"
        );
        self.sw = perm.to_vec();
    }

    /// Compute Lagrange basis polynomial L_i(x) coefficients.
    pub fn lagrange(&self, i: usize) -> Vec<Int> {
        self.interp.liofx(i, &self.list)
    }

    /// Set left input coefficients for a wire.
    pub fn set_v(&mut self, wire: usize, coefs: Vec<Int>) {
        assert!(wire < self.m, "Wire index out of bounds");
        assert_eq!(
            coefs.len(),
            self.n,
            "Coefficient count must equal n"
        );
        self.vj[wire] = coefs;
    }

    /// Set right input coefficients for a wire.
    pub fn set_w(&mut self, wire: usize, coefs: Vec<Int>) {
        assert!(wire < self.m, "Wire index out of bounds");
        assert_eq!(
            coefs.len(),
            self.n,
            "Coefficient count must equal n"
        );
        self.wj[wire] = coefs;
    }

    /// Set output coefficients for a wire.
    pub fn set_y(&mut self, wire: usize, coefs: Vec<Int>) {
        assert!(wire < self.m, "Wire index out of bounds");
        assert_eq!(
            coefs.len(),
            self.n,
            "Coefficient count must equal n"
        );
        self.yj[wire] = coefs;
    }

    /// Add two Lagrange polynomial coefficient vectors.
    pub fn add_coefs(&self, a: &[Int], b: &[Int]) -> Vec<Int> {
        let m = &self.interp.modulus;
        a.iter().zip(b.iter()).map(|(x, y)| m.add(x, y)).collect()
    }

    /// Build the final QAP with precomputed h(x) table.
    pub fn build(self) -> Qap {
        // Compute h(x) coefficient table
        let htable = self.interp.all_lilj(&self.list);

        Qap {
            n: self.n,
            m: self.m,
            l: self.l,
            sw: self.sw,
            list: self.list,
            vj: self.vj,
            wj: self.wj,
            yj: self.yj,
            htable,
        }
    }
}

/// Create the example QAP from snark_qap.c
///
/// This represents a specific arithmetic circuit with:
/// - 5 gates
/// - 10 wires
/// - 4 public inputs
pub fn create_example_qap(modulus: &Int) -> Qap {
    let mut builder = QapBuilder::new(5, 10, 4, modulus);

    // Set gate evaluation points (primes)
    builder.set_gate_points(&[
        Int::from(31),
        Int::from(37),
        Int::from(41),
        Int::from(43),
        Int::from(47),
    ]);

    // Set wire permutation for statement/witness ordering
    builder.set_permutation(&[0, 1, 2, 5, 7, 3, 4, 6, 8, 9]);

    // Set up constraint polynomials (from snark_qap.c)
    // v0 = L_0
    let l0 = builder.lagrange(0);
    builder.set_v(0, l0.clone());

    // w1 = L_0 + L_2
    let l2 = builder.lagrange(2);
    let w1 = builder.add_coefs(&l0, &l2);
    builder.set_w(1, w1);

    // v2 = L_1
    let l1 = builder.lagrange(1);
    builder.set_v(2, l1.clone());

    // w2 = L_2
    builder.set_w(2, l2.clone());

    // w4 = L_3
    let l3 = builder.lagrange(3);
    builder.set_w(4, l3.clone());

    // w3 = L_1 + L_3
    let w3 = builder.add_coefs(&l1, &l3);
    builder.set_w(3, w3);

    // v5 = L_2
    builder.set_v(5, l2.clone());

    // y5 = L_0
    builder.set_y(5, l0.clone());

    // v6 = L_3
    builder.set_v(6, l3.clone());

    // y6 = L_1
    builder.set_y(6, l1.clone());

    // v7 = L_4
    let l4 = builder.lagrange(4);
    builder.set_v(7, l4.clone());

    // y7 = L_2
    builder.set_y(7, l2.clone());

    // w8 = L_4
    builder.set_w(8, l4.clone());

    // y8 = L_3
    builder.set_y(8, l3);

    // y9 = L_4
    builder.set_y(9, l4);

    builder.build()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_qap_builder_basic() {
        let modulus = Int::from(101);
        let builder = QapBuilder::new(3, 5, 2, &modulus);

        // Just verify construction doesn't panic
        assert_eq!(builder.n, 3);
        assert_eq!(builder.m, 5);
        assert_eq!(builder.l, 2);
    }

    #[test]
    fn test_example_qap() {
        // Use a large prime (like the torsion value would be)
        let modulus = Int::from(1000003);
        let qap = create_example_qap(&modulus);

        assert_eq!(qap.n, 5);
        assert_eq!(qap.m, 10);
        assert_eq!(qap.l, 4);
        assert_eq!(qap.list.len(), 5);
        assert_eq!(qap.htable.len(), 10); // 5*4/2 = 10 pairs
    }
}
