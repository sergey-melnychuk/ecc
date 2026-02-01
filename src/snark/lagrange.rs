//! Lagrange interpolation utilities for QAP construction.
//!
//! Ported from snarkbase.c - provides functions for computing Lagrange
//! interpolation polynomials and related operations needed for SNARKs.

use rug::Integer as Int;

use crate::modulus::Modulus;

/// Lagrange interpolation utilities for SNARK QAP construction.
pub struct LagrangeInterpolator {
    /// The prime modulus for all operations
    pub modulus: Modulus,
}

impl LagrangeInterpolator {
    /// Create a new interpolator with the given modulus
    pub fn new(modulus: &Int) -> Self {
        Self {
            modulus: Modulus::new(modulus),
        }
    }

    /// Compute p_i = product of (r_i - r_j) for all j != i
    ///
    /// This is the denominator of the Lagrange basis polynomial L_i(x).
    /// Maps to C function: p_i()
    ///
    /// # Arguments
    /// * `i` - Index of the Lagrange basis polynomial
    /// * `list` - List of evaluation points (r_0, r_1, ..., r_{n-1})
    pub fn p_i(&self, i: usize, list: &[Int]) -> Int {
        let n = list.len();
        let mut result = Int::from(1);

        for j in 0..n {
            if j == i {
                continue;
            }
            let term = self.modulus.sub(&list[i], &list[j]);
            result = self.modulus.mul(&result, &term);
        }

        result
    }

    /// Create list of start and end indexes for coefficient computation.
    ///
    /// Maps to C function: startk()
    fn startk(
        &self,
        rows: usize,
        limit: usize,
    ) -> Vec<(usize, usize)> {
        let mut k = vec![(0usize, 0usize); rows];

        let mut index = 0;
        for m in 0..rows {
            k[m].0 = index;
            index += 1;
        }

        let mut index = limit - 1;
        for m in (0..rows).rev() {
            k[m].1 = index;
            if index > 0 {
                index -= 1;
            }
        }

        k
    }

    /// Compute Lagrange interpolant coefficients.
    ///
    /// If i == j: Returns coefficients for L_i(x) (degree n-1)
    /// If i != j: Returns coefficients for L_i(x) * L_j(x) / t(x) (degree n-2)
    ///
    /// Coefficients are in decreasing power: [x^deg, x^(deg-1), ..., x^0]
    ///
    /// Maps to C function: li_lj()
    ///
    /// # Arguments
    /// * `i` - First Lagrange index
    /// * `j` - Second Lagrange index (same as i for single polynomial)
    /// * `list` - List of evaluation points
    pub fn li_lj(
        &self,
        i: usize,
        j: usize,
        list: &[Int],
    ) -> Vec<Int> {
        let n = list.len();
        let cf_limit = if i == j { n - 1 } else { n - 2 };

        // Build subset list: -r_m for m != i and m != j
        let mut sublist: Vec<Int> = Vec::with_capacity(cf_limit);
        for (m, item) in list.iter().enumerate() {
            if m != i && m != j {
                sublist.push(self.modulus.neg(item));
            }
        }

        let mut coef = vec![Int::ZERO; cf_limit + 1];
        coef[0] = Int::from(1);

        for cfdex in 1..=cf_limit {
            let mut sum = Int::ZERO;
            let mut k = self.startk(cfdex, cf_limit);
            let mut done = false;

            while !done {
                // Compute product of selected terms
                let mut term = sublist[k[0].0].clone();
                for m in 1..cfdex {
                    term = self.modulus.mul(&term, &sublist[k[m].0]);
                }
                sum = self.modulus.add(&sum, &term);

                // Advance to next combination
                let mut m = cfdex - 1;
                let mut bmp = false;

                while !bmp && !done {
                    if k[m].0 != k[m].1 {
                        bmp = true;
                        k[m].0 += 1;
                    } else {
                        loop {
                            if m == 0 {
                                done = true;
                                bmp = true;
                                break;
                            }
                            m -= 1;
                            if k[m].0 != k[m].1 {
                                k[m].0 += 1;
                                while m < cfdex - 1 {
                                    m += 1;
                                    k[m].0 = k[m - 1].0 + 1;
                                }
                                bmp = true;
                                break;
                            }
                        }
                    }
                }
            }
            coef[cfdex] = sum;
        }

        coef
    }

    /// Compute coefficients for L_i(x).
    ///
    /// Returns n coefficients for degree n-1 polynomial in order [x^(n-1), ..., x^0]
    ///
    /// Maps to C function: liofx()
    pub fn liofx(&self, i: usize, list: &[Int]) -> Vec<Int> {
        let n = list.len();
        let mut coef = self.li_lj(i, i, list);

        // Divide by p_i to get actual Lagrange polynomial
        let pi = self.p_i(i, list);
        let pi_inv =
            self.modulus.inv(&pi).expect("p_i should be invertible");

        for c in coef.iter_mut().take(n) {
            *c = self.modulus.mul(c, &pi_inv);
        }

        coef
    }

    /// Compute coefficients for L_i(x) * L_j(x) / t(x).
    ///
    /// Returns n-1 coefficients for degree n-2 polynomial in order [x^(n-2), ..., x^0]
    ///
    /// Maps to C function: liljofx()
    pub fn liljofx(
        &self,
        i: usize,
        j: usize,
        list: &[Int],
    ) -> Vec<Int> {
        let n = list.len();
        let mut coef = self.li_lj(i, j, list);

        // Divide by p_i * p_j
        let pi = self.p_i(i, list);
        let pj = self.p_i(j, list);
        let pi_pj = self.modulus.mul(&pi, &pj);
        let inv = self
            .modulus
            .inv(&pi_pj)
            .expect("pi*pj should be invertible");

        for c in coef.iter_mut().take(n - 1) {
            *c = self.modulus.mul(c, &inv);
        }

        coef
    }

    /// Create table of all L_i(x) * L_j(x) / t(x) coefficients.
    ///
    /// Returns a flattened table where each entry contains n-1 coefficients.
    /// Total entries: n*(n-1)/2 cross products.
    ///
    /// Maps to C function: all_lilj()
    pub fn all_lilj(&self, list: &[Int]) -> Vec<Vec<Int>> {
        let n = list.len();
        let num_pairs = n * (n - 1) / 2;
        let mut table = Vec::with_capacity(num_pairs);

        for i in 0..(n - 1) {
            for j in (i + 1)..n {
                table.push(self.liljofx(i, j, list));
            }
        }

        table
    }

    /// Evaluate polynomial at z using coefficients.
    ///
    /// Coefficients are in decreasing power: [c_0*x^deg + c_1*x^(deg-1) + ... + c_deg]
    /// Uses Horner's method for efficiency.
    ///
    /// Maps to C function: lcalc()
    pub fn lcalc(&self, z: &Int, coef: &[Int]) -> Int {
        if coef.is_empty() {
            return Int::ZERO;
        }

        let mut result = coef[0].clone();
        for c in coef.iter().skip(1) {
            result = self.modulus.mul(&result, z);
            result = self.modulus.add(&result, c);
        }

        result
    }

    /// Compute t(z) = product of (z - r_i) for all evaluation points.
    ///
    /// This is the "target polynomial" that vanishes at all gate points.
    ///
    /// Maps to C function: tofzgrth()
    pub fn tofz(&self, z: &Int, list: &[Int]) -> Int {
        let mut result = Int::from(1);

        for item in list {
            let term = self.modulus.sub(z, item);
            result = self.modulus.mul(&result, &term);
        }

        result
    }

    /// Matrix-vector multiply and sum columns.
    ///
    /// Computes: for each column j, sum over rows i of (mat[i][j] * coef[i])
    ///
    /// This is used for combining CRS values with witness coefficients.
    ///
    /// Maps to C function: matflat()
    ///
    /// # Arguments
    /// * `mat` - Matrix as row-major vector (length rows × width)
    /// * `width` - Number of columns
    /// * `coef` - Coefficient vector (length rows)
    pub fn matflat(
        &self,
        mat: &[Int],
        width: usize,
        coef: &[Int],
    ) -> Vec<Int> {
        let length = coef.len();
        let mut result = vec![Int::ZERO; width];

        for j in 0..width {
            for i in 0..length {
                let product =
                    self.modulus.mul(&mat[i * width + j], &coef[i]);
                result[j] = self.modulus.add(&result[j], &product);
            }
        }

        result
    }

    /// Flatten nested table into single vector
    pub fn flatten_table(table: &[Vec<Int>]) -> Vec<Int> {
        table.iter().flatten().cloned().collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_p_i() {
        // Use prime 101
        let m = Int::from(101);
        let interp = LagrangeInterpolator::new(&m);

        // Simple case: list = [1, 2, 3]
        let list = vec![Int::from(1), Int::from(2), Int::from(3)];

        // p_0 = (1-2)(1-3) = (-1)(-2) = 2
        let p0 = interp.p_i(0, &list);
        assert_eq!(p0, Int::from(2));

        // p_1 = (2-1)(2-3) = (1)(-1) = -1 = 100 mod 101
        let p1 = interp.p_i(1, &list);
        assert_eq!(p1, Int::from(100));

        // p_2 = (3-1)(3-2) = (2)(1) = 2
        let p2 = interp.p_i(2, &list);
        assert_eq!(p2, Int::from(2));
    }

    #[test]
    fn test_tofz() {
        let m = Int::from(101);
        let interp = LagrangeInterpolator::new(&m);

        let list = vec![Int::from(1), Int::from(2), Int::from(3)];

        // t(4) = (4-1)(4-2)(4-3) = 3*2*1 = 6
        let t4 = interp.tofz(&Int::from(4), &list);
        assert_eq!(t4, Int::from(6));

        // t(1) = 0 (since 1 is in the list)
        let t1 = interp.tofz(&Int::from(1), &list);
        assert_eq!(t1, Int::ZERO);
    }

    #[test]
    fn test_lcalc_horner() {
        let m = Int::from(101);
        let interp = LagrangeInterpolator::new(&m);

        // p(x) = 2x^2 + 3x + 5 (coef in decreasing order: [2, 3, 5])
        let coef = vec![Int::from(2), Int::from(3), Int::from(5)];

        // p(2) = 2*4 + 3*2 + 5 = 8 + 6 + 5 = 19
        let result = interp.lcalc(&Int::from(2), &coef);
        assert_eq!(result, Int::from(19));

        // p(0) = 5
        let result = interp.lcalc(&Int::from(0), &coef);
        assert_eq!(result, Int::from(5));
    }

    #[test]
    fn test_matflat() {
        let m = Int::from(101);
        let interp = LagrangeInterpolator::new(&m);

        // 2x3 matrix (2 rows, 3 columns):
        // [1 2 3]
        // [4 5 6]
        let mat = vec![
            Int::from(1),
            Int::from(2),
            Int::from(3),
            Int::from(4),
            Int::from(5),
            Int::from(6),
        ];

        // Coefficients [2, 3]
        let coef = vec![Int::from(2), Int::from(3)];

        // Result[0] = 1*2 + 4*3 = 2 + 12 = 14
        // Result[1] = 2*2 + 5*3 = 4 + 15 = 19
        // Result[2] = 3*2 + 6*3 = 6 + 18 = 24
        let result = interp.matflat(&mat, 3, &coef);
        assert_eq!(
            result,
            vec![Int::from(14), Int::from(19), Int::from(24)]
        );
    }
}
