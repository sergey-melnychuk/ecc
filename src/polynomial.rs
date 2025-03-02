use rug::Integer as Int;

use crate::modulus::Modulus;

#[derive(Clone, Debug, PartialEq)]
pub struct Polynomial {
    coef: Vec<Int>,
}

#[allow(clippy::needless_range_loop)]
impl Polynomial {
    pub fn zeros(len: usize) -> Self {
        Self {
            coef: vec![Int::ZERO; len],
        }
    }

    #[allow(clippy::new_without_default)]
    pub fn new() -> Self {
        Self { coef: Vec::new() }
    }

    pub fn degree(&self) -> usize {
        self.coef.len() - 1
    }

    pub fn is_zero(&self) -> bool {
        self.coef.is_empty() || self.coef.iter().all(|k| k.is_zero())
    }

    pub fn trim(mut self) -> Self {
        self.coef = self
            .coef
            .into_iter()
            .skip_while(|k| k.is_zero())
            .collect::<Vec<_>>();
        self
    }

    pub fn set(&mut self, pow: usize, int: Int) {
        if pow > self.degree() {
            self.coef.resize(pow + 1, Int::ZERO);
        }
        self.coef[pow] = int;
    }

    pub fn get(&self, pow: usize) -> Int {
        self.coef.get(pow).cloned().unwrap_or(Int::ZERO)
    }

    pub fn neg(&self) -> Self {
        let mut ret = Self::new();
        for i in 0..ret.coef.len() {
            ret.coef[i] = Int::from(-&self.coef[i]);
        }
        ret
    }

    pub fn add(&self, that: &Self) -> Self {
        let mut ret = Self::new();
        let len = self.coef.len().max(that.coef.len());
        ret.coef.resize(len, Int::ZERO);
        for i in 0..len {
            let lhs = self.coef.get(i).cloned().unwrap_or(Int::ZERO);
            let rhs = that.coef.get(i).cloned().unwrap_or(Int::ZERO);
            ret.coef[i] = lhs + rhs;
        }
        ret.trim()
    }

    pub fn sub(&self, that: &Self) -> Self {
        let that = that.neg();
        self.add(&that)
    }

    // Polynomials multiplication & reduction by irreducible polynomial over prime field
    pub fn mul(
        &self,
        that: &Self,
        irreducible: &Self,
        m: &Modulus,
    ) -> Self {
        let r = self.degree() + that.degree();
        let mut ret = Self::zeros(r + 1);

        for i in 0..=r {
            for j in 0..=i {
                let x = m.mul(&self.get(i), &that.get(i - j));
                let k = m.add(&ret.get(i), &x);
                ret.set(i, k);
            }
        }

        // Listing 8.7
        if r < irreducible.degree() {
            return ret;
        }

        let mut out = Self::zeros(irreducible.degree());
        for i in 0..irreducible.degree() {
            out.set(i, irreducible.get(i))
        }

        let table = irreducible.mulprep(m);
        for i in irreducible.degree()..r {
            for j in 0..irreducible.degree() {
                let t = m.mul(&ret.get(i), &table[i][j]);
                let k = m.add(&out.get(j), &t);
                out.set(j, k);
            }
        }

        out.trim()
    }

    fn mulprep(&self, m: &Modulus) -> Vec<Vec<Int>> {
        let degree = self.degree();
        let mut table = vec![vec![Int::ZERO; degree]; degree * 2];

        // Listing 8.2
        for i in 0..degree {
            table[i][i] = Int::from(1);
        }

        // Listing 8.3
        let norm = self.normal(m);
        for j in 0..degree {
            table[degree][j] = m.neg(&norm.get(j));
        }

        // Listing 8.4
        for i in (degree + 1)..(degree * 2) {
            for j in 1..degree {
                let t = m.mul(
                    &table[degree][j],
                    &table[i - 1][degree - 1],
                );
                table[i][j] = m.add(&table[i - 1][j - 1], &t);
            }
            table[i][0] =
                m.mul(&table[degree][0], &table[i - 1][degree - 1]);
        }

        table
    }

    // Listing 8.9
    fn normal(&self, m: &Modulus) -> Self {
        let mut ret = self.clone().trim();

        let k = ret.get(ret.degree());
        if ret.is_zero() || k == 1 {
            return ret;
        }

        let d = m.inv(&k).unwrap();
        for k in ret.coef.iter_mut() {
            *k = m.mul(k, &d);
        }
        ret
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_polynomial_example() {
        // Listing 8.{11,12}
        const M: i32 = 7;

        let m = Modulus::new(&Int::from(M));
        let mut p = Polynomial::zeros(5);
        p.set(4, Int::from(1));
        p.set(3, Int::from(2));
        p.set(2, Int::from(1));
        p.set(1, Int::from(3));
        p.set(0, Int::from(5));

        let table = p.mulprep(&m);
        for row in &table {
            println!("{:?}", row);
        }

        assert_eq!(table[4], vec![M - 5, M - 3, M - 1, M - 2]);
        assert_eq!(table[7], vec![0, 6, 1, 5]);
    }
}
