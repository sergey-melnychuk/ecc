use rug::Integer as Int;

#[derive(Debug, Default, PartialEq)]
pub struct Polynomial {
    coef: Vec<Int>,
}

impl Polynomial {
    pub fn new() -> Self {
        Self { coef: Vec::new() }
    }

    pub fn degree(&self) -> usize {
        self.coef.len() - 1
    }

    pub fn is_zero(&self) -> bool {
        self.coef.is_empty() || self.coef.iter().all(|k| k.is_zero())
    }

    pub fn normalize(&mut self) {
        while !self.coef.is_empty() && self.coef[0].is_zero() {
            self.coef.remove(0);
        }
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

    pub fn add(&self, that: &Self) -> Self {
        let mut ret = Self::new();
        let len = self.coef.len().max(that.coef.len());
        ret.coef.resize(len, Int::ZERO);
        for i in 0..len {
            let lhs = self.coef.get(i).cloned().unwrap_or(Int::ZERO);
            let rhs = that.coef.get(i).cloned().unwrap_or(Int::ZERO);
            ret.coef[i] = lhs + rhs;
        }
        ret.normalize();
        ret
    }

    pub fn neg(&self) -> Self {
        let mut ret = Self::new();
        for i in 0..ret.coef.len() {
            ret.coef[i] = Int::from(-&self.coef[i]);
        }
        ret
    }

    pub fn sub(&self, that: &Self) -> Self {
        let that = that.neg();
        self.add(&that)
    }
}
