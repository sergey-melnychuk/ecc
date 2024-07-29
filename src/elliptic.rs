use rug::Integer as Int;

use crate::modulus::Modulus;

// y^2 = x^3 + ax + b mod n
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Curve {
    pub modulus: Int,
    pub order: Int,
    pub cofactor: u64,
    pub base: Point,
    pub a: Int,
    pub b: Int,
}

impl Curve {
    pub fn new(
        modulus: Int,
        order: Int,
        cofactor: u64,
        base: Point,
        a: Int,
        b: Int,
    ) -> Self {
        Self {
            modulus,
            order,
            cofactor,
            base,
            a,
            b,
        }
    }

    pub fn add(&self, p: &Point, q: &Point) -> Point {
        if p.is_inf() {
            return q.clone();
        }
        if q.is_inf() {
            return p.clone();
        }

        let m = Modulus::new(&self.modulus);

        let (num, den) = if m.add(&p.y, &q.y).is_zero() {
            let num = m.sub(&q.y, &p.y);
            let den = m.sub(&q.x, &q.y);
            (num, den)
        } else {
            let num = m.add(
                &m.add(&m.mul(&p.x, &p.x), &m.mul(&p.x, &q.x)),
                &m.add(&m.mul(&q.x, &q.x), &self.a),
            );
            let den = m.add(&p.y, &q.y);
            (num, den)
        };
        let den = m.inv(&den).expect("mod inverse");

        let lambda = m.mul(&num, &den);
        let x = m.sub(&m.sub(&m.mul(&lambda, &lambda), &p.x), &q.x);
        let y = m.sub(&m.mul(&lambda, &m.sub(&p.x, &x)), &p.y);

        let ret = Point::new(x, y);
        assert!(self.fits(&ret), "result point must fit the curve");
        ret
    }

    pub fn mul(&self, p: &Point, k: &Int) -> Point {
        if p.is_inf() {
            return p.clone();
        }
        let mut acc = p.clone();
        let mut bit = k.significant_bits() - 1;
        while bit > 0 {
            acc = self.add(&acc, &acc);
            if k.get_bit(bit - 1) {
                acc = self.add(&acc, p);
            }
            bit -= 1;
        }
        assert!(self.fits(&acc), "result point must fit the curve");
        acc
    }

    pub fn fits(&self, p: &Point) -> bool {
        let m = Modulus::new(&self.modulus);
        let lhs = m.mul(&p.y, &p.y);
        let rhs = m.add(
            &m.mul(&m.mul(&p.x, &p.x), &p.x),
            &m.add(&m.mul(&self.a, &p.x), &self.b),
        );
        lhs == rhs
    }

    pub fn apply(&self, x: &Int) -> Option<(Point, Point)> {
        let m = Modulus::new(&self.modulus);
        let rhs = m.add(
            &m.mul(&m.mul(x, x), x),
            &m.add(&m.mul(&self.a, x), &self.b),
        );
        let y = m.sqrt(&rhs)?;

        let hi = Point::new(x.clone(), y.clone());
        assert!(self.fits(&hi), "result point must fit the curve");
        let lo = Point::new(x.clone(), m.neg(&y));
        assert!(self.fits(&lo), "result point must fit the curve");
        if lo.y < hi.y {
            Some((lo, hi))
        } else {
            Some((hi, lo))
        }
    }

    pub fn find(&self, x: &Int, mut limit: usize) -> Option<Point> {
        let mut x = x.clone();
        while limit > 0 {
            if let Some((ret, _)) = self.apply(&x) {
                return Some(ret);
            }
            x += 1;
            limit -= 1;
        }
        None
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Point {
    pub x: Int,
    pub y: Int,
}

impl Point {
    pub fn new(x: Int, y: Int) -> Self {
        Self { x, y }
    }

    pub fn inf() -> Self {
        Self {
            x: Int::ZERO,
            y: Int::ZERO,
        }
    }

    pub fn is_inf(&self) -> bool {
        self.x.is_zero() && self.y.is_zero()
    }
}

pub mod curves {
    use super::*;
    use crate::{dec, hex};

    // https://github.com/drmike8888/Elliptic-curve-pairings/blob/main/Build_all/Curve_256_params.dat
    pub fn curve_256() -> Curve {
        let modulus = hex("2b000000000000000000000000000000000000000000000000000000000000001");
        let order = hex("2b0000000000000000000000000000002e7f521c85bba055a6e2161b956a47f69");
        let cofactor = 1;
        let a = hex("1");
        let b = hex("a87");

        let x = hex("2310115d283e49377820195c8e67781b6f112a625b14b747fa4cc13d06eba0919");
        let y = hex("51277aeb91946f0cb83053a10f67c5a9ef00a4f0cf2466b3bedf4fdcd774b574");
        let base = Point::new(x, y);

        let curve =
            Curve::new(modulus, order, cofactor, base.clone(), a, b);
        assert!(curve.fits(&base));
        curve
    }

    // https://hackmd.io/@jpw/bn254
    pub fn curve_bn254() -> Curve {
        let modulus = dec("21888242871839275222246405745257275088696311157297823662689037894645226208583");
        let order = dec("21888242871839275222246405745257275088548364400416034343698204186575808495617");
        let cofactor = 1;
        let a = dec("0");
        let b = dec("3");

        let x = dec("1");
        let y = dec("2");
        let base = Point::new(x, y);

        let curve =
            Curve::new(modulus, order, cofactor, base.clone(), a, b);
        assert!(curve.fits(&base));
        curve
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::elliptic::curves::{curve_256, curve_bn254};

    #[test]
    fn test_curve_256() {
        let ec = curve_256();
        assert!(ec.fits(&ec.base));
    }

    #[test]
    fn test_curve_bn254() {
        let ec = curve_bn254();
        assert!(ec.fits(&ec.base));
    }

    #[test]
    fn test_curve_bn254_math_1_plus_1_eq_2() {
        let ec = curve_bn254();

        let p = &ec.base;
        let q = ec.add(p, p);
        let k = &Int::from(2);
        let z: Point = ec.mul(p, k);
        assert_eq!(z, q);
    }

    #[test]
    fn test_curve_bn254_math_3_plus_4_eq_7() {
        let ec = curve_bn254();
        let p = &ec.base;
        let a = ec.mul(p, &Int::from(3));
        let b = ec.mul(p, &Int::from(4));
        let c = ec.add(&a, &b);
        let d = ec.mul(p, &Int::from(7));
        assert_eq!(c, d);
    }

    #[test]
    fn test_curve_bn254_math_30_plus_12_eq_42() {
        let ec = curve_bn254();
        let p = &ec.base;
        let a = ec.mul(p, &Int::from(30));
        let b = ec.mul(p, &Int::from(12));
        let c = ec.add(&a, &b);
        let d = ec.mul(p, &Int::from(42));
        assert_eq!(c, d);
    }

    #[test]
    fn test_curve_bn254_math_sqrt() {
        let ec = curve_bn254();

        for k in 1..101 {
            let p = ec.mul(&ec.base, &Int::from(k));
            if let Some((a, b)) = ec.apply(&p.x) {
                assert!(a == p || b == p, "failed: k={k}");
            } else {
                assert!(false, "failed: k={k}");
            }
        }
    }
}
