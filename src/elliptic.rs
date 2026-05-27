use std::fs;
use std::io;
use std::path::Path;

use rug::Integer as Int;

use crate::modulus::Modulus;

// y^2 = x^3 + ax + b (mod n)
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Curve {
    pub modulus: Int,
    pub order: Int,
    pub base: Point,
    pub a: Int,
    pub b: Int,
}

impl Curve {
    pub fn new(
        modulus: Int,
        order: Int,
        base: Point,
        a: Int,
        b: Int,
    ) -> Self {
        Self {
            modulus,
            order,
            base,
            a,
            b,
        }
    }

    /// Load a curve from the book's `Curve_*_params.dat` text format
    /// (Chapter 6). The format is six labeled sections of hex / decimal
    /// integers:
    ///
    /// ```text
    /// prime
    /// <hex>
    /// order
    /// <hex>
    /// cofactor
    /// <decimal>
    /// curve(a4   a6)
    /// <hex>
    /// <hex>
    /// basepoint(x   y)
    /// <hex>
    /// <hex>
    /// ```
    pub fn from_dat<P: AsRef<Path>>(path: P) -> io::Result<Self> {
        let text = fs::read_to_string(path)?;
        // Strip labels — keep only lines that are entirely [0-9a-fA-F]
        // (which covers both hex values and the decimal cofactor).
        let nums: Vec<&str> = text
            .lines()
            .map(str::trim)
            .filter(|l| !l.is_empty())
            .filter(|l| l.chars().all(|c| c.is_ascii_hexdigit()))
            .collect();
        if nums.len() != 7 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "expected 7 numeric lines, got {}",
                    nums.len()
                ),
            ));
        }
        let parse_hex = |s: &str| {
            Int::from_str_radix(s, 16).map_err(|e| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    e.to_string(),
                )
            })
        };
        let parse_dec = |s: &str| {
            Int::from_str_radix(s, 10).map_err(|e| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    e.to_string(),
                )
            })
        };
        let modulus = parse_hex(nums[0])?;
        let order = parse_hex(nums[1])?;
        let _cofactor = parse_dec(nums[2])?;
        let a = parse_hex(nums[3])?;
        let b = parse_hex(nums[4])?;
        let x = parse_hex(nums[5])?;
        let y = parse_hex(nums[6])?;
        let base = Point::new(x, y);
        let curve = Curve::new(modulus, order, base.clone(), a, b);
        if !curve.fits(&base) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "loaded base point is not on the curve",
            ));
        }
        Ok(curve)
    }

    pub fn add(&self, p: &Point, q: &Point) -> Point {
        if p.is_inf() {
            return q.clone();
        }
        if q.is_inf() {
            return p.clone();
        }

        let m = Modulus::new(&self.modulus);

        // P + (-P) = O (point at infinity).
        // Covers both p == -q (different points, y1 = -y2) and doubling
        // a 2-torsion point (p == q with y = 0).
        if p.x == q.x && m.add(&p.y, &q.y).is_zero() {
            return Point::inf();
        }

        let num = m.add(
            &m.add(&m.mul(&p.x, &p.x), &m.mul(&p.x, &q.x)),
            &m.add(&m.mul(&q.x, &q.x), &self.a),
        );
        let den = m.add(&p.y, &q.y);

        let lambda = m.div(&num, &den).expect("lambda");
        let x = m.sub(&m.sub(&m.mul(&lambda, &lambda), &p.x), &q.x);
        let y = m.sub(&m.mul(&lambda, &m.sub(&p.x, &x)), &p.y);

        let ret = Point::new(x, y);
        assert!(
            self.fits(&ret),
            "add result point must fit the curve"
        );
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
        assert!(
            self.fits(&acc),
            "mul result point must fit the curve"
        );
        acc
    }

    pub fn fits(&self, p: &Point) -> bool {
        if p.is_inf() {
            return true;
        }
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
    inf: bool,
}

impl Point {
    pub fn new(x: Int, y: Int) -> Self {
        Self { x, y, inf: false }
    }

    pub fn inf() -> Self {
        Self {
            x: Int::ZERO,
            y: Int::ZERO,
            inf: true,
        }
    }

    pub fn is_inf(&self) -> bool {
        self.inf
    }
}

pub mod curves {
    use super::*;
    use crate::{dec, hex};

    // Parameters from Chapter 6 — see aux/.../Build_all/Curve_*_params.dat.
    // Each curve is y² = x³ + a4·x + a6 over F_p, with a generator of the
    // (prime) order subgroup baked in.

    // Curve_160_params.dat
    pub fn curve_160() -> Curve {
        let modulus =
            hex("ac000000000000000000000000000000000000001");
        let order = hex("ac0000000000000000006543ba11adf8eb6345c77");
        let a = hex("1");
        let b = hex("782e");

        let x = hex("1680bbdc87647f3c382902d2f58d2754b39bca877");
        let y = hex("a08957b09764ae59da8fb3058efef9c428e497268");
        let base = Point::new(x, y);

        let curve = Curve::new(modulus, order, base.clone(), a, b);
        assert!(curve.fits(&base));
        curve
    }

    // Curve_256_params.dat
    pub fn curve_256() -> Curve {
        let modulus = hex("2b000000000000000000000000000000000000000000000000000000000000001");
        let order = hex("2b0000000000000000000000000000002e7f521c85bba055a6e2161b956a47f69");
        let a = hex("1");
        let b = hex("a87");

        let x = hex("2310115d283e49377820195c8e67781b6f112a625b14b747fa4cc13d06eba0919");
        let y = hex("51277aeb91946f0cb83053a10f67c5a9ef00a4f0cf2466b3bedf4fdcd774b574");
        let base = Point::new(x, y);

        let curve = Curve::new(modulus, order, base.clone(), a, b);
        assert!(curve.fits(&base));
        curve
    }

    // Curve_384_params.dat
    pub fn curve_384() -> Curve {
        let modulus = hex("2e00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001");
        let order = hex("2e00000000000000000000000000000000000000000000002275cc5f2f7fcc15352a2c993900a851b3a75365a9ac54733");
        let a = hex("1");
        let b = hex("310");

        let x = hex("23c0d9fcfaa3dc18b1eff7e89bf7678636580d17dd84a873b14b9c0e1680bbdc87647f3c382902d2f58d2754b39bca874");
        let y = hex("28d7205f1be0a725d2aa7c3386f2e0b0ea7c558ca19f9770cdc72f91a1cbc262687810d4c5bd536818ccfa49aae2ed0cc");
        let base = Point::new(x, y);

        let curve = Curve::new(modulus, order, base.clone(), a, b);
        assert!(curve.fits(&base));
        curve
    }

    // Curve_512_params.dat
    pub fn curve_512() -> Curve {
        let modulus = hex("e20000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001");
        let order = hex("e2000000000000000000000000000000000000000000000000000000000000007788830d091dc57e3af7d7bbd15386ee9414602d88d1e6489cd056336922bbf4d");
        let a = hex("1");
        let b = hex("41");

        let x = hex("518f204fe6846aeb6f58174d57a3372363c0d9fcfaa3dc18b1eff7e89bf7678636580d17dd84a873b14b9c0e1680bbdc87647f3c382902d2f58d2754b39bca875");
        let y = hex("c9fe7223aca476cde61f206be285898475f1dcbaefeda90057d3b8bae5146f3016ebf2139daa73f39417193e8609a4229cd4c58389e4b9095fafcd68362b310fe");
        let base = Point::new(x, y);

        let curve = Curve::new(modulus, order, base.clone(), a, b);
        assert!(curve.fits(&base));
        curve
    }

    // https://hackmd.io/@jpw/bn254
    pub fn curve_bn254() -> Curve {
        let modulus = dec("21888242871839275222246405745257275088696311157297823662689037894645226208583");
        let order = dec("21888242871839275222246405745257275088548364400416034343698204186575808495617");
        let a = dec("0");
        let b = dec("3");

        let x = dec("1");
        let y = dec("2");
        let base = Point::new(x, y);

        let curve = Curve::new(modulus, order, base.clone(), a, b);
        assert!(curve.fits(&base));
        curve
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::elliptic::curves::{
        curve_160, curve_256, curve_384, curve_512, curve_bn254,
    };

    #[test]
    fn test_curve_160() {
        let ec = curve_160();
        assert!(ec.fits(&ec.base));
        // [order] · G must be the identity (G generates the prime-order subgroup).
        assert!(ec.mul(&ec.base, &ec.order).is_inf());
    }

    #[test]
    fn test_curve_256() {
        let ec = curve_256();
        assert!(ec.fits(&ec.base));
        assert!(ec.mul(&ec.base, &ec.order).is_inf());
    }

    #[test]
    fn test_curve_384() {
        let ec = curve_384();
        assert!(ec.fits(&ec.base));
        assert!(ec.mul(&ec.base, &ec.order).is_inf());
    }

    #[test]
    fn test_curve_512() {
        let ec = curve_512();
        assert!(ec.fits(&ec.base));
        assert!(ec.mul(&ec.base, &ec.order).is_inf());
    }

    /// Every hard-coded curve must match the corresponding `.dat` file
    /// shipped in `aux/`. If the .dat ever changes upstream we want to know.
    #[test]
    fn test_from_dat_matches_hardcoded() {
        let aux = "aux/drmike8888-Elliptic-curve-pairings/Build_all";
        for (bits, expected) in [
            (160, curve_160()),
            (256, curve_256()),
            (384, curve_384()),
            (512, curve_512()),
        ] {
            let path = format!("{aux}/Curve_{bits}_params.dat");
            let loaded = Curve::from_dat(&path)
                .unwrap_or_else(|e| panic!("loading {path}: {e}"));
            assert_eq!(
                loaded.modulus, expected.modulus,
                "{bits}: modulus"
            );
            assert_eq!(loaded.order, expected.order, "{bits}: order");
            assert_eq!(loaded.a, expected.a, "{bits}: a4");
            assert_eq!(loaded.b, expected.b, "{bits}: a6");
            assert_eq!(loaded.base, expected.base, "{bits}: base");
        }
    }

    #[test]
    fn test_from_dat_rejects_off_curve_base() {
        // Build a .dat with a base point that isn't on the curve.
        let tmp = std::env::temp_dir().join("ecc_bad_curve.dat");
        std::fs::write(
            &tmp,
            "prime\n2b000000000000000000000000000000000000000000000000000000000000001\n\
             order\n1\ncofactor\n1\ncurve(a4 a6)\n1\na87\n\
             basepoint(x y)\n1\n1\n",
        )
        .unwrap();
        let err = Curve::from_dat(&tmp).unwrap_err();
        assert!(err.to_string().contains("not on the curve"));
        let _ = std::fs::remove_file(tmp);
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
    fn test_point_inf_is_distinct_from_origin() {
        // Infinity must not compare equal to an affine (0, 0) point —
        // otherwise curves with b == 0 would alias.
        let inf = Point::inf();
        let origin = Point::new(Int::ZERO, Int::ZERO);
        assert!(inf.is_inf());
        assert!(!origin.is_inf());
        assert_ne!(inf, origin);
    }

    #[test]
    fn test_curve_bn254_p_plus_neg_p_is_inf() {
        let ec = curve_bn254();
        let p = ec.base.clone();
        let modulus = Modulus::new(&ec.modulus);
        let neg_p = Point::new(p.x.clone(), modulus.neg(&p.y));
        let sum = ec.add(&p, &neg_p);
        assert!(
            sum.is_inf(),
            "P + (-P) must be the point at infinity"
        );
    }

    #[test]
    fn test_curve_bn254_inf_is_identity() {
        let ec = curve_bn254();
        let p = ec.base.clone();
        let inf = Point::inf();
        assert_eq!(ec.add(&p, &inf), p);
        assert_eq!(ec.add(&inf, &p), p);
    }

    #[test]
    fn test_curve_bn254_order_times_base_is_inf() {
        // [n]G = O for any G in the prime-order subgroup.
        let ec = curve_bn254();
        let result = ec.mul(&ec.base, &ec.order);
        assert!(result.is_inf());
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
