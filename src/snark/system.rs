//! Signature system parameters for pairing-based cryptography.
//!
//! Ported from the SIG_SYSTEM structure in pairing.h and signature.c

use rug::Integer as Int;
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Read, Write};

use crate::elliptic::{Curve, Point};
use crate::modulus::Modulus;

use super::field_ext::{Poly, PolyCurve, PolyPoint, MAX_DEGREE};

/// Signature system parameters for pairing-based cryptography.
///
/// Contains all necessary parameters for a pairing-friendly curve setup.
#[derive(Clone)]
pub struct SigSystem {
    /// Field prime p
    pub prime: Int,
    /// Base curve E: y^2 = x^3 + a4*x + a6 (mod p)
    pub e: Curve,
    /// Cardinality of E(F_p)
    pub card_e: Int,
    /// Torsion subgroup order (divides card_e)
    pub tor: Int,
    /// Cofactor for base curve: card_e / tor
    pub cobse: Int,
    /// Generator point on G1
    pub g1: Point,
    /// Irreducible polynomial for extension field
    pub irrd: Poly,
    /// Cardinality of E(F_p^k)
    pub card_ex: Int,
    /// Cofactor for extension: card_ex / tor
    pub coxtd: Int,
    /// Generator point on G2 (extension field)
    pub g2: PolyPoint,
    /// Extended curve Ex: same equation but over F_p^k
    pub ex: PolyCurve,
}

impl SigSystem {
    /// Create an empty system (for later initialization)
    pub fn new() -> Self {
        Self {
            prime: Int::ZERO,
            e: Curve::new(
                Int::ZERO,
                Int::ZERO,
                Point::inf(),
                Int::ZERO,
                Int::ZERO,
            ),
            card_e: Int::ZERO,
            tor: Int::ZERO,
            cobse: Int::ZERO,
            g1: Point::inf(),
            irrd: Poly::new(),
            card_ex: Int::ZERO,
            coxtd: Int::ZERO,
            g2: PolyPoint::new(),
            ex: PolyCurve::new(),
        }
    }

    /// Get the modulus for the base field
    pub fn modulus(&self) -> Modulus {
        Modulus::new(&self.prime)
    }

    /// Get the modulus for the torsion group
    pub fn torsion_modulus(&self) -> Modulus {
        Modulus::new(&self.tor)
    }

    /// Load system parameters from a binary file.
    ///
    /// File format matches the C implementation's get_system().
    pub fn load(filename: &str) -> io::Result<Self> {
        let file = File::open(filename)?;
        let mut reader = BufReader::new(file);
        let mut sys = Self::new();

        // Read prime
        sys.prime = read_mpz_raw(&mut reader)?;

        // Read curve coefficients (a4, a6)
        let a4 = read_mpz_raw(&mut reader)?;
        let a6 = read_mpz_raw(&mut reader)?;

        // Read cardinality, torsion, cofactor
        sys.card_e = read_mpz_raw(&mut reader)?;
        sys.tor = read_mpz_raw(&mut reader)?;
        sys.cobse = read_mpz_raw(&mut reader)?;

        // Read G1 point
        let g1_x = read_mpz_raw(&mut reader)?;
        let g1_y = read_mpz_raw(&mut reader)?;
        sys.g1 = Point::new(g1_x, g1_y);

        // Set up base curve
        sys.e = Curve::new(
            sys.prime.clone(),
            sys.card_e.clone(),
            sys.g1.clone(),
            a4.clone(),
            a6.clone(),
        );

        // Read irreducible polynomial
        sys.irrd = read_poly(&mut reader)?;

        // Read extension field parameters
        sys.card_ex = read_mpz_raw(&mut reader)?;
        sys.coxtd = read_mpz_raw(&mut reader)?;

        // Read G2 point
        sys.g2 = read_poly_point(&mut reader)?;

        // Set up extension curve (same coefficients, lifted to extension)
        sys.ex.a4 = Poly::constant(a4);
        sys.ex.a6 = Poly::constant(a6);

        Ok(sys)
    }

    /// Save system parameters to a binary file
    pub fn save(&self, filename: &str) -> io::Result<()> {
        let file = File::create(filename)?;
        let mut writer = BufWriter::new(file);

        write_mpz_raw(&mut writer, &self.prime)?;
        write_mpz_raw(&mut writer, &self.e.a)?;
        write_mpz_raw(&mut writer, &self.e.b)?;
        write_mpz_raw(&mut writer, &self.card_e)?;
        write_mpz_raw(&mut writer, &self.tor)?;
        write_mpz_raw(&mut writer, &self.cobse)?;
        write_mpz_raw(&mut writer, &self.g1.x)?;
        write_mpz_raw(&mut writer, &self.g1.y)?;

        write_poly(&mut writer, &self.irrd)?;

        write_mpz_raw(&mut writer, &self.card_ex)?;
        write_mpz_raw(&mut writer, &self.coxtd)?;

        write_poly_point(&mut writer, &self.g2)?;

        Ok(())
    }

    /// Convert a G1 point to a G2 point (lift to extension field).
    ///
    /// Maps to C function: tog2()
    pub fn to_g2(&self, g1: &Point) -> PolyPoint {
        let mut g2 = PolyPoint::new();
        g2.x.deg = 0;
        g2.x.coef[0] = g1.x.clone();
        g2.y.deg = 0;
        g2.y.coef[0] = g1.y.clone();
        g2
    }
}

impl Default for SigSystem {
    fn default() -> Self {
        Self::new()
    }
}

// Binary I/O helpers matching GMP's mpz_inp_raw / mpz_out_raw format

/// Read an integer in GMP raw format
pub fn read_mpz_raw<R: Read>(reader: &mut R) -> io::Result<Int> {
    // GMP raw format: 4-byte size (big-endian), then |size| bytes of data
    let mut size_buf = [0u8; 4];
    reader.read_exact(&mut size_buf)?;
    let size = i32::from_be_bytes(size_buf);

    let abs_size = size.unsigned_abs() as usize;
    if abs_size == 0 {
        return Ok(Int::ZERO);
    }

    let mut data = vec![0u8; abs_size];
    reader.read_exact(&mut data)?;

    // GMP stores in big-endian format
    let mut int = Int::from_digits(&data, rug::integer::Order::Msf);

    if size < 0 {
        int = -int;
    }

    Ok(int)
}

/// Write an integer in GMP raw format
pub fn write_mpz_raw<W: Write>(
    writer: &mut W,
    int: &Int,
) -> io::Result<()> {
    if int.is_zero() {
        writer.write_all(&[0u8; 4])?;
        return Ok(());
    }

    let is_neg = int < &Int::ZERO;
    let abs_int = int.clone().abs();
    let bytes = abs_int.to_digits::<u8>(rug::integer::Order::Msf);

    let size = if is_neg {
        -(bytes.len() as i32)
    } else {
        bytes.len() as i32
    };

    writer.write_all(&size.to_be_bytes())?;
    writer.write_all(&bytes)?;

    Ok(())
}

/// Read a polynomial
pub fn read_poly<R: Read>(reader: &mut R) -> io::Result<Poly> {
    // Read degree - C uses fwrite(&k, sizeof(long), 1, f) where k is int
    // On 64-bit systems, sizeof(long) = 8, so we read 8 bytes but only use lower 4
    let mut deg_buf = [0u8; 8];
    reader.read_exact(&mut deg_buf)?;
    // The actual degree is in the first 4 bytes (little-endian on most systems)
    let deg = i32::from_le_bytes([
        deg_buf[0], deg_buf[1], deg_buf[2], deg_buf[3],
    ]) as usize;

    let mut poly = Poly::new();
    poly.deg = deg;

    // Read deg+1 coefficients
    for i in 0..=deg.min(MAX_DEGREE - 1) {
        poly.coef[i] = read_mpz_raw(reader)?;
    }

    Ok(poly)
}

/// Write a polynomial
pub fn write_poly<W: Write>(
    writer: &mut W,
    poly: &Poly,
) -> io::Result<()> {
    // Match the binary format: 8-byte field with degree in lower 4 bytes (LE)
    let mut deg_buf = [0u8; 8];
    let deg_bytes = (poly.deg as i32).to_le_bytes();
    deg_buf[0..4].copy_from_slice(&deg_bytes);
    writer.write_all(&deg_buf)?;

    for i in 0..=poly.deg {
        write_mpz_raw(writer, &poly.coef[i])?;
    }

    Ok(())
}

/// Read a point on extension curve
pub fn read_poly_point<R: Read>(
    reader: &mut R,
) -> io::Result<PolyPoint> {
    let x = read_poly(reader)?;
    let y = read_poly(reader)?;
    Ok(PolyPoint::from_coords(x, y))
}

/// Write a point on extension curve
pub fn write_poly_point<W: Write>(
    writer: &mut W,
    point: &PolyPoint,
) -> io::Result<()> {
    write_poly(writer, &point.x)?;
    write_poly(writer, &point.y)?;
    Ok(())
}

/// Read a base curve point
pub fn read_point<R: Read>(reader: &mut R) -> io::Result<Point> {
    let x = read_mpz_raw(reader)?;
    let y = read_mpz_raw(reader)?;
    Ok(Point::new(x, y))
}

/// Write a base curve point
pub fn write_point<W: Write>(
    writer: &mut W,
    point: &Point,
) -> io::Result<()> {
    write_mpz_raw(writer, &point.x)?;
    write_mpz_raw(writer, &point.y)?;
    Ok(())
}

// Re-export I/O helpers for other modules
pub use read_mpz_raw as read_int;
pub use read_poly as read_polynomial;
pub use read_poly_point as read_ext_point;
pub use write_mpz_raw as write_int;
pub use write_poly as write_polynomial;
pub use write_poly_point as write_ext_point;
