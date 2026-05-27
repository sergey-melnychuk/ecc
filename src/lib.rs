#![allow(clippy::needless_range_loop)]
#![allow(clippy::explicit_counter_loop)]
#![allow(clippy::manual_memcpy)]

pub mod bls;
pub mod cm_curve;
pub mod curve_search;
pub mod elliptic;
pub mod hash;
pub mod modulus;
pub mod pairing;
pub mod poly_elliptic;
pub mod polynomial;
pub mod snark;

pub use rug::Integer as Int;

pub fn hex(hex: &str) -> Int {
    Int::from_str_radix(hex, 16).expect("hex int")
}

pub fn dec(dec: &str) -> Int {
    Int::from_str_radix(dec, 10).expect("dec int")
}
