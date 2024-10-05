pub mod elliptic;
pub mod hash;
pub mod modulus;
pub mod polynomial;

pub use rug::Integer as Int;

pub fn hex(hex: &str) -> Int {
    Int::from_str_radix(hex, 16).expect("hex int")
}

pub fn dec(dec: &str) -> Int {
    Int::from_str_radix(dec, 10).expect("dec int")
}
