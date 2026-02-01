//! SNARK (Succinct Non-interactive Argument of Knowledge) implementation
//!
//! This module provides a Groth16-style SNARK implementation based on
//! Quadratic Arithmetic Programs (QAP) and bilinear pairings.
//!
//! Ported from drmike8888's C implementation.

pub mod crs;
pub mod field_ext;
pub mod lagrange;
pub mod pairing;
pub mod prover;
pub mod qap;
pub mod system;
pub mod verifier;

pub use crs::Crs;
pub use field_ext::{Poly, PolyCurve, PolyPoint};
pub use lagrange::LagrangeInterpolator;
pub use pairing::Pairing;
pub use prover::Prover;
pub use qap::Qap;
pub use system::SigSystem;
pub use verifier::Verifier;
