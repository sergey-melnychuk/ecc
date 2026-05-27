//! Chapter 14 demo: sweep for a pairing-friendly curve at a given embedding
//! degree, then build the curve via CM (Hilbert class polynomial → j → a4/a6).
//! Closes the loop from "book code" to "reproducible curve parameters".
//!
//! Usage:
//!     cargo run --release --bin search -- <k> <lg2_r_max>
//!     cargo run --release --bin search -- 5 40
//!
//! Output is printed in the same labeled-text format the book uses for its
//! `Curve_*_params.dat` files, so you can dump it straight into a `.dat`.

use std::env;

use ecc::cm_curve::{find_cm_discriminant, get_curve, HilbertTable};
use ecc::curve_search::sweep;
use ecc::Int;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<String> = env::args().collect();
    let k: u32 = args
        .get(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(5);
    let lg2r: u32 = args
        .get(2)
        .and_then(|s| s.parse().ok())
        .unwrap_or(40);

    println!(
        "# Searching for pairing-friendly curves: k = {k}, log2(r) ≤ {lg2r}\n"
    );

    let candidates = sweep(k, lg2r);
    if candidates.is_empty() {
        eprintln!("No candidates found. Try a larger lg2_r or a different k.");
        return Ok(());
    }

    println!("Found {} candidate(s). Top by rho:\n", candidates.len());
    for (i, c) in candidates.iter().take(5).enumerate() {
        println!(
            "  [{i}] alpha = {alpha}, x = {x}, rho = {rho:.4}, |r| = {rsz}b, |q| = {qsz}b",
            alpha = c.alpha,
            x = c.x,
            rho = c.rho,
            rsz = c.r.significant_bits(),
            qsz = c.q.significant_bits(),
        );
    }

    let hilbert = HilbertTable::load(
        "aux/drmike8888-Elliptic-curve-pairings/Build_all/Hilbert_Polynomials.list",
    )?;

    // Try each candidate until one admits a CM construction we can complete.
    println!("\n# Building a curve from the first viable candidate...\n");
    for (idx, c) in candidates.iter().enumerate() {
        let Some(d) = find_cm_discriminant(&hilbert, &c.q, &c.t) else {
            println!("  [{idx}] no CM discriminant in table");
            continue;
        };
        println!("  [{idx}] CM discriminant: D = -{d}");

        let Some(curve) = get_curve(&hilbert, d, &c.q, &c.t) else {
            println!("  [{idx}] get_curve failed");
            continue;
        };

        // Print in the labeled-text format the book uses for Curve_*_params.dat.
        let card_e: Int = c.q.clone() + 1 - &c.t;
        let cofactor: Int = card_e.clone() / &c.r;
        println!("\n# Curve_k{k}_alpha{alpha}_x{x}_D{d}.dat\n",
            alpha = c.alpha,
            x = c.x,
        );
        println!("prime");
        println!("{}", curve.modulus.to_string_radix(16));
        println!("order");
        println!("{}", curve.order.to_string_radix(16));
        println!("cofactor");
        println!("{cofactor}");
        println!("curve(a4   a6)");
        println!("{}", curve.a.to_string_radix(16));
        println!("{}", curve.b.to_string_radix(16));
        println!("basepoint(x   y)");
        println!("{}", curve.base.x.to_string_radix(16));
        println!("{}", curve.base.y.to_string_radix(16));
        println!();
        println!("# Pairing metadata:");
        println!("#   embedding degree k = {k}");
        println!("#   torsion           r = {}", c.r);
        println!("#   trace             t = {}", c.t);
        println!("#   rho               = {:.4}", c.rho);
        return Ok(());
    }

    eprintln!("\nNo candidate could be turned into a curve. Try larger lg2_r.");
    Ok(())
}
