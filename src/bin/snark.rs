//! Full SNARK workflow: QAP → CRS → Prove → Verify
//!
//! Runs the complete SNARK pipeline in one binary.

use ecc::snark::{qap::create_example_qap, Crs, Prover, SigSystem, Verifier};
use ecc::Int;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let base_path = "aux/drmike8888-Elliptic-curve-pairings/Build_all";

    // Step 1: Load system parameters
    println!("=== Step 1: Load System Parameters ===");
    let sys = SigSystem::load(&format!("{}/curve_11_parameters.bin", base_path))?;
    println!("  Prime: {} bits", sys.prime.significant_bits());
    println!("  Torsion: {} bits", sys.tor.significant_bits());
    println!("  Extension degree: {}", sys.irrd.deg);

    // Step 2: Create QAP
    println!("\n=== Step 2: Create QAP ===");
    let qap = create_example_qap(&sys.tor);
    println!("  Gates (n): {}", qap.n);
    println!("  Wires (m): {}", qap.m);
    println!("  Public inputs (l): {}", qap.l);

    // Step 3: Generate CRS (trusted setup)
    println!("\n=== Step 3: Generate CRS (Trusted Setup) ===");
    let crs = Crs::generate(&qap, &sys);
    println!("  Generated {} G1 points", crs.z_g.len() + crs.theta_g.len() + crs.zt_g.len() + 3);
    println!("  Generated {} G2 points", crs.z_h.len() + 3);

    // Step 4: Create proof
    println!("\n=== Step 4: Generate Proof ===");
    // Example values from snark_proof.c
    let medicine = Int::from(2036);
    let dose = Int::from(1700000);
    let patient = Int::from(49); // witness (private)

    println!("  Statement (public):");
    println!("    medicine = {}", medicine);
    println!("    dose = {}", dose);
    println!("  Witness (private):");
    println!("    patient = {}", patient);

    let prover = Prover::new(sys.clone(), qap, crs.clone());
    let record = prover.prove(&medicine, &dose, &patient);

    println!("  Proof A: ({}, {})", 
        &record.proof.a.x.to_string()[..20.min(record.proof.a.x.to_string().len())], 
        "...");
    println!("  Proof C: ({}, {})", 
        &record.proof.c.x.to_string()[..20.min(record.proof.c.x.to_string().len())], 
        "...");

    // Step 5: Verify proof
    println!("\n=== Step 5: Verify Proof ===");
    let verifier = Verifier::new(sys, crs);
    let valid = verifier.verify(&record);

    println!("\n{}", "=".repeat(50));
    if valid {
        println!("✓ VERIFICATION PASSED");
        println!("  The proof is valid!");
    } else {
        println!("✗ VERIFICATION FAILED");
        println!("  The proof is invalid.");
    }
    println!("{}", "=".repeat(50));

    Ok(())
}
