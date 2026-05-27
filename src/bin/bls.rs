//! Chapter 18 demo: BLS pairing-based signatures on the tiny p=43 curve.
//!
//! Demonstrates:
//!   - key generation
//!   - sign / verify round trip
//!   - rejection of tampered messages
//!   - signature aggregation on a shared message
//!
//! Port of the educational parts of `signatures_11.c` (without the
//! curve_11_parameters.bin loader; we synthesize a tiny equivalent system).

use ecc::bls::tiny_system;

fn main() {
    let sys = tiny_system();
    println!("BLS demo — tiny pairing-friendly setup");
    println!("  base prime: {}", sys.e.modulus);
    println!("  torsion r : {}", sys.tor);
    println!("  cobse     : {}", sys.cobse);
    println!("  G1        : ({}, {})", sys.g1.x, sys.g1.y);
    println!("  G2 (deg-1 coords elided for brevity)");

    println!("\n=== single signer ===");
    let (sk, pk) = sys.keygen();
    println!("sk = {sk}");
    let msg = b"hello, pairings";
    let sig = sys.sign(&sk, msg).expect("sign");
    println!("sig = ({}, {})", sig.x, sig.y);

    if sys.verify(&pk, msg, &sig) {
        println!("verify(pk, msg, sig)         => OK");
    } else {
        println!("verify(pk, msg, sig)         => FAIL (bug!)");
    }

    let tampered = b"hello, tampered";
    if !sys.verify(&pk, tampered, &sig) {
        println!("verify(pk, tampered, sig)    => REJECTED (expected)");
    } else {
        println!("verify(pk, tampered, sig)    => OK (bug — should reject!)");
    }

    println!("\n=== aggregate ===");
    let signers: Vec<_> = (0..3).map(|_| sys.keygen()).collect();
    let msg = b"common payload";
    let sigs: Vec<_> = signers
        .iter()
        .map(|(sk, _)| sys.sign(sk, msg).expect("sign"))
        .collect();
    let pks: Vec<_> = signers.iter().map(|(_, pk)| pk.clone()).collect();
    let agg_sig = sys.aggregate_sigs(&sigs);
    let agg_pk = sys.aggregate_pks(&pks);
    println!(
        "{} signers; aggregate sig = ({}, {})",
        signers.len(),
        agg_sig.x,
        agg_sig.y
    );
    if sys.verify(&agg_pk, msg, &agg_sig) {
        println!("verify(agg_pk, msg, agg_sig) => OK");
    } else {
        println!("verify(agg_pk, msg, agg_sig) => FAIL (bug!)");
    }
}
