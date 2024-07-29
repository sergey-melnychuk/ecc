use ecc::hash::hash;

fn main() {
    let a = hash(b"Alice");
    let b = hash(b"Bob");

    let ec = ecc::elliptic::curves::curve_256();
    let a = ec.mul(&ec.base, &a);
    let b = ec.mul(&ec.base, &b);

    let ab = ec.add(&b, &a);
    let ba = ec.add(&a, &b);

    println!("A+B: {:x}", ab.x);
    println!("B+A: {:x}", ba.x);
    println!(" ok: {}", ab == ba);
}
