use ecc::hash::hash;

fn main() {
    let a = hash(b"x");
    let b = hash(b"Bob");

    let ec = ecc::elliptic::curves::curve_256();
    let pa = ec.mul(&ec.base, &a);
    let pb = ec.mul(&ec.base, &b);

    let ab = ec.mul(&pb, &a);
    let ba = ec.mul(&pa, &b);
    assert_eq!(ab, ba);

    println!("secret A: {a:x}");
    println!("public A: {:x}", pa.x);
    println!("---");
    println!("secret B: {b:x}");
    println!("public B: {:x}", pb.x);
    println!("---");
    println!("shared A: {:x}", ab.x);
    println!("shared B: {:x}", ba.x);
}
