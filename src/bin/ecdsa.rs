use ecc::{
    elliptic::{Curve, Point},
    hash::{hash, hash_n},
    modulus::Modulus,
    Int,
};

fn main() {
    let sk = "The quick brown fox jumps over the lazy dog";
    let sk = hash(sk.as_bytes());

    let ec: Curve = ecc::elliptic::curves::curve_256();
    let pk = ec.mul(&ec.base, &sk);

    println!("SK: {sk:x}");
    println!("PK: {:x}", pk.x);
    println!("---");

    let m = b"Attack at dawn!";
    let hash = hash(m);

    let (r, s) = sign(&ec, &hash, &sk);
    println!("---");
    let ok = check(&ec, &hash, &pk, &r, &s);
    println!("ok: {ok}");
}

fn sign(ec: &Curve, hash: &Int, sk: &Int) -> (Int, Int) {
    // r = p * k
    // s = (e + rx * sk) / k

    let h = hash.clone().modulo(&ec.order);

    let k = hash_n(&[
        sk.to_string_radix(16).as_bytes(),
        hash.to_string_radix(16).as_bytes(),
    ]);
    let k = k.modulo(&ec.order);

    let rx = ec.mul(&ec.base, &k).x.modulo(&ec.order);

    let m = &Modulus::new(&ec.order);
    let s = m.mul(
        &m.inv(&k).expect("k^(-1)"),
        &m.add(&h, &m.mul(&rx, sk)),
    );

    println!(" h: {h:x}\nrx: {rx:x}\n s: {s:x}");
    (rx, s)
}

fn check(
    ec: &Curve,
    hash: &Int,
    pk: &Point,
    rx: &Int,
    s: &Int,
) -> bool {
    // p * e / s + pk * rx / s
    // p * e / s + (sk * p) * rx / s
    // (e + rx * sk) / s * p
    // (e + rx * sk) / (e + rx * sk) * k * p
    // k * p

    let h = hash.clone().modulo(&ec.order);

    let m = &Modulus::new(&ec.order);
    let s = m.inv(s).expect("s^(-1)");
    let a = ec.mul(&ec.base, &m.mul(&h, &s));
    let b = ec.mul(pk, &m.mul(rx, &s));
    let r = ec.add(&a, &b);

    println!(" h: {h:x}\nrx: {:x}", r.x);
    &r.x == rx
}
