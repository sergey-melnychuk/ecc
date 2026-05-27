//! Chapter 17 demo: Tate pairing on the tiny p=43, k=2 curve.
//!
//! Enumerates every point, then runs four bilinearity checks at order 11
//! across G1×G1, G1×G2, G2×G1 and G2×G2.
//!
//! Port of `tate_6_bit_pairing.c`.

use rug::ops::Pow;

use ecc::modulus::Modulus;
use ecc::pairing::tate;
use ecc::poly_elliptic::{PolyCurve, PolyPoint};
use ecc::polynomial::Polynomial;
use ecc::Int;

fn lift(v: i64) -> Polynomial {
    let mut p = Polynomial::zeros(1);
    p.set(0, Int::from(v));
    p
}

fn print_poly(p: &Polynomial) {
    if p.degree() == 0 {
        print!("{}", p.get(0));
    } else {
        for i in (0..=p.degree()).rev() {
            let c = p.get(i);
            if i > 0 {
                print!("{c}*x^{i} + ");
            } else {
                print!("{c}");
            }
        }
    }
}

fn print_poly_point(label: &str, p: &PolyPoint) {
    print!("{label}: x = ");
    print_poly(&p.x);
    print!(", y = ");
    print_poly(&p.y);
    println!();
}

fn bump(x: &Polynomial, ec: &PolyCurve) -> Polynomial {
    let m = &ec.p;
    let n = ec.irrd.degree();
    let one = Int::from(1);
    let mut out = Polynomial::zeros(n);
    for i in 0..n {
        out.set(i, m.add(&x.get(i), &Int::ZERO));
    }
    let mut i = 0;
    while i < n {
        let v = m.add(&out.get(i), &one);
        let wrapped = v.is_zero();
        out.set(i, v);
        if !wrapped {
            break;
        }
        i += 1;
    }
    out.trim()
}

fn enumerate_points(
    ec: &PolyCurve,
    factors: &[Int],
) -> Vec<(PolyPoint, Int)> {
    let mut out = Vec::new();
    let total: Int = ec.p.n.clone().pow(ec.irrd.degree() as u32);
    let mut x = lift(0);
    let mut steps = Int::from(0);
    while steps < total {
        let Some((p1, p2)) = ec.embed(&x, 1) else {
            x = bump(&x, ec);
            steps += 1;
            continue;
        };
        let Some(ord) = ec.order(&p1, factors) else {
            x = bump(&p1.x, ec);
            steps += 1;
            continue;
        };
        out.push((p1.clone(), ord.clone()));
        out.push((p2, ord));
        x = bump(&p1.x, ec);
        steps += 1;
    }
    out
}

fn run_tate(
    ec: &PolyCurve,
    label: &str,
    p: &PolyPoint,
    q: &PolyPoint,
    t: &PolyPoint,
    s: &PolyPoint,
    tor: &Int,
) {
    println!("\n--- {label} ---");
    print_poly_point("P", p);
    print_poly_point("Q", q);
    print_poly_point("T", t);
    print_poly_point("S", s);

    let t_pq = tate(ec, p, q, s, tor);
    print!("  t(P, Q)         = ");
    print_poly(&t_pq);
    println!();

    let t_pt = tate(ec, p, t, s, tor);
    print!("  t(P, T)         = ");
    print_poly(&t_pt);
    println!();

    let prod = t_pq.mul(&t_pt, &ec.irrd, &ec.p);
    print!("  t(P,Q)·t(P,T)   = ");
    print_poly(&prod);
    println!();

    let tpq = ec.add(t, q);
    let t_p_tpq = tate(ec, p, &tpq, s, tor);
    print!("  t(P, T+Q)       = ");
    print_poly(&t_p_tpq);
    println!();

    if prod == t_p_tpq {
        println!("  BILINEARITY HOLDS");
    } else {
        println!("  BILINEARITY FAILED");
    }
}

fn main() {
    let m = Modulus::new(&Int::from(43));
    let irrd = Polynomial::find_irreducible(2, &m)
        .expect("irreducible deg-2");
    print!("Working in F_43[x] / <");
    print_poly(&irrd);
    println!(">");
    println!("Curve: y² = x³ + 23x + 42");

    let mut a4 = Polynomial::zeros(1);
    a4.set(0, Int::from(23));
    let mut a6 = Polynomial::zeros(1);
    a6.set(0, Int::from(42));
    let ec = PolyCurve::new(a4, a6, irrd, m);

    let factors: Vec<Int> = [3, 5, 11, 15, 33, 55, 165, 1815]
        .into_iter()
        .map(Int::from)
        .collect();
    println!("Enumerating all points on E(F_43²)...");
    let points = enumerate_points(&ec, &factors);
    println!("Found {} points.\n", points.len());

    let tor = Int::from(11);

    let g1: Vec<&PolyPoint> = points
        .iter()
        .filter(|(p, o)| o == &tor && p.g1g2() == 1)
        .map(|(p, _)| p)
        .collect();
    let g2: Vec<&PolyPoint> = points
        .iter()
        .filter(|(p, o)| o == &tor && p.g1g2() > 1)
        .map(|(p, _)| p)
        .collect();
    let aux = points
        .iter()
        .find(|(_, o)| o != &tor)
        .map(|(p, _)| p)
        .expect("aux of different order");

    println!(
        "Order-11 points: {} in G1, {} in G2",
        g1.len(),
        g2.len()
    );

    if g1.len() >= 3 {
        run_tate(&ec, "Tate G1 x G1", g1[0], g1[1], g1[2], aux, &tor);
    }
    if !g1.is_empty() && g2.len() >= 2 {
        run_tate(&ec, "Tate G1 x G2", g1[0], g2[0], g2[1], aux, &tor);
    }
    if !g1.is_empty() && g2.len() >= 2 {
        run_tate(&ec, "Tate G2 x G1", g2[0], g1[0], g1[1], aux, &tor);
    }
    if g2.len() >= 3 {
        run_tate(&ec, "Tate G2 x G2", g2[0], g2[1], g2[2], aux, &tor);
    }
}
