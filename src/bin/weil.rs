//! Chapter 16 demo: enumerate every point of E(F_{43^2}) for the tiny curve
//! y^2 = x^3 + 23x + 42, then verify Weil-pairing bilinearity on G1×G1,
//! G1×G2 and G2×G2 subgroups.
//!
//! Port of `weil_6_bit_pairing.c`.

use rug::ops::Pow;

use ecc::modulus::Modulus;
use ecc::pairing::weil;
use ecc::poly_elliptic::{PolyCurve, PolyPoint};
use ecc::polynomial::Polynomial;
use ecc::Int;

fn lift(v: i64) -> Polynomial {
    let mut p = Polynomial::zeros(1);
    p.set(0, Int::from(v));
    p
}

fn print_poly_point(label: &str, p: &PolyPoint) {
    print!("{label}: x = ");
    print_poly(&p.x);
    print!(", y = ");
    print_poly(&p.y);
    println!();
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

/// Enumerate every point on the curve by sweeping x through every value in
/// F_{p^k} (here 43² = 1849 values). For each x where f(x) is a QR we get
/// two points (P, -P). Returns the list along with each point's order.
fn enumerate_points(
    ec: &PolyCurve,
    factors: &[Int],
) -> Vec<(PolyPoint, Int)> {
    let mut out = Vec::new();
    let total = ec.p.n.clone().pow(ec.irrd.degree() as u32);
    let mut x = lift(0);
    let mut steps: Int = Int::from(0);
    while steps < total {
        let (p1, p2) = match ec.embed(&x, 1) {
            Some(pair) => pair,
            None => {
                x = bump(&x, ec);
                steps += 1;
                continue;
            }
        };
        let ord = ec
            .order(&p1, factors)
            .expect("order must be found in factor list");
        out.push((p1.clone(), ord.clone()));
        out.push((p2, ord));
        x = bump(&p1.x, ec);
        steps += 1;
    }
    out
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

fn main() {
    let m = Modulus::new(&Int::from(43));
    let irrd =
        Polynomial::find_irreducible(2, &m).expect("irreducible deg-2");
    println!("Working in F_43[x] / <");
    print_poly(&irrd);
    println!(">");
    println!("Curve: y^2 = x^3 + 23x + 42");

    let mut a4 = Polynomial::zeros(1);
    a4.set(0, Int::from(23));
    let mut a6 = Polynomial::zeros(1);
    a6.set(0, Int::from(42));
    let ec = PolyCurve::new(a4, a6, irrd, m);

    let factors: Vec<Int> = [3, 5, 11, 15, 33, 55, 165, 1815]
        .into_iter()
        .map(Int::from)
        .collect();

    println!("Enumerating all points on E(F_43^2)...");
    let points = enumerate_points(&ec, &factors);
    println!("Found {} points (expected 1815).", points.len());

    let tor = Int::from(11);

    // Pick representative points: P, Q, T of order 11 in G1, and S of
    // larger order so the Weil construction is well-defined.
    let g1_eleven: Vec<&PolyPoint> = points
        .iter()
        .filter(|(p, o)| o == &tor && p.g1g2() == 1)
        .map(|(p, _)| p)
        .collect();
    let g2_eleven: Vec<&PolyPoint> = points
        .iter()
        .filter(|(p, o)| o == &tor && p.g1g2() > 1)
        .map(|(p, _)| p)
        .collect();
    let aux: &PolyPoint = points
        .iter()
        .find(|(_, o)| o != &tor)
        .map(|(p, _)| p)
        .expect("aux point");

    println!("\nG1 x G1 test (order-11 in base subgroup) -----------");
    if g1_eleven.len() >= 3 {
        let p = g1_eleven[0];
        let q = g1_eleven[1];
        let t = g1_eleven[2];
        run_bilinearity_check(&ec, p, q, t, aux, &tor);
    } else {
        println!("Not enough G1 order-11 points to test ({} found)", g1_eleven.len());
    }

    println!("\nG1 x G2 test (mixed) --------------------------------");
    if !g1_eleven.is_empty() && g2_eleven.len() >= 2 {
        let p = g1_eleven[0];
        let q = g2_eleven[0];
        let t = g2_eleven[1];
        run_bilinearity_check(&ec, p, q, t, aux, &tor);
    } else {
        println!("Not enough mixed points to test");
    }

    println!("\nG2 x G2 test ---------------------------------------");
    if g2_eleven.len() >= 3 {
        let p = g2_eleven[0];
        let q = g2_eleven[1];
        let t = g2_eleven[2];
        run_bilinearity_check(&ec, p, q, t, aux, &tor);
    } else {
        println!("Not enough G2 order-11 points to test ({} found)", g2_eleven.len());
    }
}

fn run_bilinearity_check(
    ec: &PolyCurve,
    p: &PolyPoint,
    q: &PolyPoint,
    t: &PolyPoint,
    s: &PolyPoint,
    tor: &Int,
) {
    print_poly_point("P", p);
    print_poly_point("Q", q);
    print_poly_point("T", t);
    print_poly_point("S", s);

    let w_pq = weil(ec, p, q, s, tor);
    print!("  e(P, Q)    = ");
    print_poly(&w_pq);
    println!();

    let w_pt = weil(ec, p, t, s, tor);
    print!("  e(P, T)    = ");
    print_poly(&w_pt);
    println!();

    let tpq = ec.add(t, q);
    let w_p_tpq = weil(ec, p, &tpq, s, tor);
    print!("  e(P, T+Q)  = ");
    print_poly(&w_p_tpq);
    println!();

    let prod = w_pt.mul(&w_pq, &ec.irrd, &ec.p);
    print!("  e(P,T)*e(P,Q) = ");
    print_poly(&prod);
    println!();

    if prod == w_p_tpq {
        println!("  BILINEARITY HOLDS");
    } else {
        println!("  BILINEARITY FAILED");
    }

    // Raise the product to torsion power; must collapse to 1.
    let raised = w_p_tpq.pow(tor, &ec.irrd, &ec.p);
    print!("  e(P, T+Q)^{tor} = ");
    print_poly(&raised);
    println!();
}
