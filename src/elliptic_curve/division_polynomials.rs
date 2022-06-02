use super::super::eea::*;
use super::super::prelude::*;
use super::super::wrapper::*;
use super::super::poly::*;
use super::*;

type P<K: Ring> = WrappingRing<PolyRing<K>>;

fn duplication<K: Ring + PartialEq>(_poly_ring: &P<K>, a4: El<P<K>>, x1: El<P<K>>, y1_y: El<P<K>>, z1: El<P<K>>, f: El<P<K>>) -> (El<P<K>>, El<P<K>>, El<P<K>>) {
    let lambda_div_y_num = x1.pow(2) * 3 + z1.pow(2) * a4;
    let lambda_den = &y1_y * &z1 * 2;
    let x3_div_lambda_den = - &x1 * 2 * lambda_den.pow(2) + &z1 * lambda_div_y_num.pow(2) * &f;
    let z3 = lambda_den.pow(3) * z1;
    let y3_y = - y1_y * lambda_den.pow(3) + lambda_div_y_num * (x1 * lambda_den.pow(2) - &x3_div_lambda_den) * f;
    return (x3_div_lambda_den * lambda_den, y3_y, z3);
}

fn xy1_addition<K: Ring + PartialEq>(_poly_ring: &P<K>, x1: El<P<K>>, y1_y: El<P<K>>, z1: El<P<K>>, f: El<P<K>>) -> (El<P<K>>, El<P<K>>, El<P<K>>) {
    let x = x1.ring().from(x1.ring().wrapped_ring().unknown());
    let lambda_y_num = &y1_y - &f * &z1;
    let lambda_den = &x1 - &x * &z1;
    let x3_div_lambda_den = -&x1 * lambda_den.pow(2) - x * &z1 * lambda_den.pow(2) + &z1 * lambda_y_num.pow(2) / f;
    let z3 = lambda_den.pow(3) * &z1;
    let y3_y = - y1_y * lambda_den.pow(3) + lambda_y_num * (x1 * lambda_den.pow(2) - &x3_div_lambda_den);
    return (x3_div_lambda_den * lambda_den, y3_y, z3);
}

///
/// Computes two polynomials f and g with the property that
/// for each point `P = (x : y : z)` on the curve `E`, have that
/// the x-coordinate of `[n]P` is `f(x)/g(x)`, if `[n]P` is affine.
/// 
pub fn division_polynomials<K: Ring + PartialEq>(E: &EllipticCurve<WrappingRing<K>>, mut n: usize) 
    -> (El<WrappingRing<PolyRing<K>>>, El<WrappingRing<PolyRing<K>>>, El<WrappingRing<PolyRing<K>>>) 
{
    assert!(n >= 2);
    let poly_ring = PolyRing::adjoint(E.base_field().wrapped_ring().clone(), "x");
    let P = poly_ring.bind_ring_by_value();
    let x = P.from(P.wrapped_ring().unknown());
    let incl = embedding(E.base_field(), &P);
    let f = x.pow(3) + incl(E.a4().clone()) * &x + incl(E.a6().clone());

    let reduce = |mut x1: El<WrappingRing<PolyRing<K>>>, mut y1_y: El<WrappingRing<PolyRing<K>>>, mut z1: El<WrappingRing<PolyRing<K>>>| {
        let d1 = gcd(&x1.ring(), x1.clone(), z1.clone());
        let d = gcd(&y1_y.ring(), d1, y1_y.clone());
        x1 /= &d;
        y1_y /= &d;
        z1 /= &d;
        let nomalize_factor = x1.lc().unwrap().inv();
        x1 = x1.scaled(&nomalize_factor);
        y1_y = y1_y.scaled(&nomalize_factor);
        z1 = z1.scaled(&nomalize_factor);
        return (x1, y1_y, z1);
    };

    let (mut x1, mut y1_y, mut z1) = (x, f.clone(), P.one());
    for i in (0..(usize::BITS - n.leading_zeros() - 1)).rev() {
        (x1, y1_y, z1) = duplication(&P, incl(E.a4().clone()), x1, y1_y, z1, f.clone());
        (x1, y1_y, z1) = reduce(x1, y1_y, z1);
        if (n >> i) & 1 == 1 {
            (x1, y1_y, z1) = xy1_addition(&P, x1, y1_y, z1, f.clone());
        }
    }
    (x1, y1_y, z1) = reduce(x1, y1_y, z1);

    return (x1, y1_y, z1);
}

#[cfg(test)]
use super::super::rational::*;
#[cfg(test)]
use test::Bencher;
#[cfg(test)]
use super::super::fq::fq_small::*;

#[test]
fn test_duplication() {
    let F = r64::RING.bind_ring_by_value();
    let E = EllipticCurve::new(F, F.zero(), F.one());
    let poly_ring = PolyRing::adjoint(&r64::RING, "x");
    let P = poly_ring.bind_ring_by_value();
    let x = P.from(P.wrapped_ring().unknown());
    let f = x.pow(3) + 1;
    let a4 = P.zero();
    let (x3, y3_y, z3) = duplication(&P, a4.clone(), x, f.clone(), P.one(), f);
    let point = EllipticCurvePoint::Affine(F.from_z(2), F.from_z(3));
    let point2 = E.mul_point(&point, &BigInt::from(2), &F);
    assert_eq!(point2.x().unwrap().borrow_ring(), x3(F.from_z(2).borrow_ring()) / z3(F.from_z(2).borrow_ring()));
    assert_eq!(point2.y().unwrap().borrow_ring(), y3_y(F.from_z(2).borrow_ring()) / z3(F.from_z(2).borrow_ring()) / F.from_z(3).borrow_ring());
}

#[test]
fn test_xy1_addition() {
    let F = r64::RING.bind_ring_by_value();
    let E = EllipticCurve::new(F, F.zero(), F.one());
    let poly_ring = PolyRing::adjoint(&r64::RING, "x");
    let P = poly_ring.bind_ring_by_value();
    let x = P.from(P.wrapped_ring().unknown());
    let f = x.pow(3) + 1;
    let a4 = P.zero();
    let (x2, y2_y, z2) = duplication(&P, a4.clone(), x, f.clone(), P.one(), f.clone());
    let (x3, y3_y, z3) = xy1_addition(&P, x2, y2_y, z2, f);
    let point = EllipticCurvePoint::Affine(F.from_z(2), F.from_z(3));
    let point3 = E.mul_point(&point, &BigInt::from(3), &F);
    assert_eq!(point3.x().unwrap().borrow_ring(), x3(F.from_z(2).borrow_ring()) / z3(F.from_z(2).borrow_ring()));
    assert_eq!(point3.y().unwrap().borrow_ring(), y3_y(F.from_z(2).borrow_ring()) / z3(F.from_z(2).borrow_ring()) / F.from_z(3).borrow_ring());
}

#[test]
fn test_division_polynomials() {
    let field = r64::RING.bind_ring();
    let E = EllipticCurve::new(field, field.zero(), field.one());
    let (f, _, h) = division_polynomials(&E, 2);
    let poly_ring = PolyRing::adjoint(field.wrapped_ring().clone(), "X").bind_ring_by_value();
    let x = poly_ring.from(poly_ring.wrapped_ring().unknown());
    let f_expected = x.pow(4) - &x * 8;
    let h_expected = x.pow(3) * 4 + 4;
    let d = gcd(&f.ring(), f.clone(), h.clone());
    let f = f / &d;
    let h = h / &d;
    let base_ring = poly_ring.wrapped_ring().base_ring().bind_ring_by_value();
    let i = embedding(&base_ring, &poly_ring);
    assert_eq!(f_expected, &f / i(f.lc().unwrap()));
    assert_eq!(h_expected, &h / i(f.lc().unwrap()));
 
    let (f, _, h) = division_polynomials(&E, 3);
    let i = z_hom(&field);
    let point = EllipticCurvePoint::Affine(i(2), i(3));
    assert!(E.is_on_curve(&point));
    let point_3 = E.mul_point(&point, &BigInt::from(3), &field);
    assert_eq!(*point_3.x().unwrap(), f(i(2)) / h(i(2)));
}

#[bench]
fn bench_division_poly(b: &mut Bencher) {
    let ring = F49.bind_ring_by_value();
    let E = EllipticCurve::new(ring.clone(), ring.zero(), ring.one());
    let mut rng = oorandom::Rand32::new(3);
    let n: i64 = 53;
    b.iter(|| {
        let (f, _, h) = division_polynomials(&E, n as usize);
        let d = gcd(&f.ring(), f.clone(), h.clone());
        let P = E.random_affine_point(|| rng.rand_u32());
        assert_eq!(f(P.x().unwrap().clone()) / h(P.x().unwrap().clone()), *E.mul_point(&P, &BigInt::from(n), E.base_field()).x().unwrap());
    });
    println!("{}", super::super::fq::zn_small::operations.load(std::sync::atomic::Ordering::SeqCst));
}
