use super::super::eea::*;
use super::super::prelude::*;
use super::super::wrapper::*;
use super::super::poly::*;
use super::*;

type P<K: Ring> = WrappingRing<PolyRing<K>>;

///
/// Computes three polynomials `(f, g_y, h)` with the property that
/// for each point `P = (x : y : z)` on the curve `E`, have that
/// `[n]P = (f(x) : y g_h(x) : z(x))`.
/// 
/// Note that giving the second polynomial multiplied by y instead of
/// itself allows us to work only with univariate polynomials, which
/// simplifies the whole thing.
/// 
pub fn division_polynomials<K: Ring>(E: &EllipticCurve<WrappingRing<K>>, n: usize) 
    -> (El<P<K>>, El<P<K>>, El<P<K>>) 
{
    assert!(n > 0);
    let F = E.base_field();
    let (K, x, y) = E.function_field();
    let P = EllipticCurvePoint::Affine(x, y);
}

#[cfg(test)]
use super::super::rational::*;
#[cfg(test)]
use test::Bencher;
#[cfg(test)]
use super::super::fq::fq_small::*;

#[test]
fn test_division_polynomials() {
    let ring = r64::RING.bind_ring();
    let E = EllipticCurve::new(ring, ring.zero(), ring.one());
    let (f, _, h) = division_polynomials(&E, 2);
    let poly_ring = PolyRing::adjoint(ring.wrapped_ring().clone(), "X").bind_ring_by_value();
    let x = poly_ring.from(poly_ring.wrapped_ring().unknown());
    let f_expected = x.pow(4) - &x * 8;
    let h_expected = x.pow(3) * 4 + 4;
    assert_eq!(f_expected, f);
    assert_eq!(h_expected, h);
}

#[bench]
fn bench_division_poly(b: &mut Bencher) {
    let ring = F49.bind_ring_by_value();
    let E = EllipticCurve::new(ring.clone(), ring.zero(), ring.one());
    let mut rng = oorandom::Rand32::new(2);
    let n: i64 = 10;
    b.iter(|| {
        let (f, _, h) = division_polynomials(&E, n as usize);
        let P = E.random_affine_point(|| rng.rand_u32());
        assert_eq!(f(P.x().unwrap().clone()) / h(P.x().unwrap().clone()), *E.mul_point(&P, &BigInt::from(n), E.base_field()).x().unwrap());
    });
    println!("{}", super::super::fq::zn_small::operations.load(std::sync::atomic::Ordering::SeqCst));
}
