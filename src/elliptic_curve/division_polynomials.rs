use super::super::eea::*;
use super::super::prelude::*;
use super::super::fraction_field::*;
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
    let P: P<K> = PolyRing::adjoint(F.wrapped_ring().clone(), "X").bind_ring_by_value();
    let K = FieldOfFractions::new(P.wrapped_ring().clone()).bind_ring_by_value();
    let incl = embedding(F, &P);
    let x = P.wrapped_ring().bind_by_value(P.wrapped_ring().unknown());
    let A = incl(E.a4().clone());
    let B = incl(E.a6().clone());
    let f = x.pow(3) + &A * &x + &B;

    if n == 1 {
        return (x, f, P.one());
    }

    let mut x1 = x.pow(4) * 18 - &x * &f * 16 + x.pow(2) * &A * 12 + A.pow(2) * 2;
    let mut y1_y = -x.pow(6) * 27 - &A * x.pow(4) * 27 + x.pow(3) * &f * 28 - A.pow(2) * x.pow(2) * 9 + &A * &x * &f * 4 - A.pow(3) - &B * &f * 8;
    let mut z1 = &f * 8;

    for i in 2..n {
        let y1_squared = x1.pow(3) + &A * &x1 * z1.pow(2) + &B * z1.pow(3);
        let x3 = (y1_squared - &y1_y * z1.pow(2) * 2 + z1.pow(3) * &f) - (&x1 - &z1 * &x).pow(2) * (&z1 * &x + &x1);
        let z3 = (&x1 - &z1 * &x).pow(2) * &z1;
        (x1, y1_y, z1) = (
            &x3 * (&x1 - &z1 * &x),
            ((&z1 * &f - &y1_y) * x3) + (&z3 * (y1_y * &x - &f * &x1)),
            z3 * (x1 - z1 * &x)
        );
        let d = gcd(&P, x1.clone(), gcd(&P, y1_y.clone(), z1.clone()));
        (x1, y1_y, z1) = (x1 / &d, y1_y / &d, z1 / &d);
    }
    return (x1, y1_y, z1)
}

#[cfg(test)]
use super::super::rational::*;

#[test]
fn test_division_polynomials() {
    let ring = r64::RING.bind_ring();
    let E = EllipticCurve::new(ring, ring.zero(), ring.one());
    let (f, g_y, h) = division_polynomials(&E, 2);
    let poly_ring = PolyRing::adjoint(ring.wrapped_ring().clone(), "X").bind_ring_by_value();
    let x = poly_ring.from(poly_ring.wrapped_ring().unknown());
}