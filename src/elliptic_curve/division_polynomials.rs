use super::super::prelude::*;
use super::super::wrapper::*;
use super::super::poly::*;
use super::*;

use std::collections::HashMap;

#[allow(type_alias_bounds)]
type P<K: Ring> = WrappingRing<PolyRingImpl<K>>;

pub struct DivisionPolyArray<'a, K: Ring> {
    E: &'a EllipticCurve<WrappingRing<K>>,
    cache: HashMap<usize, El<P<K>>>,
    P: P<K>
}

impl<'a, K: Ring> DivisionPolyArray<'a, K> {

    pub fn new(E: &'a EllipticCurve<WrappingRing<K>>) -> Self {
        let F = E.base_field();
        let P = WrappingRing::new(PolyRingImpl::adjoint(F.wrapped_ring().clone(), "x"));
        let x = P.unknown();
        let i = embedding(F, &P);
    
        let mut cache = HashMap::new();
        cache.insert(0, P.zero());
        cache.insert(1, P.one());
        cache.insert(2, P.from_z(2));
        cache.insert(3, x.pow(4) * 3 + x.pow(2) * i(E.A.clone()) * 6 + &x * i(E.B.clone()) * 12 - i(E.A.pow(2)));
        cache.insert(4, x.pow(6) * 4 + x.pow(4) * i(E.A.clone()) * 20 + x.pow(3) * i(E.B.clone()) * 80 - x.pow(2) * i(E.A.pow(2)) * 20 
            - x * i(&E.A * &E.B) * 16 - i(E.B.pow(2)) * 32 - i(E.A.pow(3)) * 4);
        return DivisionPolyArray { E, P, cache };
    }

    fn x(&self) -> El<P<K>> {
        self.P.unknown()
    }

    fn i(&self) -> LiftedHom<K, PolyRingImpl<K>, StandardEmbedding<K, PolyRingImpl<K>>> {
        self.P.embedding()
    }

    fn short_division_polynomial_rec(&mut self, n: usize) {
        if self.cache.contains_key(&n) {
            return;
        }
        if n % 2 == 0 {
            let m = n / 2;
            self.short_division_polynomial_rec(m - 2);
            self.short_division_polynomial_rec(m - 1);
            self.short_division_polynomial_rec(m);
            self.short_division_polynomial_rec(m + 1);
            self.short_division_polynomial_rec(m + 2);
            let [phi0, phi1, phi2, phi3, phi4] = [
                self.cache.get(&(m - 2)).unwrap(),
                self.cache.get(&(m - 1)).unwrap(),
                self.cache.get(&(m)).unwrap(),
                self.cache.get(&(m + 1)).unwrap(),
                self.cache.get(&(m + 2)).unwrap()
            ];
            self.cache.insert(n, phi2 * (phi4 * phi1 * phi1 - phi0 * phi3 * phi3) / 2);
        } else {
            let m = (n - 1) / 2;
            self.short_division_polynomial_rec(m - 1);
            self.short_division_polynomial_rec(m);
            self.short_division_polynomial_rec(m + 1);
            self.short_division_polynomial_rec(m + 2);
            let [phi1, phi2, phi3, phi4] = [
                self.cache.get(&(m - 1)).unwrap(),
                self.cache.get(&(m)).unwrap(),
                self.cache.get(&(m + 1)).unwrap(),
                self.cache.get(&(m + 2)).unwrap()
            ];
            let x = self.x();
            let i = self.i();
            let f = x.pow(3) + x * i(self.E.A.clone()) + i(self.E.B.clone());
            if m % 2 == 0 {
                self.cache.insert(n, phi4 * phi2 * phi2 * phi2 * f.pow(2) - phi1 * phi3 * phi3 * phi3);
            } else {
                self.cache.insert(n, phi4 * phi2 * phi2 * phi2 - phi1 * phi3 * phi3 * phi3 * f.pow(2));
            }
        }
    }

    ///
    /// Returns the n-th short division polynomial, which is the polynomial
    /// ```text
    ///     n (X - x1) ... (X - xr)
    /// ```
    /// where xi runs through the x-coordinates of the points E[n] \ E[2].
    /// Note that each x-coordinate occurs only once, even though though there
    /// might be two points (x, y) and (x, -y) corresponding to that point.
    /// Hence, we see that the result has degree (n - 1)/2 if n is odd and
    /// (n - 4)/2 if n is even.
    /// 
    pub fn get(&mut self, n: usize) -> El<P<K>> {
        self.short_division_polynomial_rec(n);
        return self.cache.get(&n).cloned().unwrap();
    }
}

///
/// Computes two polynomials f and g with the property that
/// for each point `P = (x : y : z)` on the curve `E`, have that
/// the x-coordinate of `[n]P` is `f(x)/g(x)`, if `[n]P` is affine.
/// 
pub fn division_polynomials<K: Ring + PartialEq>(E: &EllipticCurve<WrappingRing<K>>, n: usize) 
    -> (El<WrappingRing<PolyRingImpl<K>>>, El<WrappingRing<PolyRingImpl<K>>>, El<WrappingRing<PolyRingImpl<K>>>) 
{
    assert!(n >= 2);
    let mut phis = DivisionPolyArray::new(E);
    let x = phis.x();
    let i = phis.i();
    let f = x.pow(3) + &x * i(E.A.clone()) + i(E.B.clone());

    if n % 2 == 0 {
        let result_x = phis.get(n).pow(2) * &x * &f - phis.get(n + 1) * phis.get(n - 1);
        let result_y = phis.get(n + 2) * phis.get(n - 1).pow(2) / 4 - phis.get(n - 2) * phis.get(n + 1).pow(2) / 4;
        return (result_x * phis.get(n), result_y, phis.get(n).pow(3) * f);
    } else {
        let result_x = phis.get(n).pow(2) * &x - phis.get(n + 1) * phis.get(n - 1) * &f;
        let result_y = phis.get(n + 2) * phis.get(n - 1).pow(2) * &f / 4 - phis.get(n - 2) * phis.get(n + 1).pow(2) * f / 4;
        return (result_x * phis.get(n), result_y, phis.get(n).pow(3));
    }
}

#[cfg(test)]
use super::super::rational::*;
#[cfg(test)]
use test::Bencher;
#[cfg(test)]
use super::super::fq::fq_small::*;
#[cfg(test)]
use super::super::eea::gcd;

#[test]
fn test_division_polynomials() {
    let field = r64::WRAPPED_RING;
    let E = EllipticCurve::new(field, field.zero(), field.one());
    let (f, _, h) = division_polynomials(&E, 2);
    let poly_ring = WrappingRing::new(PolyRingImpl::adjoint(field.wrapped_ring().clone(), "X"));
    let x = poly_ring.from(poly_ring.wrapped_ring().unknown());
    let f_expected = x.pow(4) - &x * 8;
    let h_expected = x.pow(3) * 4 + 4;
    let d = gcd(&f.ring(), f.clone(), h.clone());
    let f = f / &d;
    let h = h / &d;
    let base_ring = poly_ring.base_ring();
    let i = embedding(&base_ring, &poly_ring);
    assert_eq!(f_expected, &f / i(f.lc().unwrap()));
    assert_eq!(h_expected, &h / i(f.lc().unwrap()));
 
    let (f, _, h) = division_polynomials(&E, 3);
    let i = z_hom(&field);
    let point = EllipticCurvePoint::Affine(i(2), i(3));
    assert!(E.is_on_curve(&point, &field));
    let point_3 = E.mul_point(&point, &StdInt::from(3), &field);
    assert_eq!(*point_3.x().unwrap(), f(i(2)) / h(i(2)));
}

#[bench]
fn bench_division_poly(b: &mut Bencher) {
    let ring = WrappingRing::new(F49);
    let E = EllipticCurve::new(ring.clone(), ring.zero(), ring.one());
    let mut rng = oorandom::Rand32::new(3);
    let n: i64 = 13;
    b.iter(|| {
        let (f, _, h) = division_polynomials(&E, n as usize);
        let P = E.random_affine_point(|| rng.rand_u32());
        assert_eq!(f(P.x().unwrap().clone()) / h(P.x().unwrap().clone()), *E.mul_point(&P, &StdInt::from(n), E.base_field()).x().unwrap());
    });
}
