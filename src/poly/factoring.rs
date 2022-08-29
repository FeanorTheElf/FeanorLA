use super::super::prelude::*;
use super::super::eea::*;
use super::super::fq::*;
use super::*;
use super::super::square_multiply::abs_square_and_multiply;

use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};

use oorandom;

fn pow_mod_f<P>(poly_ring: &P, g: &El<P>, f: &El<P>, pow: &BigInt) -> El<P>
    where P: PolyRing + EuclideanInfoRing
{
    assert!(*pow >= 0);
    return abs_square_and_multiply(
        g, 
        pow, 
        BigInt::RING, 
        |x, y| poly_ring.euclidean_rem(poly_ring.mul(x, y), f), 
        |x, y| poly_ring.euclidean_rem(poly_ring.mul_ref(x, y), f), 
        poly_ring.one()
    );
}

pub fn distinct_degree_factorization<P>(poly_ring: P, p: &BigInt, mut f: El<P>) -> Vec<El<P>>
    where P: PolyRing + EuclideanInfoRing, P::BaseRing: FiniteRing
{
    assert!(poly_ring.base_ring().size() == *p);
    assert!(poly_ring.base_ring().is_field().can_use());
    assert!(!poly_ring.is_zero(&f));
    let mut result = Vec::new();
    result.push(poly_ring.one());
    let mut x_power_q_mod_f = poly_ring.unknown();
    while poly_ring.deg(&f) != Some(0) {
        // technically, we could just compute gcd(f, X^(p^i) - X), however p^i might be
        // really large and eea will be very slow. Hence, we do the first modulo operation
        // X^(p^i) mod f using square-and-multiply in the ring F[X]/(f)
        x_power_q_mod_f = pow_mod_f(&poly_ring, &x_power_q_mod_f, &f, p);
        let fq_defining_poly_mod_f = poly_ring.sub_ref_fst(&x_power_q_mod_f, poly_ring.unknown());
        let deg_i_factor = gcd(&poly_ring, f.clone(), fq_defining_poly_mod_f.clone());
        f = poly_ring.euclidean_div(f, &deg_i_factor);
        result.push(deg_i_factor);
    }
    result[0] = poly_ring.mul_ref(&result[0], &f);
    return result;
}

///
/// Uses the Cantor-Zassenhaus algorithm to find a nontrivial, factor of a polynomial f
/// over a finite field, that is squarefree and consists only of irreducible factors of 
/// degree d.
/// 
/// # Algorithm
/// 
/// The algorithm relies on the fact that for some monic polynomial T over Fp have
/// ```text
/// T^q - T = T (T^((q - 1)/2) + 1) (T^((q - 1)/2) - 1)
/// ```
/// where `q = p^d`. Furthermore, the three factors are pairwise coprime.
/// Since `X^q - X` divides `T^q - T`, and f is squarefree, we see that
/// f divides `T^q - T` and so
/// ```text
/// f = gcd(T, F) gcd(T (T^((q - 1)/2) + 1, f) gcd(T^((q - 1)/2) - 1, f)
/// ```
/// The idea is now to choose a random T and check whether `gcd(T^((q - 1)/2) - 1, f)`
/// gives a nontrivial factor of f. When f has two irreducible factors, with roots a, b
/// in Fq, then this works if exactly one of them maps to zero under the polynomial
/// `T^((q - 1)/2) - 1`. Now observe that this is the case if and only if `T(a)` resp.
/// `T(b)` is a square in Fq. Now apparently, for a polynomial chosen uniformly at random
/// among all monic polynomials of degree 2d in Fp[X], the values T(a) and T(b) are close
/// to independent and uniform on Fq, and thus the probability that one is a square and
/// the other is not is approximately 1/2.
/// 
/// ## Why is the degree of T equal to 2d ?
/// 
/// Pick an Fp-vector space basis of Fq and write a, b as dxs matrices A, B over Fp, where the
/// i-th column is the representation of `a^i` resp. `b^i` w.r.t. that basis. Then the
/// evaluation of `T(a)` resp. `T(b)` is the matrix-vector product `w^T A` resp. `w^T B` where
/// w is the coefficient vector of T (a vector over Fp). We want that `w^T A` and `w^T B` are
/// uniform and independent. Hence, we want all `w^T c` to be uniform and independent, where
/// c runs through the columns of A resp. B. There are 2d such columns in total, thus s = 2d
/// will do (note that all columns are different, as `1, a, ..., a^(d - 1)` is a basis of Fq
/// and similarly for b). 
///
#[allow(non_snake_case)]
pub fn cantor_zassenhaus<P>(
    poly_ring: P, 
    p: &BigInt, 
    f: El<P>, 
    d: usize
) -> El<P>
    where P: PolyRing + EuclideanInfoRing, P::BaseRing: FiniteRing
{
    assert!(p.is_odd());
    assert!(poly_ring.base_ring().size() == *p);
    assert!(poly_ring.deg(&f).unwrap() % d == 0);
    assert!(poly_ring.deg(&f).unwrap() > d);

    let mut hasher = DefaultHasher::new();
    p.hash(&mut hasher);
    let mut rng = oorandom::Rand32::new(hasher.finish());

    loop {
        let mut T = poly_ring.zero();
        let mut power_x = poly_ring.one();
        for _ in 0..(2 * d + 1) {
            T = poly_ring.add(T, poly_ring.mul(
                poly_ring.from(poly_ring.base_ring().random_element(|| rng.rand_u32())),
                power_x.clone()
            ));
            power_x = poly_ring.mul(power_x, poly_ring.unknown());
        }
        T = poly_ring.add(T, power_x);
        let exp = (p.clone().pow(d as u32) - 1).floor_div_small(2);
        let G = poly_ring.sub(pow_mod_f(&poly_ring, &T, &f, &exp), poly_ring.one());
        let g = eea(&poly_ring, f.clone(), G.clone()).2;
        if !poly_ring.is_unit(&g) && poly_ring.quotient(&g, &f).is_none() {
            return g;
        }
    }
}

pub fn poly_squarefree_part<P>(poly_ring: &P, poly: El<P>) -> El<P>
    where P: PolyRing + EuclideanInfoRing + DivisibilityInfoRing
{
    let derivate = poly_ring.derive(poly.clone());
    let square_part = eea(&poly_ring, poly.clone(), derivate).2;
    return poly_ring.quotient(&poly, &square_part).unwrap();
}

#[cfg(test)]
use super::super::rational::*;
#[cfg(test)]
use super::super::fq::zn_small::*;

#[test]
fn test_poly_squarefree_part() {
    let ring = WrappingRing::new(PolyRingImpl::adjoint(r64::RING, "X"));
    let x = ring.unknown();
    let a = (&x + 4) * (&x + 6).pow(2) * (&x - 2).pow(2) * (&x + 8).pow(3);
    let b = (&x + 4) * (&x + 6) * (&x - 2) * (&x + 8);
    let mut squarefree_part = ring.from(poly_squarefree_part(ring.wrapped_ring(), a.val().clone()));
    squarefree_part.normalize();
    assert_eq!(b, squarefree_part);
}

#[test]
fn test_distinct_degree_factorization() {
    let field = ZnEl::<2>::RING;
    let ring = WrappingRing::new(PolyRingImpl::adjoint(field, "X"));
    let x = ring.unknown();

    let a = &x * (&x + 1) * (&x * &x + &x + 1) * (&x * &x * &x + &x + 1) * (&x * &x * &x + &x * &x + 1);
    let a0 = ring.from_z(1);
    let a1 = &x * (&x + 1);
    let a2 = &x * &x + &x + 1;
    let a3 = (&x * &x * &x + &x + 1) * (&x * &x * &x + &x * &x + 1);
    let expected = vec![a0, a1, a2, a3];
    let distinct_degree_factorization = distinct_degree_factorization(ring.wrapped_ring(), &BigInt::from(2), a.val().clone());
    assert_eq!(expected.len(), distinct_degree_factorization.len());
    for (f, e) in distinct_degree_factorization.into_iter().zip(expected.into_iter()) {
        let mut f = ring.from(f);
        f.normalize();
        assert_eq!(e, f);
    }
}

#[test]
fn test_cantor_zassenhaus() {
    type Z7 = ZnEl<7>;
    let ring = WrappingRing::new(PolyRingImpl::adjoint(Z7::RING, "X"));
    let x = ring.unknown();
    let incl = ring.embedding();

    let f = &x * &x + 1;
    let g = &x * &x + &x + 3;
    let p = &f * &g;
    let mut factor = ring.from(cantor_zassenhaus(ring.wrapped_ring(), &BigInt::from(7), p.val().clone(), 2));
    factor.normalize();
    assert!(factor == f || factor == g);
}