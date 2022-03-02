use super::super::super::ring::*;
use super::super::super::la::mat::*;
use super::super::super::bigint::*;
use super::super::eea::*;
use super::super::fq::*;
use super::uni_var::*;
use super::ops::*;

use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};

use oorandom;

fn pow_mod_f<F>(poly_ring: &PolyRing<F>, g: &<PolyRing<F> as Ring>::El, f: &<PolyRing<F> as Ring>::El, pow: &BigInt) -> <PolyRing<F> as Ring>::El
    where F: DivisibilityInfoRing
{
    if *pow == 0 {
        return poly_ring.one();
    }
    let mut result = poly_ring.one();
    for i in (0..(pow.abs_log2_floor() + 1)).rev() {
        if pow.is_bit_set(i) {
            result = poly_ring.mul(poly_ring.mul_ref(&result, g), result);
        } else {
            result = poly_ring.mul_ref(&result, &result);
        }
        result = poly_ring.euclidean_rem(result, &f);
    }
    return result;
}

// we cannot really check if the given field is really a prime field, so just assume that it is.
// furthermore, the input polynomial must be square-free
pub fn distinct_degree_factorization<F>(prime_field: F, p: &BigInt, mut f: Vector<VectorOwned<F::El>, F::El>) -> Vec<Vector<VectorOwned<F::El>, F::El>>
    where F: FiniteRing
{
    assert!(prime_field.size() == *p);
    let poly_ring = PolyRing::adjoint(prime_field.clone(), "X");
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
/// 
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
/// among all monic polynomials of degree ? in Fp[X], the values T(a) and T(b) are close
/// to independent and uniform on Fq, and thus the probability that one is a square and
/// the other is not is approximately 1/2.
///
#[allow(non_snake_case)]
pub fn cantor_zassenhaus<F>(prime_field: F, p: &BigInt, f: Vector<VectorOwned<<F as Ring>::El>, <F as Ring>::El>, d: usize) -> Vector<VectorOwned<<F as Ring>::El>, <F as Ring>::El>
    where F: FiniteRing
{
    assert!(*p != 2);
    assert!(prime_field.size() == *p);
    assert!(poly_degree(&prime_field, f.as_ref()).unwrap() % d == 0);
    assert!(poly_degree(&prime_field, f.as_ref()).unwrap() > d);
    let poly_ring = PolyRing::adjoint(prime_field, "X");

    let mut hasher = DefaultHasher::new();
    p.hash(&mut hasher);
    let mut rng = oorandom::Rand32::new(hasher.finish());

    loop {
        let mut T = poly_ring.zero();
        let mut power_x = poly_ring.one();
        for _ in 0..(2 * d - 1) {
            T = poly_ring.add(T, poly_ring.mul(
                poly_ring.from_z_big(BigInt::get_uniformly_random_oorandom(&mut rng, p, 5)),
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

pub fn poly_squarefree_part<R>(ring: &R, poly: Vector<VectorOwned<R::El>, R::El>) -> Vector<VectorOwned<R::El>, R::El>
    where R: Ring
{
    assert!(ring.is_field().can_use());
    let poly_ring = PolyRing::adjoint(ring, "X");
    let derivate = poly_ring.derive(poly.clone());
    let square_part = eea(&poly_ring, poly.clone(), derivate).2;
    return poly_ring.div(poly, &square_part);
}

#[cfg(test)]
use super::super::rat::*;
#[cfg(test)]
use super::super::fq::zn_small::*;
#[cfg(test)]
use super::super::super::wrapper::*;
#[cfg(test)]
use super::super::super::primitive::*;
#[cfg(test)]
use super::super::super::embedding::*;

#[test]
fn test_poly_squarefree_part() {
    let ring = PolyRing::adjoint(r64::RING, "X");
    let x = ring.bind(ring.unknown());
    let a = (&x + 4) * (&x + 6).pow(2) * (&x - 2).pow(2) * (&x + 8).pow(3);
    let b = (&x + 4) * (&x + 6) * (&x - 2) * (&x + 8);
    let mut squarefree_part = ring.bind(poly_squarefree_part(ring.base_ring(), a.val().clone()));
    squarefree_part = squarefree_part.ring().quotient(&squarefree_part, &ring.bind(ring.from(*ring.lc(squarefree_part.val()).unwrap()))).unwrap();
    assert_eq!(b, squarefree_part);
}

#[test]
fn test_distinct_degree_factorization() {
    let field = ZnEl::<2>::RING;
    let ring = PolyRing::adjoint(field, "X");
    let x = ring.bind(ring.unknown());

    let a = &x * (&x + 1) * (&x * &x + &x + 1) * (&x * &x * &x + &x + 1) * (&x * &x * &x + &x * &x + 1);
    let a0 = ring.bind(ring.from_z(1));
    let a1 = &x * (&x + 1);
    let a2 = &x * &x + &x + 1;
    let a3 = (&x * &x * &x + &x + 1) * (&x * &x * &x + &x * &x + 1);
    let expected = vec![a0, a1, a2, a3];
    let distinct_degree_factorization = distinct_degree_factorization(field, &BigInt::from(2), a.val().clone());
    assert_eq!(expected.len(), distinct_degree_factorization.len());
    for (f, e) in distinct_degree_factorization.into_iter().zip(expected.into_iter()) {
        let mut f = ring.bind(f);
        f = f.ring().quotient(&f, &ring.bind(ring.from(*ring.lc(f.val()).unwrap()))).unwrap();
        assert_eq!(e, f);
    }
}

#[test]
fn test_cantor_zassenhaus() {
    type Z7 = ZnEl<7>;
    let ring = PolyRing::adjoint(Z7::RING, "X");
    let x = ring.bind(ring.unknown());
    let incl = ring.bind_ring().embedding(ring.base_ring().bind_ring());

    let f = &x * &x + 1;
    let g = &x * &x + &x + 3;
    let p = &f * &g;
    let mut factor = ring.bind(cantor_zassenhaus(Z7::RING, &BigInt::from(7), p.val().clone(), 2));
    factor = factor.clone() / incl(ring.base_ring().bind(ring.lc(factor.val()).unwrap().clone()));
    assert!(factor == f || factor == g);
}