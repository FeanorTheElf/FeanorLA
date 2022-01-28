use super::super::super::alg::*;
use super::super::super::la::mat::*;
use super::super::eea::*;
use super::super::bigint::*;
use super::super::primality::*;
use super::uni_var::*;

use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};

use oorandom;

fn pow_mod_f<F>(poly_ring: &PolyRing<F>, g: &<PolyRing<F> as Ring>::El, f: &<PolyRing<F> as Ring>::El, pow: &BigInt) -> <PolyRing<F> as Ring>::El
    where F: DivisibilityInformationRing
{
    let mut result = poly_ring.one();
    for i in (0..(pow.log2_floor() + 1)).rev() {
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
fn distinct_degree_factorization<F>(prime_field: F, p: BigInt, mut f: Vector<VectorOwned<F::El>, F::El>) -> Vec<Vector<VectorOwned<F::El>, F::El>>
    where F: DivisibilityInformationRing
{
    let poly_ring = PolyRing::adjoint(prime_field, "X");
    let mut result = Vec::new();
    let mut x_power_q_mod_f = poly_ring.unknown();
    if let Some(factor) = poly_ring.quotient(&f, &poly_ring.unknown()) {
        f = factor;
        result.push(poly_ring.unknown());
    } else {
        result.push(poly_ring.one());
    }
    while poly_ring.deg(&f) != None && poly_ring.deg(&f) != Some(0)  {
        // technically, we could just compute gcd(f, X^(p^i) - X), however p^i might be
        // really large and eea will be very slow. Hence, we do the first modulo operation
        // X^(p^i) mod f using square-and-multiply in the ring F[X]/(f)
        x_power_q_mod_f = pow_mod_f(&poly_ring, &x_power_q_mod_f, &f, &p);
        let fq_defining_poly_mod_f = poly_ring.sub_ref_fst(&x_power_q_mod_f, poly_ring.unknown());
        let deg_i_factor = eea(&poly_ring, f.clone(), fq_defining_poly_mod_f).2;
        f = poly_ring.euclidean_div(f, &deg_i_factor);
        result.push(deg_i_factor);
    }
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
/// ```
/// T^q - T = T (T^((q - 1)/2) + 1) (T^((q - 1)/2) - 1)
/// ```
/// where `q = p^d`. Furthermore, the three factors are pairwise coprime.
/// Since `X^q - X` divides `T^q - T`, and f is squarefree, we see that
/// f divides `T^q - T` and so
/// ```
/// f = gcd(T, F) gcd(T (T^((q - 1)/2) + 1, f) gcd(T^((q - 1)/2) - 1, f)
/// ```
/// The idea is now to choose a random T and check whether `gcd(T^((q - 1)/2) - 1, f)`
/// gives a nontrivial factor of f. When f has two irreducible factors, with roots a, b
/// in Fq, then this works if exactly one of them maps to zero under the polynomial
/// `T^((q - 1)/2) - 1`. Now observe that this is the case if and only if `T(a)` resp.
/// `T(b)` is a square in Fq. For a polynomial chosen uniformly at random among among
/// all monic Fp[X] of degree < 2d, T(a) and T(b) are independent and distributed uniformly
/// on Fq. Hence, we successfully separate them with probability 1/2.
/// 
#[allow(non_snake_case)]
fn cantor_zassenhaus<F>(prime_field: F, p: &BigInt, f: Vector<VectorOwned<F::El>, F::El>, d: usize) -> Vector<VectorOwned<F::El>, F::El>
    where F: DivisibilityInformationRing
{
    assert!(*p != 2);
    let poly_ring = PolyRing::adjoint(prime_field, "X");

    let mut hasher = DefaultHasher::new();
    p.hash(&mut hasher);
    let mut rng = oorandom::Rand32::new(hasher.finish());

    loop {
        let T0 = poly_ring.from_z_big(BigInt::get_uniformly_random_oorandom(&mut rng, p, 5));
        let T1 = poly_ring.mul(
            poly_ring.from_z_big(BigInt::get_uniformly_random_oorandom(&mut rng, p, 5)), 
            poly_ring.unknown()
        );
        let T = poly_ring.add(T0, T1);
        let exp = (p.clone() - 1) / 2;
        let g = eea(&poly_ring, f.clone(), pow_mod_f(&poly_ring, &T, &f, &exp)).2;
        if !poly_ring.is_unit(&g) && poly_ring.quotient(&g, &f).is_none() {
            return g;
        }
    }
}

impl<R> PolyRing<R>
    where R: Ring
{
    ///
    /// For a polynomial f, returns the square-free part of f, i.e. the
    /// product of all distinct factors of f. This is only defined up
    /// to multiplication by units.
    /// 
    /// Requires that the underlying ring is a field.
    /// 
    pub fn poly_squarefree_part(&self, poly: Vector<VectorOwned<R::El>, R::El>) -> Vector<VectorOwned<R::El>, R::El>
        where R: Ring
    {
        assert!(self.base_ring().is_field());
        let derivate = self.derive(poly.clone());
        let square_part = eea(self, poly.clone(), derivate).2;
        return self.div(poly, &square_part);
    }
}

#[cfg(test)]
use super::super::super::alg_env::*;
#[cfg(test)]
use super::super::rat::*;
#[cfg(test)]
use super::super::zn::*;

#[test]
fn test_poly_squarefree_part() {
    let ring = PolyRing::adjoint(r64::RING, "X");
    let x = ring.bind::<RingAxiomsEuclideanRing>(ring.unknown());
    let a = (&x + 4) * (&x + 6).pow(2) * (&x - 2).pow(2) * (&x + 8).pow(3);
    let b = (&x + 4) * (&x + 6) * (&x - 2) * (&x + 8);
    let mut squarefree_part = ring.bind(ring.poly_squarefree_part(a.unwrap().clone()));
    squarefree_part /= ring.bind(ring.from(*ring.lc(squarefree_part.unwrap()).unwrap()));
    assert_eq!(b, squarefree_part);
}

#[test]
fn test_distinct_degree_factorization() {
    let field = ZnEl::<2>::RING;
    let ring = PolyRing::adjoint(field, "X");
    let x = ring.bind::<RingAxiomsEuclideanRing>(ring.unknown());

    let a = &x * (&x + 1) * (&x + 2) * (&x * &x + &x + 1) * (&x * &x * &x + &x + 1) * (&x * &x * &x + &x * &x + 1);
    let a0 = x.clone();
    let a1 = (&x + 1) * (&x + 2);
    let a2 = &x * &x + &x + 1;
    let a3 = (&x * &x * &x + &x + 1) * (&x * &x * &x + &x * &x + 1);
    let expected = vec![a0, a1, a2, a3];
    let distinct_degree_factorization = distinct_degree_factorization(field, BigInt::from(2), a.unwrap().clone());
    assert_eq!(expected.len(), distinct_degree_factorization.len());
    for (f, e) in distinct_degree_factorization.into_iter().zip(expected.into_iter()) {
        let mut f = ring.bind::<RingAxiomsEuclideanRing>(f);
        f /= ring.bind(ring.from(*ring.lc(f.unwrap()).unwrap()));
        assert_eq!(e, f);
    }
}

#[ignore]
#[test]
fn test_cantor_zassenhaus() {
    type Z7 = ZnEl<7>;
    let ring = PolyRing::adjoint(Z7::RING, "X");
    let x = ring.bind::<RingAxiomsEuclideanRing>(ring.unknown());

    let f = &x * &x + 1;
    let g = &x * &x + &x + 3;
    let p = &f * &g;
    let factor = cantor_zassenhaus(Z7::RING, &BigInt::from(7), p.unwrap().clone(), 2);
    println!("{}", ring.display(&factor));
    assert!(false);
}