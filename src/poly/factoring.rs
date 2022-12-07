use super::super::prelude::*;
use super::super::eea::*;
use super::super::finite_extension::finite_extension_impl::*;
use super::super::finite_extension::*;
use super::super::fq::*;
use super::*;
use super::super::square_multiply::abs_square_and_multiply;
use super::hensel;

use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};

use oorandom;

fn pow_mod_f<P>(poly_ring: &P, g: &El<P>, f: &El<P>, pow: &StdInt) -> El<P>
    where P: PolyRing + EuclideanInfoRing
{
    assert!(*pow >= 0);
    return abs_square_and_multiply(
        g, 
        pow, 
        StdInt::RING, 
        |x, y| poly_ring.euclidean_rem(poly_ring.mul(x, y), f), 
        |x, y| poly_ring.euclidean_rem(poly_ring.mul_ref(x, y), f), 
        poly_ring.one()
    );
}

pub fn distinct_degree_factorization<P>(poly_ring: P, mut f: El<P>) -> Vec<El<P>>
    where P: PolyRing + EuclideanInfoRing, P::BaseRing: FiniteRing
{
    let p = poly_ring.base_ring().size();
    assert!(poly_ring.base_ring().is_field().can_use());
    assert!(poly_ring.is_euclidean().can_use());
    assert!(!poly_ring.is_zero(&f));

    let mut result = Vec::new();
    result.push(poly_ring.one());
    let mut x_power_q_mod_f = poly_ring.unknown();
    while poly_ring.deg(&f) != Some(0) {
        // technically, we could just compute gcd(f, X^(p^i) - X), however p^i might be
        // really large and eea will be very slow. Hence, we do the first modulo operation
        // X^(p^i) mod f using square-and-multiply in the ring F[X]/(f)
        x_power_q_mod_f = pow_mod_f(&poly_ring, &x_power_q_mod_f, &f, &p);
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
    f: El<P>, 
    d: usize
) -> El<P>
    where P: PolyRing + EuclideanInfoRing, P::BaseRing: FiniteRing
{
    let p = poly_ring.base_ring().size();
    assert!(p.is_odd());
    assert!(poly_ring.base_ring().is_field().can_use());
    assert!(poly_ring.is_euclidean().can_use());
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
        let exp = (p.pow(d as u32) - 1) / 2;
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
    assert!(poly_ring.is_euclidean().can_use());
    assert!(poly_ring.is_divisibility_computable().can_use());

    let derivate = poly_ring.derive(poly.clone());
    let square_part = eea(&poly_ring, poly.clone(), derivate).2;
    return poly_ring.quotient(&poly, &square_part).unwrap();
}

pub fn factor_complete<'a, P>(poly_ring: &'a P, mut el: El<P>) -> VecMap<RingElWrapper<&'a P>, usize> 
    where P: PolyRing + EuclideanInfoRing + DivisibilityInfoRing, P::BaseRing: FiniteRing
{
    assert!(poly_ring.base_ring().is_field().can_use());
    assert!(poly_ring.is_euclidean().can_use());
    assert!(poly_ring.is_divisibility_computable().can_use());
    assert!(!poly_ring.is_zero(&el));

    let mut result = VecMap::new();
    let mut unit = poly_ring.base_ring().one();
    let wrapped_ring = WrappingRing::new(poly_ring);

    // we repeatedly remove the square-free part
    while !poly_ring.is_unit(&el) {
        let sqrfree_part = poly_squarefree_part(poly_ring, el.clone());

        // factor the square-free part into distinct-degree factors
        for (d, el) in distinct_degree_factorization(poly_ring, sqrfree_part.clone()).into_iter().enumerate() {
            let mut stack = Vec::new();
            stack.push(el);
            
            // and finally extract each individual factor
            while let Some(el) = stack.pop() {
                let (el, scaling) = poly_ring.normalize(el);
                unit = poly_ring.base_ring().mul(unit, scaling);
                if poly_ring.is_one(&el) {
                    continue;
                } else if poly_ring.deg(&el) == Some(d) {
                    let wrapped_el = wrapped_ring.wrap(el);
                    if let Some(power) = result.get_mut(&wrapped_el) {
                        *power += 1;
                    } else {
                        debug_assert!(!poly_ring.is_unit(wrapped_el.val()));
                        result.insert(wrapped_el, 1);
                    }
                } else {
                    let factor = factoring::cantor_zassenhaus(poly_ring, el.clone(), d);
                    stack.push(poly_ring.quotient(&el, &factor).unwrap());
                    stack.push(factor);
                }
            }
        }
        el = poly_ring.quotient(&el, &sqrfree_part).unwrap();
    }
    unit = poly_ring.base_ring().mul(unit, poly_ring.coefficient_at(&el, 0));
    debug_assert!(!poly_ring.base_ring().is_zero(&unit));
    if !poly_ring.base_ring().is_one(&unit) {
        result.insert(wrapped_ring.wrap(poly_ring.from(unit)), 1);
    }
    return result;
}

///
/// Compute a factorization of `el` into irreducible elements in a polynomial ring over `Z/p^eZ`.
/// Since the base ring is not integral anymore, the polynomial ring is obviously not a UFD, and
/// this factorization is not necessarily unique. This function behaves as follows:
///  - If the leading coefficient of `el` is a unit in `Z/p^eZ`, there is a canonical factorization
///    given by lifting a factorization of `Z/pZ` via Hensel's Lemma. This canonical factorization
///    will be computed by the function.
///  - If the leading coefficient is not a unit in `Z/p^eZ`, this function just tries to compute an
///    arbitrary factorization into irreducible elements. No further guarantees are made.
///  - If the reduction of `el` to `Z/pZ` is not square-free, the factorization cannot be lifted
///    by Hensel's Lemma. In this case, the function fails and returns the factorization of the
///    reduction, which then contains square factors.  
/// 
pub fn factor_lifted<'a, 'b, P, R>(poly_ring: &'a P, el: El<P>, prime_poly_ring: &'b R) -> Result<VecMap<RingElWrapper<&'a P>, usize>, VecMap<RingElWrapper<&'b R>, usize>>
    where P: PolyRing + DivisibilityInfoRing, P::BaseRing: IntegerQuotientRing, R: PolyRing + EuclideanInfoRing, R::BaseRing: IntegerQuotientRing
{
    let is_prime_power = poly_ring.base_ring().characteristic().is_prime_power();
    assert!(is_prime_power.is_some());
    assert!(!poly_ring.is_zero(&el));
    let (p, _) = is_prime_power.unwrap();
    assert_eq!(p, prime_poly_ring.base_ring().size());
    assert!(prime_poly_ring.base_ring().is_field().can_use());

    let wrapped_ring = WrappingRing::new(poly_ring);
    let lifting_ring = poly_ring.base_ring().lifting_ring();
    let reduce = |x: &El<P>| hensel::map_poly(poly_ring, x, prime_poly_ring, |y| prime_poly_ring.base_ring().from_z_gen(poly_ring.base_ring().lift(&y, &lifting_ring), &lifting_ring));
    let factorization = factor_complete(prime_poly_ring, reduce(&el));
    let mut factorization_iter = factorization.iter();
    let mut current = el;
    let mut result = VecMap::new();
    while let Some((factor, power)) = <_ as Iterator>::next(&mut factorization_iter) {
        if *power > 1 {
            return Err(factorization);
        }
        let f = factor.clone().into_val();
        let g = factorization_iter.clone().map(|(factor, _)| factor).fold(prime_poly_ring.one(), |a, b| prime_poly_ring.mul_ref(&a, b.val()));
        let lifting = hensel::hensel_lift(prime_poly_ring, poly_ring, f, g, &current, prime_poly_ring);
        if let Ok((lifted_f, lifted_g)) = lifting {
            current = lifted_g;
            result.insert(wrapped_ring.wrap(lifted_f), 1);
        } else {
            return Err(factorization);
        }
    }
    return Ok(result);
}

pub fn is_prime<P>(poly_ring: &P, el: &El<P>) -> bool
    where P: PolyRing + DivisibilityInfoRing + EuclideanInfoRing, P::BaseRing: FiniteRing
{
    assert!(poly_ring.base_ring().is_field().can_use());
    let deg = poly_ring.deg(el);
    if deg == None || deg == Some(0) {
        return false;
    }
    if deg == Some(1) {
        return true;
    }
    let deg = deg.unwrap();
    let modulus = poly_ring.normalize(el.clone()).0;
    let ring = FiniteExtensionImpl::adjoin_element(poly_ring.base_ring().clone(), modulus, poly_ring, "X");
    let q = poly_ring.base_ring().characteristic();
    let mut current = ring.one();
    let x = ring.generator();
    let mut current_pow_x = x.clone();
    for d in 1..deg {
        current_pow_x = ring.pow_big(&current_pow_x, &q);
        let factor = ring.sub_ref_snd(current_pow_x.clone(), &x);
        ring.mul_assign(&mut current, &factor);
    }
    return poly_ring.is_unit(&gcd(poly_ring, ring.polynomial_repr(poly_ring, &current), el.clone()));
}

#[cfg(test)]
use super::super::rational::*;
#[cfg(test)]
use super::super::fq::zn_small::*;
#[cfg(test)]
use super::super::fq::zn_big::*;

#[test]
fn test_poly_squarefree_part() {
    let ring = WrappingRing::new(PolyRingImpl::adjoint(r64::RING, "X"));
    let x = ring.unknown();
    let a = (&x + 4) * (&x + 6).pow(2) * (&x - 2).pow(2) * (&x + 8).pow(3);
    let b = (&x + 4) * (&x + 6) * (&x - 2) * (&x + 8);
    let mut squarefree_part = ring.wrap(poly_squarefree_part(ring.wrapped_ring(), a.val().clone()));
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
    let distinct_degree_factorization = distinct_degree_factorization(ring.wrapped_ring(), a.val().clone());
    assert_eq!(expected.len(), distinct_degree_factorization.len());
    for (f, e) in distinct_degree_factorization.into_iter().zip(expected.into_iter()) {
        let mut f = ring.wrap(f);
        f.normalize();
        assert_eq!(e, f);
    }
}

#[test]
fn test_cantor_zassenhaus() {
    type Z7 = ZnEl<7>;
    let ring = WrappingRing::new(PolyRingImpl::adjoint(Z7::RING, "X"));
    let x = ring.unknown();

    let f = &x * &x + 1;
    let g = &x * &x + &x + 3;
    let p = &f * &g;
    let mut factor = ring.wrap(cantor_zassenhaus(ring.wrapped_ring(), p.val().clone(), 2));
    factor.normalize();
    assert!(factor == f || factor == g);
}

#[test]
fn test_is_prime() {
    let coeff_ring = Zn::new(i64::RING, i64::RING.from_z(2));
    let ring = WrappingRing::new(PolyRingImpl::adjoint(&coeff_ring, "X"));
    let x = ring.unknown();

    let p = x.pow(4) + &x + 1;
    let q = x.pow(4) + x.pow(2) + 1;
    let a = ring.from_z(1);
    assert_eq!(true, ring.is_prime(&p));
    assert_eq!(false, ring.is_prime(&q));
    assert_eq!(false, ring.is_prime(&a));

    let coeff_ring = Zn::new(i64::RING, i64::RING.from_z(2));
    let ring = WrappingRing::new(PolyRingImpl::adjoint(&coeff_ring, "X"));
    let x = ring.unknown();
}

#[test]
fn test_factor_lifted() {
    let ring= WrappingRing::new(PolyRingImpl::adjoint(ZnEl::<25>::RING, "X"));
    let prime_poly_ring = PolyRingImpl::adjoint(ZnEl::<5>::RING, "X");
    let x = ring.unknown();
    let f = x.pow(2) + &x * 5 + 2;
    let g1 = &x + 15;
    let g2 = &x + 4;
    let h = &x + 1;
    let three = ring.from_z(3);
    let factorization = factor_lifted(ring.wrapped_ring(), (&f * &g1 * &g2 * &h * 3).into_val(), &prime_poly_ring).unwrap();
    let mut expected = VecMap::new();
    expected.insert(f.borrow_ring(), 1);
    expected.insert(g1.borrow_ring(), 1);
    expected.insert(g2.borrow_ring(), 1);
    expected.insert(h.borrow_ring(), 1);
    expected.insert(three.borrow_ring(), 1);
    assert_eq!(expected, factorization);

    let f = x.pow(2) * 5 + &x + 1;
    let g = &x + 5;
    let factorization = factor_lifted(ring.wrapped_ring(), (&f * &g).into_val(), &prime_poly_ring).unwrap();
    let mut expected = VecMap::new();
    expected.insert(f.borrow_ring(), 1);
    expected.insert(g.borrow_ring(), 1);
    assert_eq!(expected, factorization);

    let f = (&x + 2).pow(2);
    let factorization = factor_lifted(ring.wrapped_ring(), f.into_val(), &prime_poly_ring);
    assert!(factorization.is_err());
}
