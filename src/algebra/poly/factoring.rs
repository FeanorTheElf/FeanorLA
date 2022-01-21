use super::super::super::alg::*;
use super::super::super::la::mat::*;
use super::super::eea::*;
use super::super::bigint::*;
use super::super::primality::*;
use super::uni_var::*;

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
            let mut tmp = poly_ring.one();
            for i in (0..(p.log2_floor() + 1)).rev() {
                if p.is_bit_set(i) {
                    tmp = poly_ring.mul(poly_ring.mul_ref(&tmp, &x_power_q_mod_f), tmp);
                } else {
                    tmp = poly_ring.mul_ref(&tmp, &tmp);
                }
                tmp = poly_ring.euclidean_rem(tmp, &f);
            }
            x_power_q_mod_f = tmp;

            let fq_defining_poly_mod_f = poly_ring.sub_ref_fst(&x_power_q_mod_f, poly_ring.unknown());
            let deg_i_factor = eea(&poly_ring, f.clone(), fq_defining_poly_mod_f).2;
            f = poly_ring.euclidean_div(f, &deg_i_factor);
            result.push(deg_i_factor);
        }
        return result;
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
    let distinct_degree_factorization = PolyRing::<StaticRing::<ZnEl::<2>>>::distinct_degree_factorization(field, BigInt::from(2), a.unwrap().clone());
    assert_eq!(expected.len(), distinct_degree_factorization.len());
    for (f, e) in distinct_degree_factorization.into_iter().zip(expected.into_iter()) {
        let mut f = ring.bind::<RingAxiomsEuclideanRing>(f);
        f /= ring.bind(ring.from(*ring.lc(f.unwrap()).unwrap()));
        assert_eq!(e, f);
    }
}