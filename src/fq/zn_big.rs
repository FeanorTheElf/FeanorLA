use super::super::prelude::*;
use super::super::eea::*;
use super::super::integer::IntegerRing;
use super::*;

use std::cell::Cell;

#[derive(Debug, Clone)]
pub struct Zn<R: IntegerRing = BigIntRing> {
    base_ring: R,
    modulus: El<R>,
    inverse_modulus: El<R>,
    inverse_modulus_bitshift: u32,
    integral: Cell<Option<bool>>
}

impl<R: IntegerRing> Zn<R> {

    pub fn new(modulus: El<R>, base_ring: R) -> Self {
        assert!(base_ring.cmp(&modulus, &base_ring.from_z(2)) != std::cmp::Ordering::Less);
        // have k such that 2^k > modulus^2
        // then (2^k / modulus) * x >> k differs at most 1 from floor(x / modulus)
        // if x < n^2, which is the case after multiplication
        let k = base_ring.abs_log2_floor(&modulus) * 2 + 2;
        let inverse_modulus = base_ring.euclidean_div(base_ring.mul_pow_2(base_ring.one(), k), &modulus);
        return Zn {
            base_ring: base_ring,
            modulus: modulus,
            inverse_modulus: inverse_modulus,
            inverse_modulus_bitshift: k,
            integral: Cell::from(None)
        };
    }

    fn project_leq_n_square(&self, n: El<R>) -> El<R> {
        let mut subtract = self.base_ring.mul_ref(&n, &self.inverse_modulus);
        subtract = self.base_ring.floor_div_pow_2(subtract, self.inverse_modulus_bitshift);
        subtract = self.base_ring.mul_ref(&subtract, &self.modulus);
        let mut m = self.base_ring.sub(n, subtract);
        if self.base_ring.cmp(&m, &self.modulus) != std::cmp::Ordering::Less {
            m = self.base_ring.sub_ref_snd(m, &self.modulus);
        }
        assert!(self.base_ring.cmp(&m, &self.modulus) == std::cmp::Ordering::Less, "The input is not smaller than {}^2", self.base_ring.display(&self.modulus));
        return m;
    }

    pub fn project(&self, n: El<R>) -> El<Self> {
        let mut red_n = n;
        let negated = self.base_ring.cmp(&red_n, &self.base_ring.zero()) == std::cmp::Ordering::Less;
        if negated {
            red_n = self.base_ring.neg(red_n);
        }
        if self.base_ring.cmp(&red_n, &self.modulus) == std::cmp::Ordering::Less {
            // already in the interval [0, self.modulus[
        } else if self.base_ring.abs_log2_floor(&red_n) + 1 < 2 * self.base_ring.abs_log2_floor(&self.modulus) {
            red_n = self.project_leq_n_square(red_n);
        } else {
            red_n = self.base_ring.sub_ref_fst(&red_n, self.base_ring.mul(self.base_ring.euclidean_div(red_n, &self.modulus), self.modulus.clone()));
        };
        if negated {
            red_n = self.base_ring.sub_ref_fst(&self.modulus.clone(), red_n);
        }
        debug_assert!(self.base_ring.cmp(&red_n, &self.modulus) == std::cmp::Ordering::Less);
        return FactorRingZEl(red_n);
    }

    ///
    /// Returns either the inverse of x (as Ok()) or a nontrivial 
    /// factor of the modulus (as Err())
    /// 
    fn invert(&self, x: El<R>) -> Result<El<R>, El<R>> {
        let (s, _, d) = eea(&self.base_ring, x.clone(), self.modulus.clone());
        if !self.base_ring.is_one(&d) && !self.base_ring.is_neg_one(&d) {
            Err(d)
        } else {
            Ok(s)
        }
    }
}

#[derive(Clone, Debug)]
pub struct ZnIterFn<R>
    where R: IntegerRing
{
    current: El<R>
}

impl<R> FiniteRingIterFn<Zn<R>> for ZnIterFn<R>
    where R: IntegerRing
{
    fn next(&mut self, ring: &Zn<R>) -> Option<FactorRingZEl<R>> {
        ring.base_ring.add_assign(&mut self.current, ring.base_ring.one());
        if ring.base_ring.cmp(&self.current, &ring.modulus) == std::cmp::Ordering::Less {
            Some(ring.project(self.current.clone()))
        } else {
            None
        }
    }
}

impl<R> FiniteRing for Zn<R>
    where R: IntegerRing
{

    type IterFn = ZnIterFn<R>;

    fn size(&self) -> BigInt {
        self.characteristic()
    }

    fn iter_fn(&self) -> Self::IterFn {
        ZnIterFn {
            current: self.base_ring.neg(self.base_ring.one())
        }
    }

    fn random_element<G>(&self, rng: G) -> El<Self> 
        where G: FnMut() -> u32
    {
        FactorRingZEl(self.base_ring.get_uniformly_random(rng, &self.modulus))
    }
}

#[derive(Debug, Clone)]
pub struct FactorRingZEl<R>(El<R>)
    where R: IntegerRing;

impl<R> RingBase for Zn<R>
    where R: IntegerRing
{
    type El = FactorRingZEl<R>;

    fn add_ref(&self, FactorRingZEl(lhs): Self::El, FactorRingZEl(rhs): &Self::El) -> Self::El {
        let mut result = self.base_ring.add_ref(lhs, &rhs);
        if self.base_ring.cmp(&result, &self.modulus) != std::cmp::Ordering::Less {
            result = self.base_ring.sub_ref_snd(result, &self.modulus);
        }

        assert!(self.base_ring.cmp(&result, &self.modulus) == std::cmp::Ordering::Less);
        return FactorRingZEl(result);
    }

    fn mul_ref(&self, FactorRingZEl(lhs): &Self::El, FactorRingZEl(rhs): &Self::El) -> Self::El {
        let result = self.project_leq_n_square(
            self.base_ring.mul_ref(lhs, rhs)
        );

        assert!(self.base_ring.cmp(&result, &self.modulus) == std::cmp::Ordering::Less);
        return FactorRingZEl(result);
    }

    fn neg(&self, FactorRingZEl(val): Self::El) -> Self::El {
        if self.base_ring.is_zero(&val) {
            self.zero()
        } else {
            FactorRingZEl(self.base_ring.sub_ref_fst(&self.modulus, val))
        }
    }

    fn zero(&self) -> Self::El {
        FactorRingZEl(self.base_ring.zero())
    }

    fn one(&self) -> Self::El {
        FactorRingZEl(self.base_ring.one())
    }

    fn eq(&self, FactorRingZEl(lhs): &Self::El, FactorRingZEl(rhs): &Self::El) -> bool {
        self.base_ring.eq(&lhs, &rhs)
    }

    fn characteristic(&self) -> BigInt {
        self.base_ring.preimage(&BigInt::RING, self.modulus)
    }

    default fn is_integral(&self) -> RingPropValue {
        RingPropValue::Unknown
    }

    fn is_noetherian(&self) -> bool {
        true
    }

    fn is_field(&self) -> RingPropValue {
        self.is_integral()
    }

    fn div(&self, lhs: Self::El, rhs: &Self::El) -> Self::El { 
        assert!(self.is_field().can_use());
        assert!(!self.is_zero(rhs));
        let FactorRingZEl(rhs_repr) = rhs;
        let inverse = self.project(self.invert(rhs_repr.clone()).unwrap());
        debug_assert!(self.is_one(&self.mul_ref(&inverse, rhs)));
        return self.mul(lhs, inverse);
    }

    fn format(&self, FactorRingZEl(el): &Self::El, f: &mut std::fmt::Formatter, _in_prod: bool) -> std::fmt::Result {
        write!(f, "[{}]_{}", self.base_ring.display(&el), self.base_ring.display(&self.modulus))
    }
}

impl<R> RingBase for Zn<R>
    where R: IntegerRing + UfdInfoRing
{
    fn is_integral(&self) -> RingPropValue {
        if self.integral.get().is_none() {
            let modulus_prime = self.base_ring.is_prime(&self.modulus);
            self.integral.set(Some(modulus_prime));
        }
        return match self.integral.get().unwrap() {
            true => RingPropValue::True,
            false => RingPropValue::False
        }
    }
}

impl<R> DivisibilityInfoRing for Zn<R>
    where R: IntegerRing + DivisibilityInfoRing
{
    fn is_divisibility_computable(&self) -> RingPropValue {
        RingPropValue::True
    }
    
    fn quotient(&self, FactorRingZEl(lhs): &Self::El, FactorRingZEl(rhs): &Self::El) -> Option<Self::El> {
        let d = gcd(&self.base_ring, lhs.clone(), rhs.clone());
        if let Ok(inv) = self.invert(self.base_ring.quotient(&rhs, &d).unwrap()) {
            return Some(self.project(self.base_ring.mul(inv, self.base_ring.quotient(lhs, &d).unwrap())));
        } else {
            return None;
        }
    }

    fn is_unit(&self, FactorRingZEl(el): &Self::El) -> bool {
        self.base_ring.is_one(&signed_gcd(el.clone(), self.modulus.clone(), &self.base_ring))
    }
}

impl<R, S> CanonicalEmbeddingInfo<Zn<S>> for Zn<R>
    where R: IntegerRing, S: IntegerRing
{
    fn has_embedding(&self, _from: &Zn<S>) -> RingPropValue {
        RingPropValue::True
    }

    fn embed(&self, from: &Zn<S>, el: El<Zn<S>>) -> Self::El {
        self.project(self.base_ring.embed(&BigInt::RING, from.base_ring.preimage(&BigInt::RING, el.0)))
    }
}

impl<R, S> CanonicalIsomorphismInfo<Zn<S>> for Zn<R>
    where R: IntegerRing, S: IntegerRing
{
    fn has_isomorphism(&self, _from: &Zn<S>) -> RingPropValue {
        RingPropValue::True
    }

    fn preimage(&self, from: &Zn<S>, el: Self::El) -> El<Zn<S>> {
        from.embed(self, el)
    }
}

#[test]
fn test_mul() {
    let z257 = Zn::new(BigInt::from_str_radix("257", 10).unwrap(), BigInt::RING);
    let x = BigInt::from_str_radix("256", 10).unwrap();
    assert!(z257.eq(&z257.one(), &z257.mul(z257.project(x.clone()), z257.project(x))));
}
