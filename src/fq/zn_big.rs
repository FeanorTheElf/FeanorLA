use super::super::prelude::*;
use super::super::eea::*;
use super::*;

use std::cell::Cell;
use std::cmp::Ordering;

#[derive(Debug, Clone)]
pub struct Zn<I: IntegerRing> {
    integer_ring: I,
    modulus: El<I>,
    inverse_modulus: El<I>,
    inverse_modulus_bitshift: u64,
    integral: Cell<Option<bool>>
}

impl<I: IntegerRing> Zn<I> {

    pub fn new(integer_ring: I, modulus: El<I>) -> Self {
        assert!(integer_ring.cmp(&modulus, &integer_ring.from_z(2)) != Ordering::Less);
        // have k such that 2^k > modulus^2
        // then (2^k / modulus) * x >> k differs at most 1 from floor(x / modulus)
        // if x < n^2, which is the case after multiplication
        let k = integer_ring.abs_log2_floor(&modulus) * 2 + 2;
        let inverse_modulus = integer_ring.euclidean_div(integer_ring.pow(&integer_ring.from_z(2), k as u32), &modulus);
        return Zn {
            integer_ring: integer_ring,
            modulus: modulus,
            inverse_modulus: inverse_modulus,
            inverse_modulus_bitshift: k,
            integral: Cell::from(None)
        };
    }

    fn project_leq_n_square(&self, n: El<I>) -> El<I> {
        assert!(self.integer_ring.cmp(&n, &self.integer_ring.zero()) != Ordering::Less);
        let mut subtract = self.integer_ring.mul_ref(&n, &self.inverse_modulus);
        subtract = self.integer_ring.euclidean_div_pow_2(subtract, self.inverse_modulus_bitshift);
        subtract = self.integer_ring.mul_ref(&subtract, &self.modulus);
        let mut m = self.integer_ring.sub(n, subtract);
        if self.integer_ring.cmp(&m, &self.modulus) != Ordering::Less {
            m = self.integer_ring.sub_ref_snd(m, &self.modulus);
        }
        assert_eq!(self.integer_ring.cmp(&m, &self.modulus), Ordering::Less, "The input is not smaller than {}^2", self.integer_ring.display(&self.modulus));
        return m;
    }

    pub fn project(&self, n: El<I>) -> El<Self> {
        let mut red_n = n;
        let negated = self.integer_ring.cmp(&red_n, &self.integer_ring.zero()) == Ordering::Less;
        if negated {
            red_n = self.integer_ring.neg(red_n);
        }
        if self.integer_ring.cmp(&red_n, &self.modulus) == Ordering::Less {
            // already in the interval [0, self.modulus[
        } else if self.integer_ring.abs_log2_floor(&red_n) + 1 < self.integer_ring.abs_log2_floor(&self.modulus) * 2 {
            red_n = self.project_leq_n_square(red_n);
        } else {
            red_n = self.integer_ring.sub(red_n.clone(), self.integer_ring.mul_ref(&self.integer_ring.euclidean_div(red_n, &self.modulus), &self.modulus));
        };
        if negated {
            red_n = self.integer_ring.sub_ref_fst(&self.modulus, red_n);
        }
        debug_assert!(self.integer_ring.cmp(&red_n, &self.modulus) == Ordering::Less);
        return FactorRingZEl(red_n);
    }

    ///
    /// Returns either the inverse of x (as Ok()) or a nontrivial 
    /// factor of the modulus (as Err())
    /// 
    pub fn invert(&self, x: El<I>) -> Result<El<I>, El<I>> {
        let (s, _, d) = eea(&self.integer_ring, x.clone(), self.modulus.clone());
        if self.integer_ring.is_neg_one(&d) || self.integer_ring.is_one(&d) {
            Ok(s)
        } else {
            Err(d)
        }
    }
}

impl<I: IntegerRing> PartialEq for Zn<I> {

    fn eq(&self, rhs: &Zn<I>) -> bool {
        self.integer_ring.is_eq(&self.modulus, &rhs.modulus)
    }
}

#[derive(Clone, Debug)]
pub struct ZnIterFn<I: IntegerRing> {
    current: El<I>
}

impl<I: IntegerRing> FiniteRingIterFn<Zn<I>> for ZnIterFn<I> {

    fn next(&mut self, ring: &Zn<I>) -> Option<FactorRingZEl<I>> {
        ring.integer_ring.add_assign(&mut self.current, ring.integer_ring.one());
        if ring.integer_ring.cmp(&self.current, &ring.modulus) == Ordering::Less {
            Some(ring.project(self.current.clone()))
        } else {
            None
        }
    }
}

impl<I: IntegerRing> FiniteRing for Zn<I> {

    type IterFn = ZnIterFn<I>;

    fn size(&self) -> BigInt {
        self.characteristic()
    }

    fn iter_fn(&self) -> Self::IterFn {
        ZnIterFn {
            current: self.integer_ring.from_z(-1)
        }
    }

    fn random_element<G>(&self, rng: G) -> El<Self> 
        where G: FnMut() -> u32
    {
        FactorRingZEl(self.integer_ring.get_uniformly_random(rng, &self.modulus))
    }
}

#[derive(Debug, Clone)]
pub struct FactorRingZEl<I: IntegerRing>(El<I>);

impl<I: IntegerRing> RingBase for Zn<I> {

    type El = FactorRingZEl<I>;

    fn add_ref(&self, FactorRingZEl(lhs): Self::El, FactorRingZEl(rhs): &Self::El) -> Self::El {
        assert!(self.integer_ring.cmp(&lhs, &self.modulus) == Ordering::Less);
        assert!(self.integer_ring.cmp(&rhs, &self.modulus) == Ordering::Less);

        let mut result = self.integer_ring.add_ref(lhs, &rhs);
        if self.integer_ring.cmp(&result, &self.modulus) != Ordering::Less {
            result = self.integer_ring.sub_ref_snd(result, &self.modulus);
        }

        assert!(self.integer_ring.cmp(&result, &self.modulus) == Ordering::Less);
        return FactorRingZEl(result);
    }

    fn mul_ref(&self, FactorRingZEl(lhs): &Self::El, FactorRingZEl(rhs): &Self::El) -> Self::El {
        assert!(self.integer_ring.cmp(&lhs, &self.modulus) == Ordering::Less);
        assert!(self.integer_ring.cmp(&rhs, &self.modulus) == Ordering::Less);

        let result = self.project_leq_n_square(
            self.integer_ring.mul_ref(lhs, rhs)
        );

        assert!(self.integer_ring.cmp(&result, &self.modulus) == Ordering::Less);
        return FactorRingZEl(result);
    }

    fn neg(&self, FactorRingZEl(val): Self::El) -> Self::El {
        assert!(self.integer_ring.cmp(&val, &self.modulus) == Ordering::Less);

        let mut result = self.integer_ring.neg(val);
        if !self.integer_ring.is_zero(&result) {
            result = self.integer_ring.add_ref(result, &self.modulus);
        }

        assert!(self.integer_ring.cmp(&result, &self.modulus) == Ordering::Less);
        return FactorRingZEl(result);
    }

    fn zero(&self) -> Self::El {
        FactorRingZEl(self.integer_ring.zero())
    }

    fn one(&self) -> Self::El {
        FactorRingZEl(self.integer_ring.one())
    }

    fn is_eq(&self, FactorRingZEl(lhs): &Self::El, FactorRingZEl(rhs): &Self::El) -> bool {
        assert!(self.integer_ring.cmp(&lhs, &self.modulus) == Ordering::Less);
        assert!(self.integer_ring.cmp(&rhs, &self.modulus) == Ordering::Less);
        self.integer_ring.is_eq(&lhs, &rhs)
    }

    fn characteristic(&self) -> BigInt {
        self.integer_ring.preimage(&BigInt::RING, self.modulus.clone())
    }

    default fn is_integral(&self) -> RingPropValue {
        return RingPropValue::Unknown;
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
        write!(f, "[{}]_{}", self.integer_ring.display(&el), self.integer_ring.display(&self.modulus))
    }
}

impl<I: IntegerRing> RingBase for Zn<I>
    where I: UfdInfoRing
{
    fn is_integral(&self) -> RingPropValue {
        if self.integral.get().is_none() {
            let modulus_prime = self.integer_ring.is_prime(&self.modulus);
            self.integral.set(Some(modulus_prime));
        }
        return match self.integral.get().unwrap() {
            true => RingPropValue::True,
            false => RingPropValue::False
        }
    }
}

impl<I: IntegerRing> DivisibilityInfoRing for Zn<I> {
    
    fn is_divisibility_computable(&self) -> RingPropValue {
        RingPropValue::True
    }
    
    fn quotient(&self, FactorRingZEl(lhs): &Self::El, FactorRingZEl(rhs): &Self::El) -> Option<Self::El> {
        let d = gcd(&self.integer_ring, lhs.clone(), rhs.clone());
        if let Ok(inv) = self.invert(self.integer_ring.quotient(rhs, &d).unwrap()) {
            return Some(self.project(self.integer_ring.mul(inv, self.integer_ring.quotient(lhs, &d).unwrap())));
        } else {
            return None;
        }
    }

    fn is_unit(&self, FactorRingZEl(el): &Self::El) -> bool {
        self.integer_ring.is_one(&signed_gcd(&self.integer_ring, el.clone(), self.modulus.clone()))
    }
}

impl<I: IntegerRing> CanonicalEmbeddingInfo<Zn<I>> for Zn<I> {

    fn has_embedding(&self, from: &Zn<I>) -> RingPropValue {
        self.integer_ring.has_embedding(&from.integer_ring) & self.integer_ring.is_eq(&self.modulus, &self.integer_ring.embed(&from.integer_ring, from.modulus.clone()))
    }

    fn embed(&self, from: &Zn<I>, FactorRingZEl(el): El<Zn<I>>) -> El<Zn<I>> {
        assert!(self.has_embedding(from).can_use());
        FactorRingZEl(self.integer_ring.embed(&from.integer_ring, el))
    }
}

impl<I: IntegerRing> CanonicalIsomorphismInfo<Zn<I>> for Zn<I> {

    fn has_isomorphism(&self, from: &Zn<I>) -> RingPropValue {
        self.integer_ring.has_isomorphism(&from.integer_ring) & self.integer_ring.is_eq(&self.modulus, &self.integer_ring.embed(&from.integer_ring, from.modulus.clone()))
    }

    fn preimage(&self, from: &Zn<I>, FactorRingZEl(el): El<Zn<I>>) -> El<Zn<I>> {
        assert!(self.has_isomorphism(from).can_use());
        FactorRingZEl(self.integer_ring.preimage(&from.integer_ring, el))
    }
}

#[test]
fn test_mul() {
    let z257 = Zn::new(BigInt::RING, BigInt::from_str_radix("257", 10).unwrap());
    let x = BigInt::from_str_radix("256", 10).unwrap();
    assert!(z257.is_eq(&z257.one(), &z257.mul(z257.project(x.clone()), z257.project(x))));
}
