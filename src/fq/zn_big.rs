use super::super::prelude::*;
use super::super::eea::*;
use super::*;

use std::cell::Cell;

#[derive(Debug, Clone)]
pub struct Zn {
    modulus: BigInt,
    inverse_modulus: BigInt,
    inverse_modulus_bitshift: usize,
    integral: Cell<Option<bool>>
}

impl Zn {

    pub fn new(modulus: BigInt) -> Self {
        assert!(modulus >= 2);
        // have k such that 2^k > modulus^2
        // then (2^k / modulus) * x >> k differs at most 1 from floor(x / modulus)
        // if x < n^2, which is the case after multiplication
        let k = modulus.abs_log2_floor() * 2 + 2;
        let inverse_modulus = BigInt::RING.euclidean_div(BigInt::power_of_two(k), &modulus);
        return Zn {
            modulus: modulus,
            inverse_modulus: inverse_modulus,
            inverse_modulus_bitshift: k,
            integral: Cell::from(None)
        };
    }

    fn project_leq_n_square(&self, n: BigInt) -> BigInt {
        let mut subtract = n.clone();
        subtract = subtract * self.inverse_modulus.clone();
        subtract = subtract >> self.inverse_modulus_bitshift;
        subtract = subtract *  self.modulus.clone();
        let mut m = n - subtract;
        if m >= self.modulus {
            m = m - self.modulus.clone();
        }
        assert!(m < self.modulus, "The input is not smaller than {}^2", self.modulus);
        return m;
    }

    pub fn project(&self, n: BigInt) -> El<Self> {
        let mut red_n = n;
        let negated = red_n < 0;
        if negated {
            red_n = -red_n;
        }
        if red_n < self.modulus {
            // already in the interval [0, self.modulus[
        } else if red_n.abs_log2_floor() + 1 < 2 * self.modulus.abs_log2_floor() {
            red_n = self.project_leq_n_square(red_n);
        } else {
            red_n = red_n.clone() - BigInt::RING.euclidean_div(red_n, &self.modulus) * self.modulus.clone();
        };
        if negated {
            red_n = self.modulus.clone() - red_n;
        }
        debug_assert!(red_n < self.modulus);
        return FactorRingZEl(red_n);
    }

    ///
    /// Returns either the inverse of x (as Ok()) or a nontrivial 
    /// factor of the modulus (as Err())
    /// 
    pub fn invert(&self, x: BigInt) -> Result<BigInt, BigInt> {
        let (s, _, d) = eea(&BigInt::RING, x.clone(), self.modulus.clone());
        if d != 1 && d != -1i64 {
            Err(d)
        } else {
            Ok(s)
        }
    }
}

impl PartialEq for Zn {

    fn eq(&self, rhs: &Zn) -> bool {
        self.modulus == rhs.modulus
    }
}

#[derive(Clone, Debug)]
pub struct ZnIterFn {
    current: BigInt
}

impl FiniteRingIterFn<Zn> for ZnIterFn {

    fn next(&mut self, ring: &Zn) -> Option<FactorRingZEl> {
        self.current += 1;
        if self.current < ring.modulus {
            Some(ring.project(self.current.clone()))
        } else {
            None
        }
    }
}

impl FiniteRing for Zn {

    type IterFn = ZnIterFn;

    fn size(&self) -> BigInt {
        self.characteristic()
    }

    fn iter_fn(&self) -> Self::IterFn {
        ZnIterFn {
            current: BigInt::from(-1)
        }
    }

    fn random_element<G>(&self, rng: G) -> El<Self> 
        where G: FnMut() -> u32
    {
        FactorRingZEl(BigInt::get_uniformly_random(rng, &self.modulus, 5))
    }
}

#[derive(Debug, Clone)]
pub struct FactorRingZEl(BigInt);

impl RingBase for Zn {

    type El = FactorRingZEl;

    fn add_ref(&self, FactorRingZEl(lhs): Self::El, FactorRingZEl(rhs): &Self::El) -> Self::El {
        assert!(lhs < self.modulus);
        assert!(rhs < &self.modulus);

        let mut result = lhs + rhs.clone();
        if result >= self.modulus {
            result = result - self.modulus.clone();
        }

        assert!(result < self.modulus);
        return FactorRingZEl(result);
    }

    fn mul_ref(&self, FactorRingZEl(lhs): &Self::El, FactorRingZEl(rhs): &Self::El) -> Self::El {
        assert!(*lhs < self.modulus);
        assert!(*rhs < self.modulus);

        let result = self.project_leq_n_square(
            BigInt::RING.mul_ref(lhs, rhs)
        );

        assert!(result < self.modulus);
        return FactorRingZEl(result);
    }

    fn neg(&self, FactorRingZEl(val): Self::El) -> Self::El {
        assert!(val < self.modulus);

        let mut result = -val;
        if result < 0 {
            result = result + self.modulus.clone();
        }

        assert!(result < self.modulus);
        return FactorRingZEl(result);
    }

    fn zero(&self) -> Self::El {
        FactorRingZEl(BigInt::ZERO)
    }

    fn one(&self) -> Self::El {
        FactorRingZEl(BigInt::RING.one())
    }

    fn is_eq(&self, FactorRingZEl(lhs): &Self::El, FactorRingZEl(rhs): &Self::El) -> bool {
        assert!(lhs < &self.modulus);
        assert!(rhs < &self.modulus);
        lhs == rhs
    }

    fn characteristic(&self) -> BigInt {
        self.modulus.clone()
    }

    fn is_integral(&self) -> RingPropValue {
        if self.integral.get().is_none() {
            let modulus_prime = BigInt::RING.is_prime(&self.modulus);
            self.integral.set(Some(modulus_prime));
        }
        return match self.integral.get().unwrap() {
            true => RingPropValue::True,
            false => RingPropValue::False
        }
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
        write!(f, "[{}]_{}", el, self.modulus)
    }
}

impl DivisibilityInfoRing for Zn {
    
    fn is_divisibility_computable(&self) -> RingPropValue {
        RingPropValue::True
    }
    
    fn quotient(&self, FactorRingZEl(lhs): &Self::El, FactorRingZEl(rhs): &Self::El) -> Option<Self::El> {
        let d = gcd(&BigInt::RING, lhs.clone(), rhs.clone());
        if let Ok(inv) = self.invert(BigInt::RING.quotient(rhs, &d).unwrap()) {
            return Some(self.project(inv * BigInt::RING.quotient(lhs, &d).unwrap()));
        } else {
            return None;
        }
    }

    fn is_unit(&self, FactorRingZEl(el): &Self::El) -> bool {
        BigInt::RING.is_one(&signed_gcd(el.clone(), self.modulus.clone(), &BigInt::RING))
    }
}

impl CanonicalEmbeddingInfo<Zn> for Zn {

    fn has_embedding(&self, from: &Zn) -> RingPropValue {
        RingPropValue::True & (self.modulus == from.modulus)
    }

    fn embed(&self, from: &Zn, el: Self::El) -> Self::El {
        assert!(self.has_embedding(from).can_use());
        el
    }
}

impl CanonicalIsomorphismInfo<Zn> for Zn {

    fn has_isomorphism(&self, from: &Zn) -> RingPropValue {
        RingPropValue::True & (self.modulus == from.modulus)
    }

    fn preimage(&self, from: &Zn, el: Self::El) -> Self::El {
        assert!(self.has_isomorphism(from).can_use());
        el
    }
}

#[test]
fn test_mul() {
    let z257 = Zn::new(BigInt::from_str_radix("257", 10).unwrap());
    let x = BigInt::from_str_radix("256", 10).unwrap();
    assert!(z257.is_eq(&z257.one(), &z257.mul(z257.project(x.clone()), z257.project(x))));
}
