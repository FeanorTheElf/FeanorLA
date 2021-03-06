use super::super::alg::*;
use super::bigint::*;
use super::eea::*;

pub struct QuotientRingZ {
    modulus: BigInt,
    inverse_modulus: BigInt,
    inverse_modulus_bitshift: usize
}

impl QuotientRingZ {

    pub fn new(modulus: BigInt) -> Self {
        assert!(modulus >= 2);
        // have k such that 2^k > modulus^2
        // then (2^k / modulus) * x >> k differs at most 1 from floor(x / modulus)
        // if x < n^2, which is the case after multiplication
        let k = modulus.log2_floor() * 2 + 2;
        let inverse_modulus = BigInt::power_of_two(k) / modulus.clone();
        return QuotientRingZ {
            modulus: modulus,
            inverse_modulus: inverse_modulus,
            inverse_modulus_bitshift: k
        };
    }

    fn project_leq_n_square(&self, mut n: BigInt) -> BigInt {
        let mut subtract = n.clone();
        subtract *= &self.inverse_modulus;
        subtract = subtract >> self.inverse_modulus_bitshift;
        subtract *= &self.modulus;
        n -= subtract;
        if n >= self.modulus {
            n -= &self.modulus;
        }
        return n;
    }

    pub fn project(&self, n: BigInt) -> <Self as Ring>::El {
        let mut red_n = if n < self.modulus {
            n
        } else if n.log2_floor() + 1 < 2 * self.modulus.log2_floor() {
            self.project_leq_n_square(n)
        } else {
            n.clone() - (n / self.modulus.clone()) * self.modulus.clone()
        };
        if red_n < BigInt::ZERO {
            red_n = -red_n;
            red_n += &self.modulus;
        }
        debug_assert!(red_n < self.modulus);
        return QuotientRingZEl(red_n);
    }

    ///
    /// Returns either the inverse of x (as Ok()) or a nontrivial 
    /// factor of the modulus (as Err())
    /// 
    pub fn invert(&self, x: BigInt) -> Result<BigInt, BigInt> {
        let (s, _t, d) = eea(&BigInt::RING, x, self.modulus.clone());
        if d != 1 && d != -1i64 {
            Err(d)
        } else {
            Ok(s)
        }
    }
}

#[derive(Debug, Clone)]
pub struct QuotientRingZEl(BigInt);

impl Ring for QuotientRingZ {

    type El = QuotientRingZEl;

    fn add_ref(&self, QuotientRingZEl(lhs): Self::El, QuotientRingZEl(rhs): &Self::El) -> Self::El {
        assert!(lhs < self.modulus);
        assert!(rhs < &self.modulus);

        let mut result = lhs + rhs;
        if result >= self.modulus {
            result -= &self.modulus;
        }

        assert!(result < self.modulus);
        return QuotientRingZEl(result);
    }

    fn mul_ref(&self, QuotientRingZEl(lhs): &Self::El, QuotientRingZEl(rhs): &Self::El) -> Self::El {
        assert!(*lhs < self.modulus);
        assert!(*rhs < self.modulus);

        let result = self.project_leq_n_square(
            BigInt::RING.mul_ref(lhs, rhs)
        );

        assert!(result < self.modulus);
        return QuotientRingZEl(result);
    }

    fn neg(&self, QuotientRingZEl(val): Self::El) -> Self::El {
        assert!(val < self.modulus);

        let mut result = -val;
        if result < 0 {
            result += &self.modulus;
        }

        assert!(result < self.modulus);
        return QuotientRingZEl(result);
    }

    fn zero(&self) -> Self::El {
        QuotientRingZEl(BigInt::zero())
    }

    fn one(&self) -> Self::El {
        QuotientRingZEl(BigInt::one())
    }

    fn eq(&self, QuotientRingZEl(lhs): &Self::El, QuotientRingZEl(rhs): &Self::El) -> bool {
        assert!(lhs < &self.modulus);
        assert!(rhs < &self.modulus);
        lhs == rhs
    }

    fn is_integral(&self) -> bool {
        unimplemented!()
    }

    fn is_euclidean(&self) -> bool {
        false
    }

    fn is_field(&self) -> bool {
        self.is_integral()
    }

    fn euclidean_div_rem(&self, _lhs: Self::El, _rhs: &Self::El) -> (Self::El, Self::El) { 
        panic!("Not a euclidean domain!");
    }

    fn div(&self, QuotientRingZEl(lhs): Self::El, QuotientRingZEl(rhs): &Self::El) -> Self::El { 
        match self.invert(rhs.clone()) {
            Err(factor) => panic!("Tried to divide in Z/{}Z, however this is not a field and the divisor is not invertible (the modulus has the nontrivial factor {})", self.modulus, factor),
            Ok(inverse) => self.mul(QuotientRingZEl(lhs), QuotientRingZEl(inverse))
        }
    }

    fn format(&self, QuotientRingZEl(el): &Self::El, f: &mut std::fmt::Formatter, _in_prod: bool) -> std::fmt::Result {
        write!(f, "{} mod {}", el, self.modulus)
    }
}

#[test]
fn test_mul() {
    let z257 = QuotientRingZ::new(BigInt::from_str_radix("257", 10).unwrap());
    let x = BigInt::from_str_radix("256", 10).unwrap();
    assert!(z257.eq(&z257.one(), &z257.mul(z257.project(x.clone()), z257.project(x))));
}