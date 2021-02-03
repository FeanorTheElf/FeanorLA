use super::bigint::BigInt;
use super::alg::*;

pub struct RingZn {
    modulus: BigInt,
    inverse_modulus: BigInt,
    inverse_modulus_bitshift: usize
}

impl RingZn {

    pub fn new(modulus: BigInt) -> Self {
        // have k such that 2^k > modulus^2
        // then (2^k / modulus) * x >> k differs at most 1 from floor(x / modulus)
        // if x < n^2, which is the case after multiplication
        let k = modulus.log2_floor() * 2 + 1;
        let inverse_modulus = BigInt::power_of_two(k) / modulus.clone();
        return RingZn {
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

    pub fn project(&self, n: BigInt) -> BigInt {
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
        return red_n;
    }
}

impl Ring for RingZn {
    type El = BigInt;

    fn add(&self, lhs: Self::El, rhs: Self::El) -> Self::El {
        debug_assert!(lhs < self.modulus);
        debug_assert!(rhs < self.modulus);

        let mut result = lhs + rhs;
        if result >= self.modulus {
            result -= &self.modulus;
        }
        return result;
    }

    fn mul(&self, lhs: Self::El, rhs: Self::El) -> Self::El {
        debug_assert!(lhs < self.modulus);
        debug_assert!(rhs < self.modulus);

        return self.project_leq_n_square(lhs * rhs);
    }

    fn neg(&self, val: Self::El) -> Self::El {
        debug_assert!(val < self.modulus);

        let mut result = -val;
        result += &self.modulus;
        return result;
    }

    fn zero(&self) -> Self::El {
        BigInt::zero()
    }

    fn one(&self) -> Self::El {
        BigInt::one()
    }

    fn eq(&self, lhs: &Self::El, rhs: &Self::El) -> bool {
        lhs == rhs
    }
}

#[test]
fn test_mul() {
    let z257 = RingZn::new(BigInt::from_str_radix("257", 10).unwrap());
    let x = BigInt::from_str_radix("256", 10).unwrap();
    assert_eq!(BigInt::one(), z257.mul(x.clone(), x));
}