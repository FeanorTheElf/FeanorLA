use super::bigint::BigInt;
use super::alg::*;

pub struct RingZn {
    modulus: BigInt,
    inverse_modulus: BigInt,
    inverse_modulus_bitshift: u64
}

impl Ring for RingZn {
    type El = BigInt;

    fn add(&self, lhs: Self::El, rhs: Self::El) -> Self::El {
        let mut result = lhs + rhs;
        if result >= self.modulus {
            result -= &self.modulus;
        }
        return result;
    }

    fn mul(&self, lhs: Self::El, rhs: Self::El) -> Self::El {
        let mut result = lhs * rhs;
        result %= &self.modulus;
        return result;
    }

    fn neg(&self, val: Self::El) -> Self::El {
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
}