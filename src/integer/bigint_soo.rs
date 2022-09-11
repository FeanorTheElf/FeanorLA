use super::*;
use super::bigint::*;
use super::super::ring::*;

///
/// 
/// # Internal contract
/// 
/// The SOO variant contains all values x with |x| <= `i128::MAX`, and the BigInt variant contains
/// all others. This has the following consequences:
///  - If the variant is BigInt, it is guaranteed that the value does not fit into i128 (exception: next note)
///  - The only value stored in BigInt that would fit into SOO is `i128::MIN` (which is `-i128::MAX - 1`)
/// 
#[derive(Clone, Debug)]
pub enum BigIntSOO {
    BigInt(BigInt), 
    SOO(i128)
}

impl BigIntSOO {

    fn operation<F, G>(f: F, g: G, lhs: BigIntSOO, rhs: BigIntSOO) -> BigIntSOO
        where F: FnOnce(i128, i128) -> Option<i128>, G: FnOnce(BigInt, BigInt) -> BigInt
    {
        let result = match (lhs, rhs) {
            (BigIntSOO::BigInt(a), BigIntSOO::BigInt(b)) => g(a, b),
            (BigIntSOO::BigInt(a), BigIntSOO::SOO(b)) => g(a, BigInt::from(b)),
            (BigIntSOO::SOO(a), BigIntSOO::BigInt(b)) => g(BigInt::from(a), b),
            (BigIntSOO::SOO(a), BigIntSOO::SOO(b)) => {
                if let Some(result) = f(a, b) {
                    return BigIntSOO::SOO(result);
                }
                g(BigInt::from(a), BigInt::from(b))
            }
        };
        match result.to_int() {
            Ok(result) => {
                return BigIntSOO::SOO(result);
            },
            Err(()) => {
                return BigIntSOO::BigInt(result);
            }
        }
    }

    fn operation_ref_fst<F, G>(f: F, g: G, lhs: &BigIntSOO, rhs: BigIntSOO) -> BigIntSOO
        where F: FnOnce(i128, i128) -> Option<i128>, G: FnOnce(&BigInt, BigInt) -> BigInt
    {
        let result = match (lhs, rhs) {
            (BigIntSOO::BigInt(a), BigIntSOO::BigInt(b)) => g(a, b),
            (BigIntSOO::BigInt(a), BigIntSOO::SOO(b)) => g(a, BigInt::from(b)),
            (BigIntSOO::SOO(a), BigIntSOO::BigInt(b)) => g(&BigInt::from(*a), b),
            (BigIntSOO::SOO(a), BigIntSOO::SOO(b)) => {
                if let Some(result) = f(*a, b) {
                    return BigIntSOO::SOO(result);
                }
                g(&BigInt::from(*a), BigInt::from(b))
            }
        };
        match result.to_int() {
            Ok(result) => {
                return BigIntSOO::SOO(result);
            },
            Err(()) => {
                return BigIntSOO::BigInt(result);
            }
        }
    }

    fn operation_ref_snd<F, G>(f: F, g: G, lhs: BigIntSOO, rhs: &BigIntSOO) -> BigIntSOO
        where F: FnOnce(i128, i128) -> Option<i128>, G: FnOnce(BigInt, &BigInt) -> BigInt
    {
        let result = match (lhs, rhs) {
            (BigIntSOO::BigInt(a), BigIntSOO::BigInt(b)) => g(a, b),
            (BigIntSOO::BigInt(a), BigIntSOO::SOO(b)) => g(a, &BigInt::from(*b)),
            (BigIntSOO::SOO(a), BigIntSOO::BigInt(b)) => g(BigInt::from(a), b),
            (BigIntSOO::SOO(a), BigIntSOO::SOO(b)) => {
                if let Some(result) = f(a, *b) {
                    return BigIntSOO::SOO(result);
                }
                g(BigInt::from(a), &BigInt::from(*b))
            }
        };
        match result.to_int() {
            Ok(result) => {
                return BigIntSOO::SOO(result);
            },
            Err(()) => {
                return BigIntSOO::BigInt(result);
            }
        }
    }

    fn operation_ref_ref<F, G>(f: F, g: G, lhs: &BigIntSOO, rhs: &BigIntSOO) -> BigIntSOO
        where F: FnOnce(i128, i128) -> Option<i128>, G: FnOnce(&BigInt, &BigInt) -> BigInt
    {
        let result = match (lhs, rhs) {
            (BigIntSOO::BigInt(a), BigIntSOO::BigInt(b)) => g(a, b),
            (BigIntSOO::BigInt(a), BigIntSOO::SOO(b)) => g(a, &BigInt::from(*b)),
            (BigIntSOO::SOO(a), BigIntSOO::BigInt(b)) => g(&BigInt::from(*a), b),
            (BigIntSOO::SOO(a), BigIntSOO::SOO(b)) => {
                if let Some(result) = f(*a, *b) {
                    return BigIntSOO::SOO(result);
                }
                g(&BigInt::from(*a), &BigInt::from(*b))
            }
        };
        match result.to_int() {
            Ok(result) => {
                return BigIntSOO::SOO(result);
            },
            Err(()) => {
                return BigIntSOO::BigInt(result);
            }
        }
    }
}

const BASE_RING: BigIntRing = BigInt::RING;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BigIntSOORing;

impl RingBase for BigIntSOORing {

    type El = BigIntSOO;

    fn add_ref(&self, lhs: Self::El, rhs: &Self::El) -> Self::El {
        BigIntSOO::operation_ref_snd(|a, b| a.checked_add(b), |a, b| BASE_RING.add_ref(a, b), lhs, rhs)
    }

    fn mul_ref(&self, lhs: &Self::El, rhs: &Self::El) -> Self::El {
        BigIntSOO::operation_ref_ref(|a, b| a.checked_mul(b), |a, b| BASE_RING.mul_ref(a, b), lhs, rhs)
    }

    fn neg(&self, val: Self::El) -> Self::El {
        match val {
            BigIntSOO::BigInt(x) => BigIntSOO::BigInt(BASE_RING.neg(x)),
            BigIntSOO::SOO(x) => BigIntSOO::SOO(-x)
        }
    }

    fn zero(&self) -> Self::El {
        BigIntSOO::SOO(0)
    }

    fn one(&self) -> Self::El {
        BigIntSOO::SOO(1)
    }

    fn is_eq(&self, lhs: &Self::El, rhs: &Self::El) -> bool {
        match (lhs, rhs) {
            (BigIntSOO::BigInt(lhs), BigIntSOO::BigInt(rhs)) => BASE_RING.is_eq(lhs, rhs),
            (BigIntSOO::SOO(lhs), BigIntSOO::SOO(rhs)) => *lhs == *rhs,
            _ => false
        }
    }

    fn sub_ref_fst(&self, lhs: &Self::El, rhs: Self::El) -> Self::El {
        BigIntSOO::operation_ref_fst(|a, b| a.checked_sub(b), |a, b| BASE_RING.sub_ref_fst(a, b), lhs, rhs)
    }

    fn sub_ref_snd(&self, lhs: Self::El, rhs: &Self::El) -> Self::El {
        BigIntSOO::operation_ref_snd(|a, b| a.checked_sub(b), |a, b| BASE_RING.sub_ref_snd(a, b), lhs, rhs)
    }

    fn add(&self, lhs: Self::El, rhs: Self::El) -> Self::El {
        BigIntSOO::operation(|a, b| a.checked_add(b), |a, b| BASE_RING.add(a, b), lhs, rhs)
    }

    fn add_assign(&self, lhs: &mut Self::El, rhs: Self::El) { 
        let value = std::mem::replace(lhs, self.unspecified_element());
        *lhs = self.add(value, rhs);
    }

    fn add_assign_ref(&self, lhs: &mut Self::El, rhs: &Self::El) { 
        let value = std::mem::replace(lhs, self.unspecified_element());
        *lhs = self.add_ref(value, rhs);
    }

    fn add_assign_int(&self, lhs: &mut Self::El, rhs: i64) {
        self.add_assign(lhs, self.from_z(rhs));
    }

    fn sum<I>(&self, data: I) -> Self::El
        where I: Iterator<Item = Self::El>
    {
        let mut sum_i128: i128 = 0;
        let mut sum_bigint: BigInt = BigInt::ZERO;
        for x in data {
            match x {
                BigIntSOO::SOO(x) => {
                    if let Some(result) = sum_i128.checked_add(x) {
                        sum_i128 = result;
                    } else {
                        BASE_RING.add_assign(&mut sum_bigint, BigInt::from(sum_i128));
                        BASE_RING.add_assign(&mut sum_bigint, BigInt::from(x));
                        sum_i128 = 0;
                    }
                },
                BigIntSOO::BigInt(x) => {
                    BASE_RING.add_assign(&mut sum_bigint, x)
                }
            }
        }
        if BASE_RING.is_zero(&sum_bigint) {
            return BigIntSOO::SOO(sum_i128);
        } else {
            return BigIntSOO::BigInt(BASE_RING.add(sum_bigint, BigInt::from(sum_i128)));
        }
    }

    fn mul(&self, lhs: Self::El, rhs: Self::El) -> Self::El {
        BigIntSOO::operation(|a, b| a.checked_mul(b), |a, b| BASE_RING.mul(a, b), lhs, rhs)
    }

    fn mul_assign(&self, lhs: &mut Self::El, rhs: Self::El) { 
        let value = std::mem::replace(lhs, self.unspecified_element());
        *lhs = self.mul(value, rhs);
    }

    fn mul_assign_int(&self, lhs: &mut Self::El, rhs: i64) {
        self.mul_assign(lhs, self.from_z(rhs));
    }

    fn product<I>(&self, data: I) -> Self::El 
        where I: Iterator<Item = Self::El>
    {
        data.fold(self.one(), |a, b| self.mul(a, b))
    }

    fn sub(&self, lhs: Self::El, rhs: Self::El) -> Self::El {
        BigIntSOO::operation(|a, b| a.checked_sub(b), |a, b| BASE_RING.sub(a, b), lhs, rhs)
    }

    fn characteristic(&self) -> BigInt {
        BASE_RING.characteristic()
    }

    fn is_integral(&self) -> RingPropValue {
        BASE_RING.is_integral()
    }

    fn is_field(&self) -> RingPropValue {
        BASE_RING.is_field()
    }

    fn is_noetherian(&self) -> bool {
        BASE_RING.is_noetherian()
    }

    fn div(&self, _: Self::El, _: &Self::El) -> Self::El {
        panic!("Not a field!")
    }

    fn from_z_big(&self, x: &BigInt) -> Self::El {
        BigIntSOO::BigInt(x.clone())
    }

    fn from_z(&self, x: i64) -> Self::El {
        BigIntSOO::SOO(x as i128)
    }

    fn format(&self, el: &Self::El, f: &mut std::fmt::Formatter, in_prod: bool) -> std::fmt::Result {
        match el {
            BigIntSOO::SOO(x) => write!(f, "{:?}", x),
            BigIntSOO::BigInt(x) => BASE_RING.format(x, f, in_prod)
        }
    }
}