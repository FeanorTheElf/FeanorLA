use super::*;
use super::bigint::*;
use super::super::ring::*;

use std::convert::TryFrom;

///
/// A big integer with small object optimization. More concretely, as long as the stored
/// integers fit into 128 bits, no dynamic memory allocation is used.
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

    pub const RING: BigIntSOORing = BigIntSOORing {};

    pub fn to_i128(&self) -> Result<i128, ()> {
        match self {
            BigIntSOO::BigInt(x) => x.to_int(),
            BigIntSOO::SOO(x) => Ok(*x)
        }
    }
    
    pub fn to_i64(&self) -> Result<i64, ()> {
        self.to_i128().and_then(|x| i64::try_from(x).map_err(|_| ()))
    }

    pub fn to_bigint(self) -> BigInt {
        match self {
            BigIntSOO::BigInt(x) => x,
            BigIntSOO::SOO(x) => BigInt::from(x)
        }
    }

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
        return Self::from(result);
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
        return Self::from(result);
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
        return Self::from(result);
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
        return Self::from(result);
    }
}

impl From<BigInt> for BigIntSOO {
    
    fn from(x: BigInt) -> BigIntSOO {
        match x.to_int() {
            Ok(result) => {
                return BigIntSOO::SOO(result);
            },
            Err(()) => {
                return BigIntSOO::BigInt(x);
            }
        }
    }
}

impl From<i128> for BigIntSOO {
    
    fn from(x: i128) -> BigIntSOO {
        if x == i128::MIN {
            return BigIntSOO::BigInt(BigInt::from(x));
        } else {
            return BigIntSOO::SOO(x)
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

    fn mul_assign_int(&self, lhs: &mut Self::El, rhs: i64) {
        self.mul_assign(lhs, &self.from_z(rhs));
    }

    fn product<I>(&self, data: I) -> Self::El 
        where I: Iterator<Item = Self::El>
    {
        data.fold(self.one(), |a, b| self.mul(a, b))
    }

    fn sub(&self, lhs: Self::El, rhs: Self::El) -> Self::El {
        BigIntSOO::operation(|a, b| a.checked_sub(b), |a, b| BASE_RING.sub(a, b), lhs, rhs)
    }

    fn characteristic(&self) -> StdInt {
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

    fn from_z_big(&self, x: &StdInt) -> Self::El {
        x.into_val()
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

impl SingletonRing for BigIntSOORing {
    
    fn singleton() -> Self {
        BigIntSOO::RING
    }
}

impl CanonicalEmbeddingInfo<BigIntSOORing> for BigIntSOORing {

    fn has_embedding(&self, _from: &BigIntSOORing) -> RingPropValue {
        RingPropValue::True
    }

    fn embed(&self, _from: &BigIntSOORing, el: BigIntSOO) -> Self::El {
        el
    }
}

impl CanonicalIsomorphismInfo<BigIntSOORing> for BigIntSOORing {

    fn has_isomorphism(&self, _from: &BigIntSOORing) -> RingPropValue {
        RingPropValue::True
    }

    fn preimage(&self, _from: &BigIntSOORing, el: BigIntSOO) -> Self::El {
        el
    }
}

impl CanonicalEmbeddingInfo<StaticRing<i64>> for BigIntSOORing {

    fn has_embedding(&self, _from: &StaticRing<i64>) -> RingPropValue {
        RingPropValue::True
    }

    fn embed(&self, _from: &StaticRing<i64>, el: i64) -> BigIntSOO {
        BigIntSOO::RING.from_z(el)
    }
}

impl CanonicalIsomorphismInfo<StaticRing<i64>> for BigIntSOORing {

    fn has_isomorphism(&self, _from: &StaticRing<i64>) -> RingPropValue {
        RingPropValue::True
    }

    fn preimage(&self, _from: &StaticRing<i64>, el: BigIntSOO) -> i64 {
        el.to_i128().expect("Overflow when embedding BigInt into i128") as i64
    }
}

impl CanonicalEmbeddingInfo<BigIntRing> for BigIntSOORing {

    fn has_embedding(&self, _from: &BigIntRing) -> RingPropValue {
        RingPropValue::True
    }

    fn embed(&self, _from: &BigIntRing, el: BigInt) -> Self::El {
        BigIntSOO::from(el)
    }
}

impl CanonicalIsomorphismInfo<BigIntRing> for BigIntSOORing {

    fn has_isomorphism(&self, _from: &BigIntRing) -> RingPropValue {
        RingPropValue::True
    }

    fn preimage(&self, _from: &BigIntRing, el: BigIntSOO) -> BigInt {
        el.to_bigint()
    }
}

impl OrderedRing for BigIntSOORing {

    fn cmp(&self, lhs: &Self::El, rhs: &Self::El) -> std::cmp::Ordering {
        match (lhs, rhs) {
            (BigIntSOO::BigInt(lhs), BigIntSOO::BigInt(rhs)) => BigInt::RING.cmp(lhs, rhs),
            (BigIntSOO::SOO(lhs), BigIntSOO::SOO(rhs)) => lhs.cmp(rhs),
            (BigIntSOO::SOO(_), BigIntSOO::BigInt(_)) => BigIntSOO::RING.cmp(rhs, lhs).reverse(),
            (BigIntSOO::BigInt(lhs), BigIntSOO::SOO(rhs)) => match (BigInt::RING.is_neg(lhs), *rhs < 0) {
                (true, true) => lhs.abs_cmp_small((-rhs) as u128).reverse(),
                (false, false) => lhs.abs_cmp_small(*rhs as u128),
                (true, false) => Ordering::Less,
                (false, true) => Ordering::Greater
            }
        }
    }
}

impl DivisibilityInfoRing for BigIntSOORing {
    
    fn is_divisibility_computable(&self) -> RingPropValue {
        RingPropValue::True
    }

    fn quotient(&self, lhs: &Self::El, rhs: &Self::El) -> Option<Self::El> {
        let (quo, rem) = self.euclidean_div_rem(lhs.clone(), rhs);
        if self.is_zero(&rem) {
            return Some(quo);
        } else {
            return None;
        }
    }
}

impl EuclideanInfoRing for BigIntSOORing {
    
    fn is_euclidean(&self) -> RingPropValue {
        RingPropValue::True
    }

    fn euclidean_div_rem(&self, lhs: Self::El, rhs: &Self::El) -> (Self::El, Self::El) {
        match (lhs, rhs.to_i64()) {
            (BigIntSOO::SOO(lhs), Ok(rhs)) => {
                let (quo, rem) = i128::RING.euclidean_div_rem(lhs, &(rhs as i128));
                return (BigIntSOO::SOO(quo), BigIntSOO::SOO(rem));
            },
            (BigIntSOO::BigInt(lhs), Ok(rhs)) => {
                let (quo, rem) = lhs.euclidean_div_rem_small(rhs);
                return (BigIntSOO::from(quo), BigIntSOO::SOO(rem as i128));
            },
            (lhs, _) => {
                let (quo, rem) = BigInt::RING.euclidean_div_rem(lhs.to_bigint(), &rhs.clone().to_bigint());
                return (BigIntSOO::from(quo), BigIntSOO::from(rem));
            }
        };
    }

    fn euclidean_deg(&self, el: Self::El) -> StdInt {
        BigIntSOO::WRAPPED_RING.from(self.abs(el))
    }
}

impl IntegerRing for BigIntSOORing {

    fn to_float_approx(&self, el: &Self::El) -> f64 {
        match el {
            BigIntSOO::BigInt(x) => BigInt::RING.to_float_approx(x),
            BigIntSOO::SOO(x) => *x as f64
        }
    }

    fn from_float_approx(&self, el: f64) -> Option<Self::El> {
        BigInt::RING.from_float_approx(el).map(BigIntSOO::from)
    }

    fn mul_pow_2(&self, el: El<Self>, power: u64) -> El<Self> {
        match el {
            BigIntSOO::BigInt(x) => BigIntSOO::from(BigInt::RING.mul_pow_2(x, power)),
            BigIntSOO::SOO(x) if power < x.leading_zeros() as u64 => BigIntSOO::SOO(x << power),
            BigIntSOO::SOO(x) => BigIntSOO::from(BigInt::RING.mul_pow_2(BigInt::from(x), power))
        }
    }

    fn euclidean_div_pow_2(&self, el: El<Self>, power: u64) -> El<Self> {
        match el {
            BigIntSOO::BigInt(x) => BigIntSOO::from(BigInt::RING.euclidean_div_pow_2(x, power)),
            BigIntSOO::SOO(x) => BigIntSOO::SOO(x.div_euclid(1 << power))
        }
    }

    fn is_odd(&self, el: &Self::El) -> bool {
        match el {
            BigIntSOO::BigInt(x) => BigInt::RING.is_odd(x),
            BigIntSOO::SOO(x) => x % 2 == 0
        }
    }

    fn abs_log2_floor(&self, el: &El<Self>) -> u64 {
        match el {
            BigIntSOO::BigInt(x) => BigInt::RING.abs_log2_floor(x),
            BigIntSOO::SOO(x) => (i128::BITS - 1 - x.abs().leading_zeros()) as u64
        }
    }

    fn abs_is_bit_set(&self, el: &El<Self>, bit: u64) -> bool {
        match el {
            BigIntSOO::BigInt(x) => BigInt::RING.abs_is_bit_set(x, bit),
            BigIntSOO::SOO(x) => (x.abs() >> bit) & 1 == 1
        }
    }

    fn get_uniformly_random<G>(
        &self,
        rng: G, 
        end_exclusive: &El<Self>
    ) -> El<Self> 
        where G: FnMut() -> u32
    {
        BigIntSOO::from(BigInt::RING.get_uniformly_random(rng, &end_exclusive.clone().to_bigint()))
    }

    fn highest_dividing_power_of_two(&self, el: &El<Self>) -> usize {
        match el {
            BigIntSOO::SOO(el) => el.trailing_zeros() as usize,
            BigIntSOO::BigInt(el) => BigInt::RING.highest_dividing_power_of_two(el)
        }
    }
}

#[test]
fn test_bigint_assumption() {
    let max = BigInt::from(i128::MAX);
    assert_eq!(Ok(i128::MAX), max.to_int());
    assert_eq!(Ok(-i128::MAX), BigInt::RING.neg(max.clone()).to_int());
    assert_eq!(Err(()), BigInt::RING.add_ref(BigInt::RING.one(), &max).to_int());
    assert_eq!(Err(()), BigInt::RING.neg(BigInt::RING.add(max, BigInt::RING.one())).to_int());
}

#[test]
fn test_sum() {
    let values = [
        BigInt::from(1 << 126),
        BigInt::from(1),
        BigInt::from(1),
        BigInt::from(1 << 126),
        BigInt::from(1),
        BigInt::from(1 << 126),
        BigInt::from(-4)
    ];
    let expected = BigIntSOO::from(BigInt::from_str_radix("BFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF", 16).unwrap());
    let actual = BigIntSOO::RING.sum(values.to_vec().into_iter().map(BigIntSOO::from));
    assert!(BigIntSOO::RING.is_eq(&expected, &actual));
}

#[test]
fn test_cmp() {
    assert_eq!(Ordering::Less, BigIntSOO::RING.cmp(
        &BigIntSOO::RING.from_z(-2), &BigIntSOO::RING.from_z(-1)
    ));
    assert_eq!(Ordering::Less, BigIntSOO::RING.cmp(
        &BigIntSOO::from(BigInt::from_str_radix("-FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF", 16).unwrap()), 
        &BigIntSOO::RING.from_z(1)
    ));
    assert_eq!(Ordering::Less, BigIntSOO::RING.cmp(
        &BigIntSOO::from(BigInt::from_str_radix("-FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF", 16).unwrap()), 
        &BigIntSOO::RING.from_z(-1)
    ));
    assert_eq!(Ordering::Less, BigIntSOO::RING.cmp(
        &BigIntSOO::from(BigInt::from_str_radix("-FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF", 16).unwrap()), 
        &BigIntSOO::from(BigInt::from_str_radix("FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF", 16).unwrap()), 
    ));
    assert_eq!(Ordering::Less, BigIntSOO::RING.cmp(
        &BigIntSOO::from(BigInt::from_str_radix("-FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF", 16).unwrap()), 
        &BigIntSOO::from(BigInt::from_str_radix("-FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFE", 16).unwrap()), 
    ));
}

#[test]
fn test_mul_pow_2() {
    assert!(BigIntSOO::RING.is_eq(
        &BigIntSOO::SOO(11 * 8), 
        &BigIntSOO::RING.mul_pow_2(BigIntSOO::RING.from_z(11), 3) 
    ));
    assert!(BigIntSOO::RING.is_eq(
        &BigIntSOO::BigInt(BigInt::from_str_radix("B00000000000000000000000000000000", 16).unwrap()), 
        &BigIntSOO::RING.mul_pow_2(BigIntSOO::RING.from_z(11), 128) 
    ));
    assert!(BigIntSOO::RING.is_eq(
        &BigIntSOO::BigInt(BigInt::from_str_radix("-B00000000000000000000000000000000", 16).unwrap()), 
        &BigIntSOO::RING.mul_pow_2(BigIntSOO::RING.from_z(-11), 128) 
    ));
}

#[test]
#[ignore]
fn test_size() {
    assert_eq!(std::mem::size_of::<BigInt>() + 8, std::mem::size_of::<BigIntSOO>());
}