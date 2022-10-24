use super::super::ring::*;
use super::super::embedding::*;
use super::super::primitive::*;
use super::super::la::vec::*;
use super::bigint_ops;

use super::primes;
use super::*;
use super::bigint_soo::*;

use std::cmp::Ordering;
use std::ops::*;
use std::hash::{Hash, Hasher};

#[derive(Debug, Clone)]
pub struct BigInt {
    /// "digits" of the (absolute value) of the number in base 2^64
    data: Vector<Vec<u64>, u64>,
    negative: bool
}

///
/// 10 to this power fits still in a u64
/// 
const BIG_POWER_TEN_ZEROS: u32 = 19;
const BIG_POWER_TEN: u64 = 10000000000000000000u64;

impl BigInt {

    const BLOCK_BITS: usize = 64;
    pub const ZERO: BigInt = BigInt { data: Vector::new(Vec::new()), negative: false };

    fn highest_set_block(&self) -> Option<usize> {
        bigint_ops::highest_set_block(self.data.as_ref())
    }

    pub fn assign(&mut self, rhs: &BigInt) {
        self.negative = rhs.negative;
        let mut data_vec = std::mem::replace(&mut self.data, Vector::new(Vec::new())).raw_data();
        data_vec.clear();
        if let Some(d) = rhs.highest_set_block() {
            if data_vec.len() < d {
                data_vec.reserve(d - self.data.len());
            }
            for i in 0..=d {
                data_vec.push(rhs.data[i]);
            }
        }
        self.data = Vector::new(data_vec);
    }

    ///
    /// calculates the number with the given representation w.r.t the given base.
    /// The passed iterator must yield the "digits" in the order starting with the
    /// highest significant digit.
    /// 
    pub fn from_radix<I, E>(data: I, base: u64) -> Result<BigInt, E> 
        where I: Iterator<Item = Result<u64, E>>
    {
        let mut result = BigInt {
            data: Vector::new(Vec::with_capacity(data.size_hint().0)),
            negative: false
        };
        for value in data {
            let val = value?;
            debug_assert!(val < base);
            bigint_ops::bigint_mul_small(&mut result.data, base);
            bigint_ops::bigint_add_small(&mut result.data, val);
        }
        return Ok(result);
    }

    pub fn from_str_radix(s: &str, base: u32) -> Result<BigInt, BigIntParseError> {
        assert!(base >= 2);
        let sign = s.chars().next().unwrap();
        let (negative, rest): (bool, &[u8]) = if sign == '+' {
            (true, &<str as AsRef<[u8]>>::as_ref(s)[1..])
        } else if sign == '-' {
            (true, &<str as AsRef<[u8]>>::as_ref(s)[1..])
        } else {
            (false, <str as AsRef<[u8]>>::as_ref(s))
        };
        // we need the -1 in Self::BLOCK_BITS to ensure that base^chunk_size is 
        // really smaller than 2^64
        let chunk_size = ((Self::BLOCK_BITS - 1) as f32 / (base as f32).log2()).floor() as usize;
        let it = rest.rchunks(chunk_size as usize).rev()
            .map(std::str::from_utf8)
            .map(|chunk| chunk.map_err(BigIntParseError::from))
            .map(|chunk| chunk.and_then(|n| 
                u64::from_str_radix(n, base).map_err(BigIntParseError::from))
            );
        let mut result = Self::from_radix(it, (base as u64).pow(chunk_size as u32));
        if let Ok(r) = &mut result {
            r.negative = negative;
        }
        return result;
    }

    fn euclidean_div_rem(&mut self, rhs: &BigInt) -> BigInt {
        let quotient = bigint_ops::bigint_div(&mut self.data, rhs.data.as_ref());
        return BigInt { data: quotient, negative: self.negative ^ rhs.negative };
    }

    ///
    /// Performs euclidean division, i.e. computes q, r such that
    /// `self = q * rhs + r` and `|r| < |rhs|`.
    /// Returns `(q, r)`
    /// 
    pub fn euclidean_div_rem_small(mut self, rhs: i64) -> (BigInt, i64) {
        let mut remainder = bigint_ops::bigint_div_small(&mut self.data, rhs.abs() as u64) as i64;
        if self.negative {
            remainder = -remainder;
        }
        self.negative = self.negative ^ (rhs < 0);
        return (self, remainder);
    }

    ///
    /// Performs floor division, i.e. computes the result of
    /// self/rhs rounded towards -inf.
    /// 
    /// # Example
    /// 
    /// ```
    /// # use feanor_la::integer::bigint::BigInt;
    /// assert_eq!(BigInt::from(-1), BigInt::from(-1).floor_div_small(2));
    /// ```
    /// 
    pub fn floor_div_small(self, rhs: i64) -> BigInt {
        let (mut quo, rem) = self.euclidean_div_rem_small(rhs);
        if rem < 0 {
            BigInt::RING.add_assign_int(&mut quo, if rhs < 0 { 1 } else { -1 });
        }
        return quo;
    }

    pub fn abs_cmp_small(&self, rhs: u128) -> Ordering {
        bigint_ops::bigint_cmp_small(self.data.as_ref(), rhs)
    }

    ///
    /// Returns `Ok(x)` if `|x| <= i128::MAX` and `Err(())` otherwise. In particular,
    /// the output is `Ok()` if and only if x fits into `i128`, except in the case 
    /// `i128::MIN` (note that `|i128::MIN| = i128::MAX + 1`).
    /// 
    /// This decision was made, as the symmetry seems more useful in applications, in 
    /// particular, we can now always compute `-x.to_int().unwrap()`.
    /// 
    pub fn to_int(&self) -> Result<i128, ()> {
        if let Some(d) = self.highest_set_block() {
            let result = if d == 1 && self.data[1] <= i64::MAX as u64 && self.negative {
                -(((self.data[1] as i128) << BigInt::BLOCK_BITS) | (self.data[0] as i128))
            } else if d == 1 && self.data[1] <= i64::MAX as u64 && !self.negative {
                ((self.data[1] as i128) << BigInt::BLOCK_BITS) | (self.data[0] as i128)
            } else if d == 0 && self.negative {
                -(self.data[0] as i128)
            } else if d == 0 && !self.negative {
                self.data[0] as i128
            } else {
                return Err(());
            };
            assert!(result.abs() <= i128::MAX);
            return Ok(result);
        } else {
            return Ok(0);
        }
    }

    ///
    /// Least significant digit first
    /// 
    pub fn base_u64_repr(self) -> Vector<Vec<u64>, u64> {
        self.data
    }

    ///
    /// Least significant digit first
    /// 
    pub fn from_base_u64_repr<V>(data: Vector<V, u64>) -> Self 
        where V: VectorView<u64>
    {
        BigInt { 
            data: data.into_owned(), 
            negative: false
        }
    }
}

impl From<i128> for BigInt {

    fn from(val: i128) -> BigInt {
        if val == i128::MIN {
            BigInt {
                negative: true,
                data: Vector::new(vec![0, 1 << 63])
            }
        } else if val < 0 {
            let val = (-val) as u128;
            BigInt {
                negative: true,
                data: Vector::new(vec![(val & ((1 << BigInt::BLOCK_BITS) - 1)) as u64, (val >> BigInt::BLOCK_BITS) as u64])
            }
        } else {
            let val = val as u128;
            BigInt {
                negative: false,
                data: Vector::new(vec![(val & ((1 << BigInt::BLOCK_BITS) - 1)) as u64, (val >> BigInt::BLOCK_BITS) as u64])
            }
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BigIntRing;

impl BigInt {
    pub const RING: BigIntRing = BigIntRing {};
    pub const WRAPPED_RING: WrappingRing<BigIntRing> = WrappingRing::new(Self::RING);
}

impl RingBase for BigIntRing {

    type El = BigInt;

    fn add_ref(&self, mut lhs: Self::El, rhs: &Self::El) -> Self::El {
        if lhs.negative == rhs.negative {
            bigint_ops::bigint_add(&mut lhs.data, rhs.data.as_ref(), 0);
        } else if bigint_ops::bigint_cmp(lhs.data.as_ref(), rhs.data.as_ref()) != Ordering::Less {
            bigint_ops::bigint_sub(&mut lhs.data, rhs.data.as_ref(), 0);
        } else {
            let mut result = rhs.clone();
            bigint_ops::bigint_sub(&mut result.data, lhs.data.as_ref(), 0);
            std::mem::swap(&mut lhs, &mut result);
        }
        return lhs;
    }

    fn mul_ref(&self, lhs: &Self::El, rhs: &Self::El) -> Self::El {
        let result = bigint_ops::bigint_mul(lhs.data.as_ref(), rhs.data.as_ref());
        return BigInt {
            negative: lhs.negative ^ rhs.negative,
            data: result
        };
    }

    fn mul(&self, lhs: Self::El, rhs: Self::El) -> Self::El {
        self.mul_ref(&lhs, &rhs)
    }

    fn add_assign_int(&self, lhs: &mut BigInt, rhs: i64) {
        if lhs.negative == (rhs < 0) {
            bigint_ops::bigint_add_small(&mut lhs.data, rhs.abs() as u64);
        } else {
            let rhs_bigint = self.from_z(rhs);
            BigInt::RING.add_assign(lhs, rhs_bigint);
        }
    }

    fn mul_assign_int(&self, lhs: &mut Self::El, rhs: i64) {
        bigint_ops::bigint_mul_small(&mut lhs.data, rhs.abs() as u64);
        lhs.negative ^= rhs < 0;
    }

    fn neg(&self, mut val: Self::El) -> Self::El {
        val.negative = !val.negative;
        return val;
    }

    fn zero(&self) -> Self::El {
        return BigInt::ZERO.clone();
    }

    fn one(&self) -> Self::El {
        return BigInt {
            negative: false,
            data: Vector::new(vec![1])
        };
    }

    fn is_eq(&self, lhs: &Self::El, rhs: &Self::El) -> bool {
        if self.is_zero(lhs) && self.is_zero(rhs) {
            return true;
        } else if lhs.negative != rhs.negative {
            return false;
        }
        let highest_block = lhs.highest_set_block();
        if highest_block != rhs.highest_set_block() {
            return false;
        }
        if let Some(d) = highest_block {
            for i in 0..=d {
                if lhs.data[i] != rhs.data[i] {
                    return false;
                }
            }
        }
        return true;
    }

    fn is_zero(&self, val: &Self::El) -> bool {
        val.highest_set_block() == None
    }

    fn characteristic(&self) -> StdInt { StdInt::zero() }
    fn is_noetherian(&self) -> bool { true }
    fn is_integral(&self) -> RingPropValue { RingPropValue::True }
    fn is_field(&self) -> RingPropValue { RingPropValue::False }

    fn div(&self, _lhs: Self::El, _rhs: &Self::El) -> Self::El {
        panic!("Not a field!")
    }

    fn from_z(&self, x: i64) -> BigInt {
        BigInt::from(x as i128)
    }

    fn from_z_big(&self, x: &StdInt) -> BigInt {
        x.clone().into_val().to_bigint()
    }

    fn format(&self, el: &BigInt, f: &mut std::fmt::Formatter, _in_prod: bool) -> std::fmt::Result {
        if el.negative {
            write!(f, "-")?;
        }
        let mut copy = el.clone();
        let mut remainders: Vec<u64> = Vec::with_capacity(
            (el.highest_set_block().unwrap_or(0) + 1) * BigInt::BLOCK_BITS / 3
        );
        while !self.is_zero(&copy) {
            let rem = bigint_ops::bigint_div_small(&mut copy.data, BIG_POWER_TEN);
            remainders.push(rem);
        }
        remainders.reverse();
        let mut it = remainders.into_iter();
        if let Some(fst) = it.next() {
            write!(f, "{}", fst)?;
            for rem in it {
                write!(f, "{:0>width$}", rem, width = BIG_POWER_TEN_ZEROS as usize)?;
            }
        } else {
            write!(f, "0")?;
        }
        return Ok(());
    }
}

impl OrderedRing for BigIntRing {

    fn cmp(&self, lhs: &Self::El, rhs: &Self::El) -> std::cmp::Ordering {
        if self.is_zero(lhs) && self.is_zero(rhs) {
            Ordering::Equal
        } else {
            match (lhs.negative, rhs.negative) {
                (true, true) => bigint_ops::bigint_cmp(rhs.data.as_ref(), lhs.data.as_ref()),
                (true, false) => Ordering::Less,
                (false, true) => Ordering::Greater,
                (false, false) => bigint_ops::bigint_cmp(lhs.data.as_ref(), rhs.data.as_ref())
            }
        }
    }
}

impl HashableElRing for BigIntRing {

    fn hash<H: std::hash::Hasher>(&self, h: &mut H, el: &BigInt) {
        el.hash(h)
    }
}

impl CanonicalEmbeddingInfo<BigIntRing> for BigIntRing {

    fn has_embedding(&self, _from: &BigIntRing) -> RingPropValue {
        RingPropValue::True
    }

    fn embed(&self, _from: &BigIntRing, el: BigInt) -> Self::El {
        el
    }
}

impl CanonicalIsomorphismInfo<BigIntRing> for BigIntRing {

    fn has_isomorphism(&self, _from: &BigIntRing) -> RingPropValue {
        RingPropValue::True
    }

    fn preimage(&self, _from: &BigIntRing, el: BigInt) -> Self::El {
        el
    }
}

impl EuclideanInfoRing for BigIntRing {

    fn is_euclidean(&self) -> RingPropValue {
        RingPropValue::True
    }

    fn euclidean_div_rem(&self, mut lhs: Self::El, rhs: &Self::El) -> (Self::El, Self::El) {
        let quo = lhs.euclidean_div_rem(rhs);
        (quo, lhs)
    }

    fn euclidean_deg(&self, el: Self::El) -> StdInt {
        StdInt::RING.wrap(BigIntSOO::from(el)).abs()
    }
}

impl PartialEq<i64> for BigInt {

    fn eq(&self, rhs: &i64) -> bool {
        if let Some(d) = self.highest_set_block() {
            (*rhs < 0 && self.negative && d == 0 && self.data[0] == (-*rhs) as u64) ||
            (*rhs > 0 && !self.negative && d == 0 && self.data[0] == *rhs as u64)
        } else {
            *rhs == 0
        }
    }
}

impl PartialEq<BigInt> for BigInt {

    fn eq(&self, rhs: &BigInt) -> bool {
        BigInt::RING.is_eq(self, rhs)
    }
}

impl PartialOrd<i64> for BigInt {

    fn partial_cmp(&self, rhs: &i64) -> Option<Ordering> {
        if BigInt::RING.is_zero(self) && *rhs == 0 {
            Some(Ordering::Equal)
        } else if self.negative && *rhs < 0 {
            Some(bigint_ops::bigint_cmp_small(self.data.as_ref(), (-*rhs) as u128).reverse())
        } else if self.negative && *rhs >= 0 {
            Some(Ordering::Less)
        } else if !self.negative && *rhs < 0 {
            Some(Ordering::Greater)
        } else {
            Some(bigint_ops::bigint_cmp_small(self.data.as_ref(), (*rhs) as u128))
        }
    }
}

impl SingletonRing for BigIntRing {
    fn singleton() -> BigIntRing {
        BigInt::RING
    }
}

impl std::fmt::Display for BigInt {

    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        BigInt::RING.format(self, f, false)
    }
}

impl Eq for BigInt {}

impl Hash for BigInt {

    fn hash<H: Hasher>(&self, hasher: &mut H) {
        if let Some(d) = self.highest_set_block() {
            self.negative.hash(hasher);
            for i in 0..=d {
                hasher.write_u64(self.data[i])
            }
        }
    }
}

#[derive(Debug, Clone)]
pub enum BigIntParseError {
    ParseIntError(std::num::ParseIntError),
    Utf8Error(std::str::Utf8Error)
}

impl From<std::num::ParseIntError> for BigIntParseError {
    fn from(e: std::num::ParseIntError) -> BigIntParseError {
        BigIntParseError::ParseIntError(e)
    }
}

impl From<std::str::Utf8Error> for BigIntParseError {
    fn from(e: std::str::Utf8Error) -> BigIntParseError {
        BigIntParseError::Utf8Error(e)
    }
}

impl std::str::FromStr for BigInt {
    type Err = BigIntParseError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        Self::from_str_radix(s, 10)
    }
}

impl DivisibilityInfoRing for BigIntRing {

    fn quotient(&self, lhs: &Self::El, rhs: &Self::El) -> Option<BigInt> {
        if self.is_zero(rhs) && !self.is_zero(lhs) {
            return None;
        } else if self.is_zero(rhs) && self.is_zero(lhs) {
            return Some(self.one());
        }
        let (quo, rem) = self.euclidean_div_rem(lhs.clone(), rhs);
        if rem == 0 {
            return Some(quo);
        } else {
            return None;
        }
    }

    fn is_divisibility_computable(&self) -> RingPropValue {
        RingPropValue::True
    }

    fn is_unit(&self, el: &Self::El) -> bool {
        self.is_one(el) || self.is_neg_one(el)
    }
}

impl UfdInfoRing for BigIntRing {

    fn is_ufd(&self) -> RingPropValue {
        RingPropValue::True
    }

    fn is_prime(&self, el: &Self::El) -> bool {
        primes::miller_rabin(&BigInt::RING, el, 10)
    }

    fn calc_factor(&self, el: &Self::El) -> Option<Self::El> {
        primes::calc_factor(&StdInt::RING.wrap(self.preimage(&BigIntSOO::RING, el.clone())))
            .map(|x| self.embed(&BigIntSOO::RING, x.into_val()))
    }
}

impl IntegerRing for BigIntRing {

    ///
    /// Returns the float that is closest to the integer. Note that 
    /// if for very big numbers (with abs() in the order of magnitude 
    /// 2^1024 or greater, this can even yield infinity)
    /// 
    fn to_float_approx(&self, el: &Self::El) -> f64 {
        if let Some(d) = el.highest_set_block() {
            let mut upper_part = el.data[d] as f64 * 
                2f64.powi(BigInt::BLOCK_BITS as i32);
            if d > 0 {
                upper_part += el.data[d - 1] as f64;
            }
            return upper_part * 2f64.powi(BigInt::BLOCK_BITS as i32).powi(d as i32 - 1)
        } else {
            return 0.;
        }
    }

    ///
    /// Returns a BigInt that has the given value, rounded. Note that
    /// for very big numbers, the float representation can be very imprecise.
    /// For Infinity and NaN, nothing is returned;
    /// 
    fn from_float_approx(&self, val: f64) -> Option<BigInt> {
        if val.is_infinite() || val.is_nan() {
            return None;
        } else if val.abs() <= 0.5 {
            return Some(BigInt::ZERO);
        } else if val.abs() <= 1.5 {
            if val.is_sign_negative() { 
                return Some(BigInt::RING.neg(BigInt::RING.one()));
            } else { 
                return Some(BigInt::RING.one());
            }
        } else {
            const MANTISSA: i32 = 52;
            let exp = std::cmp::max(val.abs().log2() as i32, MANTISSA) - MANTISSA;
            let int = (val.abs() / 2f64.powi(exp)).trunc() as u64;
            let blocks = exp as usize / BigInt::BLOCK_BITS;
            let within_block_shift = exp as usize % BigInt::BLOCK_BITS;
            let mut result = (0..blocks).map(|_| 0).collect::<Vec<_>>();
            result.push(int << within_block_shift);
            if within_block_shift != 0 {
                result.push(int >> (BigInt::BLOCK_BITS - within_block_shift));
            }
            return Some(BigInt {
                negative: val.is_sign_negative(),
                data: Vector::new(result)
            });
        }
    }

    fn mul_pow_2(&self, el: El<Self>, power: u64) -> El<Self> { 
        let add_blocks = power as usize / BigInt::BLOCK_BITS;
        let shift_amount = power as usize % BigInt::BLOCK_BITS;
        let keep_bits = BigInt::BLOCK_BITS - shift_amount;
        let keep_mask = (1u64.checked_shl(keep_bits as u32).unwrap_or(0).wrapping_sub(1)) << shift_amount;
        let carry_mask = (1 << shift_amount) - 1;
        let mut carry = 0;
        let mut el_data = el.data.raw_data();
        for i in 0..el_data.len() {
            let rotated = el_data[i].rotate_left(shift_amount as u32);
            el_data[i] = (rotated & keep_mask) | carry;
            carry = rotated & carry_mask;
        }
        if carry != 0 {
            el_data.push(carry);
        }
        let mut tmp = (0..add_blocks).map(|_| 0).collect();
        std::mem::swap(&mut tmp, &mut el_data);
        el_data.extend(tmp.into_iter());
        return BigInt {
            negative: el.negative,
            data: Vector::new(el_data)
        };
    }

    fn euclidean_div_pow_2(&self, el: El<Self>, power: u64) -> El<Self> {
        let drop_blocks = power as usize / BigInt::BLOCK_BITS;
        let shift_amount = power as usize % BigInt::BLOCK_BITS;
        let keep_bits = BigInt::BLOCK_BITS - shift_amount;
        let mut el_data = el.data.raw_data();
        if el_data.len() <= drop_blocks {
            return BigInt::ZERO;
        } else {
            el_data.drain(0..drop_blocks);
            el_data[0] >>= shift_amount;
            for i in 1..el_data.len() {
                let rotated = el_data[i].rotate_right(shift_amount as u32);
                el_data[i] = rotated & (u64::MAX >> shift_amount);
                el_data[i - 1] |= rotated & u64::MAX.wrapping_shl(keep_bits as u32);
            }
            return BigInt {
                negative: el.negative,
                data: Vector::new(el_data)
            };
        }
    }

    fn abs_log2_floor(&self, el: &BigInt) -> u64 {
        assert!(!self.is_zero(el));
        let d = el.highest_set_block().unwrap();
        return BigInt::BLOCK_BITS as u64 - el.data[d].leading_zeros() as u64 - 1 + 
            d as u64 * BigInt::BLOCK_BITS as u64;
    }

    fn abs_is_bit_set(&self, val: &BigInt, i: u64) -> bool {
        if let Some(d) = val.highest_set_block() {
            let block = i as usize / BigInt::BLOCK_BITS;
            let bit = i as usize % BigInt::BLOCK_BITS;
            if block > d {
                false
            } else {
                ((val.data[block] >> bit) & 1) == 1
            }
        } else {
            false
        }
    }

    fn get_uniformly_random<G>(
        &self,
        mut rng: G, 
        end_exclusive: &BigInt
    ) -> BigInt 
        where G: FnMut() -> u32
    {
        assert!(*end_exclusive > 0);
        let k = BigInt::RING.abs_log2_floor(&end_exclusive) as usize + 1;
        let mut rng64 = || ((rng() as u64) << u32::BITS) | (rng() as u64);
        loop {
            let random_most_significant_bits = rng64() >> (BigInt::BLOCK_BITS - (k % BigInt::BLOCK_BITS));
            let data = (0..k/BigInt::BLOCK_BITS).map(|_| rng64())
                .chain(std::iter::once(random_most_significant_bits));
            let result = BigInt {
                data: Vector::new(data.collect()),
                negative: false
            };
            if self.cmp(&result, end_exclusive) == Ordering::Less {
                return result;
            }
        }
    }

    fn highest_dividing_power_of_two(&self, n: &BigInt) -> usize {
        if let Some(d) = n.highest_set_block() {
            for i in 0..=d {
                if n.data[i] != 0 {
                    return n.data[i].trailing_zeros() as usize + 
                        i * BigInt::BLOCK_BITS;
                }
            }
            unreachable!()
        } else {
            return 0;
        }
    }
}

impl CanonicalEmbeddingInfo<StaticRing<i64>> for BigIntRing {

    fn has_embedding(&self, _from: &StaticRing<i64>) -> RingPropValue {
        RingPropValue::True
    }

    fn embed(&self, _from: &StaticRing<i64>, el: i64) -> BigInt {
        BigInt::from(el as i128)
    }
}

impl CanonicalIsomorphismInfo<StaticRing<i64>> for BigIntRing {

    fn has_isomorphism(&self, _from: &StaticRing<i64>) -> RingPropValue {
        RingPropValue::True
    }

    fn preimage(&self, _from: &StaticRing<i64>, el: BigInt) -> i64 {
        el.to_int().expect("Overflow when embedding BigInt into i128") as i64
    }
}

impl CanonicalEmbeddingInfo<BigIntSOORing> for BigIntRing {

    fn has_embedding(&self, _from: &BigIntSOORing) -> RingPropValue {
        RingPropValue::True
    }

    fn embed(&self, _from: &BigIntSOORing, el: BigIntSOO) -> BigInt {
        el.to_bigint()
    }
}

impl CanonicalIsomorphismInfo<BigIntSOORing> for BigIntRing {

    fn has_isomorphism(&self, _from: &BigIntSOORing) -> RingPropValue {
        RingPropValue::True
    }

    fn preimage(&self, _from: &BigIntSOORing, el: BigInt) -> BigIntSOO {
        BigIntSOO::from(el)
    }
}

#[cfg(test)]
use std::str::FromStr;
#[cfg(test)]
use vector_map::VecMap;

#[test]
fn test_print_power_2() {
    let x = BigInt {
        negative: false,
        data: Vector::new(vec![0, 0, 1])
    };
    assert_eq!("340282366920938463463374607431768211456", format!("{}", x));
}

#[test]
fn test_from() {
    assert_eq!(BigInt { negative: false, data: Vector::new(vec![]) }, BigInt::from(0));
    assert_eq!(BigInt { negative: false, data: Vector::new(vec![2138479]) }, BigInt::from(2138479));
    assert_eq!(BigInt { negative: true, data: Vector::new(vec![2138479]) }, BigInt::from(-2138479));
    assert_eq!(BigInt { negative: false, data: Vector::new(vec![0x38691a350bf12fca, 0x1]) }, BigInt::from(0x138691a350bf12fca));
}

#[test]
fn test_to_int() {
    assert_eq!(0, BigInt { negative: false, data: Vector::new(vec![]) }.to_int().unwrap());
    assert_eq!(2138479, BigInt { negative: false, data: Vector::new(vec![2138479]) }.to_int().unwrap());
    assert_eq!(-2138479, BigInt { negative: true, data: Vector::new(vec![2138479]) }.to_int().unwrap());
    assert_eq!(0x138691a350bf12fca, BigInt { negative: false, data: Vector::new(vec![0x38691a350bf12fca, 0x1]) }.to_int().unwrap());
    assert_eq!(Err(()), BigInt { negative: false, data: Vector::new(vec![0x38691a350bf12fca, 0x38691a350bf12fca, 0x1]) }.to_int());
    assert_eq!(i128::MAX, BigInt { negative: false, data: Vector::new(vec![(i128::MAX & ((1 << 64) - 1)) as u64, (i128::MAX >> 64) as u64]) }.to_int().unwrap());
    assert_eq!(i128::MIN + 1, BigInt { negative: true, data: Vector::new(vec![(i128::MAX & ((1 << 64) - 1)) as u64, (i128::MAX >> 64) as u64]) }.to_int().unwrap());
    // this is the possibly surprising, exceptional case
    assert_eq!(Err(()), BigInt { negative: true, data: Vector::new(vec![0, (i128::MAX >> 64) as u64 + 1]) }.to_int());
    assert_eq!(i64::MAX as i128 + 1, BigInt { negative: false, data: Vector::new(vec![i64::MAX as u64 + 1]) }.to_int().unwrap());
    assert_eq!(u64::MAX as i128, BigInt { negative: false, data: Vector::new(vec![u64::MAX]) }.to_int().unwrap());
}

#[test]
fn test_from_str_radix() {
    let x = BigInt::from_str_radix("fa3032c0ae8202135", 16).unwrap();
    assert_eq!("288447441784111374645", format!("{}", x));

    let y = BigInt::from_str_radix("1738495302390560118908327", 10).unwrap();
    assert_eq!("1738495302390560118908327", format!("{}", y));

    let u = BigInt::from_str_radix("87112285931760246650091887388390057836920", 10).unwrap();
    let v = BigInt::from_str_radix("10000000000000000BC0000000000000178", 16).unwrap();
    assert_eq!(u, v);
}

#[test]
fn test_sub_assign() {
    let mut x = "4294836225".parse::<BigInt>().unwrap();
    let y =     "4294967297".parse::<BigInt>().unwrap();
    let z =        "-131072".parse::<BigInt>().unwrap();
    x = BigInt::RING.sub_ref_fst(&x, y);
    assert_eq!(z, x);
}

#[test]
fn test_shift_right() {
    let mut x = BigInt::from_str_radix("9843a756781b34567f81394", 16).unwrap();
    let z = BigInt::from_str_radix("9843a756781b34567", 16).unwrap();
    x = BigInt::RING.euclidean_div_pow_2(x, 24);
    assert_eq!(z, x);

    let mut x = BigInt::from_str_radix("-9843a756781b34567f81394", 16).unwrap();
    let z = BigInt::from_str_radix("-9843a756781b34567", 16).unwrap();
    x = BigInt::RING.euclidean_div_pow_2(x, 24);
    assert_eq!(z, x);
}

#[test]
fn test_assumptions_integer_division() {
    assert_eq!(-1, -3 / 2);
    assert_eq!(-1, 3 / -2);
    assert_eq!(1, -3 / -2);
    assert_eq!(1, 3 / 2);

    assert_eq!(-1, -3 % 2);
    assert_eq!(1, 3 % -2);
    assert_eq!(-1, -3 % -2);
    assert_eq!(1, 3 % 2);
}

#[test]
fn test_axioms() {
    const NUMBERS: [&'static str; 10] = [
        "5444517870735015415413993718908291383295", // power of two - 1
        "5444517870735015415413993718908291383296", // power of two
        "-5444517870735015415413993718908291383295",
        "-5444517870735015415413993718908291383296",
        "3489", // the rest is random
        "891023591340178345678931246518793456983745682137459364598623489512389745698237456890239238476873429872346579",
        "172365798123602365091834765607185713205612370956192783561461248973265193754762751378496572896497125361819754136",
        "0",
        "-231780567812394562346324763251741827457123654871236548715623487612384752328164",
        "+1278367182354612381234568509783420989356938472561078564732895634928563482349872698723465"
    ];
    let ns = NUMBERS.iter().cloned().map(BigInt::from_str).map(Result::unwrap).collect::<Vec<_>>();
    let l = ns.len();
    for i in 0..l {
        assert_eq!(BigInt::ZERO, BigInt::RING.sub(ns[i].clone(), ns[i].clone()));
        if !BigInt::RING.is_zero(&ns[i]) {
            assert_eq!(BigInt::RING.one(), BigInt::RING.euclidean_div(ns[i].clone(), &ns[i]));
        }
        assert_eq!(ns[i], BigInt::RING.add_ref(BigInt::ZERO, &ns[i]));
        assert_eq!(ns[i], BigInt::RING.mul_ref(&ns[i], &BigInt::RING.one()));
    }
    for i in 0..l {
        for j in 0..l {
            if !BigInt::RING.is_zero(&ns[j]) {
                assert_eq!(ns[i], BigInt::RING.add(
                    BigInt::RING.mul_ref(
                        &BigInt::RING.euclidean_div(ns[i].clone(), &ns[j]), 
                        &ns[j]
                    ), 
                        BigInt::RING.euclidean_rem(ns[i].clone(), &ns[j])
                    )
                );
            }
            assert_eq!(BigInt::RING.add(ns[i].clone(), ns[j].clone()), BigInt::RING.add(ns[j].clone(), ns[i].clone()));
            assert_eq!(BigInt::RING.mul_ref(&ns[i], &ns[j]), BigInt::RING.mul_ref(&ns[j], &ns[i]));
        }
    }
    for i in 0..l {
        for j in 0..l {
            for k in 0..l {
                assert_eq!(
                    BigInt::RING.mul_ref(&ns[k], &BigInt::RING.add_ref(ns[i].clone(), &ns[j])), 
                    BigInt::RING.add(BigInt::RING.mul_ref(&ns[k], &ns[j]), BigInt::RING.mul_ref(&ns[k], &ns[i]))
                );assert_eq!(
                    BigInt::RING.mul_ref(&BigInt::RING.mul_ref(&ns[i], &ns[j]), &ns[k]), 
                    BigInt::RING.mul_ref(&BigInt::RING.mul_ref(&ns[k], &ns[i]), &ns[j])
                );
                assert_eq!(
                    BigInt::RING.add_ref(BigInt::RING.add_ref(ns[i].clone(), &ns[j]), &ns[k]), 
                    BigInt::RING.add_ref(BigInt::RING.add_ref(ns[k].clone(), &ns[i]), &ns[j])
                );
            }
        }
    }
}

#[bench]
fn bench_mul(bencher: &mut test::Bencher) {
    let x = BigInt::from_str_radix("2382385687561872365981723456981723456987134659834659813491964132897159283746918732563498628754", 10).unwrap();
    let y = BigInt::from_str_radix("48937502893645789234569182735646324895723409587234", 10).unwrap();
    let z = BigInt::from_str_radix("116588006478839442056346504147013274749794691549803163727888681858469844569693215953808606899770104590589390919543097259495176008551856143726436", 10).unwrap();
    bencher.iter(|| {
        let p = BigInt::RING.mul_ref(&x, &y);
        assert_eq!(z, p);
    })
}

#[test]
fn from_to_float_approx() {
    let x: f64 = 83465209236517892563478156042389675783219532497861237985328563.;
    let y = BigInt::RING.to_float_approx(&BigInt::RING.from_float_approx(x).unwrap());
    assert!(x * 0.99 < y);
    assert!(y < x * 1.01);
}

#[bench]
fn bench_div(bencher: &mut test::Bencher) {
    let x = BigInt::from_str_radix("2382385687561872365981723456981723456987134659834659813491964132897159283746918732563498628754", 10).unwrap();
    let y = BigInt::from_str_radix("48937502893645789234569182735646324895723409587234", 10).unwrap();
    let z = BigInt::from_str_radix("48682207850683149082203680872586784064678018", 10).unwrap();
    bencher.iter(|| {
        let q = BigInt::RING.euclidean_div(x.clone(), &y);
        assert_eq!(z, q);
    })
}

#[test]
fn test_eq() {
    fn calculate_hash<T: Hash>(t: &T) -> u64 {
        let mut s = std::collections::hash_map::DefaultHasher::new();
        t.hash(&mut s);
        s.finish()
    }

    let a = BigInt { negative: false, data: Vector::new(vec![98711]) };
    let b = BigInt { negative: false, data: Vector::new(vec![98711, 0]) };
    assert!(a == 98711);
    assert!(a == b);
    assert!(b == a);
    assert!(a != BigInt::RING.neg(a.clone()));
    assert!(calculate_hash(&a) == calculate_hash(&b));
    assert!(a != BigInt::RING.one());
    // the next line could theoretically fail, but it is very improbable and we definitly should test hash inequality
    assert!(calculate_hash(&a) != calculate_hash(&BigInt::RING.one()));
}

#[test]
fn test_is_zero() {
    let zero = BigInt::from(0);
    let nonzero = BigInt::RING.mul_pow_2(BigInt::RING.one(), 83124);
    assert!(BigInt::RING.is_zero(&zero));
    assert!(BigInt::RING.is_zero(&BigInt::RING.neg(zero)));
    assert!(!BigInt::RING.is_zero(&nonzero));
    assert!(!BigInt::RING.is_zero(&BigInt::RING.neg(nonzero)));
}

#[test]
fn test_cmp_small() {
    assert!("-23678".parse::<BigInt>().unwrap() < 0);
}

#[test]
fn test_factor() {
    let ring = WrappingRing::new(&BigInt::RING);
    let mut expected = VecMap::new();
    expected.insert(ring.wrap(BigInt::from(7)), 2);
    expected.insert(ring.wrap(BigInt::from(2)), 1);
    assert_eq!(expected, BigInt::RING.factor(BigInt::from(98)));
    expected = VecMap::new();
    expected.insert(ring.wrap(BigInt::from(3)), 5);
    assert_eq!(expected, BigInt::RING.factor(BigInt::from(243)));
}

#[test]
fn test_is_prime() {
    assert_eq!(false, BigInt::RING.is_prime(&BigInt::from(81)));
}

#[test]
fn test_cmp() {
    assert_eq!(true, BigInt::RING.is_lt(&BigInt::from(-1), &BigInt::from(2)));
    assert_eq!(true, BigInt::RING.is_lt(&BigInt::from(1), &BigInt::from(2)));
    assert_eq!(false, BigInt::RING.is_lt(&BigInt::from(2), &BigInt::from(2)));
    assert_eq!(false, BigInt::RING.is_lt(&BigInt::from(3), &BigInt::from(2)));
    assert_eq!(true, BigInt::RING.is_gt(&BigInt::from(-1), &BigInt::from(-2)));
}

#[test]
fn test_mul_pow_2() {
    assert_eq!(BigInt::from(2), BigInt::RING.mul_pow_2(BigInt::from(2), 0));
    assert_eq!(BigInt::from(4829192 * 8), BigInt::RING.mul_pow_2(BigInt::from(4829192), 3));
    assert_eq!(BigInt::RING.mul(BigInt::from(4829192), BigInt::RING.mul_pow_2(BigInt::RING.one(), 64)), BigInt::RING.mul_pow_2(BigInt::from(4829192), 64));
}

#[test]
fn test_get_uniformly_random() {
    let end_exclusive = BigInt::from(3);
    let mut rng = Rand32::new(0);
    let data: Vec<BigInt> = (0..100).map(|_| BigInt::RING.get_uniformly_random(|| rng.rand_u32(), &end_exclusive)).collect();
    assert!(data.iter().any(|x| *x == 0));
    assert!(data.iter().any(|x| *x == 1));
    assert!(data.iter().any(|x| *x == 2));
}

#[test]
fn test_from_overflow() {
    assert_eq!(BigInt { data: Vector::new(vec![0, 1 << 63]), negative: true }, BigInt::from(i128::MIN));
    assert_eq!(format!("{}", i128::MIN), format!("{}", BigInt::from(i128::MIN)));
}