use std::cmp::Ordering;
use std::ops::*;
use super::alg::*;

#[derive(Debug, Clone)]
pub struct BigInt {
    /// "digits" of the (absolute value) of the number in base 2^64
    data: Vec<u64>,
    negative: bool
}

fn div_rem<T>(x: T, y: T) -> (T, T) 
    where T: std::ops::Div<Output = T> + std::ops::Rem<Output = T> + Copy
{
    let quot = x / y;
    let rem = x % y;
    (quot, rem)
}

///
/// 10 to this power fits still in a u64
/// 
const BIG_POWER_TEN_ZEROS: u32 = 19;
const BIG_POWER_TEN: u64 = 10000000000000000000u64;

impl BigInt {

    const BLOCK_BITS: usize = 64;
    pub const ZERO: BigInt = BigInt { data: Vec::new(), negative: false };

    pub fn power_of_two(n: usize) -> BigInt {
        let blocks = n / Self::BLOCK_BITS;
        let shift = n % Self::BLOCK_BITS;
        let mut result = Vec::with_capacity(blocks + 1);
        result.resize(blocks, 0);
        result.push(1 << shift);
        return BigInt {
            negative: false,
            data: result
        };
    }

    ///
    /// Calculate abs(self) += rhs * (1 << BLOCK_BITS)^block_offset
    /// 
    /// the sign bit will be left unchanged.
    /// 
    fn abs_addition(&mut self, rhs: &BigInt, block_offset: usize) {
        let mut buffer: bool = false;
        let mut i = 0;
        while i < rhs.data.len() || buffer {
            let rhs_val = *rhs.data.get(i).unwrap_or(&0);
            let j = i + block_offset;
            while j >= self.data.len() {
                self.data.push(0);
            }
            let (sum, overflow) = self.data[j].overflowing_add(rhs_val);
            if buffer {
                let (carry_sum, carry_overflow) = sum.overflowing_add(1);
                self.data[j] = carry_sum;
                buffer = overflow || carry_overflow;
            } else {
                self.data[j] = sum;
                buffer = overflow;
            }
            i += 1;
        }
    }

    pub fn is_zero(&self) -> bool {
        self.highest_set_block().is_none()
    }
    
    fn highest_set_block(&self) -> Option<usize> {
        for i in (0..self.data.len()).rev() {
            if self.data[i] != 0 {
                return Some(i);
            }
        }
        return None;
    }

    ///
    /// Computes abs(self) <=> abs(rhs)
    /// 
    fn abs_compare(&self, rhs: &BigInt) -> Ordering {
        match (self.highest_set_block(), rhs.highest_set_block()) {
           (None, None) => return Ordering::Equal,
           (Some(_), None) => return Ordering::Greater,
           (None, Some(_)) => return Ordering::Less,
           (Some(x), Some(y)) => match x.cmp(&y) {
                Ordering::Less => return Ordering::Less,
                Ordering::Greater => return Ordering::Greater,
                Ordering::Equal => {
                    for i in (0..=x).rev() {
                        match self.data[i].cmp(&rhs.data[i]) {
                            Ordering::Less => return Ordering::Less,
                            Ordering::Greater => return Ordering::Greater,
                            _ => {}
                        }
                    }
                    return Ordering::Equal;
                }
           }
        };
    }

    ///
    /// Computes abs(self) <=> abs(rhs)
    ///
    fn abs_compare_small(&self, rhs: u64) -> Ordering {
        match self.highest_set_block() {
           None => 0.cmp(&rhs),
           Some(0) => self.data[0].cmp(&rhs),
           Some(_) => Ordering::Greater,
        }
    }

    ///
    /// Calculate abs(self) -= rhs * (1 << BLOCK_BITS)^block_offset
    /// 
    /// the sign bit will be left unchanged.
    /// 
    fn abs_subtraction(&mut self, rhs: &BigInt, block_offset: usize) {
        debug_assert!(self.abs_compare(rhs) != Ordering::Less);

        let rhs_high = if let Some(d) = rhs.highest_set_block() {
            d
        } else {
            return;
        };
        
        let mut buffer: bool = false;
        let mut i = 0;
        while i <= rhs_high || buffer {
            let rhs_val = *rhs.data.get(i).unwrap_or(&0);
            let j = i + block_offset;
            debug_assert!(j < self.data.len());
            let (difference, overflow) = self.data[j].overflowing_sub(rhs_val);
            if buffer {
                let (carry_difference, carry_overflow) = difference.overflowing_sub(1);
                self.data[j] = carry_difference;
                buffer = overflow || carry_overflow;
            } else {
                self.data[j] = difference;
                buffer = overflow;
            }
            i += 1;
        }
    }

    fn abs_multiplication(&self, rhs: &BigInt) -> BigInt {
        let mut result = BigInt {
            negative: false,
            data: Vec::with_capacity(self.highest_set_block().unwrap_or(0) + rhs.highest_set_block().unwrap_or(0) + 2)
        };
        if let Some(d) = rhs.highest_set_block() {
            let mut val = BigInt::zero();
            for i in 0..=d {
                val.assign(self);
                val.abs_multiplication_small(rhs.data[i]);
                result.abs_addition(&val, i);
            }
        }
        return result;
    }

    ///
    /// Same as division_step, but for self_high == rhs_high == d
    /// 
    fn division_step_last(&mut self, rhs: &BigInt, d: usize, tmp: &mut BigInt) -> u64 {
        assert!(self.data[d] != 0);
        assert!(rhs.data[d] != 0);

        let self_high_blocks: u128 = ((self.data[d] as u128) << Self::BLOCK_BITS) | (self.data[d - 1] as u128);
        let rhs_high_blocks: u128 = ((rhs.data[d] as u128) << Self::BLOCK_BITS) | (rhs.data[d - 1] as u128);

        if rhs_high_blocks == u128::MAX {
            if self.abs_compare(&rhs) != Ordering::Less {
                self.abs_subtraction(rhs, 0);
                return 1;
            } else {
                return 0;
            }
        } else {
            let mut quotient = (self_high_blocks / (rhs_high_blocks + 1)) as u64;
            tmp.assign(rhs);
            tmp.abs_multiplication_small(quotient);
            self.abs_subtraction(&tmp, 0);
            if self.abs_compare(&rhs) != Ordering::Less {
                self.abs_subtraction(&rhs, 0);
                quotient += 1;
            }
            if self.abs_compare(&rhs) != Ordering::Less {
                self.abs_subtraction(&rhs, 0);
                quotient += 1;
            }
            // we have been at most 2 away from the real quotient, so here it must be done so far
            debug_assert!(self.abs_compare(&rhs) == Ordering::Less);
            return quotient;
        }
    }

    ///
    /// Finds some integer d such that subtracting d * rhs from self clears the top
    /// block of self. self will be assigned the value after the subtraction and d
    /// will be returned as d = (u * 2 ^ block_bits + l) * 2 ^ (k * block_bits) 
    /// where the return value is (u, l, k)
    /// 
    fn division_step(&mut self, rhs: &BigInt, self_high: usize, rhs_high: usize, tmp: &mut BigInt) -> (u64, u64, usize) {
        assert!(self_high > rhs_high);
        assert!(self.data[self_high] != 0);
        assert!(rhs.data[rhs_high] != 0);

        // the best approximation of the quotient we get through self_high_blocks / (rhs_high_blocks + 1)
        // the +1 is required to ensure that the quotient is smaller than the real quotient
        // one can prove that subtracting this is at most 2 * rhs away from the actual remainder
        
        // Therefore, the uppermost block may not be completely cleared. Therefore, perform the division again
        // with the one block shifted rhs. Here, we only use the upper 64 bits of rhs as otherwise, the truncating
        // division will only yield 0 and we will get smaller than rhs * shift, but may still have upper bits
        // uncleared (as rhs may have upper bits uncleared)

        let mut result_upper = 0;
        let mut result_lower = 0;

        {
            let self_high_blocks: u128 = ((self.data[self_high] as u128) << Self::BLOCK_BITS) | (self.data[self_high - 1] as u128);
            let rhs_high_blocks: u128 = ((rhs.data[rhs_high] as u128) << Self::BLOCK_BITS) | (rhs.data[rhs_high - 1] as u128);

            if rhs_high_blocks != u128::MAX && self_high_blocks >= (rhs_high_blocks + 1) {
                let mut quotient = (self_high_blocks / (rhs_high_blocks + 1)) as u64;
                debug_assert!(quotient != 0);
                tmp.assign(rhs);
                tmp.abs_multiplication_small(quotient);
                self.abs_subtraction(&tmp, self_high - rhs_high);

                // we might be up to 2 away from the real quotient
                if self.data[self_high] > rhs.data[rhs_high] {
                    self.abs_subtraction(rhs, self_high - rhs_high - 1);
                    quotient += 1;
                }
                if self.data[self_high] > rhs.data[rhs_high] {
                    self.abs_subtraction(rhs, self_high - rhs_high - 1);
                    quotient += 1;
                }
                result_upper = quotient;
            }
        }

        {
            let self_high_blocks: u128 = ((self.data[self_high] as u128) << Self::BLOCK_BITS) | (self.data[self_high - 1] as u128);

            if self.data[self_high] != 0 {
                let  quotient = (self_high_blocks / (rhs.data[rhs_high] as u128 + 1)) as u64;
                tmp.assign(rhs);
                tmp.abs_multiplication_small(quotient);
                self.abs_subtraction(&tmp, self_high - rhs_high - 1);
                
                result_lower = quotient;
            }
            debug_assert!(self.data[self_high] == 0);
            return (result_upper, result_lower, self_high - rhs_high - 1);
        }
    }

    fn abs_division(&mut self, rhs: &BigInt) -> BigInt {
        assert!(!rhs.is_zero());

        if let Some(mut d) = self.highest_set_block() {
            let mut tmp = BigInt::zero();
            let k = rhs.highest_set_block().expect("Division by zero");
            if d < k {
                return Self::ZERO.clone();
            } else if k == 0 {
                let rem = self.abs_division_small(rhs.data[0]);
                let rem_data = vec![rem];
                let div_data = std::mem::replace(&mut self.data, rem_data);
                return BigInt {
                    negative: false,
                    data: div_data
                };
            } else {
                let mut result_data = Vec::new();
                result_data.resize(d + 1 - k, 0);
                while d > k {
                    if self.data[d] != 0 {
                        let (quo_upper, quo_lower, quo_power) = self.division_step(&rhs, d, k, &mut tmp);
                        result_data[quo_power] += quo_lower;
                        result_data[quo_power + 1] += quo_upper;
                        debug_assert!(self.data[d] == 0);
                    }
                    d -= 1;
                }
                let quo = self.division_step_last(&rhs, d, &mut tmp);
                result_data[0] += quo;
                return BigInt {
                    negative: false,
                    data: result_data
                };
            }
        } else {
            return Self::ZERO.clone();
        }
    }

    ///
    /// Calculates self /= divisor and returns the remainder of the division.
    /// This only works for positive numbers, as for negative numbers, as the remainder
    /// must be returned as a u64 to avoid overflow. Instead of throwing, this function
    /// therefore works with abs(self) instead of self.
    /// 
    /// the sign bit will be left unchanged.
    /// 
    fn abs_division_small(&mut self, divisor: u64) -> u64 {
        assert!(divisor != 0);
        let highest_block_opt = self.highest_set_block();
        if highest_block_opt == Some(0) {
            let (quo, rem) = div_rem(self.data[0], divisor);
            self.data[0] = quo;
            return rem;
        } else if let Some(highest_block) = highest_block_opt {
            let (quo, rem) = div_rem(self.data[highest_block], divisor);
            let mut buffer = rem as u128;
            self.data[highest_block] = quo;
            for i in (0..highest_block).rev() {
                buffer = (buffer << Self::BLOCK_BITS) | (self.data[i] as u128);
                let (quo, rem) = div_rem(buffer, divisor as u128);
                debug_assert!(quo <= u64::MAX as u128);
                self.data[i] = quo as u64;
                buffer = rem;
            }
            return buffer as u64;
        } else {
            return 0;
        }
    }

    pub fn assign(&mut self, rhs: &BigInt) {
        self.negative = rhs.negative;
        self.data.clear();
        if let Some(d) = rhs.highest_set_block() {
            if self.data.len() < d {
                self.data.reserve(d - self.data.len());
            }
            for i in 0..=d {
                self.data.push(rhs.data[i]);
            }
        }
    }

    fn abs_multiplication_small(&mut self, factor: u64) {
        if let Some(d) = self.highest_set_block() {
            let mut buffer: u64 = 0;
            for i in 0..=d {
                let prod = self.data[i] as u128 * factor as u128 + buffer as u128;
                self.data[i] = (prod & ((1u128 << Self::BLOCK_BITS) - 1)) as u64;
                buffer = (prod >> Self::BLOCK_BITS) as u64;
            }
            if d + 1 < self.data.len() {
                self.data[d + 1] = buffer;
            } else {
                self.data.push(buffer);
            }
        }
    }

    ///
    /// Calculates abs(self) += summand * (1 << BLOCK_BITS)^block_offset
    /// 
    /// the sign bit will be left unchanged.
    /// 
    fn abs_addition_small(&mut self, summand: u64) {
        if self.data.len() > 0 {
            let (sum, mut buffer) = self.data[0].overflowing_add(summand);
            self.data[0] = sum;
            let mut i = 1;
            while buffer {
                if self.data.len() <= i {
                    self.data.push(0);
                }
                let (sum, overflow) = self.data[i].overflowing_add(1);
                buffer = overflow;
                self.data[i] = sum;
                i += 1;
            }
        } else {
            self.data.push(summand);
        }
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
            data: Vec::with_capacity(data.size_hint().0),
            negative: false
        };
        for value in data {
            let val = value?;
            debug_assert!(val < base);
            result.abs_multiplication_small(base);
            result.abs_addition_small(val);
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
            .map(|chunk| chunk.and_then(|n| u64::from_str_radix(n, base).map_err(BigIntParseError::from)));
        let mut result = Self::from_radix(it, (base as u64).pow(chunk_size as u32));
        if let Ok(r) = &mut result {
            r.negative = negative;
        }
        return result;
    }

    pub fn div_rem(&mut self, rhs: &BigInt) -> BigInt {
        let mut quotient = self.abs_division(rhs);
        quotient.negative = self.negative ^ rhs.negative;
        return quotient;
    }

    pub fn div_rem_small(mut self, rhs: i64) -> (BigInt, i64) {
        let mut remainder = self.abs_division_small(rhs.abs() as u64) as i64;
        if self.negative {
            remainder = -remainder;
        }
        self.negative = self.negative ^ (rhs < 0);
        return (self, remainder);
    }

    pub fn normalize(&mut self) {
        if let Some(d) = self.highest_set_block() {
            self.data.truncate(d + 1);
        } else {
            self.data.truncate(0);
        }
    }

    pub fn log2_floor(&self) -> usize {
        if let Some(d) = self.highest_set_block() {
            return Self::BLOCK_BITS - self.data[d].leading_zeros() as usize - 1 + d * Self::BLOCK_BITS;
        } else {
            // the number is zero, so the result would be -inf
            panic!("log2 is undefined for 0");
        }
    }

    ///
    /// Returns the float that is closest to the integer. Note that if for very big numbers
    /// (with abs() in the order of magnitude 2^1024 or greater, this can even yield infinity)
    /// 
    pub fn to_float_approx(&self) -> f64 {
        if let Some(d) = self.highest_set_block() {
            let mut upper_part = self.data[d] as f64 * 2f64.powi(Self::BLOCK_BITS as i32);
            if d > 0 {
                upper_part += self.data[d - 1] as f64;
            }
            return upper_part * 2f64.powi(Self::BLOCK_BITS as i32).powi(d as i32 - 1)
        } else {
            return 0.;
        }
    }

    ///
    /// Returns a BigInt that has the given value, rounded. Note that
    /// for very big numbers, the float representation can be very imprecise.
    /// 
    pub fn from_float_approx(val: f64) -> BigInt {
        if val.abs() <= 0.5 {
            return BigInt::zero()
        } else if val.abs() <= 1.5 {
            return if val.is_sign_negative() { -BigInt::one() } else { BigInt::one() };
        } else {
            const MANTISSA: i32 = 52;
            let exp = std::cmp::max(val.abs().log2() as i32, MANTISSA) - MANTISSA;
            let int = (val.abs() / 2f64.powi(exp)).trunc() as u64;
            let blocks = exp as usize / Self::BLOCK_BITS;
            let within_block_shift = exp as usize % Self::BLOCK_BITS;
            let mut result = (0..blocks).map(|_| 0).collect::<Vec<_>>();
            result.push(int << within_block_shift);
            if within_block_shift != 0 {
                result.push(int >> (Self::BLOCK_BITS - within_block_shift));
            }
            return BigInt {
                negative: val.is_sign_negative(),
                data: result
            };
        }
    }

    pub fn pow(self, power: u64) -> BigInt {
        StaticRing::<RingAxiomsEuclideanRing, BigInt>::RING.pow(self, power)
    }

    pub fn abs(mut self) -> BigInt {
        self.negative = false;
        return self;
    }
}

impl From<i64> for BigInt {
    fn from(val: i64) -> BigInt {
        if val < 0 {
            BigInt {
                negative: true,
                data: vec![(-val) as u64]
            }
        } else {
            BigInt {
                negative: false,
                data: vec![val as u64]
            }
        }
    }
}

impl From<u64> for BigInt {
    fn from(val: u64) -> BigInt {
        BigInt {
            negative: false,
            data: vec![val]
        }
    }
}

impl PartialEq for BigInt {

    fn eq(&self, rhs: &BigInt) -> bool {
        let highest_block = self.highest_set_block();
        if highest_block != rhs.highest_set_block() {
            return false;
        }
        if let Some(d) = highest_block {
            for i in 0..=d {
                if self.data[i] != rhs.data[i] {
                    return false;
                }
            }
        }
        return true;
    }
}

impl PartialEq<u64> for BigInt {

    fn eq(&self, rhs: &u64) -> bool {
        (self.highest_set_block() == Some(0) && self.data[0] == *rhs && !self.negative) ||
        (self.is_zero() && *rhs == 0)
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

impl Eq for BigInt {}

impl PartialOrd<u64> for BigInt {

    fn partial_cmp(&self, rhs: &u64) -> Option<Ordering> {
        if self.negative && !self.is_zero() {
            Some(Ordering::Less)
        } else {
            Some(self.abs_compare_small(*rhs))
        }
    }
}

impl PartialOrd for BigInt {

    fn partial_cmp(&self, rhs: &BigInt) -> Option<Ordering> {
        Some(self.cmp(rhs))
    }
}

impl Ord for BigInt {

    fn cmp(&self, rhs: &BigInt) -> Ordering {
        if self.is_zero() && rhs.is_zero() {
            Ordering::Equal
        } else {
            match (self.negative, rhs.negative) {
                (true, true) => rhs.abs_compare(self),
                (true, false) => Ordering::Greater,
                (false, true) => Ordering::Less,
                (false, false) => self.abs_compare(rhs)
            }
        }
    }
}

impl<T> Add<T> for BigInt 
    where BigInt: AddAssign<T>
{
    type Output = BigInt;

    fn add(mut self, rhs: T) -> Self::Output {
        self += rhs;
        return self;
    }
}

impl AddAssign<BigInt> for BigInt {

    fn add_assign(&mut self, mut rhs: BigInt) {
        if self.negative == rhs.negative {
            self.abs_addition(&rhs, 0);
        } else if self.abs_compare(&rhs) != Ordering::Less {
            self.abs_subtraction(&rhs, 0);
        } else {
            rhs.abs_subtraction(&self, 0);
            std::mem::swap(self, &mut rhs);
        }
    }
}

impl AddAssign<&BigInt> for BigInt {

    fn add_assign(&mut self, rhs: &BigInt) {
        if self.negative == rhs.negative {
            self.abs_addition(rhs, 0);
        } else if self.abs_compare(&rhs) != Ordering::Less {
            self.abs_subtraction(rhs, 0);
        } else {
            let mut result = rhs.clone();
            result.abs_subtraction(&self, 0);
            std::mem::swap(self, &mut result);
        }
    }
}

impl AddAssign<i64> for BigInt {

    fn add_assign(&mut self, rhs: i64) {
        if self.negative == (rhs < 0) {
            self.abs_addition_small(rhs.abs() as u64);
        } else {
            let rhs_bigint = BigInt {
                negative: rhs < 0,
                data: vec![rhs.abs() as u64]
            };
            *self += rhs_bigint;
        }
    }
}

impl<T> Sub<T> for BigInt 
    where BigInt: SubAssign<T>
{
    type Output = BigInt;

    fn sub(mut self, rhs: T) -> Self::Output {
        self -= rhs;
        return self;
    }
}

impl SubAssign for BigInt {

    fn sub_assign(&mut self, rhs: BigInt) {
        self.negative = !self.negative;
        self.add_assign(rhs);
        self.negative = !self.negative;
        self.normalize();
    }
}

impl SubAssign<&BigInt> for BigInt {

    fn sub_assign(&mut self, rhs: &BigInt) {
        self.negative = !self.negative;
        *self += rhs;
        self.negative = !self.negative;
        self.normalize();
    }
}

impl SubAssign<i64> for BigInt {

    fn sub_assign(&mut self, rhs: i64) {
        self.negative = !self.negative;
        *self += rhs;
        self.negative = !self.negative;
    }
}

impl<T> Mul<T> for BigInt 
    where BigInt: MulAssign<T>
{
    type Output = BigInt;

    fn mul(mut self, rhs: T) -> Self::Output {
        self *= rhs;
        return self;
    }
}

impl MulAssign for BigInt {

    fn mul_assign(&mut self, rhs: BigInt) {
        let sign = self.negative ^ rhs.negative;
        *self = self.abs_multiplication(&rhs);
        self.negative = sign;
    }
}

impl MulAssign<&BigInt> for BigInt {

    fn mul_assign(&mut self, rhs: &BigInt) {
        let sign = self.negative ^ rhs.negative;
        *self = self.abs_multiplication(rhs);
        self.negative = sign;
    }
}

impl MulAssign<i64> for BigInt {

    fn mul_assign(&mut self, rhs: i64) {
        self.abs_multiplication_small(rhs.abs() as u64);
        self.negative ^= rhs < 0;
    }
}

impl<T> Div<T> for BigInt
    where BigInt: DivAssign<T>
{
    type Output = BigInt;

    fn div(mut self, rhs: T) -> BigInt {
        self /= rhs;
        return self;
    }
}

impl DivAssign for BigInt {

    fn div_assign(&mut self, rhs: BigInt) {
        *self = self.div_rem(&rhs);
    }
}

impl DivAssign<&BigInt> for BigInt {

    fn div_assign(&mut self, rhs: &BigInt) {
        *self = self.div_rem(rhs);
    }
}

impl DivAssign<u64> for BigInt {

    fn div_assign(&mut self, rhs: u64) {
        self.abs_division_small(rhs);
    }
}

impl<T> Rem<T> for BigInt 
    where BigInt: RemAssign<T>
{
    type Output = BigInt;

    fn rem(mut self, rhs: T) -> BigInt {
        self %= rhs;
        return self;
    }
}

impl RemAssign for BigInt {

    fn rem_assign(&mut self, rhs: BigInt) {
        self.div_rem(&rhs);
        self.normalize();
    }
}

impl RemAssign<&BigInt> for BigInt {

    fn rem_assign(&mut self, rhs: &BigInt) {
        self.div_rem(&rhs);
        self.normalize();
    }
}

impl Rem<i64> for BigInt {

    type Output = i64;

    fn rem(self, rhs: i64) -> Self::Output {
        self.div_rem_small(rhs).1
    }
}

impl Shr<usize> for BigInt {

    type Output = BigInt;

    fn shr(mut self, rhs: usize) -> Self::Output {
        let drop_blocks = rhs / Self::BLOCK_BITS;
        let shift_amount = rhs % Self::BLOCK_BITS;
        let keep_bits = Self::BLOCK_BITS - shift_amount;
        if self.data.len() <= drop_blocks {
            self = Self::ZERO.clone();
        } else {
            self.data.drain(0..drop_blocks);
            self.data[0] >>= shift_amount;
            for i in 1..self.data.len() {
                let rotated = self.data[i].rotate_right(shift_amount as u32);
                self.data[i] = rotated & (u64::MAX >> shift_amount);
                self.data[i - 1] |= rotated & (u64::MAX << keep_bits);
            }
        }
        return self;
    }
}

impl One for BigInt {

    fn one() -> BigInt {
        BigInt {
            negative: false,
            data: vec![1]
        }
    }
}

impl Zero for BigInt {

    fn zero() -> BigInt {
        BigInt::ZERO.clone()
    }
}

impl Neg for BigInt {

    type Output = BigInt;

    fn neg(mut self) -> BigInt {
        self.negative = !self.negative;
        self
    }
}

impl std::fmt::Display for BigInt {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        if self.negative {
            write!(f, "-")?;
        }
        let mut copy = self.clone();
        let mut remainders: Vec<u64> = Vec::with_capacity((self.highest_set_block().unwrap_or(0) + 1) * Self::BLOCK_BITS / 3);
        while !copy.is_zero() {
            let rem = copy.abs_division_small(BIG_POWER_TEN);
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

impl RingEl for BigInt {

    type Axioms = RingAxiomsEuclideanRing;
}

impl EuclideanRingEl for BigInt {

    fn div_rem(&mut self, rhs: Self) -> Self
    {
        self.div_rem(&rhs)
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

#[cfg(test)]
use std::str::FromStr;

#[test]
fn test_print_power_2() {
    let x = BigInt {
        negative: false,
        data: vec![0, 0, 1]
    };
    assert_eq!("340282366920938463463374607431768211456", format!("{}", x));
}

#[test]
fn test_from_str_radix() {
    let x = BigInt::from_str_radix("Fa3032c0ae8202135", 16).unwrap();
    assert_eq!("288447441784111374645", format!("{}", x));

    let y = BigInt::from_str_radix("1738495302390560118908327", 10).unwrap();
    assert_eq!("1738495302390560118908327", format!("{}", y));

    let u = BigInt::from_str_radix("87112285931760246650091887388390057836920", 10).unwrap();
    let v = BigInt::from_str_radix("10000000000000000BC0000000000000178", 16).unwrap();
    assert_eq!(u, v);
}

#[test]
fn test_abs_subtraction() {
    let mut x = "923645871236598172365987287530543".parse::<BigInt>().unwrap();
    let y =      "58430657823473456743684735863478".parse::<BigInt>().unwrap();
    let z =     "865215213413124715622302551667065".parse::<BigInt>().unwrap();
    x.abs_subtraction(&y, 0);
    assert_eq!(z, x);
}

#[test]
fn test_abs_subtraction_with_carry() {
    let mut x = BigInt::from_str_radix("1000000000000000000", 16).unwrap();
    let y =      BigInt::from_str_radix("FFFFFFFFFFFFFFFF00", 16).unwrap();
    x.abs_subtraction(&y, 0);
    assert_eq!(x, 256u64);
}

#[test]
fn test_sub_assign() {
    let mut x = "4294836225".parse::<BigInt>().unwrap();
    let y =     "4294967297".parse::<BigInt>().unwrap();
    let z =        "-131072".parse::<BigInt>().unwrap();
    x -= &y;
    assert_eq!(z, x);
}

#[test]
fn test_abs_addition() {
    let mut x = "923645871236598172365987287530543".parse::<BigInt>().unwrap();
    let y =      "58430657823473456743684735863478".parse::<BigInt>().unwrap();
    let z =     "982076529060071629109672023394021".parse::<BigInt>().unwrap();
    x.abs_addition(&y, 0);
    assert_eq!(z, x);
}

#[test]
fn test_abs_addition_with_carry() {
    let mut x =             BigInt::from_str_radix("1BC00000000000000BC", 16).unwrap();
    let y =  BigInt::from_str_radix("FFFFFFFFFFFFFFFF0000000000000000BC", 16).unwrap();
    let z = BigInt::from_str_radix("10000000000000000BC0000000000000178", 16).unwrap();
    x.abs_addition(&y, 0);
    assert_eq!(z, x);
}

#[test]
fn test_abs_multiplication() {
    let x = BigInt::from_str_radix("57873674586797895671345345", 10).unwrap();
    let y = BigInt::from_str_radix("21308561789045691782534873921650342768903561413264128756389247568729346542359871235465", 10).unwrap();
    let z = BigInt::from_str_radix("1233204770891906354921751949503652431220138020953161094405729272872607166072371117664593787957056214903826660425", 10).unwrap();
    assert_eq!(z, x.abs_multiplication(&y));
}

#[test]
fn test_abs_division_no_remainder() {
    let mut x = BigInt::from_str_radix("578435387FF0582367863200000000000000000000", 16).unwrap();
    let y =                          BigInt::from_str_radix("200000000000000000000", 16).unwrap();
    let z = BigInt::from_str_radix("2BC21A9C3FF82C11B3C319", 16).unwrap();
    let quotient = x.abs_division(&y);
    assert_eq!(BigInt::ZERO, x);
    assert_eq!(z, quotient);
}

#[test]
fn test_abs_division_remainder() {
    let mut x = BigInt::from_str_radix("578435387FF0582367863200000000007651437856", 16).unwrap();
    let y =                          BigInt::from_str_radix("200000000000000000000", 16).unwrap();
    let z = BigInt::from_str_radix("2BC21A9C3FF82C11B3C319", 16).unwrap();
    let r = BigInt::from_str_radix("7651437856", 16).unwrap();
    let quotient = x.abs_division(&y);
    assert_eq!(r, x);
    assert_eq!(z, quotient);
}

#[test]
fn test_abs_division_big() {
    let mut x = BigInt::from_str_radix("581239456149785691238569872349872348569871269871234657986123987237865847935698734296434575367565723846982523852347", 10).unwrap();
    let y = BigInt::from_str_radix("903852718907268716125180964783634518356783568793426834569872365791233387356325", 10).unwrap();
    let q = BigInt::from_str_radix("643068769934649368349591185247155725", 10).unwrap();
    let r = BigInt::from_str_radix("265234469040774335115597728873888165088018116561138613092906563355599185141722", 10).unwrap();
    let quotient = x.abs_division(&y);
    assert_eq!(r, x);
    assert_eq!(q, quotient);
}

#[test]
fn test_abs_division_small() {
    let mut x = BigInt::from_str_radix("891023591340178345678931246518793456983745682137459364598623489512389745698237456890239238476873429872346579", 10).unwrap();
    let q = BigInt::from_str_radix("255380794307875708133829534685810678413226048190730686328066348384175908769916152734376393945793473738133", 10).unwrap();
    x.abs_division_small(3489);
    assert_eq!(q, x);
}

#[test]
fn test_shift_right() {
    let mut x = BigInt::from_str_radix("9843a756781b34567f81394", 16).unwrap();
    let z = BigInt::from_str_radix("9843a756781b34567", 16).unwrap();
    x = x >> 24;
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
        assert_eq!(BigInt::ZERO, ns[i].clone() - ns[i].clone());
        if !ns[i].is_zero() {
            assert_eq!(BigInt::one(), ns[i].clone() / ns[i].clone());
        }
        assert_eq!(ns[i], ns[i].clone() + BigInt::zero());
        assert_eq!(ns[i], ns[i].clone() * BigInt::one());
    }
    for i in 0..l {
        for j in 0..l {
            if !ns[j].is_zero() {
                assert_eq!(ns[i], (ns[i].clone() / ns[j].clone()) * ns[j].clone() + (ns[i].clone() % ns[j].clone()));
            }
            assert_eq!(ns[i].clone() + ns[j].clone(), ns[j].clone() + ns[i].clone());
            assert_eq!(ns[i].clone() * ns[j].clone(), ns[j].clone() * ns[i].clone());
        }
    }
    for i in 0..l {
        for j in 0..l {
            for k in 0..l {
                assert_eq!(ns[k].clone() * (ns[i].clone() + ns[j].clone()), ns[k].clone() * ns[j].clone() + ns[k].clone() * ns[i].clone());
                assert_eq!(ns[k].clone() * (ns[i].clone() * ns[j].clone()), (ns[k].clone() * ns[i].clone()) * ns[j].clone());
                assert_eq!(ns[k].clone() + (ns[i].clone() + ns[j].clone()), (ns[k].clone() + ns[i].clone()) + ns[j].clone());
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
        let p = x.clone() * y.clone();
        assert_eq!(z, p);
    })
}

#[test]
fn from_to_float_approx() {
    let x: f64 = 83465209236517892563478156042389675783219532497861237985328563.;
    let y = BigInt::from_float_approx(x).to_float_approx();
    assert!(x * 0.99 < y);
    assert!(y < x * 1.01);
}

#[bench]
fn bench_div(bencher: &mut test::Bencher) {
    let x = BigInt::from_str_radix("2382385687561872365981723456981723456987134659834659813491964132897159283746918732563498628754", 10).unwrap();
    let y = BigInt::from_str_radix("48937502893645789234569182735646324895723409587234", 10).unwrap();
    let z = BigInt::from_str_radix("48682207850683149082203680872586784064678018", 10).unwrap();
    bencher.iter(|| {
        let q = x.clone() / y.clone();
        assert_eq!(z, q);
    })
}

#[test]
fn test_eq() {
    assert!(BigInt::from_str_radix("98711", 10).unwrap() != BigInt::one());
}

#[test]
fn test_cmp_small() {
    assert!("-23678".parse::<BigInt>().unwrap() < 0);
}