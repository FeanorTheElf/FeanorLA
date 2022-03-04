use super::ring::*;
use super::algebra::primes;

use std::cmp::Ordering;
use std::ops::*;
use std::hash::{Hash, Hasher};

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

    pub fn signum(&self) -> i32 {
        if self.is_zero() {
            return 0;
        } else if self.negative {
            return -1;
        } else {
            return 1;
        }
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
    pub fn abs_compare(&self, rhs: &BigInt) -> Ordering {
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
            data: Vec::with_capacity(
                self.highest_set_block().unwrap_or(0) + 
                rhs.highest_set_block().unwrap_or(0) + 2
            )
        };
        if let Some(d) = rhs.highest_set_block() {
            let mut val = BigInt::ZERO;
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

        let self_high_blocks: u128 = ((self.data[d] as u128) << Self::BLOCK_BITS) | 
            (self.data[d - 1] as u128);
        let rhs_high_blocks: u128 = ((rhs.data[d] as u128) << Self::BLOCK_BITS) | 
            (rhs.data[d - 1] as u128);

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
            // we have been at most 2 away from the real quotient, so here 
            // it must be done so far
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
    /// Complexity O(log(n))
    /// 
    fn division_step(
        &mut self, 
        rhs: &BigInt, 
        self_high: usize, 
        rhs_high: usize, 
        tmp: &mut BigInt
    ) -> (u64, u64, usize) 
    {
        assert!(self_high > rhs_high);
        assert!(self.data[self_high] != 0);
        assert!(rhs.data[rhs_high] != 0);

        // the best approximation of the quotient we get through 
        // self_high_blocks / (rhs_high_blocks + 1)
        // the +1 is required to ensure that the quotient is smaller 
        // than the real quotient
        // one can prove that subtracting this is at most 2 * rhs away 
        // from the actual remainder
        
        // Therefore, the uppermost block may not be completely cleared. 
        // Therefore, perform the division again with the one block shifted 
        // rhs. Here, we only use the upper 64 bits of rhs as otherwise, 
        // the truncating division will only yield 0 and we will get smaller 
        // than rhs * shift, but may still have upper bits
        // uncleared (as rhs may have upper bits uncleared)

        let mut result_upper = 0;
        let mut result_lower = 0;

        {
            let self_high_blocks: u128 = ((self.data[self_high] as u128) << Self::BLOCK_BITS) | 
                (self.data[self_high - 1] as u128);
            let rhs_high_blocks: u128 = ((rhs.data[rhs_high] as u128) << Self::BLOCK_BITS) | 
                (rhs.data[rhs_high - 1] as u128);

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
            let self_high_blocks: u128 = ((self.data[self_high] as u128) << Self::BLOCK_BITS) | 
                (self.data[self_high - 1] as u128);

            if self.data[self_high] != 0 {
                let quotient = (self_high_blocks / (rhs.data[rhs_high] as u128 + 1)) as u64;
                tmp.assign(rhs);
                tmp.abs_multiplication_small(quotient);
                self.abs_subtraction(&tmp, self_high - rhs_high - 1);
                
                result_lower = quotient;
            }
            debug_assert!(self.data[self_high] == 0);
            return (result_upper, result_lower, self_high - rhs_high - 1);
        }
    }

    ///
    /// Calculates abs(self) = abs(self) % abs(rhs) and returns the quotient
    /// of the division abs(self) / abs(rhs). The sign bit of self is ignored
    /// and left unchanged.
    /// 
    /// Complexity O(log(n)^2)
    /// 
    fn abs_division(&mut self, rhs: &BigInt) -> BigInt {
        assert!(!rhs.is_zero());

        if let Some(mut d) = self.highest_set_block() {
            let mut tmp = BigInt::ZERO;
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
                        result_data[quo_power] = quo_lower;
                        let (new_upper_part, overflow) = 
                            result_data[quo_power + 1].overflowing_add(quo_upper);
                        result_data[quo_power + 1] = new_upper_part;
                        if overflow {
                            result_data[quo_power + 2] += 1;
                        }
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
    /// This only works for positive numbers, as for negative numbers, 
    /// as the remainder must be returned as a u64 to avoid overflow. 
    /// Instead of throwing, this function therefore works with abs(self) 
    /// instead of self.
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

    ///
    /// Complexity O(log(n))
    /// 
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
    /// Amortized complexity O(1)
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
            .map(|chunk| chunk.and_then(|n| 
                u64::from_str_radix(n, base).map_err(BigIntParseError::from))
            );
        let mut result = Self::from_radix(it, (base as u64).pow(chunk_size as u32));
        if let Ok(r) = &mut result {
            r.negative = negative;
        }
        return result;
    }

    pub fn euclidean_div_rem(&mut self, rhs: &BigInt) -> BigInt {
        let mut quotient = self.abs_division(rhs);
        quotient.negative = self.negative ^ rhs.negative;
        return quotient;
    }

    ///
    /// Performs euclidean division, i.e. computes q, r such that
    /// `self = q * rhs + r` and `|r| < |rhs|`.
    /// Returns `(q, r)`
    /// 
    pub fn euclidean_div_rem_small(mut self, rhs: i64) -> (BigInt, i64) {
        let mut remainder = self.abs_division_small(rhs.abs() as u64) as i64;
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
    /// # use feanor_la::bigint::BigInt;
    /// assert_eq!(BigInt::from(-1), BigInt::from(-1).floor_div_small(2));
    /// ```
    /// 
    pub fn floor_div_small(self, rhs: i64) -> BigInt {
        let (quo, rem) = self.euclidean_div_rem_small(rhs);
        if rem < 0 {
            quo + if rhs < 0 { 1 } else { -1 }
        } else {
            quo
        }
    }

    pub fn normalize(&mut self) {
        if let Some(d) = self.highest_set_block() {
            self.data.truncate(d + 1);
        } else {
            self.data.truncate(0);
        }
    }

    pub fn abs_log2_floor(&self) -> usize {
        if let Some(d) = self.highest_set_block() {
            return Self::BLOCK_BITS - self.data[d].leading_zeros() as usize - 1 + 
                d * Self::BLOCK_BITS;
        } else {
            // the number is zero, so the result would be -inf
            panic!("log2 is undefined for 0");
        }
    }

    ///
    /// Returns the float that is closest to the integer. Note that 
    /// if for very big numbers (with abs() in the order of magnitude 
    /// 2^1024 or greater, this can even yield infinity)
    /// 
    pub fn to_float_approx(&self) -> f64 {
        if let Some(d) = self.highest_set_block() {
            let mut upper_part = self.data[d] as f64 * 
                2f64.powi(Self::BLOCK_BITS as i32);
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
    /// For Infinity and NaN, nothing is returned;
    /// 
    pub fn from_float_approx(val: f64) -> Option<BigInt> {
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
            let blocks = exp as usize / Self::BLOCK_BITS;
            let within_block_shift = exp as usize % Self::BLOCK_BITS;
            let mut result = (0..blocks).map(|_| 0).collect::<Vec<_>>();
            result.push(int << within_block_shift);
            if within_block_shift != 0 {
                result.push(int >> (Self::BLOCK_BITS - within_block_shift));
            }
            return Some(BigInt {
                negative: val.is_sign_negative(),
                data: result
            });
        }
    }

    pub fn pow(&self, power: u32) -> BigInt {
        Self::RING.pow(self, power)
    }

    pub fn pow_big(&self, power: &BigInt) -> BigInt {
        Self::RING.pow_big(self, power)
    }

    ///
    /// Given an increasing, continuous function f: R -> R that is negative for some x1 and 
    /// positive for some x2, finds the floor of some root of f (if f is strictly increasing, 
    /// this is unique).
    /// 
    /// # General case
    /// 
    /// This function also works in a slightly more general context. Assume that
    /// f(x) is negative for all sufficiently small x and positive for all suffiently 
    /// large x. Then this function will return the floor of some root of f. Note that
    /// this might not be a root of f, even if f has integral roots.
    /// 
    /// # Complexity
    /// 
    /// This function runs in O((T + log(d)) * log(d)) where d is the error made in 
    /// approx (i.e. the difference between the found root x and approx) and T is the
    /// time required for computing f on a value between x - d and x + d.
    /// 
    pub fn find_zero_floor<F>(mut f: F, approx: BigInt) -> BigInt
        where F: FnMut(&BigInt) -> BigInt
    {
        let mut begin = approx.clone();
        let mut step = BigInt::RING.one();
        while f(&begin).signum() > 0 {
            begin = BigInt::RING.sub_ref_snd(begin, &step);
            step *= 2;
        }
        let mut end = approx;
        step = BigInt::RING.one();
        while f(&end).signum() < 0 {
            end = BigInt::RING.add_ref(end, &step);
            step *= 2;
        }
        return Self::bisect(f, begin, end);
    }

    ///
    /// Given a continuous function f: R -> R that is negative on `begin` and 
    /// positive on `end`, finds the floor of some root of f. Note that even
    /// if f has integral roots, the returned value does not have to be a root
    /// of f.
    /// 
    /// # Complexity
    /// 
    /// This function runs in O((T + log(d)) * log(d)) where d is the difference between
    /// begin and end and T is the time required for computing f on a value between 
    /// begin and end. 
    /// 
    pub fn bisect<F>(mut f: F, mut start: BigInt, mut end: BigInt) -> BigInt
        where F: FnMut(&BigInt) -> BigInt
    {
        assert!(f(&start).signum() <= 0);
        assert!(f(&end).signum() >= 0);
        if f(&end) == 0 {
            return end;
        }
        loop {
            let mid = BigInt::RING.add_ref(start.clone(), &end).floor_div_small(2);
            if mid == start {
                return start;
            }
            match f(&mid).signum() {
                -1 => {
                    start = mid;
                },
                1 => {
                    end = mid;
                },
                _ => {
                    return mid;
                }
            }
        }
    }

    ///
    /// Computes the n-th root of this number.
    /// 
    /// # Complexity
    /// 
    /// The asymptotic worst-case complexity is O(log(n)^2), however it
    /// will be quite fast on most inputs due to internal use of floating
    /// point approximations.
    /// 
    pub fn root_floor(self, n: usize) -> BigInt {
        assert!(n > 0);
        let root_approx = self.to_float_approx().powf(1. / n as f64);
        if n % 2 == 0 {
            return BigInt::find_zero_floor(
                |x| BigInt::RING.sub_ref_snd(BigInt::RING.mul(x.clone().abs(), BigInt::RING.pow(x, (n - 1) as u32)), &self), 
                BigInt::from_float_approx(root_approx).unwrap_or(BigInt::ZERO)
            );
        } else {
            return BigInt::find_zero_floor(
                |x| BigInt::RING.sub_ref_snd(BigInt::RING.pow(x, n as u32), &self), 
                BigInt::from_float_approx(root_approx).unwrap_or(BigInt::ZERO)
            );
        }
    }

    pub fn log_floor(self, base: BigInt) -> BigInt {
        let log_approx = self.to_float_approx().log(base.to_float_approx());
        return BigInt::find_zero_floor(
            |x| base.clone().pow_big(x), 
            BigInt::from_float_approx(log_approx).unwrap_or(BigInt::ZERO)
        );
    }

    pub fn abs(mut self) -> BigInt {
        self.negative = false;
        return self;
    }

    ///
    /// Generates a uniformly random number from the range 0 to end_exclusive, using
    /// entropy from the given rng.
    /// 
    /// The distribution may not be perfectly uniform, but within 
    /// l1-statistical distance 2^(-statistical_distance_bound) of a true 
    /// uniformly random distribution
    /// 
    pub fn get_uniformly_random_oorandom(
        rng: &mut oorandom::Rand32,
        end_exclusive: &BigInt,
        statistical_distance_bound: usize
    ) -> BigInt {
        Self::get_uniformly_random(|| ((rng.rand_u32() as u64) << 32) | (rng.rand_u32() as u64), end_exclusive, statistical_distance_bound)
    }

    ///
    /// Generates a uniformly random number from the range 0 to end_exclusive, using
    /// entropy from the given rng.
    /// 
    /// The distribution may not be perfectly uniform, but within 
    /// l1-statistical distance 2^(-statistical_distance_bound) of a true 
    /// uniformly random distribution
    /// 
    pub fn get_uniformly_random<G>(
        mut rng: G, 
        end_exclusive: &BigInt, 
        statistical_distance_bound: usize
    ) -> BigInt 
        where G: FnMut() -> u64
    {
        assert!(*end_exclusive > 0);
        let k = statistical_distance_bound + end_exclusive.abs_log2_floor();
        let random_blocks = k / Self::BLOCK_BITS + 1;
        // generate a uniform random number in the range from 0 to 2^k 
        // and take it modulo end_exclusive the number of bigints 
        // between 0 and 2^k that give a fixed x differs at most by one. 
        // Therefore the probability difference to get any to distinct 
        // numbers is at most 1/2^k. The l1-distance
        // between the distributions is therefore bounded 
        // by n/2^k <= 2^(-statistical_distance_bound)
        let mut result = BigInt {
            data: (0..random_blocks).map(|_| rng()).collect(),
            negative: false
        };
        result.abs_division(&end_exclusive);
        return result;
    }

    pub fn highest_dividing_power_of_two(&self) -> usize {
        if let Some(d) = self.highest_set_block() {
            for i in 0..=d {
                if self.data[i] != 0 {
                    return self.data[i].trailing_zeros() as usize + 
                        i * Self::BLOCK_BITS;
                }
            }
            unreachable!()
        } else {
            return 0;
        }
    }

    pub fn is_bit_set(&self, i: usize) -> bool {
        if let Some(d) = self.highest_set_block() {
            let block = i / Self::BLOCK_BITS;
            let bit = i % Self::BLOCK_BITS;
            if block > d {
                false
            } else {
                ((self.data[block] >> bit) & 1) == 1
            }
        } else {
            false
        }
    }

    pub fn to_int(&self) -> Result<i64, ()> {
        if let Some(d) = self.highest_set_block() {
            if d == 0 && self.data[0] <= i64::MAX as u64 && self.negative {
                Ok(-(self.data[0] as i64))
            } else if d == 0 && self.data[0] <= i64::MAX as u64 && !self.negative {
                Ok(self.data[0] as i64)
            } else {
                Err(())
            }
        } else {
            Ok(0)
        }
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

#[derive(Debug, Clone, Copy)]
pub struct BigIntRing;

impl BigInt {
    pub const RING: BigIntRing = BigIntRing {};
}

impl Ring for BigIntRing {

    type El = BigInt;

    fn add_ref(&self, mut lhs: Self::El, rhs: &Self::El) -> Self::El {
        if lhs.negative == rhs.negative {
            lhs.abs_addition(rhs, 0);
        } else if lhs.abs_compare(&rhs) != Ordering::Less {
            lhs.abs_subtraction(rhs, 0);
        } else {
            let mut result = rhs.clone();
            result.abs_subtraction(&lhs, 0);
            std::mem::swap(&mut lhs, &mut result);
        }
        return lhs;
    }

    fn mul_ref(&self, lhs: &Self::El, rhs: &Self::El) -> Self::El {
        let mut result = lhs.abs_multiplication(rhs);
        result.negative = lhs.negative ^ rhs.negative;
        return result;
    }

    fn mul(&self, lhs: Self::El, rhs: Self::El) -> Self::El {
        let mut result = lhs.abs_multiplication(&rhs);
        result.negative = lhs.negative ^ rhs.negative;
        return result;
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
            data: vec![1]
        };
    }

    fn eq(&self, lhs: &Self::El, rhs: &Self::El) -> bool {
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

    fn is_noetherian(&self) -> bool { true }
    fn is_integral(&self) -> RingPropValue { RingPropValue::True }
    fn is_field(&self) -> RingPropValue { RingPropValue::False }

    fn div(&self, _lhs: Self::El, _rhs: &Self::El) -> Self::El {
        panic!("Not a field!")
    }

    fn from_z(&self, x: i64) -> BigInt {
        BigInt::from(x)
    }

    fn from_z_big(&self, x: &BigInt) -> BigInt {
        x.clone()
    }

    fn format(&self, el: &BigInt, f: &mut std::fmt::Formatter, _in_prod: bool) -> std::fmt::Result {
        if el.negative {
            write!(f, "-")?;
        }
        let mut copy = el.clone();
        let mut remainders: Vec<u64> = Vec::with_capacity(
            (el.highest_set_block().unwrap_or(0) + 1) * BigInt::BLOCK_BITS / 3
        );
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

impl EuclideanInfoRing for BigIntRing {

    fn is_euclidean(&self) -> RingPropValue {
        RingPropValue::True
    }

    fn euclidean_div_rem(&self, mut lhs: Self::El, rhs: &Self::El) -> (Self::El, Self::El) {
        let quo = lhs.euclidean_div_rem(rhs);
        (quo, lhs)
    }

    fn euclidean_deg(&self, el: Self::El) -> BigInt {
        el.abs()
    }
}

impl BigIntRing {

    pub fn z_embedding<'a, R>(&self, target: &'a R) -> impl 'a + FnMut(BigInt) -> R::El
        where R: Ring
    {
        move |x| target.from_z_big(&x)
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
        BigInt::RING.eq(self, rhs)
    }
}

impl PartialOrd<i64> for BigInt {

    fn partial_cmp(&self, rhs: &i64) -> Option<Ordering> {
        if self.is_zero() && *rhs == 0 {
            Some(Ordering::Equal)
        } else if self.negative && *rhs < 0 {
            Some(self.abs_compare_small((-*rhs) as u64).reverse())
        } else if self.negative && *rhs >= 0 {
            Some(Ordering::Less)
        } else if !self.negative && *rhs < 0 {
            Some(Ordering::Greater)
        } else {
            Some(self.abs_compare_small(*rhs as u64))
        }
    }
}

impl PartialOrd for BigInt {

    fn partial_cmp(&self, rhs: &BigInt) -> Option<Ordering> {
        Some(self.cmp(rhs))
    }
}

impl Neg for BigInt {

    type Output = BigInt;

    fn neg(self) -> BigInt {
        BigInt::RING.neg(self)
    }
}

impl Ord for BigInt {

    fn cmp(&self, rhs: &BigInt) -> Ordering {
        if self.is_zero() && rhs.is_zero() {
            Ordering::Equal
        } else {
            match (self.negative, rhs.negative) {
                (true, true) => rhs.abs_compare(self),
                (true, false) => Ordering::Less,
                (false, true) => Ordering::Greater,
                (false, false) => self.abs_compare(rhs)
            }
        }
    }
}

impl SingletonRing for BigIntRing {
    fn singleton() -> BigIntRing {
        BigInt::RING
    }
}

impl Add<i64> for BigInt {
    type Output = BigInt;

    fn add(mut self, rhs: i64) -> Self::Output {
        self += rhs;
        return self;
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
            take_mut::take_or_recover(self, || BigInt::ZERO, |x| BigInt::RING.add(x, rhs_bigint));
        }
    }
}

impl Add<BigInt> for BigInt {

    type Output = BigInt;

    fn add(self, rhs: BigInt) -> Self::Output {
        BigInt::RING.add(self, rhs)
    }
}

impl Sub<i64> for BigInt {
    type Output = BigInt;

    fn sub(mut self, rhs: i64) -> Self::Output {
        self -= rhs;
        return self;
    }
}

impl SubAssign<i64> for BigInt {

    fn sub_assign(&mut self, rhs: i64) {
        self.negative = !self.negative;
        *self += rhs;
        self.negative = !self.negative;
    }
}

impl Sub<BigInt> for BigInt {
    type Output = BigInt;

    fn sub(self, rhs: BigInt) -> Self::Output {
        BigInt::RING.sub(self, rhs)
    }
}

impl Mul<i64> for BigInt {
    type Output = BigInt;

    fn mul(mut self, rhs: i64) -> Self::Output {
        self *= rhs;
        return self;
    }
}

impl Mul<BigInt> for BigInt {
    type Output = BigInt;

    fn mul(self, rhs: BigInt) -> Self::Output {
        BigInt::RING.mul(self, rhs)
    }
}

impl MulAssign<i64> for BigInt {

    fn mul_assign(&mut self, rhs: i64) {
        self.abs_multiplication_small(rhs.abs() as u64);
        self.negative ^= rhs < 0;
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

    fn is_divisibility_computable(&self) -> bool {
        true
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
        primes::miller_rabin(el, 10)
    }

    fn calc_factor(&self, el: &Self::El) -> Option<Self::El> {
        primes::calc_factor(el)
    }
}

#[cfg(test)]
use std::str::FromStr;
#[cfg(test)]
use super::wrapper::*;
#[cfg(test)]
use vector_map::VecMap;

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
    assert_eq!(x, 256);
}

#[test]
fn test_sub_assign() {
    let mut x = "4294836225".parse::<BigInt>().unwrap();
    let y =     "4294967297".parse::<BigInt>().unwrap();
    let z =        "-131072".parse::<BigInt>().unwrap();
    x = x - y;
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
            assert_eq!(BigInt::RING.one(), BigInt::RING.euclidean_div(ns[i].clone(), &ns[i]));
        }
        assert_eq!(ns[i], ns[i].clone() + BigInt::ZERO);
        assert_eq!(ns[i], ns[i].clone() * BigInt::RING.one());
    }
    for i in 0..l {
        for j in 0..l {
            if !ns[j].is_zero() {
                assert_eq!(ns[i], BigInt::RING.euclidean_div(ns[i].clone(), &ns[j]) * ns[j].clone() + BigInt::RING.euclidean_rem(ns[i].clone(), &ns[j]));
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
    let y = BigInt::from_float_approx(x).unwrap().to_float_approx();
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

    let a = "98711".parse::<BigInt>().unwrap();
    let mut b = a.clone();
    b.data.push(0);
    assert!(a == 98711);
    assert!(a == b);
    assert!(b == a);
    assert!(a != -a.clone());
    assert!(calculate_hash(&a) == calculate_hash(&b));
    assert!(a != BigInt::RING.one());
    // the next line could theoretically fail, but it is very improbable and we definitly should test hash inequality
    assert!(calculate_hash(&a) != calculate_hash(&BigInt::RING.one()));
}

#[test]
fn test_is_zero() {
    let zero = BigInt::from(0);
    let nonzero = BigInt::power_of_two(83124);
    assert!(BigInt::RING.is_zero(&zero));
    assert!(BigInt::RING.is_zero(&(-zero)));
    assert!(!BigInt::RING.is_zero(&nonzero));
    assert!(!BigInt::RING.is_zero(&(-nonzero)));
}

#[test]
fn test_cmp_small() {
    assert!("-23678".parse::<BigInt>().unwrap() < 0);
}

#[test]
fn test_find_zero_floor() {
    let f = |x: &BigInt| BigInt::RING.mul_ref(x, x) - 234867;
    assert_eq!(BigInt::from(484), BigInt::find_zero_floor(f, BigInt::ZERO));

    let f = |x: &BigInt| x.clone();
    assert_eq!(BigInt::ZERO, BigInt::find_zero_floor(f, BigInt::ZERO));
}

#[test]
fn test_root_floor() {
    let n = BigInt::from(7681).pow(32);
    assert_eq!(BigInt::from(7681), n.root_floor(32));
}

#[test]
fn test_factor() {
    let mut expected = VecMap::new();
    expected.insert(BigInt::RING.bind(BigInt::from(7)), 2);
    expected.insert(BigInt::RING.bind(BigInt::from(2)), 1);
    assert_eq!(expected, BigInt::RING.factor(BigInt::from(98)));
    expected = VecMap::new();
    expected.insert(BigInt::RING.bind(BigInt::from(3)), 5);
    assert_eq!(expected, BigInt::RING.factor(BigInt::from(243)));
}

#[test]
fn test_is_prime() {
    assert_eq!(false, BigInt::RING.is_prime(&BigInt::from(81)));
}

#[test]
fn test_cmp() {
    assert_eq!(true, BigInt::from(-1) < BigInt::from(2));
    assert_eq!(true, BigInt::from(1) < BigInt::from(2));
    assert_eq!(false, BigInt::from(2) < BigInt::from(2));
    assert_eq!(false, BigInt::from(3) < BigInt::from(2));
    assert_eq!(true, BigInt::from(-1) > BigInt::from(-2));
}
