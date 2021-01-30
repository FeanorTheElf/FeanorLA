use std::cmp::Ordering;

#[derive(Debug, Clone)]
struct BigInt {
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

const BIG_POWER_TEN_ZEROS: u32 = 19;
const BIG_POWER_TEN: u64 = 10000000000000000000u64;

impl BigInt {

    const BLOCK_BITS: usize = 64;
    const ZERO: BigInt = BigInt { data: Vec::new(), negative: false };

    ///
    /// Calculate self += rhs * (1 << BLOCK_BITS)^block_offset
    /// 
    fn do_addition(&mut self, rhs: &BigInt, block_offset: usize) {
        assert_eq!(self.negative, rhs.negative);

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

    fn is_zero(&self) -> bool {
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

    fn abs_cmp(&self, rhs: &BigInt) -> Ordering {
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
    /// Calculate self -= rhs * (1 << BLOCK_BITS)^block_offset
    /// 
    fn do_subtraction(&mut self, rhs: &BigInt, block_offset: usize) {
        assert_eq!(self.negative, rhs.negative);
        debug_assert!(self.abs_cmp(rhs) != Ordering::Less);

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

    fn do_multiplication(&self, rhs: &BigInt) -> BigInt {
        let mut result = BigInt {
            negative: self.negative ^ rhs.negative,
            data: Vec::with_capacity(self.highest_set_block().unwrap_or(0) + rhs.highest_set_block().unwrap_or(0) + 2)
        };
        if let Some(d) = rhs.highest_set_block() {
            for i in 0..=d {
                let mut val = self.clone();
                val.do_multiplication_small(rhs.data[i]);
                result.do_addition(&val, i);
            }
        }
        return result;
    }

    ///
    /// Same as division_step, but for self_high == rhs_high == d
    /// 
    fn division_step_last(&mut self, rhs: &BigInt, d: usize) -> u64 {
        assert!(self.data[d] != 0);
        assert!(rhs.data[d] != 0);

        let self_high_blocks: u128 = ((self.data[d] as u128) << Self::BLOCK_BITS) | (self.data[d - 1] as u128);
        let rhs_high_blocks: u128 = ((rhs.data[d] as u128) << Self::BLOCK_BITS) | (rhs.data[d - 1] as u128);

        if rhs_high_blocks == u128::MAX {
            return 0;
        } else {
            let mut quotient = (self_high_blocks / (rhs_high_blocks + 1)) as u64;
            let mut subtract = rhs.clone();
            subtract.do_multiplication_small(quotient);
            self.do_subtraction(&subtract, 0);
            if self.abs_cmp(&rhs) != Ordering::Less {
                self.do_subtraction(&rhs, 0);
                quotient += 1;
            }
            if self.abs_cmp(&rhs) != Ordering::Less {
                self.do_subtraction(&rhs, 0);
                quotient += 1;
            }
            // we have been at most 2 away from the real quotient, so here it must be done so far
            debug_assert!(self.abs_cmp(&rhs) == Ordering::Less);
            return quotient;
        }
    }

    ///
    /// Finds some integer d such that subtracting d * rhs from self clears the top
    /// block of self. self will be assigned the value after the subtraction and d
    /// will be returned as d = (u * 2^block_bits + l) * 2^(k log_bits) where the return value is (u, l, k)
    /// 
    fn division_step(&mut self, rhs: &BigInt, self_high: usize, rhs_high: usize) -> (u64, u64, usize) {
        assert!(self_high > rhs_high);
        assert!(self.data[self_high] != 0);
        assert!(rhs.data[rhs_high] != 0);

            // the best approximation of the quotient we get through self_high_blocks / (rhs_high_blocks + 1)
            // the +1 is required to ensure that the quotient is smaller than the real quotient
            // one can prove that this quotient is at most 2 away from the actual quotient
            
        let mut result_upper = 0;
        let mut result_lower = 0;

        {
            let self_high_blocks: u128 = ((self.data[self_high] as u128) << Self::BLOCK_BITS) | (self.data[self_high - 1] as u128);
            let rhs_high_blocks: u128 = ((rhs.data[rhs_high] as u128) << Self::BLOCK_BITS) | (rhs.data[rhs_high - 1] as u128);

            if rhs_high_blocks != u128::MAX && self_high_blocks >= (rhs_high_blocks + 1) {
                let quotient = (self_high_blocks / (rhs_high_blocks + 1)) as u64;
                debug_assert!(quotient != 0);
                let mut subtract = rhs.clone();
                subtract.do_multiplication_small(quotient);
                self.do_subtraction(&subtract, self_high - rhs_high);
                result_upper = quotient;
            }
        }

        {
            let self_high_blocks: u128 = ((self.data[self_high] as u128) << Self::BLOCK_BITS) | (self.data[self_high - 1] as u128);

            if self.data[self_high] != 0 {
                let mut quotient = (self_high_blocks / (rhs.data[rhs_high] as u128 + 1)) as u64;
                let mut subtract = rhs.clone();
                subtract.do_multiplication_small(quotient);
                self.do_subtraction(&subtract, self_high - rhs_high - 1);
                
                // we might be up to 2 away from the real quotient
                if self.data[self_high] != 0 {
                    self.do_subtraction(rhs, self_high - rhs_high - 1);
                    quotient += 1;
                }
                if self.data[self_high] != 0 {
                    self.do_subtraction(rhs, self_high - rhs_high - 1);
                    quotient += 1;
                }
                result_lower = quotient;
            }
            debug_assert!(self.data[self_high] == 0);
            return (result_upper, result_lower, self_high - rhs_high - 1);
        }
    }

    fn do_division(&mut self, rhs: &BigInt) -> BigInt {
        assert!(!self.negative);
        assert!(!rhs.negative);

        if let Some(mut d) = self.highest_set_block() {
            let k = rhs.highest_set_block().expect("Division by zero");
            if d < k {
                return Self::ZERO.clone();
            } else if k == 0 {
                let rem = self.truncating_div_small(rhs.data[0]);
                let rem_data = vec![rem];
                let div_data = std::mem::replace(&mut self.data, rem_data);
                return BigInt {
                    negative: false,
                    data: div_data
                };
            }
            let mut result_data = Vec::new();
            result_data.resize(d + 1 - k, 0);
            while d > k {
                if self.data[d] != 0 {
                    let (quo_upper, quo_lower, quo_power) = self.division_step(&rhs, d, k);
                    result_data[quo_power] += quo_lower;
                    result_data[quo_power + 1] += quo_upper;
                    debug_assert!(self.data[d] == 0);
                }
                d -= 1;
            }
            let quo = self.division_step_last(&rhs, d);
            result_data[0] += quo;
            return BigInt {
                negative: false,
                data: result_data
            };
        } else {
            return Self::ZERO.clone();
        }
    }

    ///
    /// Calculates self /= divisor and returns the remainder of the division.
    /// This only works for positive numbers, as for negative numbers, as the remainder
    /// must be returned as a u64 to avoid overflow
    /// 
    fn truncating_div_small(&mut self, divisor: u64) -> u64 {
        assert!(!self.negative);
        assert!(divisor != 0);
        let highest_block_opt = self.highest_set_block();
        if highest_block_opt == Some(0) {
            let (quo, rem) = div_rem(self.data[0], divisor);
            self.data[0] = quo;
            return rem;
        } else if let Some(highest_block) = highest_block_opt {
            let mut buffer: u128 = self.data[highest_block] as u128;
            self.data[highest_block] = 0;
            for i in (0..highest_block).rev() {
                buffer = (buffer << Self::BLOCK_BITS) | (self.data[i] as u128);
                let (quo, rem) = div_rem(buffer, divisor as u128);
                self.data[i] = quo as u64;
                buffer = rem;
            }
            return buffer as u64;
        } else {
            return 0;
        }
    }

    fn do_multiplication_small(&mut self, factor: u64) {
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
    /// Calculates self += summand * (1 << BLOCK_BITS)^block_offset
    /// 
    fn do_addition_small(&mut self, summand: u64, block_offset: usize) {
        assert!(!self.negative);

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
    fn from_radix<I, E>(data: I, base: u64) -> Result<BigInt, E> 
        where I: Iterator<Item = Result<u64, E>>
    {
        let mut result = BigInt {
            data: Vec::with_capacity(data.size_hint().0),
            negative: false
        };
        for value in data {
            let val = value?;
            debug_assert!(val < base);
            result.do_multiplication_small(base);
            result.do_addition_small(val, 0);
        }
        return Ok(result);
    }

    fn from_str_radix(s: &str, base: u32) -> Result<BigInt, BigIntParseError> {
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
}

impl PartialEq for BigInt {

    fn eq(&self, rhs: &BigInt) -> bool {
        let highest_block = self.highest_set_block();
        if highest_block != rhs.highest_set_block() {
            return false;
        }
        if let Some(d) = highest_block {
            for i in 0..d {
                if self.data[i] != rhs.data[i] {
                    return false;
                }
            }
        }
        return true;
    }
}

impl Eq for BigInt {}

impl PartialEq<u64> for BigInt {

    fn eq(&self, rhs: &u64) -> bool {
        self.highest_set_block() == Some(0) && self.data[0] == *rhs
    }
}

impl std::fmt::Display for BigInt {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let mut copy = self.clone();
        let mut remainders: Vec<u64> = Vec::with_capacity((self.highest_set_block().unwrap_or(0) + 1) * Self::BLOCK_BITS / 3);
        while !copy.is_zero() {
            let rem = copy.truncating_div_small(BIG_POWER_TEN);
            remainders.push(rem);
        }
        remainders.reverse();
        let mut it = remainders.into_iter();
        if let Some(fst) = it.next() {
            write!(f, "{}", fst);
            for rem in it {
                write!(f, "{:0>width$}", rem, width = BIG_POWER_TEN_ZEROS as usize)?;
            }
        } else {
            write!(f, "0")?;
        }
        return Ok(());
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
fn test_do_subtraction() {
    let mut x = "923645871236598172365987287530543".parse::<BigInt>().unwrap();
    let y =      "58430657823473456743684735863478".parse::<BigInt>().unwrap();
    let z =     "865215213413124715622302551667065".parse::<BigInt>().unwrap();
    x.do_subtraction(&y, 0);
    assert_eq!(z, x);
}

#[test]
fn test_do_subtraction_carry() {
    let mut x = BigInt::from_str_radix("1000000000000000000", 16).unwrap();
    let y =      BigInt::from_str_radix("FFFFFFFFFFFFFFFF00", 16).unwrap();
    x.do_subtraction(&y, 0);
    assert_eq!(x, 256);
}

#[test]
fn test_do_addition() {
    let mut x = "923645871236598172365987287530543".parse::<BigInt>().unwrap();
    let y =      "58430657823473456743684735863478".parse::<BigInt>().unwrap();
    let z =     "982076529060071629109672023394021".parse::<BigInt>().unwrap();
    x.do_addition(&y, 0);
    assert_eq!(z, x);
}

#[test]
fn test_do_addition_carry() {
    let mut x =             BigInt::from_str_radix("1BC00000000000000BC", 16).unwrap();
    let y =  BigInt::from_str_radix("FFFFFFFFFFFFFFFF0000000000000000BC", 16).unwrap();
    let z = BigInt::from_str_radix("10000000000000000BC0000000000000178", 16).unwrap();
    x.do_addition(&y, 0);
    assert_eq!(z, x);
}

#[test]
fn test_do_multiplication() {
    let x = BigInt::from_str_radix("57873674586797895671345345", 10).unwrap();
    let y = BigInt::from_str_radix("21308561789045691782534873921650342768903561413264128756389247568729346542359871235465", 10).unwrap();
    let z = BigInt::from_str_radix("1233204770891906354921751949503652431220138020953161094405729272872607166072371117664593787957056214903826660425", 10).unwrap();
    assert_eq!(z, x.do_multiplication(&y));
}

#[test]
fn test_do_division_no_remainder() {
    let mut x = BigInt::from_str_radix("578435387FF0582367863200000000000000000000", 16).unwrap();
    let y =                          BigInt::from_str_radix("200000000000000000000", 16).unwrap();
    let z = BigInt::from_str_radix("2BC21A9C3FF82C11B3C319", 16).unwrap();
    let quotient = x.do_division(&y);
    assert_eq!(BigInt::ZERO, x);
    assert_eq!(z, quotient);
}