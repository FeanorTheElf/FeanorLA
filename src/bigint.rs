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

impl BigInt {

    const ENTRY_BITS: usize = 64;

    fn do_addition(&mut self, rhs: &BigInt) {
        assert_eq!(self.negative, rhs.negative);

        let mut buffer: bool = false;
        let mut i = 0;
        while i < rhs.data.len() || buffer {
            let rhs_val = *rhs.data.get(i).unwrap_or(&0);
            if i == self.data.len() {
                self.data.push(0);
            }
            let (sum, overflow) = self.data[i].overflowing_add(rhs_val);
            if buffer {
                let (carry_sum, carry_overflow) = sum.overflowing_add(1);
                self.data[i] = carry_sum;
                buffer = overflow || carry_overflow;
            } else {
                self.data[i] = sum;
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

    fn do_subtraction(&mut self, rhs: &BigInt) {
        assert_eq!(self.negative, rhs.negative);
        debug_assert!(self.abs_cmp(rhs) != Ordering::Less);
        
        let mut buffer: bool = false;
        let mut i = 0;
        while i < rhs.data.len() || buffer {
            let rhs_val = *rhs.data.get(i).unwrap_or(&0);
            debug_assert!(i < self.data.len());
            let (difference, overflow) = self.data[i].overflowing_sub(rhs_val);
            if buffer {
                let (carry_difference, carry_overflow) = difference.overflowing_sub(1);
                self.data[i] = carry_difference;
                buffer = overflow || carry_overflow;
            } else {
                self.data[i] = difference;
                buffer = overflow;
            }
            i += 1;
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
                buffer = (buffer << Self::ENTRY_BITS) | (self.data[i] as u128);
                let (quo, rem) = div_rem(buffer, divisor as u128);
                self.data[i] = quo as u64;
                buffer = rem;
            }
            return buffer as u64;
        } else {
            return 0;
        }
    }
}

impl std::fmt::Display for BigInt {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let mut copy = self.clone();
        const BIG_POWER_TEN_ZEROS: u32 = 19;
        const BIG_POWER_TEN: u64 = 10000000000000000000u64;
        let mut remainders: Vec<u64> = Vec::with_capacity((self.highest_set_block().unwrap_or(0) + 1) * Self::ENTRY_BITS / 3);
        while !copy.is_zero() {
            let rem = copy.truncating_div_small(BIG_POWER_TEN);
            remainders.push(rem);
        }
        remainders.reverse();
        let mut it = remainders.into_iter();
        if let Some(fst) = it.next() {
            write!(f, "{}", fst);
            for rem in it {
                write!(f, "{0:>width$}", rem, width = BIG_POWER_TEN_ZEROS as usize)?;
            }
        } else {
            write!(f, "0")?;
        }
        return Ok(());
    }
}

#[test]
fn test_print_power_2() {
    let x = BigInt {
        negative: false,
        data: vec![0, 0, 1]
    };
    assert_eq!("340282366920938463463374607431768211456".to_owned(), format!("{}", x));
}