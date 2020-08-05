#![allow(non_camel_case_types)]

use super::arith::{One, Zero};

use std::cmp::{Ord, Ordering, PartialEq, PartialOrd};
use std::convert::From;
use std::fmt::{Debug, Display, Formatter};
use std::ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign};

///
/// Overflow contract: r64 may overflow, if the naive formulas for the
/// computations would cause an overflow, given all participating fractions
/// are completly reduced and i64 is used. Only in this case, the operation
/// may panic, otherwise it has to yield the correct result.
/// The naive formulas are:
///  - Add: (a.num * b.den + b.num * a.den) / (a.den * b.den)
///  - Mul: (a.num * b.num) / (a.den * b.den)
///  - Eq: never overflows
///
/// Fractions may not be reduced as long as possible for performance reasons,
/// and a debug output may print a non-reduced fraction. Display output will
/// always print a correctly reduced form of this fraction.
///
#[derive(Clone, Copy)]
pub struct r64 {
    numerator: i64,
    denominator: i64,
}

pub const NAN: r64 = r64 {
    numerator: 0,
    denominator: 0,
};
pub const INFINITY: r64 = r64 {
    numerator: 1,
    denominator: 0,
};
pub const ZERO: r64 = r64 {
    numerator: 0,
    denominator: 1,
};
pub const ONE: r64 = r64 {
    numerator: 1,
    denominator: 1,
};

impl r64 {
    pub fn new(numerator: i64, denominator: i64) -> r64 {
        r64 {
            numerator: numerator,
            denominator: denominator,
        }
    }

    // a < 0 => ggT(a, b) < 0
    // a > 0 => ggT(a, b) > 0
    // sign of b is irrelevant
    // ggT(0, 0) = 1
    fn ggT(mut a: i64, mut b: i64) -> i64 {
        while b != 0 {
            let c = a % b;
            if c == 0 {
                return if (a ^ b) & i64::min_value() != 0 {
                    -b
                } else {
                    b
                };
            }
            b = b % c;
            a = c;
        }
        return if a == 0 { 1 } else { a };
    }

    pub fn reduce(&mut self) {
        let ggT: i64 = r64::ggT(self.denominator, self.numerator);
        self.denominator /= ggT;
        self.numerator /= ggT;
    }

    pub fn compare_zero(&self) -> Option<Ordering> {
        if self.denominator == 0 && self.numerator == 0 {
            None
        } else {
            Some((self.numerator * self.denominator.signum()).cmp(&0))
        }
    }

    pub fn is_nan(&self) -> bool {
        self.denominator == 0 && self.numerator == 0
    }

    pub fn num(&self) -> i64 {
        self.numerator
    }

    pub fn den(&self) -> i64 {
        self.denominator
    }
}

impl From<i64> for r64 {
    fn from(value: i64) -> Self {
        r64::new(value, 1)
    }
}

impl<'a> From<&'a i32> for r64 {
    fn from(value: &'a i32) -> Self {
        r64::new(*value as i64, 1)
    }
}

impl Zero for r64 {
    fn zero() -> Self {
        ZERO
    }
}

impl One for r64 {
    fn one() -> Self {
        ONE
    }
}

macro_rules! assign_or_reduce_on_failure {
    ($e:expr => $target:expr; $reduce_fst:stmt; $reduce_snd:stmt; $reduce_trd:stmt; $reduce_fth:stmt) => {
        if let Some(value) = $e {
            $target = value;
        } else {
            $reduce_fst
            if let Some(value) = $e {
                $target = value;
            } else {
                $reduce_snd
                if let Some(value) = $e {
                    $target = value;
                } else {
                    $reduce_trd
                    if let Some(value) = $e {
                        $target = value;
                    } else {
                        $reduce_fth
                        $target = ($e).unwrap();
                    }
                }
            }
        }
    };
    ($e:expr => $target:expr; $reduce_fst:stmt; $reduce_snd:stmt) => {
        if let Some(value) = $e {
            $target = value;
        } else {
            $reduce_fst
            if let Some(value) = $e {
                $target = value;
            } else {
                $reduce_snd
                $target = ($e).unwrap();
            }
        }
    };
    ($e:expr => $target:expr; $reduce_fst:expr) => {
        if let Some(value) = $e {
            $target = value;
        } else {
            $reduce_fst
            $target = ($e).unwrap();
        }
    };
}

impl AddAssign<r64> for r64 {
    fn add_assign(&mut self, mut rhs: r64) {
        // inf + inf or -inf + -inf
        if self.denominator == 0
            && rhs.denominator == 0
            && (self.numerator > 0 && rhs.numerator > 0 || self.numerator < 0 && rhs.numerator < 0)
        {
            return;
        }
        assign_or_reduce_on_failure!(self.numerator.checked_mul(rhs.denominator) => self.numerator; self.reduce(); rhs.reduce());
        assign_or_reduce_on_failure!(self.denominator.checked_mul(rhs.numerator) => rhs.numerator; self.reduce(); {
            let ggT = r64::ggT(rhs.denominator, rhs.numerator);
            rhs.denominator /= ggT;
            rhs.numerator /= ggT;
            self.numerator /= ggT;
        });
        assign_or_reduce_on_failure!(self.numerator.checked_add(rhs.numerator) => self.numerator; {
            let ggT = r64::ggT(self.denominator, self.numerator);
            self.denominator /= ggT;
            self.numerator /= ggT;
            rhs.numerator /= ggT;
        }; {
            let ggT = r64::ggT(rhs.denominator, rhs.numerator);
            rhs.denominator /= ggT;
            rhs.numerator /= ggT;
            self.numerator /= ggT;
        });
        assign_or_reduce_on_failure!(self.denominator.checked_mul(rhs.denominator) => self.denominator; self.reduce(); {
            let ggT = r64::ggT(self.numerator, rhs.denominator);
            self.numerator /= ggT;
            rhs.denominator /= ggT;
        });
    }
}

impl MulAssign<r64> for r64 {
    fn mul_assign(&mut self, mut rhs: r64) {
        assign_or_reduce_on_failure!(self.numerator.checked_mul(rhs.numerator) => self.numerator; self.reduce(); rhs.reduce());
        assign_or_reduce_on_failure!(self.denominator.checked_mul(rhs.denominator) => self.denominator; self.reduce(); {
            let ggT = r64::ggT(rhs.numerator, rhs.denominator);
            rhs.denominator /= ggT;
            self.numerator /= ggT;
        });
    }
}

impl DivAssign<r64> for r64 {
    fn div_assign(&mut self, rhs: r64) {
        self.mul_assign(r64::new(rhs.denominator, rhs.numerator));
    }
}

impl SubAssign<r64> for r64 {
    fn sub_assign(&mut self, rhs: r64) {
        self.add_assign(-rhs)
    }
}

impl Neg for r64 {
    type Output = r64;

    fn neg(mut self) -> r64 {
        self.numerator = -self.numerator;
        return self;
    }
}

impl PartialEq for r64 {
    fn eq(&self, rhs: &r64) -> bool {
        if self.denominator == 0 || rhs.denominator == 0 {
            self.numerator != 0
                && rhs.numerator != 0
                && (self.numerator > 0) == (rhs.numerator > 0)
                && self.denominator == rhs.denominator
        } else {
            self.numerator as i128 * rhs.denominator as i128
                == self.denominator as i128 * rhs.numerator as i128
        }
    }
}

impl PartialOrd for r64 {
    fn partial_cmp(&self, rhs: &r64) -> Option<Ordering> {
        if self.is_nan() || rhs.is_nan() {
            return None;
        } else if self.denominator == 0 && rhs.denominator == 0 {
            return Some(self.numerator.cmp(&rhs.numerator));
        } else {
            let ord = (self.numerator as i128 * rhs.denominator as i128)
                .cmp(&(self.denominator as i128 * rhs.numerator as i128));
            return Some(if (self.denominator < 0) ^ (rhs.denominator < 0) {
                ord.reverse()
            } else {
                ord
            });
        }
    }
}

impl Add<Self> for r64 {
    type Output = r64;

    fn add(mut self, rhs: r64) -> r64 {
        self += rhs;
        return self;
    }
}

impl Mul<Self> for r64 {
    type Output = r64;

    fn mul(mut self, rhs: r64) -> r64 {
        self *= rhs;
        return self;
    }
}

impl Sub<Self> for r64 {
    type Output = r64;

    fn sub(mut self, rhs: r64) -> r64 {
        self -= rhs;
        return self;
    }
}

impl Div<Self> for r64 {
    type Output = r64;

    fn div(mut self, rhs: r64) -> r64 {
        self /= rhs;
        return self;
    }
}

impl Debug for r64 {
    fn fmt(&self, f: &mut Formatter) -> std::fmt::Result {
        write!(f, "({}/{})", self.numerator, self.denominator)
    }
}

impl Display for r64 {
    fn fmt(&self, f: &mut Formatter) -> std::fmt::Result {
        let mut reduced: r64 = *self;
        reduced.reduce();
        return write!(f, "({}/{})", reduced.numerator, reduced.denominator);
    }
}

#[test]
fn test_ggT() {
    assert_eq!(3, r64::ggT(15, 6));
    assert_eq!(3, r64::ggT(6, 15));

    assert_eq!(7, r64::ggT(0, 7));
    assert_eq!(7, r64::ggT(7, 0));
    assert_eq!(1, r64::ggT(0, 0));

    assert_eq!(1, r64::ggT(9, 1));
    assert_eq!(1, r64::ggT(1, 9));

    assert_eq!(1, r64::ggT(13, 300));
    assert_eq!(1, r64::ggT(300, 13));

    assert_eq!(-3, r64::ggT(-15, 6));
    assert_eq!(3, r64::ggT(6, -15));
    assert_eq!(-3, r64::ggT(-6, -15));
}

#[test]
fn test_add_assign() {
    let mut a: r64 = r64::from(1);
    a += a;
    assert_eq!(r64::from(2), a);
    let b = r64::new(1, 2);
    a -= b;
    assert_eq!(r64::new(3, 2), a);
}

#[test]
fn test_add_assign_overflow() {
    let mut a = r64::new(-(1 << 60), 81);
    let mut b = r64::new(1 << 60, -(1 << 55));
    a += b;
    assert_eq!(r64::new((1 << 5) * 81 + (1 << 60), -81), a);

    a = r64::new((1 << 62) / 50 /* divisible by 3 */, -81);
    b = r64::new((1 << 62) / 81, 50);
    a -= b;
    assert_eq!(
        r64::new(((1 << 62) / 150) * 50 + ((1 << 62) / 81) * 27, -27 * 50),
        a
    );
}

#[test]
fn test_mul_assign_overflow() {
    let mut a = r64::new(1 << 52, 1 << 62);
    let b = r64::new(1 << 50, -(1 << 60));
    a *= b;
    assert_eq!(r64::new(-2, 1 << 21), a);

    let c = r64::new(1 << 60, 1 << 20);
    a /= c;
    assert_eq!(r64::new(-1, 1 << 60), a);
}

#[bench]
fn benchmark_combined_add_mult_eq(bencher: &mut test::Bencher) {
    let numerator_a = (1 << 62) / 50; // divisible by 3
    let numerator_b = (1 << 62) / 81;
    let result_numerator = 44211 * 69540552025927 + 1350;
    let not_optimized: i64 = std::hint::black_box(0);
    bencher.iter(|| {
        let mut a = r64::new(numerator_a + not_optimized, 81);
        let b = r64::new(numerator_b + not_optimized, 50);
        a += b;
        let c = r64::new(1 + not_optimized, 44211);
        a *= c;
        assert_eq!(r64::new(69540552025927 + not_optimized, 1350), a);
        if b == a {
            a /= r64::from(100);
        } else {
            a /= c;
        }
        let mut d = r64::new(32 + not_optimized, 1024);
        d *= d;
        d *= d;
        d *= d;
        d *= c;
        if r64::new(1 + not_optimized, (1 << 40) * 44211) == d {
            a += r64::from(1);
        }
        assert_eq!(r64::new(result_numerator, 1350), a);
    });
}

#[test]
#[should_panic]
fn test_real_overflow() {
    let mut a = r64::new(1, 2 << 33);
    let b = r64::new(1, (2 << 32) + 1);
    a += b;
}

#[test]
fn test_cmp_NaN_Inf() {
    assert_eq!(INFINITY, r64::from(1) / r64::from(0));
    assert_eq!(-INFINITY, r64::from(-2) / r64::from(0));
    assert_eq!(INFINITY, INFINITY);
    assert_eq!(-INFINITY, -INFINITY);
    assert_ne!(INFINITY, -INFINITY);
    assert_ne!(-INFINITY, INFINITY);
    assert_ne!(r64::from(0), INFINITY);
    assert_ne!(r64::from(1), INFINITY);
    assert_ne!(r64::from(0), -INFINITY);
    assert_ne!(r64::from(-1), -INFINITY);
    assert_ne!(NAN, NAN);
    assert_ne!(INFINITY, NAN);
    assert_ne!(NAN, -INFINITY);
    assert_eq!(r64::from(0), r64::from(0));
    assert_ne!(r64::from(0), r64::from(1));
}

#[test]
fn test_calculate_NaN_Inf() {
    assert_eq!(INFINITY, INFINITY + INFINITY);
    assert!((INFINITY - INFINITY).is_nan());
    assert_eq!(-INFINITY, r64::from(1) - INFINITY);
    assert_eq!(-INFINITY, -INFINITY - INFINITY);
    assert_eq!(INFINITY, r64::from(-100) * -INFINITY);
    assert_eq!(-INFINITY, r64::from(-3) * INFINITY);
    assert!((INFINITY * r64::from(0)).is_nan());
    assert!((-INFINITY * r64::from(0)).is_nan());
}

#[test]
fn test_ord_NaN_Inf() {
    assert!(-INFINITY < INFINITY);
    assert!(r64::from(10) < INFINITY);
    assert!(r64::from(0) < INFINITY);
    assert!(r64::from(-10) < INFINITY);
    assert!(r64::from(10) > -INFINITY);
    assert!(r64::from(0) > -INFINITY);
    assert!(r64::from(-10) > -INFINITY);
    assert!(!(r64::from(0) < NAN));
    assert!(!(r64::from(0) == NAN));
    assert!(!(r64::from(0) > NAN));
}

#[test]
fn test_ord() {
    assert!(r64::from(-1) < r64::from(0));
    assert!(r64::from(1) > r64::from(0));
    assert!(r64::new(1, 2) < r64::from(1));
    assert!(r64::new(1, -2) < r64::new(-1, 3));
}
