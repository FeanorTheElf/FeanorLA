use super::bigint::BigInt;
use super::alg::*;
use super::eea::signed_eea;

use std::ops::{AddAssign, MulAssign, SubAssign, DivAssign, Add, Mul, Sub, Div, Neg};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GeneralZnEl<const N: u64, const IS_FIELD: bool> {
    repr: u64
}

const fn is_prime(n: u64) -> bool {
    assert!(n >= 2);
    const fn is_b_smooth(b: u64, n: u64) -> bool {
        if b <= 1 {
            true
        } else if n % b == 0 {
            false
        } else {
            is_b_smooth(b - 1, n)
        }
    }
    return is_b_smooth(n - 1, n);
}

impl<const N: u64, const IS_FIELD: bool> GeneralZnEl<N, IS_FIELD> {

    pub const fn project(v: i64) -> GeneralZnEl<N, IS_FIELD> {
        assert!(N <= u64::MAX / 2);
        assert!(N <= i64::MAX as u64);
        assert!(IS_FIELD == is_prime(N));
        let mut result = v % N as i64;
        if result < 0 {
            result += N as i64;
        }
        debug_assert!((result as u64) < N);
        GeneralZnEl {
            repr: result as u64
        }
    }

    pub const ZERO: GeneralZnEl<N, IS_FIELD> = GeneralZnEl::project(0);
    pub const ONE: GeneralZnEl<N, IS_FIELD> = GeneralZnEl::project(1);
}

impl<const N: u64, const IS_FIELD: bool> AddAssign for GeneralZnEl<N, IS_FIELD> {

    fn add_assign(&mut self, rhs: GeneralZnEl<N, IS_FIELD>) {
        debug_assert!(self.repr < N);
        debug_assert!(rhs.repr < N);
        self.repr += rhs.repr;
        if self.repr >= N {
            self.repr -= N;
        }
        debug_assert!(self.repr < N);
    }
}

impl<const N: u64, const IS_FIELD: bool> SubAssign for GeneralZnEl<N, IS_FIELD> {

    fn sub_assign(&mut self, rhs: GeneralZnEl<N, IS_FIELD>) {
        debug_assert!(self.repr < N);
        debug_assert!(rhs.repr < N);
        self.repr += N;
        self.repr -= rhs.repr;
        if self.repr >= N {
            self.repr -= N;
        }
        debug_assert!(self.repr < N);
    }
}

impl<const N: u64, const IS_FIELD: bool> MulAssign for GeneralZnEl<N, IS_FIELD> {

    fn mul_assign(&mut self, rhs: GeneralZnEl<N, IS_FIELD>) {
        debug_assert!(self.repr < N);
        debug_assert!(rhs.repr < N);
        self.repr = ((self.repr as u128 * rhs.repr as u128) % N as u128) as u64;
        debug_assert!(self.repr < N);
    }
}

impl<const N: u64> DivAssign for GeneralZnEl<N, true> {

    fn div_assign(&mut self, rhs: GeneralZnEl<N, true>) {
        assert!(is_prime(N));
        debug_assert!(self.repr < N);
        debug_assert!(rhs.repr < N);
        let (s, _t, _d) = signed_eea(rhs.repr as i64, N as i64);
        let mut result = ((s as i128 * self.repr as i128) % N as i128) as i64;
        if result < 0 {
            result += N as i64;
        }
        self.repr = result as u64;
        debug_assert!(self.repr < N);
    }
}

impl<const N: u64, const IS_FIELD: bool> Add for GeneralZnEl<N, IS_FIELD> {

    type Output = GeneralZnEl<N, IS_FIELD>;

    fn add(mut self, rhs: GeneralZnEl<N, IS_FIELD>) -> Self::Output {
        self += rhs;
        self
    }
}

impl<const N: u64, const IS_FIELD: bool> Sub for GeneralZnEl<N, IS_FIELD> {

    type Output = GeneralZnEl<N, IS_FIELD>;

    fn sub(mut self, rhs: GeneralZnEl<N, IS_FIELD>) -> Self::Output {
        self -= rhs;
        self
    }
}

impl<const N: u64, const IS_FIELD: bool> Mul for GeneralZnEl<N, IS_FIELD> {

    type Output = GeneralZnEl<N, IS_FIELD>;

    fn mul(mut self, rhs: GeneralZnEl<N, IS_FIELD>) -> Self::Output {
        self *= rhs;
        self
    }
}

impl<const N: u64> Div for GeneralZnEl<N, true> {

    type Output = GeneralZnEl<N, true>;

    fn div(mut self, rhs: GeneralZnEl<N, true>) -> Self::Output {
        self /= rhs;
        self
    }
}

impl<const N: u64, const IS_FIELD: bool> Neg for GeneralZnEl<N, IS_FIELD> {

    type Output = GeneralZnEl<N, IS_FIELD>;

    fn neg(mut self) -> Self::Output {
        assert!(self.repr < N);
        if self.repr != 0 {
            self.repr = N - self.repr;
        }
        debug_assert!(self.repr < N);
        return self;
    }
}

impl<const N: u64, const IS_FIELD: bool> std::fmt::Display for GeneralZnEl<N, IS_FIELD> {

    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "[{}]", self.repr)
    }
}

impl<const N: u64, const IS_FIELD: bool> Zero for GeneralZnEl<N, IS_FIELD> {

    fn zero() -> GeneralZnEl<N, IS_FIELD> {
        Self::project(0)
    }
}

impl<const N: u64, const IS_FIELD: bool> One for GeneralZnEl<N, IS_FIELD> {

    fn one() -> GeneralZnEl<N, IS_FIELD> {
        Self::project(1)
    }
}

pub type ZnEl<const N: u64> = GeneralZnEl<N, {is_prime(N)}>;


impl<const N: u64> RingEl for GeneralZnEl<N, true> {
    type Axioms = RingAxiomsField;
}

impl<const N: u64> RingEl for GeneralZnEl<N, false> {
    type Axioms = RingAxiomsZeroDivisorRing;
}

impl<const N: u64> FieldEl for GeneralZnEl<N, true> {}

#[test]
fn test_is_prime() {
    assert_eq!(true, is_prime(17));
    assert_eq!(false, is_prime(49));
}

#[cfg(test)]
type Z17 = ZnEl<17>;

#[test]
fn test_zn_el_add() {
    let a = Z17::project(6);
    let b = Z17::project(12);
    assert_eq!(Z17::project(1), a + b);
}

#[test]
fn test_zn_el_sub() {
    let a = Z17::project(6);
    let b = Z17::project(12);
    assert_eq!(Z17::project(11), a - b);
}

#[test]
fn test_zn_el_mul() {
    let a = Z17::project(6);
    let b = Z17::project(12);
    assert_eq!(Z17::project(4), a * b);
}

#[test]
fn test_zn_el_div() {
    let a = Z17::project(6);
    let b = Z17::project(12);
    assert_eq!(Z17::project(9), a / b);
}

#[test]
fn fn_test_div_impossible() {
    let _a = ZnEl::<18>::project(4);
    // the following line should give a compiler error
    // _a / _a;
}