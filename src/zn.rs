use super::bigint::BigInt;
use super::alg::*;
use super::eea::signed_eea;

use std::ops::{AddAssign, MulAssign, SubAssign, DivAssign, Add, Mul, Sub, Div, Neg};

pub struct RingZn {
    modulus: BigInt,
    inverse_modulus: BigInt,
    inverse_modulus_bitshift: usize
}

impl RingZn {

    pub fn new(modulus: BigInt) -> Self {
        // have k such that 2^k > modulus^2
        // then (2^k / modulus) * x >> k differs at most 1 from floor(x / modulus)
        // if x < n^2, which is the case after multiplication
        let k = modulus.log2_floor() * 2 + 1;
        let inverse_modulus = BigInt::power_of_two(k) / modulus.clone();
        return RingZn {
            modulus: modulus,
            inverse_modulus: inverse_modulus,
            inverse_modulus_bitshift: k
        };
    }

    fn project_leq_n_square(&self, mut n: BigInt) -> BigInt {
        let mut subtract = n.clone();
        subtract *= &self.inverse_modulus;
        subtract = subtract >> self.inverse_modulus_bitshift;
        subtract *= &self.modulus;
        n -= subtract;
        if n >= self.modulus {
            n -= &self.modulus;
        }
        return n;
    }

    pub fn project(&self, n: BigInt) -> BigInt {
        let mut red_n = if n < self.modulus {
            n
        } else if n.log2_floor() + 1 < 2 * self.modulus.log2_floor() {
            self.project_leq_n_square(n)
        } else {
            n.clone() - (n / self.modulus.clone()) * self.modulus.clone()
        };
        if red_n < BigInt::ZERO {
            red_n = -red_n;
            red_n += &self.modulus;
        }
        debug_assert!(red_n < self.modulus);
        return red_n;
    }
}

impl Ring for RingZn {
    type El = BigInt;

    fn add_ref(&self, lhs: Self::El, rhs: &Self::El) -> Self::El {
        assert!(lhs < self.modulus);
        assert!(rhs < &self.modulus);

        let mut result = lhs + rhs;
        if result >= self.modulus {
            result -= &self.modulus;
        }
        return result;
    }

    fn mul_ref(&self, lhs: Self::El, rhs: &Self::El) -> Self::El {
        assert!(lhs < self.modulus);
        assert!(rhs < &self.modulus);

        return self.project_leq_n_square(lhs * rhs);
    }

    fn neg(&self, val: Self::El) -> Self::El {
        assert!(val < self.modulus);

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

    fn eq(&self, lhs: &Self::El, rhs: &Self::El) -> bool {
        assert!(lhs < &self.modulus);
        assert!(rhs < &self.modulus);
        lhs == rhs
    }
}

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
        GeneralZnEl {
            repr: result as u64
        }
    }
}

impl<const N: u64, const IS_FIELD: bool> AddAssign for GeneralZnEl<N, IS_FIELD> {

    fn add_assign(&mut self, rhs: GeneralZnEl<N, IS_FIELD>) {
        debug_assert!(self.repr < N);
        debug_assert!(rhs.repr < N);
        self.repr += rhs.repr;
        if self.repr > N {
            self.repr -= N;
        }
    }
}

impl<const N: u64, const IS_FIELD: bool> SubAssign for GeneralZnEl<N, IS_FIELD> {

    fn sub_assign(&mut self, rhs: GeneralZnEl<N, IS_FIELD>) {
        debug_assert!(self.repr < N);
        debug_assert!(rhs.repr < N);
        self.repr += N;
        self.repr -= rhs.repr;
        if self.repr > N {
            self.repr -= N;
        }
    }
}

impl<const N: u64, const IS_FIELD: bool> MulAssign for GeneralZnEl<N, IS_FIELD> {

    fn mul_assign(&mut self, rhs: GeneralZnEl<N, IS_FIELD>) {
        debug_assert!(self.repr < N);
        debug_assert!(rhs.repr < N);
        self.repr = ((self.repr as u128 * rhs.repr as u128) % N as u128) as u64;
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
        self.repr = N - self.repr;
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

impl<const N: u64, const IS_FIELD: bool> RingEl for GeneralZnEl<N, IS_FIELD> {}
impl<const N: u64> IntegralRingEl for GeneralZnEl<N, true> {}
impl<const N: u64> FieldEl for GeneralZnEl<N, true> {}

#[test]
fn test_is_prime() {
    assert_eq!(true, is_prime(17));
    assert_eq!(false, is_prime(49));
}

#[test]
fn test_mul() {
    let z257 = RingZn::new(BigInt::from_str_radix("257", 10).unwrap());
    let x = BigInt::from_str_radix("256", 10).unwrap();
    assert_eq!(BigInt::one(), z257.mul(x.clone(), x));
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