use super::super::eea::*;
use super::*;
use super::super::super::primitive::*;
use super::super::super::la::constant::*;

use std::ops::{AddAssign, MulAssign, SubAssign, DivAssign, Add, Mul, Sub, Div, Neg};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ZnElImpl<const N: u64, const IS_FIELD: bool> {
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

impl<const N: u64, const IS_FIELD: bool> ZnElImpl<N, IS_FIELD> {

    pub const IS_FIELD: bool = IS_FIELD;

    pub const fn project(v: i64) -> ZnElImpl<N, IS_FIELD> {
        assert!(N <= u64::MAX / 2);
        assert!(N <= i64::MAX as u64);
        assert!(!IS_FIELD || is_prime(N));
        let mut result = v % N as i64;
        if result < 0 {
            result += N as i64;
        }
        debug_assert!((result as u64) < N);
        ZnElImpl {
            repr: result as u64
        }
    }

    pub const ZERO: ZnElImpl<N, IS_FIELD> = ZnElImpl::project(0);
    pub const ONE: ZnElImpl<N, IS_FIELD> = ZnElImpl::project(1);
}

impl<const N: u64, const IS_FIELD: bool> From<i8> for ZnElImpl<N, IS_FIELD> {

    fn from(x: i8) -> Self {
        Self::project(x as i64)
    }
}

impl<const N: u64, const IS_FIELD: bool> AddAssign for ZnElImpl<N, IS_FIELD> {

    fn add_assign(&mut self, rhs: ZnElImpl<N, IS_FIELD>) {
        debug_assert!(self.repr < N);
        debug_assert!(rhs.repr < N);
        self.repr += rhs.repr;
        if self.repr >= N {
            self.repr -= N;
        }
        debug_assert!(self.repr < N);
    }
}

impl<const N: u64, const IS_FIELD: bool> SubAssign for ZnElImpl<N, IS_FIELD> {

    fn sub_assign(&mut self, rhs: ZnElImpl<N, IS_FIELD>) {
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

impl<const N: u64, const IS_FIELD: bool> MulAssign for ZnElImpl<N, IS_FIELD> {

    fn mul_assign(&mut self, rhs: ZnElImpl<N, IS_FIELD>) {
        debug_assert!(self.repr < N);
        debug_assert!(rhs.repr < N);
        self.repr = ((self.repr as u128 * rhs.repr as u128) % N as u128) as u64;
        debug_assert!(self.repr < N);
    }
}

impl<const N: u64> DivAssign for ZnElImpl<N, true> {

    fn div_assign(&mut self, rhs: ZnElImpl<N, true>) {
        assert!(is_prime(N));
        debug_assert!(self.repr < N);
        debug_assert!(rhs.repr < N);
        let (s, _t, _d) = signed_eea(rhs.repr as i64, N as i64, &i64::RING);
        let mut result = ((s as i128 * self.repr as i128) % N as i128) as i64;
        if result < 0 {
            result += N as i64;
        }
        self.repr = result as u64;
        debug_assert!(self.repr < N);
    }
}

impl<const N: u64, const IS_FIELD: bool> Add for ZnElImpl<N, IS_FIELD> {

    type Output = ZnElImpl<N, IS_FIELD>;

    fn add(mut self, rhs: ZnElImpl<N, IS_FIELD>) -> Self::Output {
        self += rhs;
        self
    }
}

impl<const N: u64, const IS_FIELD: bool> Sub for ZnElImpl<N, IS_FIELD> {

    type Output = ZnElImpl<N, IS_FIELD>;

    fn sub(mut self, rhs: ZnElImpl<N, IS_FIELD>) -> Self::Output {
        self -= rhs;
        self
    }
}

impl<const N: u64, const IS_FIELD: bool> Mul for ZnElImpl<N, IS_FIELD> {

    type Output = ZnElImpl<N, IS_FIELD>;

    fn mul(mut self, rhs: ZnElImpl<N, IS_FIELD>) -> Self::Output {
        self *= rhs;
        self
    }
}

impl<const N: u64> Div for ZnElImpl<N, true> {

    type Output = ZnElImpl<N, true>;

    fn div(mut self, rhs: ZnElImpl<N, true>) -> Self::Output {
        self /= rhs;
        self
    }
}

impl<const N: u64, const IS_FIELD: bool> Neg for ZnElImpl<N, IS_FIELD> {

    type Output = ZnElImpl<N, IS_FIELD>;

    fn neg(mut self) -> Self::Output {
        assert!(self.repr < N);
        if self.repr != 0 {
            self.repr = N - self.repr;
        }
        debug_assert!(self.repr < N);
        return self;
    }
}

impl<const N: u64, const IS_FIELD: bool> std::fmt::Display for ZnElImpl<N, IS_FIELD> {

    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "[{}]_{}", self.repr, N)
    }
}

impl<const N: u64, const IS_FIELD: bool> Zero for ZnElImpl<N, IS_FIELD> {

    fn zero() -> ZnElImpl<N, IS_FIELD> {
        Self::project(0)
    }
}

impl<const N: u64, const IS_FIELD: bool> One for ZnElImpl<N, IS_FIELD> {

    fn one() -> ZnElImpl<N, IS_FIELD> {
        Self::project(1)
    }
}

impl<const N: u64> RingEl for ZnElImpl<N, true> {
    type Axioms = RingAxiomsField;
    type RingType = StaticRing<Self>;
    const RING: Self::RingType = StaticRing::<Self>::RING;
}

impl<const N: u64> FieldEl for ZnElImpl<N, true> {}

impl<const N: u64> RingEl for ZnElImpl<N, false> {
    type Axioms = RingAxiomsGeneral;
    type RingType = StaticRing<Self>;
    const RING: Self::RingType = StaticRing::<Self>::RING;
}

#[derive(Clone, Debug)]
pub struct StaticZnIterFn {
    current: i64
}

impl<const N: u64> FiniteRingIterFn<StaticRing<ZnElImpl<N, true>>> for StaticZnIterFn {

    fn next(&mut self, _: &StaticRing<ZnElImpl<N, true>>) -> Option<ZnElImpl<N, true>> {
        self.current += 1;
        if self.current as u64 <= N {
            Some(ZnElImpl::project(self.current as i64))
        } else {
            None
        }
    }
}

impl<const N: u64> FiniteRingIterFn<StaticRing<ZnElImpl<N, false>>> for StaticZnIterFn {

    fn next(&mut self, _: &StaticRing<ZnElImpl<N, false>>) -> Option<ZnElImpl<N, false>> {
        self.current += 1;
        if self.current as u64 <= N {
            Some(ZnElImpl::project(self.current as i64))
        } else {
            None
        }
    }
}

impl<const N: u64> FiniteRing for StaticRing<ZnElImpl<N, true>> {

    type VectorBasisType = VectorConstant<Self::El>;
    type IterFn = StaticZnIterFn;

    fn characteristic(&self) -> BigInt {
        BigInt::from(N as i64)
    }

    fn size(&self) -> BigInt {
        BigInt::from(N as i64)
    }
    
    fn iter_fn(&self) -> Self::IterFn {
        StaticZnIterFn {
            current: -1
        }
    }
}

impl<const N: u64> DivisibilityInfoRing for StaticRing<ZnElImpl<N, false>> {
    
    fn is_divisibility_computable(&self) -> bool {
        true
    }

    fn quotient(&self, lhs: &Self::El, rhs: &Self::El) -> Option<Self::El> {
        let (s, _, d) = signed_eea(rhs.repr as i64, N as i64, &i64::RING);
        if lhs.repr as i64 % d == 0 {
            Some(ZnElImpl::project(lhs.repr as i64 / d * s))
        } else {
            None
        }
    }

    fn is_unit(&self, el: &Self::El) -> bool {
        signed_gcd(N as i64, el.repr as i64, &i64::RING) == 1
    }
}

impl<const N: u64> FiniteRing for StaticRing<ZnElImpl<N, false>> {

    type VectorBasisType = VectorConstant<Self::El>;
    type IterFn = StaticZnIterFn;

    fn characteristic(&self) -> BigInt {
        BigInt::from(N as i64)
    }

    fn size(&self) -> BigInt {
        BigInt::from(N as i64)
    }
    
    fn iter_fn(&self) -> Self::IterFn {
        StaticZnIterFn {
            current: -1
        }
    }
}

pub type ZnEl<const N: u64> = ZnElImpl<N, {is_prime(N)}>;

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
