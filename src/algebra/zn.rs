use super::super::alg::*;
use super::bigint::*;
use super::eea::*;

use std::ops::{AddAssign, MulAssign, SubAssign, DivAssign, Add, Mul, Sub, Div, Neg};

#[derive(Debug, Clone)]
pub struct FactorRingZ {
    modulus: BigInt,
    inverse_modulus: BigInt,
    inverse_modulus_bitshift: usize
}

impl FactorRingZ {

    pub fn new(modulus: BigInt) -> Self {
        assert!(modulus >= 2);
        // have k such that 2^k > modulus^2
        // then (2^k / modulus) * x >> k differs at most 1 from floor(x / modulus)
        // if x < n^2, which is the case after multiplication
        let k = modulus.log2_floor() * 2 + 2;
        let inverse_modulus = BigInt::power_of_two(k) / modulus.clone();
        return FactorRingZ {
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

    pub fn project(&self, n: BigInt) -> <Self as Ring>::El {
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
        return QuotientRingZEl(red_n);
    }

    ///
    /// Returns either the inverse of x (as Ok()) or a nontrivial 
    /// factor of the modulus (as Err())
    /// 
    pub fn invert(&self, x: BigInt) -> Result<BigInt, BigInt> {
        let (s, _t, d) = eea(&BigInt::RING, x, self.modulus.clone());
        if d != 1 && d != -1i64 {
            Err(d)
        } else {
            Ok(s)
        }
    }
}

#[derive(Debug, Clone)]
pub struct QuotientRingZEl(BigInt);

impl Ring for FactorRingZ {

    type El = QuotientRingZEl;

    fn add_ref(&self, QuotientRingZEl(lhs): Self::El, QuotientRingZEl(rhs): &Self::El) -> Self::El {
        assert!(lhs < self.modulus);
        assert!(rhs < &self.modulus);

        let mut result = lhs + rhs;
        if result >= self.modulus {
            result -= &self.modulus;
        }

        assert!(result < self.modulus);
        return QuotientRingZEl(result);
    }

    fn mul_ref(&self, QuotientRingZEl(lhs): &Self::El, QuotientRingZEl(rhs): &Self::El) -> Self::El {
        assert!(*lhs < self.modulus);
        assert!(*rhs < self.modulus);

        let result = self.project_leq_n_square(
            BigInt::RING.mul_ref(lhs, rhs)
        );

        assert!(result < self.modulus);
        return QuotientRingZEl(result);
    }

    fn neg(&self, QuotientRingZEl(val): Self::El) -> Self::El {
        assert!(val < self.modulus);

        let mut result = -val;
        if result < 0 {
            result += &self.modulus;
        }

        assert!(result < self.modulus);
        return QuotientRingZEl(result);
    }

    fn zero(&self) -> Self::El {
        QuotientRingZEl(BigInt::zero())
    }

    fn one(&self) -> Self::El {
        QuotientRingZEl(BigInt::one())
    }

    fn eq(&self, QuotientRingZEl(lhs): &Self::El, QuotientRingZEl(rhs): &Self::El) -> bool {
        assert!(lhs < &self.modulus);
        assert!(rhs < &self.modulus);
        lhs == rhs
    }

    fn is_integral(&self) -> bool {
        unimplemented!()
    }

    fn is_euclidean(&self) -> bool {
        false
    }

    fn is_field(&self) -> bool {
        self.is_integral()
    }

    fn euclidean_div_rem(&self, _lhs: Self::El, _rhs: &Self::El) -> (Self::El, Self::El) { 
        panic!("Not a euclidean domain!");
    }

    fn div(&self, QuotientRingZEl(lhs): Self::El, QuotientRingZEl(rhs): &Self::El) -> Self::El { 
        match self.invert(rhs.clone()) {
            Err(factor) => panic!("Tried to divide in Z/{}Z, however this is not a field and the divisor is not invertible (the modulus has the nontrivial factor {})", self.modulus, factor),
            Ok(inverse) => self.mul(QuotientRingZEl(lhs), QuotientRingZEl(inverse))
        }
    }

    fn format(&self, QuotientRingZEl(el): &Self::El, f: &mut std::fmt::Formatter, _in_prod: bool) -> std::fmt::Result {
        write!(f, "{} mod {}", el, self.modulus)
    }
}

#[test]
fn test_mul() {
    let z257 = FactorRingZ::new(BigInt::from_str_radix("257", 10).unwrap());
    let x = BigInt::from_str_radix("256", 10).unwrap();
    assert!(z257.eq(&z257.one(), &z257.mul(z257.project(x.clone()), z257.project(x))));
}



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
        let (s, _t, _d) = signed_eea(rhs.repr as i64, N as i64);
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
        write!(f, "[{}]", self.repr)
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

impl<const N: u64> RingEl for ZnElImpl<N, false> {
    type Axioms = RingAxiomsRing;
    type RingType = StaticRing<Self>;
    const RING: Self::RingType = StaticRing::<Self>::RING;
}

impl<const N: u64> FieldEl for ZnElImpl<N, true> {}

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
