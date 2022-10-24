use std::fmt::{Debug};
use libc;
use super::prelude::*;
use super::la::vec::*;
use super::integer::bigint::*;
use super::integer::bigint_soo::*;
use super::wrapper::*;

mod mpz_bindings;

pub struct MPZ {
    integer: mpz_bindings::__mpz_struct
}

impl MPZ {

    fn new() -> Self {
        unsafe {
            let mut result = mpz_bindings::UNINIT_MPZ;
            mpz_bindings::__gmpz_init(&mut result as mpz_bindings::mpz_ptr);
            return MPZ { integer: result };
        }
    }

    pub fn assign(&mut self, rhs: &MPZ) {
        unsafe {
            mpz_bindings::__gmpz_set(&mut self.integer as mpz_bindings::mpz_ptr, &rhs.integer as mpz_bindings::mpz_srcptr);
        }
    }
}

impl Clone for MPZ {

    fn clone(&self) -> Self {
        let mut result = MPZ::new();
        result.assign(self);
        return result;
    }
}

impl Debug for MPZ {

    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        MPZ::RING.format(self, f, false)
    }
}

impl Drop for MPZ {

    fn drop(&mut self) {
        unsafe {
            mpz_bindings::__gmpz_clear(&mut self.integer as mpz_bindings::mpz_ptr)
        }
    }
}

///
/// Note that this consumes the memory from src, but still
/// allocates memory for dst. In this sense, this is suboptimal,
/// as one could also reuse the memory. However, this would require
/// going into the internals and break encapsulation.
/// 
fn stdint_to_mpz(dst: &mut MPZ, src: StdInt) {
    unsafe {
        // for performance reasons, we distinguish the three cases that the value
        //  - fits into i64
        //  - fits into i128
        //  - needs BigInt
        match src.val().to_i128() {
            Ok(x) => {
                if x.abs() <= i64::MAX as i128 {
                    mpz_bindings::__gmpz_set_si(&mut dst.integer as mpz_bindings::mpz_ptr, x as i64);
                } else {
                    let upper = (x.abs() as u128 >> u64::BITS) as u64;
                    let lower = (x.abs() as u128 & ((1 << u64::BITS) - 1)) as u64;
                    mpz_bindings::__gmpz_set_ui(&mut dst.integer as mpz_bindings::mpz_ptr, upper);
                    mpz_bindings::__gmpz_mul_2exp(&mut dst.integer as mpz_bindings::mpz_ptr, &dst.integer as mpz_bindings::mpz_srcptr, u64::BITS as u64);
                    mpz_bindings::__gmpz_add_ui(&mut dst.integer as mpz_bindings::mpz_ptr, &dst.integer as mpz_bindings::mpz_srcptr, lower);
                    if x < 0 {
                        mpz_bindings::__gmpz_neg(&mut dst.integer as mpz_bindings::mpz_ptr, &dst.integer as mpz_bindings::mpz_srcptr);
                    }
                }
            },
            Err(()) => {
                let is_neg = src < 0;
                let data: Vec<mpz_bindings::mpir_ui> = src.into_val().to_bigint().base_u64_repr().into_owned().raw_data();
                // should never happen (we would be in the first case), however, with c it pays off to be paranoid
                assert!(data.len() > 0);
                mpz_bindings::__gmpz_import(
                    &mut dst.integer as mpz_bindings::mpz_ptr, 
                    data.len(), 
                    -1i32,
                    (u64::BITS / 8) as libc::size_t,
                    0, 
                    0, 
                    (data.as_ptr() as *const mpz_bindings::mpir_ui) as *const libc::c_void
                );
                if is_neg {
                    mpz_bindings::__gmpz_neg(&mut dst.integer as mpz_bindings::mpz_ptr, &dst.integer as mpz_bindings::mpz_srcptr);
                }
            }
        }
    }
}

fn mpz_to_stdint(src: &MPZ) -> StdInt {
    unsafe {
        // for performance reasons, we distinguish the three cases that the value
        //  - fits into i64
        //  - fits into i128
        //  - needs BigInt
        let size = mpz_bindings::__gmpz_sizeinbase(&src.integer as mpz_bindings::mpz_srcptr, 2);
        if size < i64::BITS as usize - 1 {
            return StdInt::from(mpz_bindings::__gmpz_get_si(&src.integer as mpz_bindings::mpz_srcptr));
        } else if size < i128::BITS as usize - 1 {
            let mut data = [0u64, 0u64];
            let mut size = 0;
            mpz_bindings::__gmpz_export(
                (data.as_mut_ptr() as *mut mpz_bindings::mpir_ui) as *mut libc::c_void, 
                &mut size, 
                -1, 
                (u64::BITS / 8) as libc::size_t, 
                0, 
                0,
                &src.integer as mpz_bindings::mpz_srcptr
            );
            assert_eq!(2, size);
            let abs_result = data[0] as i128 | ((data[1] as i128) << u64::BITS as i128);
            let is_neg: bool = mpz_bindings::mpz_is_neg(&src.integer as mpz_bindings::mpz_srcptr);
            if is_neg {
                return StdInt::RING.from_z_gen(
                    -abs_result,
                    &i128::RING
                );
            } else {
                return StdInt::RING.from_z_gen(
                    abs_result,
                    &i128::RING
                );
            }
        } else {
            let len = (size + u64::BITS as usize - 1) / u64::BITS as usize;
            let mut data = Vec::new();
            data.resize(len, 0u64);
            let mut size = 0;

            mpz_bindings::__gmpz_export(
                (data.as_mut_ptr() as *mut mpz_bindings::mpir_ui) as *mut libc::c_void,
                &mut size,
                -1,
                (u64::BITS / 8) as libc::size_t,
                0,
                0,
                &src.integer as mpz_bindings::mpz_srcptr
            );
            assert_eq!(len, size);
            let mut result = BigInt::from_base_u64_repr(Vector::new(data));
            if mpz_bindings::mpz_is_neg(&src.integer as mpz_bindings::mpz_srcptr) {
                result = BigInt::RING.neg(result);
            }
            return StdInt::RING.from_z_gen(
                result,
                &BigInt::RING
            );
        }
    }
}

pub type MPZInt = RingElWrapper<MPZRing>;

impl RingEl for MPZInt {

    type Axioms = RingAxiomsEuclidean;
    type RingType = WrappingRing<MPZRing>;
    const RING: Self::RingType = WrappingRing::new(MPZ::RING);
    const WRAPPED_RING: WrappingRing<Self::RingType> = WrappingRing::new(Self::RING);

    fn characteristic() -> StdInt {
        StdInt::zero()
    }
}

impl MPZ {
    pub const RING: MPZRing = MPZRing;
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MPZRing;

impl RingBase for MPZRing {
    
    type El = MPZ;

    fn add_ref(&self, mut lhs: Self::El, rhs: &Self::El) -> Self::El {
        self.add_assign_ref(&mut lhs, rhs);
        return lhs;
    }

    fn mul_ref(&self, lhs: &Self::El, rhs: &Self::El) -> Self::El {
        unsafe {
            let mut result = MPZ::new();
            mpz_bindings::__gmpz_mul(&mut result.integer as mpz_bindings::mpz_ptr, &lhs.integer as mpz_bindings::mpz_srcptr, &rhs.integer as mpz_bindings::mpz_srcptr);
            return result;
        }
    }

    fn neg(&self, mut val: Self::El) -> Self::El {
        unsafe {
            mpz_bindings::__gmpz_neg(&mut val.integer as mpz_bindings::mpz_ptr, &val.integer as mpz_bindings::mpz_srcptr);
            return val;
        }
    }

    fn zero(&self) -> Self::El {
        MPZ::new()
    }

    fn one(&self) -> Self::El {
        self.from_z(1)
    }

    fn is_eq(&self, lhs: &Self::El, rhs: &Self::El) -> bool {
        self.cmp(lhs, rhs) == std::cmp::Ordering::Equal
    }

    fn unspecified_element(&self) -> Self::El {
        MPZ::new()
    }

    fn sub_ref_fst(&self, lhs: &Self::El, mut rhs: Self::El) -> Self::El {
        unsafe {
            mpz_bindings::__gmpz_sub(&mut rhs.integer as mpz_bindings::mpz_ptr, &lhs.integer as mpz_bindings::mpz_srcptr, &rhs.integer as mpz_bindings::mpz_srcptr);
            return rhs;
        }
    }

    fn sub_ref_snd(&self, mut lhs: Self::El, rhs: &Self::El) -> Self::El {
        unsafe {
            mpz_bindings::__gmpz_sub(&mut lhs.integer as mpz_bindings::mpz_ptr, &lhs.integer as mpz_bindings::mpz_srcptr, &rhs.integer as mpz_bindings::mpz_srcptr);
            return lhs;
        }
    }

    fn add_assign(&self, lhs: &mut Self::El, rhs: Self::El) { 
        self.add_assign_ref(lhs, &rhs);
    }

    fn add_assign_ref(&self, lhs: &mut Self::El, rhs: &Self::El) { 
        unsafe {
            mpz_bindings::__gmpz_add(&mut lhs.integer as mpz_bindings::mpz_ptr, &lhs.integer as mpz_bindings::mpz_srcptr, &rhs.integer as mpz_bindings::mpz_srcptr);
        }
    }

    fn add_assign_int(&self, lhs: &mut Self::El, rhs: i64) {
        unsafe {
            if rhs < 0 {
                mpz_bindings::__gmpz_add_ui(&mut lhs.integer as mpz_bindings::mpz_ptr, &lhs.integer as mpz_bindings::mpz_srcptr, rhs as u64);
            } else {
                mpz_bindings::__gmpz_sub_ui(&mut lhs.integer as mpz_bindings::mpz_ptr, &lhs.integer as mpz_bindings::mpz_srcptr, (-rhs) as u64);
            }  
        }
    }

    fn mul_assign(&self, lhs: &mut Self::El, rhs: &Self::El) { 
        unsafe {
            mpz_bindings::__gmpz_mul(&mut lhs.integer as mpz_bindings::mpz_ptr, &lhs.integer as mpz_bindings::mpz_srcptr, &rhs.integer as mpz_bindings::mpz_srcptr);
        }
    }

    fn mul_assign_int(&self, lhs: &mut Self::El, rhs: i64) {
        unsafe {
            mpz_bindings::__gmpz_mul_ui(&mut lhs.integer as mpz_bindings::mpz_ptr, &lhs.integer as mpz_bindings::mpz_srcptr, rhs.abs() as u64);
            if rhs < 0 {
                mpz_bindings::__gmpz_neg(&mut lhs.integer as mpz_bindings::mpz_ptr, &lhs.integer as mpz_bindings::mpz_srcptr)
            }
        }
    }

    fn sub(&self, lhs: Self::El, rhs: Self::El) -> Self::El {
        self.sub_ref_snd(lhs, &rhs)
    }

    fn is_zero(&self, val: &Self::El) -> bool { 
        unsafe {
            mpz_bindings::__gmpz_cmp_si(&val.integer as mpz_bindings::mpz_srcptr, 0) == 0
        }
    }
    
    fn is_one(&self, val: &Self::El) -> bool { 
        unsafe {
            mpz_bindings::__gmpz_cmp_si(&val.integer as mpz_bindings::mpz_srcptr, 1) == 0
        }
    }

    fn is_neg_one(&self, val: &Self::El) -> bool {
        unsafe {
            mpz_bindings::__gmpz_cmp_si(&val.integer as mpz_bindings::mpz_srcptr, -1) == 0
        }
    }

    fn characteristic(&self) -> StdInt {
        StdInt::RING.zero()
    }

    fn is_integral(&self) -> RingPropValue {
        RingPropValue::True
    }

    fn is_field(&self) -> RingPropValue {
        RingPropValue::False
    }

    fn is_noetherian(&self) -> bool {
        true
    }

    fn div(&self, _lhs: Self::El, _rhs: &Self::El) -> Self::El {
        panic!("Not a field")
    }

    fn from_z_big(&self, x: &StdInt) -> Self::El {
        let mut result = MPZ::new();
        stdint_to_mpz(&mut result, x.clone());
        return result;
    }

    fn from_z(&self, x: i64) -> Self::El {
        unsafe {
            let mut result = MPZ::new();
            mpz_bindings::__gmpz_set_si(&mut result.integer as mpz_bindings::mpz_ptr, x);
            return result;
        }
    }

    fn format(&self, el: &Self::El, f: &mut std::fmt::Formatter, in_prod: bool) -> std::fmt::Result {
        let x = mpz_to_stdint(el);
        StdInt::RING.format(&x, f, in_prod)
    }
}

impl CanonicalEmbeddingInfo<MPZRing> for MPZRing {

    fn has_embedding(&self, _from: &MPZRing) -> RingPropValue {
        RingPropValue::True
    }

    fn embed(&self, _from: &MPZRing, el: MPZ) -> MPZ {
        el
    }
}

impl CanonicalIsomorphismInfo<MPZRing> for MPZRing {

    fn has_isomorphism(&self, _from: &MPZRing) -> RingPropValue {
        RingPropValue::True
    }

    fn preimage(&self, _from: &MPZRing, el: MPZ) -> MPZ {
        el
    }
}

impl SingletonRing for MPZRing {

    fn singleton() -> Self {
        MPZ::RING
    }
}

impl CanonicalEmbeddingInfo<StaticRing<i64>> for MPZRing {

    fn has_embedding(&self, _from: &StaticRing<i64>) -> RingPropValue {
        RingPropValue::True
    }

    fn embed(&self, _from: &StaticRing<i64>, el: i64) -> MPZ {
        MPZ::RING.from_z(el)
    }
}

impl CanonicalIsomorphismInfo<StaticRing<i64>> for MPZRing {

    fn has_isomorphism(&self, _from: &StaticRing<i64>) -> RingPropValue {
        RingPropValue::True
    }

    fn preimage(&self, _from: &StaticRing<i64>, el: MPZ) -> i64 {
        mpz_to_stdint(&el).to_i64()
    }
}

impl CanonicalEmbeddingInfo<BigIntSOORing> for MPZRing {

    fn has_embedding(&self, _from: &BigIntSOORing) -> RingPropValue {
        RingPropValue::True
    }

    fn embed(&self, _from: &BigIntSOORing, el: BigIntSOO) -> MPZ {
        let mut result = MPZ::new();
        stdint_to_mpz(&mut result, StdInt::RING.wrap(el));
        return result;
    }
}

impl CanonicalIsomorphismInfo<BigIntSOORing> for MPZRing {

    fn has_isomorphism(&self, _from: &BigIntSOORing) -> RingPropValue {
        RingPropValue::True
    }

    fn preimage(&self, _from: &BigIntSOORing, el: MPZ) -> BigIntSOO {
        mpz_to_stdint(&el).into_val()
    }
}


impl DivisibilityInfoRing for MPZRing {

    fn is_divisibility_computable(&self) -> RingPropValue {
        RingPropValue::True
    }

    fn is_divisible_by(&self, lhs: &Self::El, rhs: &Self::El) -> bool {
        self.is_zero(&self.euclidean_rem(lhs.clone(), rhs))
    }

    fn quotient(&self, lhs: &Self::El, rhs: &Self::El) -> Option<Self::El> {
        let (q, r) = self.euclidean_div_rem(lhs.clone(), rhs);
        if self.is_zero(&r) {
            return Some(q);
        } else {
            return None;
        }
    }

    fn is_unit(&self, el: &Self::El) -> bool {
        self.is_one(el) || self.is_neg_one(el)
    }
}

impl OrderedRing for MPZRing {

    fn cmp(&self, lhs: &Self::El, rhs: &Self::El) -> std::cmp::Ordering {
        unsafe {
            let res = mpz_bindings::__gmpz_cmp(
                &lhs.integer as mpz_bindings::mpz_srcptr,
                &rhs.integer as mpz_bindings::mpz_srcptr
            );
            if res < 0 {
                return std::cmp::Ordering::Less;
            } else if res > 0 {
                return std::cmp::Ordering::Greater;
            } else {
                return std::cmp::Ordering::Equal;
            }
        }
    }
}

impl EuclideanInfoRing for MPZRing {
    
    fn is_euclidean(&self) -> RingPropValue {
        RingPropValue::True
    }

    fn euclidean_div_rem(&self, mut lhs: Self::El, rhs: &Self::El) -> (Self::El, Self::El) {
        unsafe {
            let mut quo = MPZ::new();
            mpz_bindings::__gmpz_tdiv_qr(
                &mut quo.integer as mpz_bindings::mpz_ptr, 
                &mut lhs.integer as mpz_bindings::mpz_ptr, 
                &lhs.integer as mpz_bindings::mpz_srcptr, 
                &rhs.integer as mpz_bindings::mpz_srcptr
            );
            return (quo, lhs);
        }
    }

    fn euclidean_rem(&self, lhs: Self::El, rhs: &Self::El) -> Self::El { 
        unsafe {
            let mut rem = MPZ::new();
            mpz_bindings::__gmpz_tdiv_r(
                &mut rem.integer as mpz_bindings::mpz_ptr, 
                &lhs.integer as mpz_bindings::mpz_srcptr, 
                &rhs.integer as mpz_bindings::mpz_srcptr
            );
            return rem;
        }
    }

    fn euclidean_div(&self, lhs: Self::El, rhs: &Self::El) -> Self::El {
        unsafe {
            let mut rem = MPZ::new();
            mpz_bindings::__gmpz_tdiv_q(
                &mut rem.integer as mpz_bindings::mpz_ptr, 
                &lhs.integer as mpz_bindings::mpz_srcptr, 
                &rhs.integer as mpz_bindings::mpz_srcptr
            );
            return rem;
        }
    }

    fn euclidean_deg(&self, el: Self::El) -> StdInt {
        StdInt::RING.wrap(self.preimage(&BigIntSOO::RING, self.abs(el)))
    }
}

impl IntegerRing for MPZRing {
    
    fn to_float_approx(&self, el: &Self::El) -> f64 {
        unsafe {
            mpz_bindings::__gmpz_get_d(&el.integer as mpz_bindings::mpz_srcptr)
        }
    }

    fn from_float_approx(&self, el: f64) -> Option<Self::El> {
        unsafe {
            let mut result = MPZ::new();
            mpz_bindings::__gmpz_set_d(&mut result.integer as mpz_bindings::mpz_ptr, el);
            return Some(result);
        }
    }

    fn mul_pow_2(&self, mut el: El<Self>, power: u64) -> El<Self> {
        unsafe {
            mpz_bindings::__gmpz_mul_2exp(&mut el.integer as mpz_bindings::mpz_ptr, &el.integer as mpz_bindings::mpz_srcptr, power);
            return el;
        }
    }

    fn euclidean_div_pow_2(&self, mut el: El<Self>, power: u64) -> El<Self> {
        unsafe {
            mpz_bindings::__gmpz_tdiv_q_2exp(&mut el.integer as mpz_bindings::mpz_ptr, &el.integer as mpz_bindings::mpz_srcptr, power);
            return el;
        }
    }

    fn is_odd(&self, el: &Self::El) -> bool {
        unsafe {
            mpz_bindings::__gmpz_get_ui(&el.integer as mpz_bindings::mpz_srcptr) % 2 == 1
        }
    }

    fn abs_log2_floor(&self, el: &El<Self>) -> u64 {
        unsafe {
            mpz_bindings::__gmpz_sizeinbase(&el.integer as mpz_bindings::mpz_srcptr, 2) as u64
        }
    }

    fn abs_is_bit_set(&self, el: &El<Self>, bit: u64) -> bool {
        unsafe {
            if mpz_bindings::mpz_is_neg(&el.integer as mpz_bindings::mpz_srcptr) {
                let value = mpz_bindings::__gmpz_tstbit(&el.integer as mpz_bindings::mpz_srcptr, bit) == 1;
                let least_significant_zero = mpz_bindings::__gmpz_scan1(&el.integer as mpz_bindings::mpz_srcptr, 0);
                if bit <= least_significant_zero {
                    value
                } else {
                    !value
                }
            } else {
                mpz_bindings::__gmpz_tstbit(&el.integer as mpz_bindings::mpz_srcptr, bit) == 1
            }
        }
    }

    fn floor_div(&self, mut a: Self::El, b: &Self::El) -> Self::El {
        unsafe {
            mpz_bindings::__gmpz_fdiv_q(&mut a.integer as mpz_bindings::mpz_ptr, &a.integer as mpz_bindings::mpz_srcptr, &b.integer as mpz_bindings::mpz_srcptr);
            return a;
        }
    }

    fn abs_cmp(&self, lhs: &Self::El, rhs: &Self::El) -> std::cmp::Ordering {
        unsafe {
            let res = mpz_bindings::__gmpz_cmpabs(&lhs.integer as mpz_bindings::mpz_srcptr, &rhs.integer as mpz_bindings::mpz_srcptr);
            if res < 0 {
                return std::cmp::Ordering::Less;
            } else if res > 0 {
                return std::cmp::Ordering::Greater;
            } else {
                return std::cmp::Ordering::Equal;
            }
        }
    }

    fn root_floor(&self, el: &Self::El, n: u64) -> Self::El {
        assert!(n > 0);
        unsafe {
            assert!(!mpz_bindings::mpz_is_neg(&el.integer as mpz_bindings::mpz_srcptr));
            let mut result = MPZ::new();
            mpz_bindings::__gmpz_nthroot(&mut result.integer as mpz_bindings::mpz_ptr, &el.integer as mpz_bindings::mpz_srcptr, n);
            return result;
        }
    }

    fn highest_dividing_power_of_two(&self, el: &El<Self>) -> usize {
        unsafe {
            mpz_bindings::__gmpz_scan1(&el.integer as mpz_bindings::mpz_srcptr, 0) as usize
        }
    }
}

#[test]
fn test_mpz_to_stdint() {
    // testing for i64 value
    let a = StdInt::from(3684);
    let mut b = MPZ::new();
    stdint_to_mpz(&mut b, a.clone());
    assert_eq!(a, mpz_to_stdint(&b));

    // testing for i128 value
    let a = StdInt::from(3684).mul_pow_2(64);
    let mut b = MPZ::new();
    stdint_to_mpz(&mut b, a.clone());
    assert_eq!(a, mpz_to_stdint(&b));
    
    // testing for negative i128 value
    let a = StdInt::from(-3684).mul_pow_2(64);
    let mut b = MPZ::new();
    stdint_to_mpz(&mut b, a.clone());
    assert_eq!(a, mpz_to_stdint(&b));

    // testing for larger values
    let a = StdInt::from(3684).mul_pow_2(128);
    let mut b = MPZ::new();
    stdint_to_mpz(&mut b, a.clone());
    assert_eq!(a, mpz_to_stdint(&b));

    // testing for the special case i128::MIN
    let a = StdInt::RING.from_z_gen(i128::MIN, &i128::RING);
    let mut b = MPZ::new();
    stdint_to_mpz(&mut b, a.clone());
    assert_eq!(a, mpz_to_stdint(&b));

    // testing for compatibility
    let a1 = StdInt::from(384).mul_pow_2(64);
    let a2 = StdInt::from(-1).mul_pow_2(128);
    let mut b1 = MPZ::new();
    let mut b2 = MPZ::new();
    stdint_to_mpz(&mut b1, a1.clone());
    stdint_to_mpz(&mut b2, a2.clone());
    let b = MPZ::RING.add(b1, b2);
    let a = a1 + a2;
    assert_eq!(a, mpz_to_stdint(&b));
}

#[test]
fn test_abs_is_bit_set() {
    let a = MPZ::RING.from_z(1 << 15);
    assert_eq!(true, MPZ::RING.abs_is_bit_set(&a, 15));
    assert_eq!(false, MPZ::RING.abs_is_bit_set(&a, 16));
    assert_eq!(false, MPZ::RING.abs_is_bit_set(&a, 14));

    let a = MPZ::RING.from_z(-7);
    assert_eq!(true, MPZ::RING.abs_is_bit_set(&a, 0));
    assert_eq!(true, MPZ::RING.abs_is_bit_set(&a, 1));
    assert_eq!(true, MPZ::RING.abs_is_bit_set(&a, 2));
    assert_eq!(false, MPZ::RING.abs_is_bit_set(&a, 3));

    let a = MPZ::RING.from_z(-1 << 15);
    assert_eq!(true, MPZ::RING.abs_is_bit_set(&a, 15));
    assert_eq!(false, MPZ::RING.abs_is_bit_set(&a, 16));
    assert_eq!(false, MPZ::RING.abs_is_bit_set(&a, 14));
}

#[test]
fn test_highest_dividing_power_of_two() {
    let a = MPZ::RING.from_z(1);
    assert_eq!(0, MPZ::RING.highest_dividing_power_of_two(&a));

    let a = MPZ::RING.from_z(83489 << 15);
    assert_eq!(15, MPZ::RING.highest_dividing_power_of_two(&a));

    let a = MPZ::RING.from_z(-83489 << 15);
    assert_eq!(15, MPZ::RING.highest_dividing_power_of_two(&a));

    let a = MPZ::RING.from_z(-1);
    assert_eq!(0, MPZ::RING.highest_dividing_power_of_two(&a));
}

#[bench]
fn bench_mul(bencher: &mut test::Bencher) {
    let x = MPZ::RING.from_z_gen(BigInt::from_str_radix("2382385687561872365981723456981723456987134659834659813491964132897159283746918732563498628754", 10).unwrap(), &BigInt::RING);
    let y = MPZ::RING.from_z_gen(BigInt::from_str_radix("48937502893645789234569182735646324895723409587234", 10).unwrap(), &BigInt::RING);
    let z = MPZ::RING.from_z_gen(BigInt::from_str_radix("116588006478839442056346504147013274749794691549803163727888681858469844569693215953808606899770104590589390919543097259495176008551856143726436", 10).unwrap(), &BigInt::RING);
    bencher.iter(|| {
        let p = MPZ::RING.mul_ref(&x, &y);
        assert!(MPZ::RING.is_eq(&z, &p));
    })
}

#[bench]
fn bench_div(bencher: &mut test::Bencher) {
    let x = MPZ::RING.from_z_gen(BigInt::from_str_radix("2382385687561872365981723456981723456987134659834659813491964132897159283746918732563498628754", 10).unwrap(), &BigInt::RING);
    let y = MPZ::RING.from_z_gen(BigInt::from_str_radix("48937502893645789234569182735646324895723409587234", 10).unwrap(), &BigInt::RING);
    let z = MPZ::RING.from_z_gen(BigInt::from_str_radix("48682207850683149082203680872586784064678018", 10).unwrap(), &BigInt::RING);
    bencher.iter(|| {
        let q = MPZ::RING.euclidean_div(x.clone(), &y);
        assert!(MPZ::RING.is_eq(&z, &q));
    })
}
