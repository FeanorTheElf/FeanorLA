use std::fmt::{Debug};
use libc;
use super::prelude::*;
use super::la::vec::*;
use super::integer::bigint::*;

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
        unimplemented!()
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

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MPZRing;

impl MPZRing {
    pub const RING: MPZRing = MPZRing;
}

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
        unsafe {
            mpz_bindings::__gmpz_cmp(&lhs.integer as mpz_bindings::mpz_srcptr, &rhs.integer as mpz_bindings::mpz_srcptr) == 0
        }
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

    fn div(&self, lhs: Self::El, rhs: &Self::El) -> Self::El {
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
    let b = MPZRing::RING.add(b1, b2);
    let a = a1 + a2;
    assert_eq!(a, mpz_to_stdint(&b));
}