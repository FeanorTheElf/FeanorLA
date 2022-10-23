use std::fmt::{Debug};
use libc;
use super::prelude::*;

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
                }
            },
            Err(()) => {
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
                )
            }
        }
    }
}

fn mpz_to_stdint(src: &MPZ) -> StdInt {
    unsafe {
        let size = mpz_bindings::__gmpz_sizeinbase(&src.integer as mpz_bindings::mpz_srcptr, 2);
        if size < i64::BITS as usize - 1 {
            return StdInt::from(mpz_bindings::__gmpz_get_si(&src.integer as mpz_bindings::mpz_srcptr));
        } else if size < i128::BITS as usize - 1 {
            let mut data = [0u64, 0u64];
            let mut size = 0;
            mpz_bindings::__gmpz_export(
                (&mut data as *mut mpz_bindings::mpir_ui) as *mut libc::c_void, 
                &mut size, 
                -1, 
                (u64::BITS / 8) as libc::size_t, 
                0, 
                0,
                &src.integer as mpz_bindings::mpz_srcptr
            );
            let mut abs_result = data[0] as i128 | ((data[1] as i128) << u64::BITS as i128);
            let sign = mpz_bindings::__gmpz_sgn(&src.integer as mpz_bindings::mpz_srcptr);
            return StdInt::RING.from_z_gen(
                abs_result * sign as i128,
                &i128::RING
            );
        } else {
            let len = (size + u64::BITS as usize - 1) / u64::BITS as usize;
            let mut data = Vec::new();

            unimplemented!()
        }
    }
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
        unimplemented!()
    }
    
    fn is_one(&self, val: &Self::El) -> bool { 
        unimplemented!()
    }

    fn is_neg_one(&self, val: &Self::El) -> bool {
        unimplemented!()
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
        unimplemented!()
    }

    fn from_z(&self, x: i64) -> Self::El {
        unsafe {
            let mut result = MPZ::new();
            mpz_bindings::__gmpz_set_si(&mut result.integer as mpz_bindings::mpz_ptr, x);
            return result;
        }
    }

    fn format(&self, el: &Self::El, f: &mut std::fmt::Formatter, _in_prod: bool) -> std::fmt::Result {
        unimplemented!()
    }

}

