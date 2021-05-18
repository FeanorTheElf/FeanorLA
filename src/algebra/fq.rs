use super::super::alg::*;
use super::eea::signed_eea;
use super::zn::*;

use std::ops::{AddAssign, MulAssign, SubAssign, DivAssign, Add, Mul, Sub, Div, Neg};

pub const MAX_DEGREE: usize = 64;

pub struct ZnExtEl<const N: u64, const D: usize, const AD: [u64; MAX_DEGREE]> {
    data: [u64; D]
}

#[test]
fn test_add() {
    const ONES: [u64; MAX_DEGREE] = [1; MAX_DEGREE];
    type F2 = ZnExtEl<2, 2, ONES>;
    let a: F2 = ZnExtEl { data: [0, 1] };
}