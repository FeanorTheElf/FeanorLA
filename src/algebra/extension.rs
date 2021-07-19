use super::super::alg::*;
use super::eea::signed_eea;
use super::poly::*;

use std::ops::{AddAssign, MulAssign, SubAssign, DivAssign, Add, Mul, Sub, Div, Neg};

pub const MAX_DEGREE: usize = 64;

///
/// Represents elements of the ring T[X]/(X^d - a_(d-1) X^(d-1) - ... - a_0)
/// where the values a_0, ..., a_(d-1) are contained in POLY.
/// 
#[derive(Copy, Clone, PartialEq, Eq, Debug)]
pub struct RingExtEl<const D: usize, const POLY: [i8; MAX_DEGREE], T: RingEl> {
    data: [T; D]
}

const fn copy<const N: usize, const M: usize, T>(src: [T; N], mut dst: [T; M], i: usize) -> [T; M] 
    where T: Copy
{
    if i == N {
        return dst;
    } else {
        dst[i] = src[i];
        return copy(src, dst, i + 1);
    }
}

pub const fn ext_poly_array<const D: usize>(data: [i8; D]) -> [i8; MAX_DEGREE] {
    assert!(D <= MAX_DEGREE);
    return copy::<D, MAX_DEGREE, i8>(data, [0; MAX_DEGREE], 0);
}

impl<const D: usize, const POLY: [i8; MAX_DEGREE], T: RingEl> RingExtEl<D, POLY, T> {

    pub const fn new(coefficients: [T; D]) -> Self {
        assert!(D > 0);
        RingExtEl {
            data: coefficients
        }
    }
}

pub struct RingExtRing<const D: usize, const POLY: [i8; MAX_DEGREE], T: RingEl> {
    phantom: std::marker::PhantomData<T>
}

impl<const D: usize, const POLY: [i8; MAX_DEGREE], T: RingEl> Ring for RingExtRing<D, POLY, T> 
    where T: Clone
{
    type El = RingExtEl<D, POLY, T>;

    fn add_ref(&self, mut lhs: Self::El, rhs: &Self::El) -> Self::El {
        for i in 0..D {
            lhs.data[i] += rhs.data[i];
        }
        return lhs;
    }

    fn mul_ref(&self, lhs: &Self::El, rhs: &Self::El) -> Self::El {
        let mut result: [[T; D]; 2] = [[T::zero(); D]; 2];
        for i in 0..D {
            for j in 0..D {
                result[(i + j) / D][(i + j) % D] += lhs.data[i] * rhs.data[j];
            }
        }
        for i in (0..D).rev() {
            for j in 0..D {
                result[(i + j) / D][(i + j) % D] += result[1][i] * T::from(POLY[j]);
            }
        }
        return RingExtEl {
            data: result[0]
        };
    }

    fn neg(&self, mut val: Self::El) -> Self::El {
        for i in 0..D {
            val.data[i] = -val.data[i];
        }
        return val;
    }

    fn zero(&self) -> Self::El {
        RingExtEl::new([T::zero(); D])
    }

    fn one(&self) -> Self::El {
        let mut result = self.zero();
        result.data[0] = T::one();
        return result;
    }

    fn eq(&self, lhs: &Self::El, rhs: &Self::El) -> bool {
        lhs.data == rhs.data
    }

    fn add(&self, mut lhs: Self::El, rhs: Self::El) -> Self::El {
        for (i, x) in std::array::IntoIter::new(rhs.data).enumerate() {
            lhs.data[i] += x;
        }
        return lhs;
    }

    fn is_integral(&self) -> bool {
        T::Axioms::is_integral()
    }

    fn is_euclidean(&self) -> bool {
        false
    }
    
    fn is_field(&self) -> bool {
        T::Axioms::is_field()
    }
    
    fn euclidean_div_rem(&self, lhs: Self::El, rhs: &Self::El) -> (Self::El, Self::El) {
        panic!("not a euclidean ring")
    }

    fn div(&self, lhs: Self::El, rhs: &Self::El) -> Self::El {
        assert!(T::Axioms::is_field());
        let poly_ring = PolyRing::adjoint(T::RING, "X");
    }

    fn from_z(&self, x: i64) -> Self::El {
        let mut result = self.zero();
        for _ in 0..x.abs() {
            result = self.add(self.one(), result);
        }
        if x < 0 {
            return self.neg(result);
        } else {
            return result;
        }
    }

    fn format(&self, el: &Self::El, f: &mut std::fmt::Formatter, _in_prod: bool) -> std::fmt::Result {
        write!(f, "{:?}", el)
    }

    fn format_in_brackets(&self, el: &Self::El, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "(")?;
        self.format(el, f, false)?;
        write!(f, ")")?;
        return Ok(());
    }
}

impl<const D: usize, const POLY: [u64; MAX_DEGREE], T: RingEl> AddAssign for RingExtEl<D, POLY, T> {

    fn add_assign(&mut self, rhs: Self) {
        for i in 0..D {
            self.data[i] += rhs.data[i];
        }
    }
}

impl<const D: usize, const POLY: [u64; MAX_DEGREE], T: RingEl> MulAssign for RingExtEl<D, POLY, T> {

    fn mul_assign(&mut self, rhs: Self) {
        
    }
}

#[cfg(test)]
use super::zn::*;

#[test]
fn test_add() {
    type F2 = RingExtEl<2, {ext_poly_array([1, 1])}, ZnEl::<2>>;
    let a: F2 = RingExtEl { data: [ ZnEl::ZERO, ZnEl::ONE ] };
}