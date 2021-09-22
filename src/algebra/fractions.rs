use super::super::alg::*;
use super::eea::*;

#[derive(Debug, Clone)]
pub struct FieldOfFractions<R>
    where R: Ring
{
    base_ring: R
}

impl<R> FieldOfFractions<R>
    where R: Ring
{
    pub fn new(base_ring: R) -> Self {
        assert!(base_ring.is_integral());
        Self { base_ring }
    }

    pub fn from(&self, el: R::El) -> <Self as Ring>::El {
        (el, self.base_ring.one())
    }

    fn format_base(&self, num: &R::El, den: &R::El, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        if self.base_ring.is_zero(&num) || self.base_ring.is_one(&den) {
            self.base_ring.format(&num, f, true)?;
        } else if self.base_ring.is_neg_one(&den) {
            write!(f, "-")?;
            self.base_ring.format(&num, f, true)?;
        } else {
            self.base_ring.format(num, f, true)?;
            write!(f, " / ")?;
            self.base_ring.format(den, f, true)?;
        }
        return Ok(());
    }

    pub fn soft_reduce(&self, el: <Self as Ring>::El) -> <Self as Ring>::El {
        if self.base_ring.is_euclidean() {
            let d = gcd(&self.base_ring, el.0.clone(), el.1.clone());
            return (self.base_ring.euclidean_div(el.0, &d), self.base_ring.euclidean_div(el.1, &d));
        } else {
            return el;
        }
    }
}

impl<R> Ring for FieldOfFractions<R>
    where R: Ring
{
    type El = (R::El, R::El);

    fn add_ref(&self, lhs: Self::El, rhs: &Self::El) -> Self::El {
        (self.base_ring.add(
            self.base_ring.mul_ref(&lhs.0, &rhs.1),
            self.base_ring.mul_ref(&lhs.1, &rhs.0)
        ), self.base_ring.mul_ref(&lhs.1, &rhs.1))
    }

    fn add(&self, lhs: Self::El, rhs: Self::El) -> Self::El {
        (self.base_ring.add(
            self.base_ring.mul_ref(&lhs.0, &rhs.1),
            self.base_ring.mul_ref(&lhs.1, &rhs.0)
        ), self.base_ring.mul(lhs.1, rhs.1))
    }

    fn mul_ref(&self, lhs: &Self::El, rhs: &Self::El) -> Self::El {
        (self.base_ring.mul_ref(&lhs.0, &rhs.0), self.base_ring.mul_ref(&lhs.1, &rhs.1))
    }

    fn neg(&self, val: Self::El) -> Self::El {
        (self.base_ring.neg(val.0), val.1)
    }

    fn zero(&self) -> Self::El {
        (self.base_ring.zero(), self.base_ring.one())
    }

    fn one(&self) -> Self::El {
        (self.base_ring.one(), self.base_ring.one())
    }

    fn eq(&self, lhs: &Self::El, rhs: &Self::El) -> bool {
        self.base_ring.eq(&self.base_ring.mul_ref(&lhs.0, &rhs.1), &self.base_ring.mul_ref(&lhs.1, &rhs.0))
    }

    fn is_zero(&self, val: &Self::El) -> bool {
        self.base_ring.is_zero(&val.0)
    }

    fn is_one(&self, val: &Self::El) -> bool {
        self.base_ring.eq(&val.0, &val.1)
    }

    fn unspecified_element(&self) -> Self::El {
        (self.base_ring.unspecified_element(), self.base_ring.unspecified_element())
    }

    fn is_integral(&self) -> bool {
        true
    }

    fn is_euclidean(&self) -> bool {
        false
    }

    fn is_field(&self) -> bool {
        true
    }
    
    fn euclidean_div_rem(&self, _: Self::El, _: &Self::El) -> (Self::El, Self::El) {
        panic!("Not a euclidean domain!")
    }

    fn div(&self, lhs: Self::El, rhs: &Self::El) -> Self::El {
        if self.is_zero(rhs) {
            panic!("division by zero!")
        } else {
            (self.base_ring.mul(lhs.0, rhs.1.clone()), self.base_ring.mul(lhs.1, rhs.0.clone()))
        }
    }

    fn from_z(&self, x: i64) -> Self::El {
        (self.base_ring.from_z(x), self.base_ring.one())
    }

    fn format(&self, el: &Self::El, f: &mut std::fmt::Formatter, _in_prod: bool) -> std::fmt::Result {
        let to_format = self.soft_reduce(el.clone());
        self.format_base(&to_format.0, &to_format.1, f)
    }
}

#[cfg(test)]
use super::bigint::BigInt;
#[cfg(test)]
use super::super::alg_env::*;
#[cfg(test)]
use super::super::alg_macros::*;

#[test]
fn test_add() {
    let rats = FieldOfFractions::new(BigInt::RING);
    let two = rats.bind::<RingAxiomsField>(rats.from(BigInt::from(2)));
    let three = rats.bind(rats.from(BigInt::from(3)));
    let two_thirds = two.clone() / three.clone();
    assert_eq!(two, rats.bind(rats.from_z(2)));
    let one_half = rats.bind(rats.one()) / two;
    assert_eq!(rats.bind::<RingAxiomsField>(rats.from_z(7)) / rats.bind(rats.from_z(6)), two_thirds + one_half);
}

#[test]
fn test_mul() {
    let rats = FieldOfFractions::new(BigInt::RING);
    let two = rats.bind::<RingAxiomsField>(rats.from(BigInt::from(2)));
    let three = rats.bind(rats.from(BigInt::from(3)));
    let two_thirds = two.clone() / three.clone();
    let one_half = rats.bind(rats.one()) / two;
    assert_eq!(rats.bind::<RingAxiomsField>(rats.from_z(1)) / rats.bind(rats.from_z(3)), two_thirds * one_half);
}