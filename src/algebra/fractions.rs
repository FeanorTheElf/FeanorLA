use super::super::ring::*;
use super::eea::*;

pub trait FieldOfFractionsInfoRing: Ring {
    type FieldOfFractions: Ring;
    type CanonicalEmbedding: FnMut(Self::El) -> <Self::FieldOfFractions as Ring>::El;

    fn field_of_fractions(&self) -> (Self::FieldOfFractions, Self::CanonicalEmbedding);
}

impl<R> FieldOfFractionsInfoRing for R
    where R: Ring
{
    type FieldOfFractions = FieldOfFractions<R>;
    type CanonicalEmbedding = FieldOfFractionsEmbedding<R>;

    fn field_of_fractions(&self) -> (Self::FieldOfFractions, Self::CanonicalEmbedding) {
        assert!(self.is_integral().can_use());
        let field = FieldOfFractions { base_ring: self.clone() };
        (field.clone(), FieldOfFractionsEmbedding { field_of_fractions: field })
    }
}

pub struct FieldOfFractionsEmbedding<R>
    where R: Ring
{
    field_of_fractions: FieldOfFractions<R>
}

impl<R> FnOnce<(R::El, )> for FieldOfFractionsEmbedding<R> 
    where R: Ring
{
    type Output = <FieldOfFractions<R> as Ring>::El;

    extern "rust-call" fn call_once(
        mut self, 
        (x, ): (R::El, )
    ) -> Self::Output {
        self.call_mut((x, ))
    }
}

impl<R> FnMut<(R::El, )> for FieldOfFractionsEmbedding<R> 
    where R: Ring
{
    extern "rust-call" fn call_mut(
        &mut self, 
        (x, ): (R::El, )
    ) -> Self::Output {
        self.field_of_fractions.from(x)
    }
}

#[derive(Debug, Clone)]
pub struct FieldOfFractions<R>
    where R: Ring
{
    base_ring: R
}

impl<R> FieldOfFractions<R>
    where R: Ring
{
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
}

trait SoftReducable: Ring {

    fn soft_reduce(&self, el: <Self as Ring>::El) -> <Self as Ring>::El;
}

impl<R> SoftReducable for FieldOfFractions<R>
    where R: Ring
{
    default fn soft_reduce(&self, el: <Self as Ring>::El) -> <Self as Ring>::El {
        el
    }
}

impl<R> SoftReducable for FieldOfFractions<R>
    where R: EuclideanInfoRing
{
    fn soft_reduce(&self, el: <Self as Ring>::El) -> <Self as Ring>::El {
        if self.base_ring.is_euclidean().can_use() {
            let d = gcd(&self.base_ring, el.0.clone(), el.1.clone());
            return (self.base_ring.euclidean_div(el.0, &d), self.base_ring.euclidean_div(el.1, &d));
        } else {
            return el;
        }
    }
}

impl<R> FieldOfFractions<R>
    where R: DivisibilityInfoRing
{
    pub fn in_base_ring(&self, (num, den): &<Self as Ring>::El) -> Option<R::El> {
        assert!(self.base_ring.is_divisibility_computable());
        self.base_ring.quotient(num, den)
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

    fn is_integral(&self) -> RingPropValue {
        RingPropValue::True
    }

    fn is_noetherian(&self) -> bool {
        true
    }

    fn is_field(&self) -> RingPropValue {
        RingPropValue::True
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

impl<R> DivisibilityInfoRing for FieldOfFractions<R>
    where R: Ring
{
    fn is_divisibility_computable(&self) -> bool {
        true
    }

    fn quotient(&self, lhs: &Self::El, rhs: &Self::El) -> Option<Self::El> {
        if self.is_zero(rhs) {
            None
        } else {
            Some(self.div(lhs.clone(), rhs))
        }
    }
}

#[cfg(test)]
use super::super::bigint::BigInt;
#[cfg(test)]
use super::super::wrapper::*;

#[test]
fn test_add() {
    let (rats, _) = BigInt::RING.field_of_fractions();
    let two = rats.bind(rats.from(BigInt::from(2)));
    let three = rats.bind(rats.from(BigInt::from(3)));
    let two_thirds = two.clone() / three.clone();
    assert_eq!(two, rats.bind(rats.from_z(2)));
    let one_half = rats.bind(rats.one()) / two;
    assert_eq!(rats.bind(rats.from_z(7)) / rats.bind(rats.from_z(6)), two_thirds + one_half);
}

#[test]
fn test_mul() {
    let (rats, _) = BigInt::RING.field_of_fractions();
    let two = rats.bind(rats.from(BigInt::from(2)));
    let three = rats.bind(rats.from(BigInt::from(3)));
    let two_thirds = two.clone() / three.clone();
    let one_half = rats.bind(rats.one()) / two;
    assert_eq!(rats.bind(rats.from_z(1)) / rats.bind(rats.from_z(3)), two_thirds * one_half);
}