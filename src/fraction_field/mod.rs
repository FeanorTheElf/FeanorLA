use super::prelude::*;
use super::rational::*;
use super::wrapper::*;
use super::eea::*;

#[derive(Debug, Clone, Copy)]
pub struct FieldOfFractions<R>
    where R: Ring
{
    base_ring: R
}

impl<R> FieldOfFractions<R>
    where R: Ring
{
    pub fn new(base_ring: R) -> Self {
        assert!(base_ring.is_integral().can_use());
        FieldOfFractions { base_ring }
    }

    ///
    /// Creates the fraction field of the given ring, but skips the check
    /// that the base ring is integral. Use with care!
    /// 
    pub fn new_no_check(base_ring: R) -> Self {
        FieldOfFractions { base_ring }
    }

    pub fn from(&self, el: R::El) -> El<Self> {
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

    pub fn base_ring(&self) -> &R {
        &self.base_ring
    }
}

impl<R> RingElWrapper<FieldOfFractions<R>>
    where R: Ring
{
    pub fn num(self) -> RingElWrapper<R> {
        let ((num, den), ring) = self.destruct();
        return ring.base_ring.bind_by_value(num);
    }

    pub fn den(self) -> RingElWrapper<R> {
        let ((num, den), ring) = self.destruct();
        return ring.base_ring.bind_by_value(den);
    }
}

impl<R> FieldOfFractions<R>
    where R: EuclideanInfoRing
{
    pub fn reduce(&self, el: El<Self>) -> El<Self> {
        assert!(self.base_ring().is_euclidean().can_use());
        let d = gcd(&self.base_ring, el.0.clone(), el.1.clone());
        return (self.base_ring.euclidean_div(el.0, &d), self.base_ring.euclidean_div(el.1, &d));
    }
}

impl<R> SingletonRing for FieldOfFractions<R>
    where R: SingletonRing
{
    fn singleton() -> Self {
        FieldOfFractions::new(R::singleton())
    }
}

trait SoftReducable: Ring {

    fn soft_reduce(&self, el: El<Self>) -> El<Self>;
}

impl<R> SoftReducable for FieldOfFractions<R>
    where R: Ring
{
    default fn soft_reduce(&self, el: El<Self>) -> El<Self> {
        el
    }
}

impl<R> SoftReducable for FieldOfFractions<R>
    where R: EuclideanInfoRing
{
    fn soft_reduce(&self, el: El<Self>) -> El<Self> {
        if self.base_ring.is_euclidean().can_use() {
            self.reduce(el)
        } else {
            el
        }
    }
}

impl<R> FieldOfFractions<R>
    where R: DivisibilityInfoRing
{
    pub fn in_base_ring(&self, (num, den): &El<Self>) -> Option<R::El> {
        assert!(self.base_ring.is_divisibility_computable().can_use());
        self.base_ring.quotient(num, den)
    }
}

impl<R> RingBase for FieldOfFractions<R>
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

    fn is_eq(&self, lhs: &Self::El, rhs: &Self::El) -> bool {
        self.base_ring.is_eq(&self.base_ring.mul_ref(&lhs.0, &rhs.1), &self.base_ring.mul_ref(&lhs.1, &rhs.0))
    }

    fn is_zero(&self, val: &Self::El) -> bool {
        self.base_ring.is_zero(&val.0)
    }

    fn is_one(&self, val: &Self::El) -> bool {
        self.base_ring.is_eq(&val.0, &val.1)
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
    
    fn characteristic(&self) -> BigInt {
        self.base_ring().characteristic()
    }

    fn is_field(&self) -> RingPropValue {
        RingPropValue::True
    }

    fn div(&self, lhs: Self::El, rhs: &Self::El) -> Self::El {
        assert!(!self.is_zero(rhs));
        (self.base_ring.mul(lhs.0, rhs.1.clone()), self.base_ring.mul(lhs.1, rhs.0.clone()))
    }

    fn from_z(&self, x: i64) -> Self::El {
        (self.base_ring.from_z(x), self.base_ring.one())
    }

    fn format(&self, el: &Self::El, f: &mut std::fmt::Formatter, _in_prod: bool) -> std::fmt::Result {
        let to_format = self.soft_reduce(el.clone());
        self.format_base(&to_format.0, &to_format.1, f)
    }
}

///
/// Finding a normal form for the fraction is not so easy, since the (reduced) 
/// denominator is only unique up to multiplication by a unit. In other words,
/// we have to assign a canonical element to each principal ideal of the ring.
/// 
/// That is why for now, we only support it for integers (the only units are +/-1,
/// so it is easy).
/// 
impl<R> HashableElRing for FieldOfFractions<R>
    where R: IntegerRing + HashableElRing + OrderedRing + EuclideanInfoRing
{
    fn hash<H: std::hash::Hasher>(&self, h: &mut H, el: &Self::El) {
        let (mut num, mut den) = self.reduce(el.clone());
        if self.base_ring().cmp(&den, &self.base_ring().zero()) == std::cmp::Ordering::Less {
            num = self.base_ring().neg(num);
            den = self.base_ring().neg(den);
        }
        self.base_ring().hash(h, &num);
        self.base_ring().hash(h, &den);
    }
}

impl<R> CanonicalEmbeddingInfo<R> for FieldOfFractions<R> 
    where R: Ring
{
    fn has_embedding(&self, from: &R) -> RingPropValue {
        self.base_ring().has_embedding(from)
    }

    fn embed(&self, from: &R, el: R::El) -> Self::El {
        assert!(self.has_embedding(from).can_use());
        self.from(self.base_ring().embed(from, el))
    }
}

impl<R> CanonicalEmbeddingInfo<R> for FieldOfFractions<&R> 
    where R: Ring
{
    fn has_embedding(&self, from: &R) -> RingPropValue {
        self.base_ring().has_embedding(&from)
    }

    fn embed(&self, from: &R, el: R::El) -> Self::El {
        assert!(self.has_embedding(from).can_use());
        self.from(self.base_ring().embed(&from, el))
    }
}

impl<R> CanonicalEmbeddingInfo<&R> for FieldOfFractions<R> 
    where R: Ring
{
    fn has_embedding(&self, from: &&R) -> RingPropValue {
        self.base_ring().has_embedding(*from)
    }

    fn embed(&self, from: &&R, el: R::El) -> Self::El {
        assert!(self.has_embedding(from).can_use());
        self.from(self.base_ring().embed(*from, el))
    }
}

impl<R> CanonicalEmbeddingInfo<FieldOfFractions<R>> for FieldOfFractions<R> 
    where R: Ring
{
    fn has_embedding(&self, from: &FieldOfFractions<R>) -> RingPropValue {
        self.base_ring().has_embedding(from.base_ring())
    }

    fn embed(&self, from: &FieldOfFractions<R>, el: Self::El) -> Self::El {
        assert!(self.has_embedding(from).can_use());
        (self.base_ring().embed(from.base_ring(), el.0), self.base_ring().embed(from.base_ring(), el.1))
    }
}

impl<R> CanonicalIsomorphismInfo<FieldOfFractions<R>> for FieldOfFractions<R> 
    where R: Ring
{
    fn has_isomorphism(&self, from: &FieldOfFractions<R>) -> RingPropValue {
        self.base_ring().has_isomorphism(from.base_ring())
    }

    fn preimage(&self, from: &FieldOfFractions<R>, el: Self::El) -> Self::El {
        assert!(self.has_isomorphism(from).can_use());
        (self.base_ring().preimage(from.base_ring(), el.0), self.base_ring().preimage(from.base_ring(), el.1))
    }
}

impl<R: IntegerRing> CanonicalEmbeddingInfo<StaticRing<r64>> for FieldOfFractions<R> {

    fn has_embedding(&self, _from: &StaticRing<r64>) -> RingPropValue {
        RingPropValue::True
    }

    fn embed(&self, _from: &StaticRing<r64>, el: r64) -> Self::El {
        self.div(self.from(self.base_ring().from_z(el.num())), &self.from(self.base_ring().from_z(el.den())))
    }
}

impl<R: IntegerRing> CanonicalIsomorphismInfo<StaticRing<r64>> for FieldOfFractions<R> {

    fn has_isomorphism(&self, _from: &StaticRing<r64>) -> RingPropValue {
        RingPropValue::True
    }

    fn preimage(&self, _from: &StaticRing<r64>, el: Self::El) -> r64 {
        let (num, den) = self.soft_reduce(el);
        r64::new(self.base_ring().preimage(&StaticRing::<i64>::RING, num), self.base_ring().preimage(&StaticRing::<i64>::RING, den))
    }
}

impl<R: IntegerRing> RationalField for FieldOfFractions<R> {

    type UnderlyingIntegers = R;

    fn num(&self, el: &Self::El) -> El<Self::UnderlyingIntegers> {
        el.0.clone()
    }

    fn den(&self, el: &Self::El) -> El<Self::UnderlyingIntegers> {
        el.1.clone()
    }

    fn underlying_integers(&self) -> Self::UnderlyingIntegers {
        self.base_ring().clone()
    }
}

impl<R> DivisibilityInfoRing for FieldOfFractions<R>
    where R: Ring
{
    fn is_divisibility_computable(&self) -> RingPropValue {
        RingPropValue::True
    }

    fn quotient(&self, lhs: &Self::El, rhs: &Self::El) -> Option<Self::El> {
        if self.is_zero(rhs) && !self.is_zero(lhs) {
            return None;
        } else if self.is_zero(rhs) && self.is_zero(lhs) {
            return Some(self.one());
        } else {
            Some(self.div(lhs.clone(), rhs))
        }
    }
}

impl<R> PartialEq for FieldOfFractions<R>
    where R: Ring + PartialEq
{
    fn eq(&self, rhs: &Self) -> bool {
        self.base_ring() == rhs.base_ring()
    }
}

#[cfg(test)]
use super::wrapper::*;

#[test]
fn test_add() {
    let rats = FieldOfFractions::new(BigInt::RING);
    let two = rats.bind(rats.from(BigInt::from(2)));
    let three = rats.bind(rats.from(BigInt::from(3)));
    let two_thirds = two.clone() / three.clone();
    assert_eq!(two, rats.bind(rats.from_z(2)));
    let one_half = rats.bind(rats.one()) / two;
    assert_eq!(rats.bind(rats.from_z(7)) / rats.bind(rats.from_z(6)), two_thirds + one_half);
}

#[test]
fn test_mul() {
    let rats = FieldOfFractions::new(BigInt::RING);
    let two = rats.bind(rats.from(BigInt::from(2)));
    let three = rats.bind(rats.from(BigInt::from(3)));
    let two_thirds = two.clone() / three.clone();
    let one_half = rats.bind(rats.one()) / two;
    assert_eq!(rats.bind(rats.from_z(1)) / rats.bind(rats.from_z(3)), two_thirds * one_half);
}

#[test]
fn test_size_zero() {
    assert_eq!(0, std::mem::size_of::<FieldOfFractions<BigIntRing>>())
}