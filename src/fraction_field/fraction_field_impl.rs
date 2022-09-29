use super::super::prelude::*;
use super::super::rational::*;
use super::*;

#[derive(Debug, Clone, Copy)]
pub struct FractionFieldImpl<R>
    where R: Ring
{
    base_ring: R
}

impl<R> FractionFieldImpl<R>
    where R: Ring
{
    pub fn new(base_ring: R) -> Self {
        assert!(base_ring.is_integral().can_use());
        FractionFieldImpl { base_ring }
    }

    ///
    /// Creates the fraction field of the given ring, but skips the check
    /// that the base ring is integral. Use with care!
    /// 
    pub fn new_no_check(base_ring: R) -> Self {
        FractionFieldImpl { base_ring }
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

impl<R> SingletonRing for FractionFieldImpl<R>
    where R: SingletonRing
{
    fn singleton() -> Self {
        FractionFieldImpl::new(R::singleton())
    }
}

impl<R> FractionFieldImpl<R>
    where R: DivisibilityInfoRing
{
    pub fn in_base_ring(&self, (num, den): &El<Self>) -> Option<R::El> {
        assert!(self.base_ring.is_divisibility_computable().can_use());
        self.base_ring.quotient(num, den)
    }
}

impl<R> RingBase for FractionFieldImpl<R>
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

    default fn format(&self, el: &Self::El, f: &mut std::fmt::Formatter, _in_prod: bool) -> std::fmt::Result {
        self.format_base(&el.0, &el.1, f)
    }
}

impl<R> RingBase for FractionFieldImpl<R>
    where R: DivisibilityInfoRing
{
    default fn format(&self, el: &Self::El, f: &mut std::fmt::Formatter, _in_prod: bool) -> std::fmt::Result {
        if self.base_ring().is_divisibility_computable().can_use() {
            let data = [&el.0, &el.1];
            let d = self.base_ring().reduce_divisor(data.iter().copied().cloned());
            let num = self.base_ring().quotient(&el.0, &d).unwrap();
            let den = self.base_ring().quotient(&el.1, &d).unwrap();
            self.format_base(&num, &den, f)
        } else {
            self.format_base(&el.0, &el.1, f)
        }
    }
}

impl<R: DivisibilityInfoRing> ReducableElementRing for FractionFieldImpl<R> {

    fn reduce_divisor<I: Iterator<Item = El<Self>>>(&self, elements: I) -> El<Self> {
        let values: Vec<El<Self>> = elements.collect();
        let total_den = self.base_ring().product(values.iter().map(|(_, d)| d.clone()));
        let num_red = self.base_ring().reduce_divisor(
            values.iter()
                .map(|(n, d)| self.base_ring().quotient(&self.base_ring().mul_ref(n, &total_den), d).unwrap())
            );
        return (num_red, total_den);
    }
}

impl<R: Ring> RingExtension for FractionFieldImpl<R> {
    
    type BaseRing = R;
    type Embedding = StandardEmbedding<R, FractionFieldImpl<R>>;

    fn is_extension(&self) -> RingPropValue {
        RingPropValue::True
    }

    fn base_ring(&self) -> &R {
        &self.base_ring
    }
    
    fn embedding(&self) -> Self::Embedding {
        embedding(self.base_ring().clone(), self.clone())
    }

    fn from(&self, el: El<R>) -> El<Self> {
        (el, self.base_ring().one())
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
impl<R> HashableElRing for FractionFieldImpl<R>
    where R: IntegerRing + HashableElRing + OrderedRing + EuclideanInfoRing
{
    fn hash<H: std::hash::Hasher>(&self, h: &mut H, el: &Self::El) {
        let data = [&el.0, &el.1];
        let d = self.base_ring().reduce_divisor(data.iter().copied().cloned());
        let mut num = self.base_ring().quotient(&el.0, &d).unwrap();
        let mut den = self.base_ring().quotient(&el.1, &d).unwrap();
        if self.base_ring().cmp(&den, &self.base_ring().zero()) == std::cmp::Ordering::Less {
            num = self.base_ring().neg(num);
            den = self.base_ring().neg(den);
        }
        self.base_ring().hash(h, &num);
        self.base_ring().hash(h, &den);
    }
}

impl<R> CanonicalEmbeddingInfo<R> for FractionFieldImpl<R> 
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

impl<R> CanonicalEmbeddingInfo<R> for FractionFieldImpl<&R> 
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

impl<R> CanonicalEmbeddingInfo<&R> for FractionFieldImpl<R> 
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

impl<R> CanonicalEmbeddingInfo<FractionFieldImpl<R>> for FractionFieldImpl<R> 
    where R: Ring
{
    fn has_embedding(&self, from: &FractionFieldImpl<R>) -> RingPropValue {
        self.base_ring().has_embedding(from.base_ring())
    }

    fn embed(&self, from: &FractionFieldImpl<R>, el: Self::El) -> Self::El {
        assert!(self.has_embedding(from).can_use());
        (self.base_ring().embed(from.base_ring(), el.0), self.base_ring().embed(from.base_ring(), el.1))
    }
}

impl<R> CanonicalIsomorphismInfo<FractionFieldImpl<R>> for FractionFieldImpl<R> 
    where R: Ring
{
    fn has_isomorphism(&self, from: &FractionFieldImpl<R>) -> RingPropValue {
        self.base_ring().has_isomorphism(from.base_ring())
    }

    fn preimage(&self, from: &FractionFieldImpl<R>, el: Self::El) -> Self::El {
        assert!(self.has_isomorphism(from).can_use());
        (self.base_ring().preimage(from.base_ring(), el.0), self.base_ring().preimage(from.base_ring(), el.1))
    }
}

impl<R: IntegerRing> CanonicalEmbeddingInfo<StaticRing<r64>> for FractionFieldImpl<R> {

    fn has_embedding(&self, _from: &StaticRing<r64>) -> RingPropValue {
        RingPropValue::True
    }

    fn embed(&self, _from: &StaticRing<r64>, el: r64) -> Self::El {
        self.div(self.from(self.base_ring().from_z(el.num())), &self.from(self.base_ring().from_z(el.den())))
    }
}

impl<R: IntegerRing> CanonicalIsomorphismInfo<StaticRing<r64>> for FractionFieldImpl<R> {

    fn has_isomorphism(&self, _from: &StaticRing<r64>) -> RingPropValue {
        RingPropValue::True
    }

    fn preimage(&self, _from: &StaticRing<r64>, el: Self::El) -> r64 {
        let data = [&el.0, &el.1];
        let d = self.base_ring().reduce_divisor(data.iter().copied().cloned());
        let num = self.base_ring().quotient(&el.0, &d).unwrap();
        let den = self.base_ring().quotient(&el.1, &d).unwrap();
        r64::new(self.base_ring().preimage(&StaticRing::<i64>::RING, num), self.base_ring().preimage(&StaticRing::<i64>::RING, den))
    }
}

impl<R: IntegerRing> FractionField for FractionFieldImpl<R> {

    fn num<'a>(&self, el: &'a El<Self>) -> &'a El<Self::BaseRing> {
        &el.0
    }

    fn den<'a>(&self, el: &'a El<Self>) -> &'a El<Self::BaseRing> {
        &el.1
    }
}

impl<R: IntegerRing> RationalField for FractionFieldImpl<R> {}

impl<R> DivisibilityInfoRing for FractionFieldImpl<R>
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

impl<R> PartialEq for FractionFieldImpl<R>
    where R: Ring + PartialEq
{
    fn eq(&self, rhs: &Self) -> bool {
        self.base_ring() == rhs.base_ring()
    }
}

#[test]
fn test_add() {
    let rats = WrappingRing::new(FractionFieldImpl::new(BigInt::RING));
    let two = rats.from_z(2);
    let three = rats.from_z(3);
    let two_thirds = two.clone() / three.clone();
    assert_eq!(two, rats.from_z(2));
    let one_half = rats.one() / two;
    assert_eq!(rats.from_z(7) / rats.from_z(6), two_thirds + one_half);
}

#[test]
fn test_mul() {
    let rats = WrappingRing::new(FractionFieldImpl::new(BigInt::RING));
    let two = rats.from_z(2);
    let three = rats.from_z(3);
    let two_thirds = two.clone() / three.clone();
    let one_half = rats.one() / two;
    assert_eq!(rats.from_z(1) / rats.from_z(3), two_thirds * one_half);
}

#[test]
fn test_size_zero() {
    assert_eq!(0, std::mem::size_of::<FractionFieldImpl<BigIntRing>>())
}