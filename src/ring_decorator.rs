use super::prelude::*;
use super::wrapper::*;

use vector_map::VecMap;

///
/// Trait for rings that are decorated versions of other rings. Implementing this
/// trait will automatically provide a default implementation of the ring trait for
/// the object, which delegates all calls to the decorated ring. 
/// It is thus discouraged to use this trait for rings that change the mathematical
/// behavior of the decorated ring. In this case, prefer `RingExtension` where applicable.
/// 
pub trait RingDecorator: Clone + std::fmt::Debug {
    type DecoratedRing: RingBase;
    fn decorated_ring(&self) -> &Self::DecoratedRing;
}

impl<R: RingDecorator> RingBase for R {
    type El = El<R::DecoratedRing>;

    fn add_ref(&self, lhs: Self::El, rhs: &Self::El) -> Self::El { self.decorated_ring().add_ref(lhs, rhs) }
    fn mul_ref(&self, lhs: &Self::El, rhs: &Self::El) -> Self::El { self.decorated_ring().mul_ref(lhs, rhs) }
    fn add_assign(&self, lhs: &mut Self::El, rhs: Self::El) { self.decorated_ring().add_assign(lhs, rhs) }
    fn add_assign_ref(&self, lhs: &mut Self::El, rhs: &Self::El) { self.decorated_ring().add_assign_ref(lhs, rhs) }
    fn mul_assign(&self, lhs: &mut Self::El, rhs: Self::El) { self.decorated_ring().mul_assign(lhs, rhs) }
    fn neg(&self, val: Self::El) -> Self::El { self.decorated_ring().neg(val) }
    fn zero(&self) -> Self::El { self.decorated_ring().zero() }
    fn one(&self) -> Self::El { self.decorated_ring().one() }
    fn is_eq(&self, lhs: &Self::El, rhs: &Self::El) -> bool { self.decorated_ring().is_eq(lhs, rhs) }
    fn unspecified_element(&self) -> Self::El { self.decorated_ring().unspecified_element() }
    fn sub_ref_fst(&self, lhs: &Self::El, rhs: Self::El) -> Self::El { self.decorated_ring().sub_ref_fst(lhs, rhs) }
    fn sub_ref_snd(&self, lhs: Self::El, rhs: &Self::El) -> Self::El { self.decorated_ring().sub_ref_snd(lhs, rhs) }
    fn add(&self, lhs: Self::El, rhs: Self::El) -> Self::El { self.decorated_ring().add(lhs, rhs) }
    fn mul(&self, lhs: Self::El, rhs: Self::El) -> Self::El { self.decorated_ring().mul(lhs, rhs) }
    fn sub(&self, lhs: Self::El, rhs: Self::El) -> Self::El { self.decorated_ring().sub(lhs, rhs) }
    fn pow(&self, basis: &Self::El, exp: u32) -> Self::El { self.decorated_ring().pow(basis, exp) }
    fn pow_big(&self, basis: &Self::El, exp: &BigInt) -> Self::El { self.decorated_ring().pow_big(basis, exp) }
    fn from_z(&self, x: i64) -> Self::El { self.decorated_ring().from_z(x) }
    fn from_z_big(&self, x: &BigInt) -> Self::El { self.decorated_ring().from_z_big(x) }
    fn is_zero(&self, val: &Self::El) -> bool { self.decorated_ring().is_zero(val) }
    fn is_one(&self, val: &Self::El) -> bool { self.decorated_ring().is_one(val) }
    fn is_neg_one(&self, val: &Self::El) -> bool { self.decorated_ring().is_neg_one(val) }
    fn is_integral(&self) -> RingPropValue { self.decorated_ring().is_integral() }
    fn characteristic(&self) -> BigInt { self.decorated_ring().characteristic() }
    fn is_field(&self) -> RingPropValue { self.decorated_ring().is_field() }
    fn is_noetherian(&self) -> bool { self.decorated_ring().is_noetherian() }
    fn div(&self, lhs: Self::El, rhs: &Self::El) -> Self::El { self.decorated_ring().div(lhs, rhs) }
    fn format(&self, el: &Self::El, f: &mut std::fmt::Formatter, in_prod: bool) -> std::fmt::Result { self.decorated_ring().format(el, f, in_prod) }
    fn format_in_brackets(&self, el: &Self::El, f: &mut std::fmt::Formatter) -> std::fmt::Result { self.decorated_ring().format_in_brackets(el, f) }
    fn add_assign_int(&self, lhs: &mut Self::El, rhs: i64) { self.decorated_ring().add_assign_int(lhs, rhs) }
    fn mul_assign_int(&self, lhs: &mut Self::El, rhs: i64) { self.decorated_ring().mul_assign_int(lhs, rhs) }
    
    fn sum<I>(&self, data: I) -> Self::El
        where I: Iterator<Item = Self::El>
    { 
        self.decorated_ring().sum(data)
    }

    fn product<I>(&self, data: I) -> Self::El 
        where I: Iterator<Item = Self::El>
    {
        self.decorated_ring().product(data)
    }
} 

impl<'a, R: RingBase> RingDecorator for &'a R {
    type DecoratedRing = R;
    fn decorated_ring(&self) -> &Self::DecoratedRing { *self }
}


impl<R> DivisibilityInfoRing for R
    where R: RingDecorator, R::DecoratedRing: DivisibilityInfoRing
{
    fn is_divisibility_computable(&self) -> RingPropValue { self.decorated_ring().is_divisibility_computable() }
    fn is_divisible_by(&self, lhs: &Self::El, rhs: &Self::El) -> bool { self.decorated_ring().is_divisible_by(lhs, rhs) }
    fn quotient(&self, lhs: &Self::El, rhs: &Self::El) -> Option<Self::El> { self.decorated_ring().quotient(lhs, rhs) }
    fn is_unit(&self, el: &Self::El) -> bool { self.decorated_ring().is_unit(el) }
}

impl<R> UfdInfoRing for R
    where R: RingDecorator, R::DecoratedRing: UfdInfoRing
{
    fn is_ufd(&self) -> RingPropValue { self.decorated_ring().is_ufd() }
    fn is_prime(&self, el: &Self::El) -> bool { self.decorated_ring().is_prime(el) }
    fn calc_factor(&self, el: &Self::El) -> Option<Self::El> { self.decorated_ring().calc_factor(el) }

    fn factor<'b>(&'b self, el: Self::El) -> VecMap<RingElWrapper<&'b Self>, usize> { 
        let result = self.decorated_ring().factor(el);
        let wrapping_ring = WrappingRing::new(self);
        return result.into_iter().map(|(el, power)| (wrapping_ring.from(el.into_val()), power)).collect();
    }
}

impl<R> HashableElRing for R
    where R: RingDecorator, R::DecoratedRing: HashableElRing
{
    fn hash<H: std::hash::Hasher>(&self, h: &mut H, el: &Self::El) { self.decorated_ring().hash(h, el) }
}

impl<R> OrderedRing for R
    where R: RingDecorator, R::DecoratedRing: OrderedRing
{
    fn cmp(&self, lhs: &Self::El, rhs: &Self::El) -> std::cmp::Ordering { self.decorated_ring().cmp(lhs, rhs) }
}

impl<'a, R> EuclideanInfoRing for R 
    where R: RingDecorator, R::DecoratedRing: EuclideanInfoRing
{
    
    fn is_euclidean(&self) -> RingPropValue { self.decorated_ring().is_euclidean() }
    fn euclidean_div_rem(&self, lhs: Self::El, rhs: &Self::El) -> (Self::El, Self::El) { self.decorated_ring().euclidean_div_rem(lhs, rhs) }
    fn euclidean_rem(&self, lhs: Self::El, rhs: &Self::El) -> Self::El { self.decorated_ring().euclidean_rem(lhs, rhs) }
    fn euclidean_div(&self, lhs: Self::El, rhs: &Self::El) -> Self::El { self.decorated_ring().euclidean_div(lhs, rhs) }
    fn euclidean_deg(&self, el: Self::El) -> BigInt { self.decorated_ring().euclidean_deg(el) }
}

impl<R> RingExtension for R
    where R: RingDecorator, R::DecoratedRing: RingExtension
{
    type BaseRing = <R::DecoratedRing as RingExtension>::BaseRing;
    type Embedding = <R::DecoratedRing as RingExtension>::Embedding;

    fn is_extension(&self) -> RingPropValue { self.decorated_ring().is_extension() }
    fn base_ring(&self) -> &Self::BaseRing { self.decorated_ring().base_ring() }
    fn embedding(&self) -> Self::Embedding { self.decorated_ring().embedding() }
    fn from(&self, el: El<Self::BaseRing>) -> El<Self> { self.decorated_ring().from(el) }
}
