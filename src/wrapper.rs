use super::ring::*;
use super::primitive::*;
use super::embedding::*;
use super::bigint::*;

use std::ops::{ 
    Add, Mul, Sub, Neg, Div,
    AddAssign, MulAssign, SubAssign, DivAssign
};

#[derive(Debug)]
pub struct RingElWrapper<R>
    where R: Ring
{
    ring: R,
    el: R::El
}

impl<R> RingElWrapper<R>
    where R: Ring
{
    pub fn ring(&self) -> WrappingRing<R> {
        WrappingRing {
            ring: self.ring.clone()
        }
    }

    pub fn val(&self) -> &R::El {
        &self.el
    }

    pub fn val_mut(&mut self) -> &mut R::El {
        &mut self.el
    } 

    pub fn into_val(self) -> R::El {
        self.el
    }

    pub fn pow(&self, exp: u32) -> RingElWrapper<R> {
        RingElWrapper {
            ring: self.ring.clone(),
            el: self.ring.pow(&self.el, exp)
        }
    }

    pub fn base_ring(&self) -> &R {
        &self.ring
    }
}

impl<R> Clone for RingElWrapper<R>
    where R: Ring
{
    fn clone(&self) -> Self {
        RingElWrapper {
            ring: self.ring.clone(),
            el: self.el.clone()
        }
    }
}

impl<R> Copy for RingElWrapper<R>
    where R: Ring + Copy, R::El: Copy
{}

impl<R> PartialEq<RingElWrapper<R>> for RingElWrapper<R>
    where R: Ring
{
    fn eq(&self, rhs: &RingElWrapper<R>) -> bool {
        self.ring.eq(&self.el, &rhs.el)
    }
}

impl<R> Eq for RingElWrapper<R>
    where R: Ring
{}

impl<R> Add<RingElWrapper<R>> for RingElWrapper<R>
    where R: Ring
{
    type Output = RingElWrapper<R>;

    fn add(self, rhs: RingElWrapper<R>) -> Self::Output {
        RingElWrapper {
            el: self.ring.add(self.el, rhs.el),
            ring: self.ring
        }
    }
}

impl<'a, R> Add<&'a RingElWrapper<R>> for RingElWrapper<R>
    where R: Ring
{
    type Output = RingElWrapper<R>;

    fn add(self, rhs: &'a RingElWrapper<R>) -> Self::Output {
        RingElWrapper {
            el: self.ring.add_ref(self.el, &rhs.el),
            ring: self.ring
        }
    }
}

impl<'b, R> Add<RingElWrapper<R>> for &'b RingElWrapper<R>
    where R: Ring
{
    type Output = RingElWrapper<R>;

    fn add(self, rhs: RingElWrapper<R>) -> Self::Output {
        RingElWrapper {
            el: self.ring.add_ref(rhs.el, &self.el),
            ring: rhs.ring
        }
    }
}

impl<'a, 'b, R> Add<&'a RingElWrapper<R>> for &'b RingElWrapper<R>
    where R: Ring
{
    type Output = RingElWrapper<R>;

    fn add(self, rhs: &'a RingElWrapper<R>) -> Self::Output {
        RingElWrapper {
            el: self.ring.add_ref(self.el.clone(), &rhs.el),
            ring: self.ring.clone()
        }
    }
}

impl<R> AddAssign<RingElWrapper<R>> for RingElWrapper<R>
    where R: Ring
{
    fn add_assign(&mut self, rhs: RingElWrapper<R>) {
        let recover_el = self.ring().unspecified_element();
        take_mut::take_or_recover(self, || recover_el, move |x| x + rhs);
    }
}

impl<R> Sub<RingElWrapper<R>> for RingElWrapper<R>
    where R: Ring
{
    type Output = RingElWrapper<R>;

    fn sub(self, rhs: RingElWrapper<R>) -> Self::Output {
        RingElWrapper {
            el: self.ring.sub(self.el, rhs.el),
            ring: self.ring
        }
    }
}

impl<'b, R> Sub<RingElWrapper<R>> for &'b RingElWrapper<R>
    where R: Ring
{
    type Output = RingElWrapper<R>;

    fn sub(self, rhs: RingElWrapper<R>) -> Self::Output {
        RingElWrapper {
            el: self.ring.sub_ref_fst(&self.el, rhs.el),
            ring: rhs.ring
        }
    }
}

impl<'a, R> Sub<&'a RingElWrapper<R>> for RingElWrapper<R>
    where R: Ring
{
    type Output = RingElWrapper<R>;

    fn sub(self, rhs: &'a RingElWrapper<R>) -> Self::Output {
        RingElWrapper {
            el: self.ring.sub_ref_snd(self.el, &rhs.el),
            ring: self.ring
        }
    }
}

impl<'a, 'b, R> Sub<&'a RingElWrapper<R>> for &'b RingElWrapper<R>
    where R: Ring
{
    type Output = RingElWrapper<R>;

    fn sub(self, rhs: &'a RingElWrapper<R>) -> Self::Output {
        RingElWrapper {
            el: self.ring.sub_ref_fst(&self.el, rhs.el.clone()),
            ring: self.ring.clone()
        }
    }
}

impl<R> SubAssign<RingElWrapper<R>> for RingElWrapper<R>
    where R: Ring
{
    fn sub_assign(&mut self, rhs: RingElWrapper<R>) {
        let recover_el = self.ring().unspecified_element();
        take_mut::take_or_recover(self, || recover_el, move |x| x - rhs);
    }
}

impl<R> Mul<RingElWrapper<R>> for RingElWrapper<R>
    where R: Ring
{
    type Output = RingElWrapper<R>;

    fn mul(self, rhs: RingElWrapper<R>) -> Self::Output {
        RingElWrapper {
            el: self.ring.mul(self.el, rhs.el),
            ring: self.ring
        }
    }
}

impl<'a, 'b, R> Mul<&'a RingElWrapper<R>> for &'b RingElWrapper<R>
    where R: Ring
{
    type Output = RingElWrapper<R>;

    fn mul(self, rhs: &'a RingElWrapper<R>) -> Self::Output {
        RingElWrapper {
            el: self.ring.mul_ref(&self.el, &rhs.el),
            ring: self.ring.clone()
        }
    }
}

impl<'a, R> Mul<&'a RingElWrapper<R>> for RingElWrapper<R>
    where R: Ring
{
    type Output = RingElWrapper<R>;

    fn mul(self, rhs: &'a RingElWrapper<R>) -> Self::Output {
        RingElWrapper {
            el: self.ring.mul_ref(&self.el, &rhs.el),
            ring: self.ring
        }
    }
}

impl<'b, R> Mul<RingElWrapper<R>> for &'b RingElWrapper<R>
    where R: Ring
{
    type Output = RingElWrapper<R>;

    fn mul(self, rhs: RingElWrapper<R>) -> Self::Output {
        RingElWrapper {
            el: self.ring.mul_ref(&self.el, &rhs.el),
            ring: rhs.ring
        }
    }
}

impl<R> MulAssign<RingElWrapper<R>> for RingElWrapper<R>
    where R: Ring
{
    fn mul_assign(&mut self, rhs: RingElWrapper<R>) {
        let recover_el = self.ring().unspecified_element();
        take_mut::take_or_recover(self, || recover_el, move |x| x * rhs);
    }
}

impl<R> Div<RingElWrapper<R>> for RingElWrapper<R>
    where R: Ring
{
    type Output = RingElWrapper<R>;

    default fn div(self, rhs: RingElWrapper<R>) -> Self::Output {
        assert!(self.ring.is_field().can_use());
        RingElWrapper {
            el: self.ring.div(self.el, &rhs.el),
            ring: self.ring
        }
    }
}

impl<R> Div<RingElWrapper<R>> for RingElWrapper<R>
    where R: DivisibilityInfoRing
{
    fn div(self, rhs: RingElWrapper<R>) -> Self::Output {
        assert!(self.ring.is_divisibility_computable());
        RingElWrapper {
            el: self.ring.quotient(&self.el, &rhs.el).unwrap(),
            ring: self.ring
        }
    }
}

impl<R> DivAssign<RingElWrapper<R>> for RingElWrapper<R>
    where R: Ring
{
    fn div_assign(&mut self, rhs: RingElWrapper<R>) {
        assert!(self.ring.is_field().can_use());
        let recover_el = self.ring().unspecified_element();
        take_mut::take_or_recover(self, || recover_el, move |x| x / rhs);
    }
}

impl<R> Neg for RingElWrapper<R>
    where R: Ring
{
    type Output = RingElWrapper<R>;

    fn neg(self) -> Self::Output {
        RingElWrapper {
            el: self.ring.neg(self.el),
            ring: self.ring
        }
    }
}

impl<R> Add<i64> for RingElWrapper<R>
    where R: Ring
{
    type Output = RingElWrapper<R>;

    fn add(self, rhs: i64) -> Self::Output {
        RingElWrapper {
            el: self.ring.add(self.el, (&self.ring).z_embedding(i64::RING)(rhs)),
            ring: self.ring
        }
    }
}

impl<'b, R> Add<i64> for &'b RingElWrapper<R>
    where R: Ring
{
    type Output = RingElWrapper<R>;

    fn add(self, rhs: i64) -> Self::Output {
        RingElWrapper {
            el: self.ring.add_ref((&self.ring).z_embedding(i64::RING)(rhs), &self.el),
            ring: self.ring.clone()
        }
    }
}

impl<R> Sub<i64> for RingElWrapper<R>
    where R: Ring
{
    type Output = RingElWrapper<R>;

    fn sub(self, rhs: i64) -> Self::Output {
        RingElWrapper {
            el: self.ring.sub(self.el, (&self.ring).z_embedding(i64::RING)(rhs)),
            ring: self.ring
        }
    }
}

impl<'b, R> Sub<i64> for &'b RingElWrapper<R>
    where R: Ring
{
    type Output = RingElWrapper<R>;

    fn sub(self, rhs: i64) -> Self::Output {
        RingElWrapper {
            el: self.ring.sub_ref_fst(&self.el, (&self.ring).z_embedding(i64::RING)(rhs)),
            ring: self.ring.clone()
        }
    }
}

impl<R> Mul<i64> for RingElWrapper<R>
    where R: Ring
{
    type Output = RingElWrapper<R>;

    fn mul(self, rhs: i64) -> Self::Output {
        RingElWrapper {
            el: self.ring.mul(self.el, (&self.ring).z_embedding(i64::RING)(rhs)),
            ring: self.ring
        }
    }
}

impl<'b, R> Mul<i64> for &'b RingElWrapper<R>
    where R: Ring
{
    type Output = RingElWrapper<R>;

    fn mul(self, rhs: i64) -> Self::Output {
        RingElWrapper {
            el: self.ring.mul_ref(&self.el, &(&self.ring).z_embedding(i64::RING)(rhs)),
            ring: self.ring.clone()
        }
    }
}

impl<R> Div<i64> for RingElWrapper<R>
    where R: Ring
{
    type Output = RingElWrapper<R>;

    default fn div(self, rhs: i64) -> Self::Output {
        assert!(self.ring.is_field().can_use());
        RingElWrapper {
            el: self.ring.div(self.el, &(&self.ring).z_embedding(i64::RING)(rhs)),
            ring: self.ring
        }
    }
}

impl<R> Div<i64> for RingElWrapper<R>
    where R: DivisibilityInfoRing
{
    fn div(self, rhs: i64) -> Self::Output {
        assert!(self.ring.is_divisibility_computable());
        RingElWrapper {
            el: self.ring.quotient(&self.el, &(&self.ring).z_embedding(i64::RING)(rhs)).unwrap(),
            ring: self.ring
        }
    }
}

impl<R> std::fmt::Display for RingElWrapper<R>
    where R: Ring
{
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        self.ring.format(&self.el, f, false)
    }
}

pub trait BindableElementRing: Ring {

    fn bind<'a>(&'a self, el: Self::El) -> RingElWrapper<&'a Self>;
    fn bind_ring<'a>(&'a self) -> WrappingRing<&'a Self>;
}

impl<R: Ring> BindableElementRing for R {
    
    fn bind<'a>(&'a self, el: Self::El) -> RingElWrapper<&'a Self> {
        RingElWrapper {
            ring: self,
            el: el
        }
    }
    
    fn bind_ring<'a>(&'a self) -> WrappingRing<&'a Self> {
        WrappingRing {
            ring: self
        }
    }
}

#[derive(Debug, Clone)]
pub struct WrappingRing<R>
    where R: Ring
{
    ring: R
}

impl<R> Copy for WrappingRing<R>
    where R: Ring + Copy
{}

impl<R> Ring for WrappingRing<R> 
    where R: Ring
{
    type El = RingElWrapper<R>;

    fn add_ref(&self, lhs: Self::El, rhs: &Self::El) -> Self::El {
        RingElWrapper {
            ring: self.ring.clone(),
            el: self.ring.add_ref(lhs.el, &rhs.el)
        }
    }

    fn mul_ref(&self, lhs: &Self::El, rhs: &Self::El) -> Self::El {
        RingElWrapper {
            ring: self.ring.clone(),
            el: self.ring.mul_ref(&lhs.el, &rhs.el)
        }
    }

    fn neg(&self, val: Self::El) -> Self::El {
        RingElWrapper {
            ring: self.ring.clone(),
            el: self.ring.neg(val.el)
        }
    }

    fn zero(&self) -> Self::El {
        RingElWrapper {
            ring: self.ring.clone(),
            el: self.ring.zero()
        }
    }

    fn one(&self) -> Self::El {
        RingElWrapper {
            ring: self.ring.clone(),
            el: self.ring.one()
        }
    }

    fn eq(&self, lhs: &Self::El, rhs: &Self::El) -> bool {
        self.ring.eq(&lhs.el, &rhs.el)
    }

    fn unspecified_element(&self) -> Self::El {
        RingElWrapper {
            ring: self.ring.clone(),
            el: self.ring.unspecified_element()
        }
    }

    fn sub_ref_fst(&self, lhs: &Self::El, rhs: Self::El) -> Self::El {
        RingElWrapper {
            ring: self.ring.clone(),
            el: self.ring.sub_ref_fst(&lhs.el, rhs.el)
        }
    }

    fn sub_ref_snd(&self, lhs: Self::El, rhs: &Self::El) -> Self::El {
        RingElWrapper {
            ring: self.ring.clone(),
            el: self.ring.sub_ref_snd(lhs.el, &rhs.el)
        }
    }

    fn add(&self, lhs: Self::El, rhs: Self::El) -> Self::El { 
        RingElWrapper {
            ring: self.ring.clone(),
            el: self.ring.add(lhs.el, rhs.el)
        }
    }

    fn mul(&self, lhs: Self::El, rhs: Self::El) -> Self::El { 
        RingElWrapper {
            ring: self.ring.clone(),
            el: self.ring.mul(lhs.el, rhs.el)
        }
    }

    fn sub(&self, lhs: Self::El, rhs: Self::El) -> Self::El {
        RingElWrapper {
            ring: self.ring.clone(),
            el: self.ring.sub(lhs.el, rhs.el)
        }
    }

    fn pow(&self, basis: &Self::El, exp: u32) -> Self::El 
        where Self::El: Clone
    {
        RingElWrapper {
            ring: self.ring.clone(),
            el: self.ring.pow(&basis.el, exp)
        }
    }

    fn pow_big(&self, basis: &Self::El, exp: &BigInt) -> Self::El 
        where Self::El: Clone
    {
        RingElWrapper {
            ring: self.ring.clone(),
            el: self.ring.pow_big(&basis.el, exp)
        }
    }

    fn is_zero(&self, val: &Self::El) -> bool { 
        self.ring.is_zero(&val.el)
    }

    fn is_one(&self, val: &Self::El) -> bool { 
        self.ring.is_one(&val.el)    
    }

    fn is_neg_one(&self, val: &Self::El) -> bool { 
        self.ring.is_neg_one(&val.el)
    }

    fn is_integral(&self) -> RingPropValue {
        self.ring.is_integral()
    }

    fn is_field(&self) -> RingPropValue {
        self.ring.is_field()
    }

    fn is_noetherian(&self) -> bool {
        self.ring.is_noetherian()
    }

    fn div(&self, lhs: Self::El, rhs: &Self::El) -> Self::El {
        assert!(self.ring.is_field().can_use());
        RingElWrapper {
            ring: self.ring.clone(),
            el: self.ring.div(lhs.el, &rhs.el)
        }
    }

    fn from_z_big(&self, x: BigInt) -> Self::El {
        RingElWrapper {
            ring: self.ring.clone(),
            el: self.ring.from_z_big(x)
        }
    }

    fn from_z(&self, x: i64) -> Self::El {
        RingElWrapper {
            ring: self.ring.clone(),
            el: self.ring.from_z(x)
        }
    }

    fn format(&self, el: &Self::El, f: &mut std::fmt::Formatter, in_prod: bool) -> std::fmt::Result {
        self.ring.format(&el.el, f, in_prod)
    }
}

impl<R> EuclideanInfoRing for WrappingRing<R>
    where R: EuclideanInfoRing
{
    fn is_euclidean(&self) -> RingPropValue {
        self.ring.is_euclidean()
    }

    fn euclidean_div_rem(&self, lhs: Self::El, rhs: &Self::El) -> (Self::El, Self::El) {
        let (quo, rem) = self.ring.euclidean_div_rem(lhs.el, &rhs.el);
        (RingElWrapper {
            ring: self.ring.clone(),
            el: quo
        }, RingElWrapper {
            ring: lhs.ring,
            el: rem
        })
    }

    fn euclidean_rem(&self, lhs: Self::El, rhs: &Self::El) -> Self::El { 
        RingElWrapper {
            el: self.ring.euclidean_rem(lhs.el, &rhs.el),
            ring: lhs.ring
        }
    }

    fn euclidean_div(&self, lhs: Self::El, rhs: &Self::El) -> Self::El {
        RingElWrapper {
            el: self.ring.euclidean_div(lhs.el, &rhs.el),
            ring: lhs.ring
        }
    }

    fn euclidean_deg(&self, el: Self::El) -> BigInt {
        self.ring.euclidean_deg(el.el)
    }
}

impl<R> DivisibilityInfoRing for WrappingRing<R>
    where R: DivisibilityInfoRing
{
    fn is_divisibility_computable(&self) -> bool {
        self.ring.is_divisibility_computable()
    }

    fn quotient(&self, lhs: &Self::El, rhs: &Self::El) -> Option<Self::El> { 
        self.ring.quotient(&lhs.el, &rhs.el).map(|x| RingElWrapper {
            el:x ,
            ring: self.ring.clone()
        })
    }

    fn is_divisible_by(&self, lhs: &Self::El, rhs: &Self::El) -> bool {
        self.ring.is_divisible_by(&lhs.el, &rhs.el)
    }

    fn is_unit(&self, x: &Self::El) -> bool {
        self.ring.is_unit(&x.el)
    }
}

pub struct WrappingEmbedding<R, S>
    where R: CanonicalEmbeddingInfo<S>, S: Ring
{
    base_embedding: R::Embedding,
    target_ring: R
}

impl<R, S> FnOnce<(<WrappingRing<S> as Ring>::El, )> for WrappingEmbedding<R, S> 
    where R: CanonicalEmbeddingInfo<S>, S: Ring
{
    type Output = <WrappingRing<R> as Ring>::El;

    extern "rust-call" fn call_once(
        mut self, 
        (x, ): (<WrappingRing<S> as Ring>::El, )
    ) -> Self::Output {
        self.call_mut((x, ))
    }
}

impl<R, S> FnMut<(<WrappingRing<S> as Ring>::El, )> for WrappingEmbedding<R, S> 
    where R: CanonicalEmbeddingInfo<S>, S: Ring
{
    extern "rust-call" fn call_mut(
        &mut self, 
        (x, ): (<WrappingRing<S> as Ring>::El, )
    ) -> Self::Output {
        self.call((x, ))
    }
}

impl<R, S> Fn<(<WrappingRing<S> as Ring>::El, )> for WrappingEmbedding<R, S> 
    where R: CanonicalEmbeddingInfo<S>, S: Ring
{
    extern "rust-call" fn call(
        &self, 
        (x, ): (<WrappingRing<S> as Ring>::El, )
    ) -> Self::Output {
        RingElWrapper {
            el: (self.base_embedding)(x.el), 
            ring: self.target_ring.clone()
        }
    }
}

impl<R, S> CanonicalEmbeddingInfo<WrappingRing<S>> for WrappingRing<R>
    where R: CanonicalEmbeddingInfo<S>, S: Ring
{
    type Embedding = WrappingEmbedding<R, S>;

    fn has_embedding(&self, from: &WrappingRing<S>) -> RingPropValue {
        self.ring.has_embedding(&from.ring)
    }

    fn embedding(self, from: WrappingRing<S>) -> Self::Embedding {
        WrappingEmbedding {
            target_ring: self.ring.clone(),
            base_embedding: self.ring.embedding(from.ring)
        }
    }
}