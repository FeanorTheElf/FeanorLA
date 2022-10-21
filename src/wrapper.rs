use super::ring::*;
use super::embedding::*;
use super::integer::*;
use super::primitive::*;
use super::ring_property::*;

use vector_map::VecMap;
use std::ops::{ 
    Add, Mul, Sub, Neg, Div, Rem,
    AddAssign, MulAssign, SubAssign, DivAssign
};
use std::iter::Step;
use std::marker::PhantomData;

pub struct RingElWrapper<R>
    where R: Ring
{
    ring: WrappingRing<R>,
    el: R::El
}

impl<R> RingElWrapper<R>
    where R: Ring
{
    pub fn new(el: El<R>, ring: R) -> Self {
        let ring = WrappingRing::new(ring);
        RingElWrapper { el, ring }
    }

    pub fn ring(&self) -> &WrappingRing<R> {
        &self.ring
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

    pub fn into_ring(self) -> R {
        self.ring.ring
    }

    pub fn pow(&self, exp: u32) -> RingElWrapper<R> {
        self.ring().from(self.parent_ring().pow(&self.el, exp))
    }

    pub fn pow_big(&self, exp: &BigInt) -> RingElWrapper<R> {
        self.ring().from(self.parent_ring().pow_big(&self.el, exp))
    }

    pub fn parent_ring(&self) -> &R {
        self.ring.wrapped_ring()
    }

    pub fn destruct(self) -> (R::El, R) {
        (self.el, self.ring.ring)
    }

    pub fn borrow_ring<'a>(&'a self) -> RingElWrapper<&'a R> {
        RingElWrapper::new(self.val().clone(), self.parent_ring())
    }

    pub fn is_zero(&self) -> bool {
        self.ring().is_zero(self)
    }

    pub fn is_one(&self) -> bool {
        self.ring().is_one(self)
    }

    pub fn ref_ring(&self) -> RingElWrapper<&R> {
        RingElWrapper::new(self.val().clone(), self.parent_ring())
    }
}

impl<R> RingElWrapper<R>
    where R: UfdInfoRing
{
    pub fn factor(self) -> VecMap<RingElWrapper<R>, usize> {
        let ring = self.ring.wrapped_ring();
        ring.factor(self.el).into_iter().map(|(k, v)| (RingElWrapper::new(k.into_val(), ring.clone()), v)).collect()
    }
}

impl<R> RingElWrapper<R>
    where R: DivisibilityInfoRing
{
    pub fn divides(&self, rhs: &RingElWrapper<R>) -> bool {
        self.parent_ring().is_divisible_by(rhs.val(), self.val())
    }
}

impl<R> RingElWrapper<&R>
    where R: Ring
{
    pub fn clone_ring(self) -> RingElWrapper<R> {
        RingElWrapper::new(self.el, (**self.ring.wrapped_ring()).clone())
    }
}

impl<R> RingElWrapper<R>
    where R: OrderedRing
{
    pub fn abs(self) -> Self {
        if self < 0 {
            -self
        } else {
            self
        }
    }

    pub fn sgn(&self) -> Self {
        if *self < 0 {
            -self.ring().one()
        } else if *self == 0 {
            self.ring().zero()
        } else {
            self.ring().one()
        }
    }
}

pub trait Invertible {
    fn inv(self) -> Self;
}

impl<R> Invertible for RingElWrapper<R> 
    where R: Ring
{
    default fn inv(self) -> Self {
        assert!(self.parent_ring().is_field().can_use());
        RingElWrapper {
            el: self.parent_ring().div(self.parent_ring().one(), &self.el),
            ring: self.ring
        }
    }
}

impl<R> Invertible for RingElWrapper<R> 
    where R: DivisibilityInfoRing
{
    fn inv(self) -> Self {
        assert!(self.parent_ring().is_unit(self.val()));
        RingElWrapper {
            el: self.parent_ring().quotient(&self.parent_ring().one(), &self.el).unwrap(),
            ring: self.ring
        }
    }
}

impl<R> Clone for RingElWrapper<R>
    where R: Ring
{
    fn clone(&self) -> Self {
        self.ring().from(self.el.clone())
    }
}

impl<R> Copy for RingElWrapper<R>
    where R: Ring + Copy, R::El: Copy
{}

impl<R> PartialEq<RingElWrapper<R>> for RingElWrapper<R>
    where R: Ring
{
    fn eq(&self, rhs: &RingElWrapper<R>) -> bool {
        self.parent_ring().is_eq(&self.el, &rhs.el)
    }
}

impl<R> Eq for RingElWrapper<R>
    where R: Ring
{}

impl<R> PartialEq<i64> for RingElWrapper<R>
    where R: Ring
{
    fn eq(&self, rhs: &i64) -> bool {
        self.parent_ring().is_eq(self.val(), &self.parent_ring().from_z(*rhs))
    }
}

impl<R> Add<RingElWrapper<R>> for RingElWrapper<R>
    where R: Ring
{
    type Output = RingElWrapper<R>;

    fn add(self, rhs: RingElWrapper<R>) -> Self::Output {
        RingElWrapper {
            el: self.ring.wrapped_ring().add(self.el, rhs.el),
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
            el: self.ring.wrapped_ring().add_ref(self.el, &rhs.el),
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
            el: self.parent_ring().add_ref(rhs.el, &self.el),
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
            el: self.parent_ring().add_ref(self.el.clone(), &rhs.el),
            ring: self.ring.clone()
        }
    }
}

impl<R> AddAssign<RingElWrapper<R>> for RingElWrapper<R>
    where R: Ring
{
    fn add_assign(&mut self, rhs: RingElWrapper<R>) {
        rhs.ring.wrapped_ring().add_assign(self.val_mut(), rhs.el);
    }
}

impl<R> AddAssign<i64> for RingElWrapper<R>
    where R: Ring
{
    fn add_assign(&mut self, rhs: i64) {
        self.ring.wrapped_ring().add_assign_int(&mut self.el, rhs);
    }
}

impl<'a, R> AddAssign<&'a RingElWrapper<R>> for RingElWrapper<R>
    where R: Ring
{
    fn add_assign(&mut self, rhs: &'a RingElWrapper<R>) {
        rhs.parent_ring().add_assign_ref(self.val_mut(), rhs.val());
    }
}

impl<R> Sub<RingElWrapper<R>> for RingElWrapper<R>
    where R: Ring
{
    type Output = RingElWrapper<R>;

    fn sub(self, rhs: RingElWrapper<R>) -> Self::Output {
        RingElWrapper {
            el: self.ring.wrapped_ring().sub(self.el, rhs.into_val()),
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
            el: self.parent_ring().sub_ref_fst(self.val(), rhs.el), 
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
            el: self.ring.wrapped_ring().sub_ref_snd(self.el, rhs.val()),
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
            el: self.parent_ring().sub_ref_fst(self.val(), rhs.val().clone()),
            ring: self.ring.clone()
        }
    }
}

impl<R> SubAssign<RingElWrapper<R>> for RingElWrapper<R>
    where R: Ring
{
    fn sub_assign(&mut self, rhs: RingElWrapper<R>) {
        self.ring.wrapped_ring().add_assign(&mut self.el, self.ring.wrapped_ring().neg(rhs.into_val()));
    }
}

impl<R> SubAssign<i64> for RingElWrapper<R>
    where R: Ring
{
    fn sub_assign(&mut self, rhs: i64) {
        self.ring.wrapped_ring().add_assign_int(&mut self.el, -rhs);
    }
}

impl<R> Mul<RingElWrapper<R>> for RingElWrapper<R>
    where R: Ring
{
    type Output = RingElWrapper<R>;

    fn mul(self, rhs: RingElWrapper<R>) -> Self::Output {
        RingElWrapper {
            el: self.ring.wrapped_ring().mul(self.el, rhs.into_val()),
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
            el: self.parent_ring().mul_ref(self.val(), rhs.val()),
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
            el: self.parent_ring().mul_ref(self.val(), rhs.val()),
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
            el: self.parent_ring().mul_ref(self.val(), rhs.val()),
            ring: rhs.ring
        }
    }
}

impl<R> MulAssign<RingElWrapper<R>> for RingElWrapper<R>
    where R: Ring
{
    fn mul_assign(&mut self, rhs: RingElWrapper<R>) {
        self.ring.wrapped_ring().mul_assign(&mut self.el, rhs.val());
    }
}

impl<'a, R> MulAssign<&'a RingElWrapper<R>> for RingElWrapper<R>
    where R: Ring
{
    fn mul_assign(&mut self, rhs: &'a RingElWrapper<R>) {
        self.ring.wrapped_ring().mul_assign(&mut self.el, rhs.val());
    }
}

impl<R> MulAssign<i64> for RingElWrapper<R>
    where R: Ring
{
    fn mul_assign(&mut self, rhs: i64) {
        self.ring.wrapped_ring().mul_assign_int(&mut self.el, rhs);
    }
}

impl<R> Div<RingElWrapper<R>> for RingElWrapper<R>
    where R: Ring
{
    type Output = RingElWrapper<R>;

    default fn div(self, rhs: RingElWrapper<R>) -> Self::Output {
        assert!(self.parent_ring().is_field().can_use());
        assert!(!self.parent_ring().is_zero(rhs.val()));
        RingElWrapper {
            el: self.ring.wrapped_ring().div(self.el, rhs.val()),
            ring: self.ring
        }
    }
}

impl<R> Div<RingElWrapper<R>> for RingElWrapper<R>
    where R: DivisibilityInfoRing
{
    fn div(self, rhs: RingElWrapper<R>) -> Self::Output {
        assert!(self.parent_ring().is_divisibility_computable().can_use());
        assert!(!self.parent_ring().is_zero(rhs.val()));
        RingElWrapper {
            el: self.parent_ring().quotient(self.val(), rhs.val()).unwrap(),
            ring: self.ring
        }
    }
}

impl<R> DivAssign<RingElWrapper<R>> for RingElWrapper<R>
    where R: Ring
{
    default fn div_assign(&mut self, rhs: RingElWrapper<R>) {
        assert!(self.parent_ring().is_field().can_use());
        assert!(!self.parent_ring().is_zero(rhs.val()));
        let value = std::mem::replace(&mut self.el, self.ring.wrapped_ring().unspecified_element());
        self.el = self.parent_ring().div(value, rhs.val());
    }
}

impl<R> DivAssign<RingElWrapper<R>> for RingElWrapper<R>
    where R: DivisibilityInfoRing
{
    fn div_assign(&mut self, rhs: RingElWrapper<R>) {
        assert!(self.parent_ring().is_divisibility_computable().can_use());
        assert!(!self.parent_ring().is_zero(rhs.val()));
        let value = std::mem::replace(&mut self.el, self.ring.wrapped_ring().unspecified_element());
        *self.val_mut() = self.parent_ring().quotient(&value, rhs.val()).unwrap();
    }
}

impl<'a, R> Div<&'a RingElWrapper<R>> for RingElWrapper<R>
    where R: Ring
{
    type Output = RingElWrapper<R>;

    default fn div(self, rhs: &'a RingElWrapper<R>) -> Self::Output {
        assert!(self.parent_ring().is_field().can_use());
        assert!(!self.parent_ring().is_zero(rhs.val()));
        RingElWrapper {
            el: self.ring.wrapped_ring().div(self.el, rhs.val()),
            ring: self.ring
        }
    }
}

impl<'a, R> Div<&'a RingElWrapper<R>> for RingElWrapper<R>
    where R: DivisibilityInfoRing
{
    fn div(self, rhs: &'a RingElWrapper<R>) -> Self::Output {
        assert!(self.parent_ring().is_divisibility_computable().can_use());
        assert!(!self.parent_ring().is_zero(rhs.val()));
        RingElWrapper {
            el: self.parent_ring().quotient(self.val(), rhs.val()).unwrap(),
            ring: self.ring
        }
    }
}

impl<'a, R> Div<RingElWrapper<R>> for &'a RingElWrapper<R>
    where R: Ring
{
    type Output = RingElWrapper<R>;

    default fn div(self, rhs: RingElWrapper<R>) -> Self::Output {
        assert!(self.parent_ring().is_field().can_use());
        assert!(!self.parent_ring().is_zero(rhs.val()));
        RingElWrapper {
            el: self.parent_ring().div(self.val().clone(), rhs.val()),
            ring: self.ring.clone()
        }
    }
}

impl<'a, R> Div<RingElWrapper<R>> for &'a RingElWrapper<R>
    where R: DivisibilityInfoRing
{
    fn div(self, rhs: RingElWrapper<R>) -> Self::Output {
        assert!(self.parent_ring().is_divisibility_computable().can_use());
        assert!(!self.parent_ring().is_zero(rhs.val()));
        RingElWrapper {
            el: self.parent_ring().quotient(self.val(), rhs.val()).unwrap(),
            ring: self.ring.clone()
        }
    }
}

impl<'a, R> Div<&'a RingElWrapper<R>> for &'a RingElWrapper<R>
    where R: Ring
{
    type Output = RingElWrapper<R>;

    default fn div(self, rhs: &'a RingElWrapper<R>) -> Self::Output {
        assert!(self.parent_ring().is_field().can_use());
        assert!(!self.parent_ring().is_zero(rhs.val()));
        RingElWrapper {
            el: self.parent_ring().div(self.val().clone(), rhs.val()),
            ring: self.ring.clone()
        }
    }
}

impl<'a, R> Div<&'a RingElWrapper<R>> for &'a RingElWrapper<R>
    where R: DivisibilityInfoRing
{
    fn div(self, rhs: &'a RingElWrapper<R>) -> Self::Output {
        assert!(self.parent_ring().is_divisibility_computable().can_use());
        assert!(!self.parent_ring().is_zero(rhs.val()));
        RingElWrapper {
            el: self.parent_ring().quotient(self.val(), rhs.val()).unwrap(),
            ring: self.ring.clone()
        }
    }
}

impl<'a, R> DivAssign<&'a RingElWrapper<R>> for RingElWrapper<R>
    where R: Ring
{
    default fn div_assign(&mut self, rhs: &'a RingElWrapper<R>) {
        assert!(self.parent_ring().is_field().can_use());
        assert!(!self.parent_ring().is_zero(rhs.val()));
        let value = std::mem::replace(&mut self.el, self.ring.wrapped_ring().unspecified_element());
        self.el = self.parent_ring().div(value, rhs.val());
    }
}

impl<'a, R> DivAssign<&'a RingElWrapper<R>> for RingElWrapper<R>
    where R: DivisibilityInfoRing
{
    fn div_assign(&mut self, rhs: &'a RingElWrapper<R>) {
        assert!(self.parent_ring().is_divisibility_computable().can_use());
        assert!(!self.parent_ring().is_zero(rhs.val()));
        let value = std::mem::replace(&mut self.el, self.ring.wrapped_ring().unspecified_element());
        self.el = self.parent_ring().quotient(&value, rhs.val()).unwrap();
    }
}

impl<R> Neg for RingElWrapper<R>
    where R: Ring
{
    type Output = RingElWrapper<R>;

    fn neg(self) -> Self::Output {
        RingElWrapper {
            el: self.ring.wrapped_ring().neg(self.el),
            ring: self.ring
        }
    }
}

impl<'a, R> Neg for &'a RingElWrapper<R>
    where R: Ring
{
    type Output = RingElWrapper<R>;

    fn neg(self) -> Self::Output {
        -self.clone()
    }
}

impl<R> Add<i64> for RingElWrapper<R>
    where R: Ring
{
    type Output = RingElWrapper<R>;

    fn add(mut self, rhs: i64) -> Self::Output {
        self.ring.wrapped_ring().add_assign_int(&mut self.el, rhs);
        return self;
    }
}

impl<'b, R> Add<i64> for &'b RingElWrapper<R>
    where R: Ring
{
    type Output = RingElWrapper<R>;

    fn add(self, rhs: i64) -> Self::Output {
        let mut result = self.clone();
        result += rhs;
        return result;
    }
}

impl<R> Sub<i64> for RingElWrapper<R>
    where R: Ring
{
    type Output = RingElWrapper<R>;

    fn sub(mut self, rhs: i64) -> Self::Output {
        self += -rhs;
        return self;
    }
}

impl<'b, R> Sub<i64> for &'b RingElWrapper<R>
    where R: Ring
{
    type Output = RingElWrapper<R>;

    fn sub(self, rhs: i64) -> Self::Output {
        let mut result = self.clone();
        result -= rhs;
        return result;
    }
}

impl<R> Mul<i64> for RingElWrapper<R>
    where R: Ring
{
    type Output = RingElWrapper<R>;

    fn mul(mut self, rhs: i64) -> Self::Output {
        self.ring.wrapped_ring().mul_assign_int(&mut self.el, rhs);
        return self;
    }
}

impl<'b, R> Mul<i64> for &'b RingElWrapper<R>
    where R: Ring
{
    type Output = RingElWrapper<R>;

    fn mul(self, rhs: i64) -> Self::Output {
        let mut result = self.clone();
        result *= rhs;
        return result;
    }
}

impl<R> Div<i64> for RingElWrapper<R>
    where R: Ring
{
    type Output = RingElWrapper<R>;

    default fn div(self, rhs: i64) -> Self::Output {
        assert!(self.parent_ring().is_field().can_use());
        assert!(rhs != 0);
        RingElWrapper {
            el: self.ring.wrapped_ring().div(self.el, &self.ring.wrapped_ring().from_z(rhs)),
            ring: self.ring
        }
    }
}

impl<R> Div<i64> for RingElWrapper<R>
    where R: DivisibilityInfoRing
{
    fn div(self, rhs: i64) -> Self::Output {
        assert!(self.parent_ring().is_divisibility_computable().can_use());
        assert!(rhs != 0);
        RingElWrapper {
            el: self.parent_ring().quotient(self.val(), &self.parent_ring().from_z(rhs)).unwrap(),
            ring: self.ring
        }
    }
}

impl<R> Rem<RingElWrapper<R>> for RingElWrapper<R>
    where R: EuclideanInfoRing
{
    type Output = RingElWrapper<R>;

    fn rem(self, rhs: RingElWrapper<R>) -> Self::Output {
        assert!(self.parent_ring().is_euclidean().can_use());
        RingElWrapper {
            el: self.ring.wrapped_ring().euclidean_rem(self.el, rhs.val()),
            ring: self.ring
        }
    }
}

impl<'b, R> Rem<&'b RingElWrapper<R>> for RingElWrapper<R>
    where R: EuclideanInfoRing
{
    type Output = RingElWrapper<R>;

    fn rem(self, rhs: &'b RingElWrapper<R>) -> Self::Output {
        assert!(self.parent_ring().is_euclidean().can_use());
        RingElWrapper {
            el: self.ring.wrapped_ring().euclidean_rem(self.el, rhs.val()),
            ring: self.ring
        }
    }
}

impl<R> std::cmp::PartialOrd for RingElWrapper<R>
    where R: OrderedRing
{
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl<R> std::cmp::PartialOrd<i64> for RingElWrapper<R>
    where R: OrderedRing
{
    fn partial_cmp(&self, other: &i64) -> Option<std::cmp::Ordering> {
        Some(self.parent_ring().cmp(self.val(), &self.parent_ring().from_z(*other)))
    }
}

impl<R> std::cmp::Ord for RingElWrapper<R>
    where R: OrderedRing
{
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.parent_ring().cmp(&self.el, &other.el)
    }
}

impl<R> std::fmt::Display for RingElWrapper<R>
    where R: Ring
{
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        self.parent_ring().format(&self.el, f, false)
    }
}

impl<R> std::fmt::Debug for RingElWrapper<R>
    where R: Ring
{
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        <Self as std::fmt::Display>::fmt(self, f)
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct WrappingRing<R>
    where R: Ring
{
    ring: R
}

pub struct LiftedHom<S, R, F> 
    where S: Ring, R: Ring, F: Fn(El<S>) -> El<R>
{
    base: F,
    from: PhantomData<S>,
    to: WrappingRing<R>
}

impl<R> WrappingRing<R>
    where R: Ring
{
    pub const fn new(ring: R) -> Self {
        WrappingRing { ring }
    }

    pub fn from(&self, x: El<R>) -> El<Self> {
        RingElWrapper {
            ring: self.clone(),
            el: x
        }
    }

    pub fn borrow_ring<'a>(&'a self) -> WrappingRing<&'a R> {
        WrappingRing {
            ring: &self.ring
        }
    }

    pub fn lift_hom<S, F>(&self, e: F) -> LiftedHom<S, R, F>
        where S: Ring, F: Fn(El<S>) -> El<R>
    {
        LiftedHom {
            base: e,
            from: PhantomData,
            to: WrappingRing::new(self.ring.clone())
        }
    }

    pub fn wrapped_ring(&self) -> &R {
        &self.ring
    }

    pub fn wrapped_ring_mut(&mut self) -> &mut R {
        &mut self.ring
    }

    pub fn wrapping_embedding<'a>(&'a self) -> impl 'a + Clone + Fn(R::El) -> <Self as RingBase>::El {
        move |x| self.from(x)
    }

    pub fn map_ring<F, S>(self, f: F) -> WrappingRing<S>
        where F: FnOnce(R) -> S, S: Ring
    {
        WrappingRing {
            ring: f(self.ring)
        }
    }
}

impl<R> WrappingRing<R>
    where R: RingExtension
{
    pub fn base_ring(&self) -> WrappingRing<&R::BaseRing> {
        WrappingRing::new(self.wrapped_ring().base_ring())
    }

    pub fn embedding(&self) -> LiftedHom<R::BaseRing, R, R::Embedding> {
        self.lift_hom(self.wrapped_ring().embedding())
    }

    pub fn embedding_ref<'a, 'b>(&'a self) -> LiftedHom<&'b R::BaseRing, R, R::Embedding> {
        self.lift_hom(self.wrapped_ring().embedding())
    }
}

impl<S, R, F> PartialEq for LiftedHom<S, R, F>
    where S: Ring, R: Ring, F: PartialEq + Fn(El<S>) -> El<R>
{
    fn eq(&self, rhs: &LiftedHom<S, R, F>) -> bool {
        self.base == rhs.base
    }
}

impl<S, R, F> Clone for LiftedHom<S, R, F>
    where S: Ring, R: Ring, F: Clone + Fn(El<S>) -> El<R>
{
    fn clone(&self) -> LiftedHom<S, R, F> {
        LiftedHom {
            base: self.base.clone(),
            from: PhantomData,
            to: self.to.clone()
        }
    }
}

impl<S, R, F> Copy for LiftedHom<S, R, F>
    where S: Ring, R: Ring + Copy, F: Copy + Fn(El<S>) -> El<R>
{}

impl<S, R, F> FnOnce<(RingElWrapper<S>,)> for LiftedHom<S, R, F>
    where S: Ring, R: Ring, F: Fn(El<S>) -> El<R>
{
    type Output = RingElWrapper<R>;

    extern "rust-call" fn call_once(
        mut self, 
        (x, ): (RingElWrapper<S>,)
    ) -> Self::Output {
        self.call_mut((x, ))
    }
}

impl<S, R, F> FnMut<(RingElWrapper<S>,)> for LiftedHom<S, R, F>
    where S: Ring, R: Ring, F: Fn(El<S>) -> El<R>
{
    extern "rust-call" fn call_mut(
        &mut self, 
        (x, ): (RingElWrapper<S>,)
    ) -> Self::Output {
        self.call((x, ))
    }
}

impl<S, R, F> Fn<(RingElWrapper<S>,)> for LiftedHom<S, R, F>
    where S: Ring, R: Ring, F: Fn(El<S>) -> El<R>
{
    extern "rust-call" fn call(
        &self, 
        (x, ): (RingElWrapper<S>,)
    ) -> Self::Output {
        self.to.from((self.base)(x.el))
    }
}

impl<'a, R> WrappingRing<&'a R>
    where R: Ring
{
    pub fn clone_ring(&self) -> WrappingRing<R> {
        WrappingRing {
            ring: self.ring.clone()
        }
    }
}

impl<R> RingBase for WrappingRing<R> 
    where R: Ring
{
    type El = RingElWrapper<R>;

    fn add_ref(&self, lhs: Self::El, rhs: &Self::El) -> Self::El {
        self.from(self.ring.add_ref(lhs.el, &rhs.el))
    }

    fn mul_ref(&self, lhs: &Self::El, rhs: &Self::El) -> Self::El {
        self.from(self.ring.mul_ref(&lhs.el, &rhs.el))
    }

    fn neg(&self, val: Self::El) -> Self::El {
        self.from(self.ring.neg(val.el))
    }

    fn zero(&self) -> Self::El {
        self.from(self.ring.zero())
    }

    fn one(&self) -> Self::El {
        self.from(self.ring.one())
    }

    fn is_eq(&self, lhs: &Self::El, rhs: &Self::El) -> bool {
        self.ring.is_eq(&lhs.el, &rhs.el)
    }

    fn unspecified_element(&self) -> Self::El {
        self.from(self.ring.unspecified_element())
    }

    fn add_assign(&self, lhs: &mut Self::El, rhs: Self::El) { 
        self.wrapped_ring().add_assign(lhs.val_mut(), rhs.el)
    }

    fn add_assign_ref(&self, lhs: &mut Self::El, rhs: &Self::El) { 
        self.wrapped_ring().add_assign_ref(lhs.val_mut(), &rhs.el)
    }

    fn sub_ref_fst(&self, lhs: &Self::El, rhs: Self::El) -> Self::El {
        self.from(self.ring.sub_ref_fst(&lhs.el, rhs.el))
    }

    fn sub_ref_snd(&self, lhs: Self::El, rhs: &Self::El) -> Self::El {
        self.from(self.ring.sub_ref_snd(lhs.el, &rhs.el))
    }

    fn add(&self, lhs: Self::El, rhs: Self::El) -> Self::El { 
        self.from(self.ring.add(lhs.el, rhs.el))
    }

    fn mul(&self, lhs: Self::El, rhs: Self::El) -> Self::El { 
        self.from(self.ring.mul(lhs.el, rhs.el))
    }

    fn sub(&self, lhs: Self::El, rhs: Self::El) -> Self::El {
        self.from(self.ring.sub(lhs.el, rhs.el))
    }

    fn pow(&self, basis: &Self::El, exp: u32) -> Self::El 
        where Self::El: Clone
    {
        self.from(self.ring.pow(&basis.el, exp))
    }

    fn pow_big(&self, basis: &Self::El, exp: &BigInt) -> Self::El 
        where Self::El: Clone
    {
        self.from(self.ring.pow_big(&basis.el, exp))
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

    fn characteristic(&self) -> BigInt {
        self.wrapped_ring().characteristic()
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
        assert!(!self.is_zero(rhs));
        self.from(self.ring.div(lhs.el, &rhs.el))
    }

    fn from_z_big(&self, x: &BigInt) -> Self::El {
        self.from(self.ring.from_z_big(x))
    }

    fn from_z(&self, x: i64) -> Self::El {
        self.from(self.ring.from_z(x))
    }

    fn format(&self, el: &Self::El, f: &mut std::fmt::Formatter, in_prod: bool) -> std::fmt::Result {
        self.ring.format(&el.el, f, in_prod)
    }
}

impl<R> SingletonRing for WrappingRing<R>
    where R: SingletonRing
{
    fn singleton() -> Self {
        WrappingRing::new(R::singleton())
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
        (self.from(quo), RingElWrapper {
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

impl<R> HashableElRing for WrappingRing<R>
    where R: HashableElRing
{
    fn hash<H: std::hash::Hasher>(&self, h: &mut H, el: &Self::El) {
        el.parent_ring().hash(h, &el.el)
    }
}

impl<R> OrderedRing for WrappingRing<R>
    where R: OrderedRing
{
    fn cmp(&self, lhs: &Self::El, rhs: &Self::El) -> std::cmp::Ordering {
        self.wrapped_ring().cmp(lhs.val(), rhs.val())
    }
}

impl<R> DivisibilityInfoRing for WrappingRing<R>
    where R: DivisibilityInfoRing
{
    fn is_divisibility_computable(&self) -> RingPropValue {
        self.ring.is_divisibility_computable()
    }

    fn quotient(&self, lhs: &Self::El, rhs: &Self::El) -> Option<Self::El> { 
        self.ring.quotient(&lhs.el, &rhs.el).map(|x| RingElWrapper::new(x, self.ring.clone()))
    }

    fn is_divisible_by(&self, lhs: &Self::El, rhs: &Self::El) -> bool {
        self.ring.is_divisible_by(&lhs.el, &rhs.el)
    }

    fn is_unit(&self, x: &Self::El) -> bool {
        self.ring.is_unit(&x.el)
    }
}

impl<R> UfdInfoRing for WrappingRing<R>
    where R: UfdInfoRing
{
    fn is_ufd(&self) -> RingPropValue {
        self.wrapped_ring().is_ufd()
    }

    fn is_prime(&self, el: &Self::El) -> bool {
        self.wrapped_ring().is_prime(el.val())
    }

    fn calc_factor(&self, el: &Self::El) -> Option<Self::El> {
        self.wrapped_ring().calc_factor(el.val()).map(|x| self.from(x))
    }

    fn factor<'a>(&'a self, el: Self::El) -> VecMap<RingElWrapper<&'a Self>, usize> {
        let wrapping_ring = WrappingRing::new(self);
        self.wrapped_ring().factor(el.into_val()).into_iter().map(|(key, val)| 
            (wrapping_ring.from(self.from(key.into_val())), val)
        ).collect()
    }
}

impl<R, S> CanonicalEmbeddingInfo<WrappingRing<S>> for WrappingRing<R>
    where R: Ring + CanonicalEmbeddingInfo<S>, S: Ring
{
    fn has_embedding(&self, from: &WrappingRing<S>) -> RingPropValue {
        self.ring.has_embedding(&from.ring)
    }

    fn embed(&self, from: &WrappingRing<S>, el: <WrappingRing<S> as RingBase>::El) -> Self::El {
        RingElWrapper::new(
            self.wrapped_ring().embed(from.wrapped_ring(), el.el),
            self.wrapped_ring().clone()
        )
    }
}

impl<R, T: RingEl> CanonicalEmbeddingInfo<StaticRing<T>> for WrappingRing<R>
    where R: Ring + CanonicalEmbeddingInfo<StaticRing<T>>
{
    fn has_embedding(&self, from: &StaticRing<T>) -> RingPropValue {
        self.ring.has_embedding(from)
    }

    fn embed(&self, from: &StaticRing<T>, el: T) -> Self::El {
        RingElWrapper::new(
            self.wrapped_ring().embed(from, el),
            self.wrapped_ring().clone()
        )
    }
}

impl<R> CanonicalEmbeddingInfo<BigIntRing> for WrappingRing<R>
    where R: Ring + CanonicalEmbeddingInfo<BigIntRing>
{
    fn has_embedding(&self, from: &BigIntRing) -> RingPropValue {
        self.ring.has_embedding(from)
    }

    fn embed(&self, from: &BigIntRing, el: BigInt) -> Self::El {
        RingElWrapper::new(
            self.wrapped_ring().embed(from, el),
            self.wrapped_ring().clone()
        )
    }
}

impl<R, S> CanonicalIsomorphismInfo<WrappingRing<S>> for WrappingRing<R>
    where R: Ring + CanonicalIsomorphismInfo<S>, S: Ring
{
    fn has_isomorphism(&self, from: &WrappingRing<S>) -> RingPropValue {
        self.ring.has_isomorphism(&from.ring)
    }

    fn preimage(&self, from: &WrappingRing<S>, el: Self::El) -> <WrappingRing<S> as RingBase>::El {
        RingElWrapper::new(
            self.wrapped_ring().preimage(from.wrapped_ring(), el.el),
            from.wrapped_ring().clone()
        )
    }
}

impl<R, T: RingEl> CanonicalIsomorphismInfo<StaticRing<T>> for WrappingRing<R>
    where R: Ring + CanonicalIsomorphismInfo<StaticRing<T>>
{
    fn has_isomorphism(&self, from: &StaticRing<T>) -> RingPropValue {
        self.ring.has_isomorphism(from)
    }

    fn preimage(&self, from: &StaticRing<T>, el: Self::El) -> T {
        self.wrapped_ring().preimage(from, el.into_val())
    }
}

impl<R> CanonicalIsomorphismInfo<BigIntRing> for WrappingRing<R>
    where R: Ring + CanonicalIsomorphismInfo<BigIntRing>
{
    fn has_isomorphism(&self, from: &BigIntRing) -> RingPropValue {
        self.ring.has_isomorphism(from)
    }

    fn preimage(&self, from: &BigIntRing, el: Self::El) -> BigInt {
        self.wrapped_ring().preimage(from, el.into_val())
    }
}

impl<R> std::iter::Sum for RingElWrapper<R>
    where R: Ring
{
    fn sum<I>(mut iter: I) -> Self
        where I: Iterator<Item = Self>
    {
        let el = iter.next().unwrap();
        return RingElWrapper {
            el: el.ring.wrapped_ring().add(el.ring.wrapped_ring().sum(iter.map(|x| x.el)), el.el),
            ring: el.ring
        }
    }
}

impl<R> std::iter::Product for RingElWrapper<R>
    where R: Ring
{
    fn product<I>(mut iter: I) -> Self
        where I: Iterator<Item = Self>
    {
        let el = iter.next().unwrap();
        return RingElWrapper {
            el: el.ring.wrapped_ring().mul(el.ring.wrapped_ring().product(iter.map(|x| x.el)), el.el),
            ring: el.ring
        }
    }
}

impl<R> std::hash::Hash for RingElWrapper<R>
    where R: HashableElRing
{
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.parent_ring().hash(state, &self.el)
    }
}

impl<I> Step for RingElWrapper<I>
    where I: IntegerRing
{
    fn steps_between(start: &Self, end: &Self) -> Option<usize> {
        let difference = end - start;
        assert!(difference >= 0);
        if difference < i64::MAX {
            let difference = start.parent_ring().preimage(&i64::RING, difference.el);
            return Some(difference as usize)
        } else {
            return None;
        }
    }

    fn forward_checked(start: Self, count: usize) -> Option<Self> {
        Some(start + (count as i64))
    }

    fn backward_checked(start: Self, count: usize) -> Option<Self> {
        Some(start - (count as i64))
    }
}