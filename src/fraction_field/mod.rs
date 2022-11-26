use self::fraction_field_impl::FractionFieldImpl;

use super::prelude::*;
use super::eea::gcd;
use super::wrapper::*;

pub mod fraction_field_impl;

pub trait ReducableElementRing: DivisibilityInfoRing {

    ///
    /// Finds a non zero divisor in the current ring such that the division
    /// by this element makes all given elements as "small" as possible.
    /// As an exception to the above, the result may be 0 if the given iterator
    /// is empty. 
    /// Naturally, the division by the resulting element must be well-defined
    /// in the ring, i.e. `self.is_divisible_by(x, self.reduce_divisor(iter))`
    /// must be true if x is an element of iter.
    /// Here "small" is an implementation-defined concept, but the general
    /// idea is that storing and computing with "small" values should be
    /// fairly efficient.
    /// 
    /// This is mainly a performance tool, to keep e.g. projective coordinates,
    /// fractions, polynomials up to units small while doing computations.
    /// 
    /// Note that there is a default implementation that returns 1.
    /// 
    fn reduce_divisor<I: Iterator<Item = El<Self>>>(&self, elements: I) -> El<Self>;
}

impl<R: DivisibilityInfoRing> ReducableElementRing for R {

    default fn reduce_divisor<I: Iterator<Item = El<Self>>>(&self, _: I) -> El<Self> {
        self.one()
    }
}

impl<R: EuclideanInfoRing> ReducableElementRing for R {

    default fn reduce_divisor<I: Iterator<Item = El<Self>>>(&self, elements: I) -> El<Self> {
        if self.is_euclidean().can_use() {
            elements.fold(self.zero(), |a, b| gcd(self, a, b))
        } else {
            self.one()
        }
    }
}

pub trait FractionField: RingExtension {

    fn num<'a>(&self, el: &'a El<Self>) -> &'a El<Self::BaseRing>;
    fn den<'a>(&self, el: &'a El<Self>) -> &'a El<Self::BaseRing>;
}

impl<'a, R: FractionField> FractionField for &'a R {
    
    fn num<'b>(&self, el: &'b El<Self>) -> &'b El<Self::BaseRing> { (**self).num(el) }
    fn den<'b>(&self, el: &'b El<Self>) -> &'b El<Self::BaseRing> { (**self).den(el) }
}

impl<R: FractionField> RingElWrapper<R> {

    pub fn num(&self) -> RingElWrapper<&R::BaseRing> {
        RingElWrapper::new(self.parent_ring().num(self.val()).clone(), self.parent_ring().base_ring())
    }

    pub fn den(&self) -> RingElWrapper<&R::BaseRing> {
        RingElWrapper::new(self.parent_ring().den(self.val()).clone(), self.parent_ring().base_ring())
    }
}

impl<R: Ring> WrappingRing<R> {

    pub fn std_fraction_field(&self) -> WrappingRing<FractionFieldImpl<&R>> {
        WrappingRing::new(FractionFieldImpl::new(self.wrapped_ring()))
    }
}