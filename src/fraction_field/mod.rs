use super::prelude::*;
use super::eea::gcd;
use super::wrapper::*;

pub mod fraction_field_impl;

pub trait ReducableElementRing: Ring {

    fn reduce(&self, num: El<Self>, den: El<Self>) -> (El<Self>, El<Self>);
}

impl<R: Ring> ReducableElementRing for R {

    default fn reduce(&self, num: El<Self>, den: El<Self>) -> (El<Self>, El<Self>) {
        (num, den)
    }
}

impl<R: EuclideanInfoRing> ReducableElementRing for R {

    fn reduce(&self, num: El<Self>, den: El<Self>) -> (El<Self>, El<Self>) {
        assert!(self.is_euclidean().can_use());
        let d = gcd(&self, num.clone(), den.clone());
        return (self.euclidean_div(num, &d), self.euclidean_div(den, &d));
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

impl<R: FractionField> RingElWrapper<R> 
{
    pub fn num(&self) -> RingElWrapper<&R::BaseRing> {
        self.parent_ring().base_ring().bind(self.parent_ring().num(self.val()).clone())
    }

    pub fn den(&self) -> RingElWrapper<&R::BaseRing> {
        self.parent_ring().base_ring().bind(self.parent_ring().den(self.val()).clone())
    }
}