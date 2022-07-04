use super::ring::*;
use super::embedding::*;
use super::ring_property::*;
use super::wrapper::*;

pub mod ops;
pub mod factoring;
pub mod uni_var;
pub mod multi_var;
pub mod sumation;

pub use uni_var::*;
pub use multi_var::*;

pub trait UnivarPolyRing: Ring + CanonicalIsomorphismInfo<PolyRing<Self::BaseRing>> + RingExtension {

    fn lc(&self, x: &El<Self>) -> Option<El<Self::BaseRing>>;
    fn deg(&self, x: &El<Self>) -> Option<usize>;
    fn coefficient_at(&self, x: &El<Self>, i: usize) -> El<Self::BaseRing>;
    fn scale(&self, x: &mut El<Self>, coeff: &El<Self::BaseRing>);
}

impl<'a, R, P> CanonicalEmbeddingInfo<PolyRing<R>> for &'a P
    where R: Ring, P: Ring + CanonicalEmbeddingInfo<PolyRing<R>>
{
    fn has_embedding(&self, from: &PolyRing<R>) -> RingPropValue { (**self).has_embedding(from) }
    fn embed(&self, from: &PolyRing<R>, el: El<PolyRing<R>>) -> Self::El  { (**self).embed(from, el) }
}

impl<'a, R, P> CanonicalIsomorphismInfo<PolyRing<R>> for &'a P
    where R: Ring, P: Ring + CanonicalIsomorphismInfo<PolyRing<R>>
{
    fn has_isomorphism(&self, from: &PolyRing<R>) -> RingPropValue { (**self).has_isomorphism(from) }
    fn preimage(&self, from: &PolyRing<R>, el: El<Self>) -> El<PolyRing<R>>  { (**self).preimage(from, el) }
}

impl<'a, P> UnivarPolyRing for &'a P
    where P: UnivarPolyRing
{
    fn lc(&self, x: &El<Self>) -> Option<El<Self::BaseRing>> { (**self).lc(x) }
    fn deg(&self, x: &El<Self>) -> Option<usize>  { (**self).deg(x) }
    fn coefficient_at(&self, x: &El<Self>, i: usize) -> El<Self::BaseRing>  { (**self).coefficient_at(x, i) }
    fn scale(&self, x: &mut El<Self>, coeff: &El<Self::BaseRing>) { (**self).scale(x, coeff) }
}