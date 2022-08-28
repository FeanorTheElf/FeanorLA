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

pub trait PolyRing: Ring + CanonicalIsomorphismInfo<PolyRingImpl<Self::BaseRing>> + RingExtension {

    fn lc(&self, x: &El<Self>) -> Option<El<Self::BaseRing>>;
    fn deg(&self, x: &El<Self>) -> Option<usize>;
    fn coefficient_at(&self, x: &El<Self>, i: usize) -> El<Self::BaseRing>;
    fn scale(&self, x: &mut El<Self>, coeff: &El<Self::BaseRing>);
    fn evaluate_at<S: Ring + CanonicalEmbeddingInfo<Self::BaseRing>>(&self, f: &El<Self>, x: El<S>, ring: &S) -> El<S>;
    fn derive(&self, el: El<Self>) -> El<Self>;
    fn unknown(&self) -> El<Self>;

    fn normalize(&self, mut f: El<Self>) -> (El<Self>, El<Self::BaseRing>) {
        assert!(self.base_ring().is_field().can_use());
        let lc = self.lc(&f).unwrap().clone();
        let lc_inv = self.base_ring().div(self.base_ring().one(), &lc);
        self.scale(&mut f, &lc_inv);
        return (f, lc);
    }
}

impl<'a, R, P> CanonicalEmbeddingInfo<PolyRingImpl<R>> for &'a P
    where R: Ring, P: Ring + CanonicalEmbeddingInfo<PolyRingImpl<R>>
{
    fn has_embedding(&self, from: &PolyRingImpl<R>) -> RingPropValue { (**self).has_embedding(from) }
    fn embed(&self, from: &PolyRingImpl<R>, el: El<PolyRingImpl<R>>) -> Self::El  { (**self).embed(from, el) }
}

impl<'a, R, P> CanonicalIsomorphismInfo<PolyRingImpl<R>> for &'a P
    where R: Ring, P: Ring + CanonicalIsomorphismInfo<PolyRingImpl<R>>
{
    fn has_isomorphism(&self, from: &PolyRingImpl<R>) -> RingPropValue { (**self).has_isomorphism(from) }
    fn preimage(&self, from: &PolyRingImpl<R>, el: El<Self>) -> El<PolyRingImpl<R>>  { (**self).preimage(from, el) }
}

impl<'a, P> PolyRing for &'a P
    where P: PolyRing
{
    fn lc(&self, x: &El<Self>) -> Option<El<Self::BaseRing>> { (**self).lc(x) }
    fn deg(&self, x: &El<Self>) -> Option<usize>  { (**self).deg(x) }
    fn coefficient_at(&self, x: &El<Self>, i: usize) -> El<Self::BaseRing>  { (**self).coefficient_at(x, i) }
    fn scale(&self, x: &mut El<Self>, coeff: &El<Self::BaseRing>) { (**self).scale(x, coeff) }
    fn evaluate_at<S: Ring + CanonicalEmbeddingInfo<Self::BaseRing>>(&self, f: &El<Self>, x: El<S>, ring: &S) -> El<S> { (**self).evaluate_at(f, x, ring) }
    fn derive(&self, el: El<Self>) -> El<Self> { (**self).derive(el) }
    fn normalize(&self, f: El<Self>) -> (El<Self>, El<Self::BaseRing>) { (**self).normalize(f) }
    fn unknown(&self) -> El<Self> { (**self).unknown() }
}

impl<S, P> FnOnce<(RingElWrapper<S>, )> for RingElWrapper<P>
    where S: Ring + CanonicalEmbeddingInfo<P::BaseRing>, P: PolyRing
{
    type Output = RingElWrapper<S>;

    extern "rust-call" fn call_once(
        mut self, 
        (x, ): (RingElWrapper<S>, )
    ) -> Self::Output {
        self.call_mut((x, ))
    }
}

impl<S, P> FnMut<(RingElWrapper<S>, )> for RingElWrapper<P>
    where S: Ring + CanonicalEmbeddingInfo<P::BaseRing>, P: PolyRing
{
    extern "rust-call" fn call_mut(
        &mut self, 
        (x, ): (RingElWrapper<S>, )
    ) -> Self::Output {
        self.call((x, ))
    }
}

impl<S, P> Fn<(RingElWrapper<S>, )> for RingElWrapper<P>
    where S: Ring + CanonicalEmbeddingInfo<P::BaseRing>, P: PolyRing
{
    extern "rust-call" fn call(
        &self, 
        (x, ): (RingElWrapper<S>, )
    ) -> Self::Output {
        let (x, ring) = x.destruct();
        ring.bind_by_value(self.parent_ring().evaluate_at(self.val(), x, &ring))
    }
}
