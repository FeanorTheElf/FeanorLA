use crate::la::vec::VectorView;

use super::prelude::*;
use super::la::vec::*;
use super::wrapper::*;

pub mod finite_extension_impl;

///
/// Traits for finite degree ring extensions that are generated by a single element.
/// Note that for implementation purposes, each such ring comes with a dedicated, fixed
/// generator.
/// 
pub trait FiniteExtension: RingExtension {

    type ModuleVectorView: VectorView<El<Self::BaseRing>>;

    fn generator(&self) -> El<Self>;
    fn rank_as_module(&self) -> usize;
    fn as_module_el(&self, x: El<Self>) -> Vector<Self::ModuleVectorView, El<Self::BaseRing>>;
}

impl<R> FiniteExtension for R
    where R: RingDecorator, R::DecoratedRing: FiniteExtension
{
    type ModuleVectorView = <R::DecoratedRing as FiniteExtension>::ModuleVectorView;

    fn generator(&self) -> El<Self> {
        self.decorated_ring().generator()
    }

    fn rank_as_module(&self) -> usize {
        self.decorated_ring().rank_as_module()
    }

    fn as_module_el(&self, x: El<Self>) -> Vector<Self::ModuleVectorView, El<Self::BaseRing>> {
        self.decorated_ring().as_module_el(x)
    }
}

impl<R> WrappingRing<R>
    where  R: FiniteExtension
{
    pub fn generator(&self) -> El<Self> {
        self.from(self.wrapped_ring().generator())
    }

    pub fn degree(&self) -> usize {
        assert!(self.is_field().can_use());
        self.wrapped_ring().rank_as_module()
    }

    pub fn rank_as_module(&self) -> usize {
        self.wrapped_ring().rank_as_module()
    }

    pub fn from_module_el<'a, V>(&self, data: Vector<V, RingElWrapper<&'a R::BaseRing>>) -> El<Self> 
        where V: VectorView<RingElWrapper<&'a R::BaseRing>>
    {
        let x = self.generator();
        (0..self.rank_as_module()).map(|i| self.embedding_ref()(data.at(i).clone()) * x.pow(i as u32)).sum()
    }
}

impl<R> RingElWrapper<R>
    where R: FiniteExtension
{
    pub fn into_module_el<'a>(&'a self) -> Vector<Vec<RingElWrapper<&'a R::BaseRing>>, RingElWrapper<&'a R::BaseRing>> {
        let base_ring = self.parent_ring().base_ring();
        let mut as_module_els = self.parent_ring().as_module_el(self.val().clone()).into_owned().raw_data().into_iter();
        let result = Vector::from_fn(self.parent_ring().rank_as_module(), |_| WrappingRing::new(base_ring).from(
            as_module_els.next().unwrap()
        ));
        assert!(as_module_els.next().is_none());
        return result;
    }
}