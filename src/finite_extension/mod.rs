use crate::la::vec::VectorView;

use super::prelude::*;
use super::la::vec::*;
use super::la::mat::*;
use super::poly::*;
use super::la::inversion::*;
use super::la::determinant::MatrixDeterminant;
use super::wrapper::*;
use super::eea::gcd;

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

    ///
    /// Computes the minimal polynomial of an element in the extension, that
    /// is the normalized polynomial of smallest degree such that f(a) = 0.
    /// 
    fn mipo<P>(&self, x: &El<Self>, poly_ring: &P) -> El<P>
        where P: PolyRing + EuclideanInfoRing, P::BaseRing: CanonicalEmbeddingInfo<Self::BaseRing>
    {
        compute_mipo(self, x, poly_ring)
    }

    fn multiplication_matrix(&self, x: El<Self>) -> Matrix<MatrixOwned<El<Self::BaseRing>>, El<Self::BaseRing>> {
        let mut m = Matrix::zero_ring(self.rank_as_module(), self.rank_as_module(), self.base_ring()).into_owned();
        let mut current = x;
        for i in 0..self.rank_as_module() {
            m.col_mut(i).assign(self.as_module_el(current.clone()));
            self.mul_assign(&mut current, &self.generator());
        }
        return m;
    }

    fn norm(&self, x: El<Self>) -> El<Self::BaseRing> {
        self.base_ring().matrix_determinant(self.multiplication_matrix(x))
    }

    fn trace(&self, x: El<Self>) -> El<Self::BaseRing> {
        self.base_ring().sum(self.multiplication_matrix(x).diag().iter().cloned())
    }

    fn from_module_el<V>(&self, x: Vector<V, El<Self::BaseRing>>) -> El<Self>
        where V: VectorView<El<Self::BaseRing>> 
    {
        let a = self.generator();
        self.sum(x.iter().enumerate().map(|(i, x)| self.mul(self.from(x.clone()), self.pow(&a, i as u32))))
    }
}

fn compute_mipo<R, P>(extension_ring: &R, x: &El<R>, poly_ring: &P) -> El<P>
    where R: FiniteExtension, P: PolyRing + EuclideanInfoRing, P::BaseRing: CanonicalEmbeddingInfo<R::BaseRing>
{
    assert!(poly_ring.base_ring().has_embedding(extension_ring.base_ring()).can_use());
    let mut m = Matrix::zero_ring(extension_ring.rank_as_module(), extension_ring.rank_as_module() + 1, extension_ring.base_ring()).into_owned();
    let mut current = extension_ring.one();
    for i in 0..=extension_ring.rank_as_module() {
        m.col_mut(i).assign(extension_ring.as_module_el(current.clone()));
        extension_ring.mul_assign(&mut current, x);
    }
    let x = poly_ring.unknown();
    let mut result = None;
    let embed = embedding(extension_ring.base_ring(), poly_ring.base_ring());
    for poly_vec in extension_ring.base_ring().calc_matrix_kernel_space(m).unwrap().cols() {
        let poly = poly_ring.sum(poly_vec.iter().enumerate().map(
            |(i, a)| poly_ring.mul(poly_ring.from(embed(a.clone())), poly_ring.pow(&x, i as u32))
        ));
        if let Some(current) = result {
            result = Some(gcd(poly_ring, current, poly));
        } else {
            result = Some(poly);
        }
    }
    
    return poly_ring.normalize(result.unwrap()).0;
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

    fn from_module_el<V>(&self, x: Vector<V, El<Self::BaseRing>>) -> El<Self>
        where V: VectorView<El<Self::BaseRing>> 
    {
        self.decorated_ring().from_module_el(x)
    }

    fn mipo<P>(&self, x: &El<Self>, poly_ring: &P) -> El<P>
        where P: PolyRing + EuclideanInfoRing, P::BaseRing: CanonicalEmbeddingInfo<Self::BaseRing>
    {
        self.decorated_ring().mipo(x, poly_ring)
    }

    fn multiplication_matrix(&self, x: El<Self>) -> Matrix<MatrixOwned<El<Self::BaseRing>>, El<Self::BaseRing>> {
        self.decorated_ring().multiplication_matrix(x)
    }

    fn norm(&self, x: El<Self>) -> El<Self::BaseRing> {
        self.decorated_ring().norm(x)
    }

    fn trace(&self, x: El<Self>) -> El<Self::BaseRing> {
        self.decorated_ring().trace(x)
    }
}

impl<R> WrappingRing<R>
    where R: FiniteExtension
{
    pub fn generator(&self) -> El<Self> {
        self.wrap(self.wrapped_ring().generator())
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
        self.wrap(self.wrapped_ring().from_module_el(data.access_by(|x| x.val())))
    }
}

impl<R> RingElWrapper<R>
    where R: FiniteExtension
{
    pub fn as_module_el<'a>(&'a self) -> Vector<Vec<RingElWrapper<&'a R::BaseRing>>, RingElWrapper<&'a R::BaseRing>> {
        let base_ring = self.parent_ring().base_ring();
        let mut as_module_els = self.parent_ring().as_module_el(self.val().clone()).into_owned().raw_data().into_iter();
        let result = Vector::from_fn(self.parent_ring().rank_as_module(), |_| WrappingRing::new(base_ring).wrap(
            as_module_els.next().unwrap()
        ));
        assert!(as_module_els.next().is_none());
        return result;
    }

    pub fn mipo(&self) -> RingElWrapper<PolyRingImpl<R::BaseRing>> {
        let ring = PolyRingImpl::adjoint(self.parent_ring().base_ring().clone(), "X");
        let result = self.parent_ring().mipo(self.val(), &ring);
        return RingElWrapper::new(result, ring);
    }

    pub fn trace<'a>(&'a self) -> RingElWrapper<&'a R::BaseRing> {
        RingElWrapper::new(self.parent_ring().trace(self.val().clone()), self.parent_ring().base_ring())
    }

    pub fn norm<'a>(&'a self) -> RingElWrapper<&'a R::BaseRing> {
        RingElWrapper::new(self.parent_ring().norm(self.val().clone()), self.parent_ring().base_ring())
    }
}

#[cfg(test)]
use super::rational::primitive_rational::*;
#[cfg(test)]
use finite_extension_impl::*;
#[cfg(test)]
use super::fq::fq_small::*;

#[test]
fn test_default_mipo() {
    let base_field = r64::RING;
    let extension_field: FiniteExtensionImpl<_> = FiniteExtensionImpl::new(base_field, Vector::new(vec![r64::RING.from_z(-1), r64::RING.from_z(0)]), "x");
    let poly_ring = PolyRingImpl::adjoint(base_field, "x");
    let i = extension_field.generator();
    let mipo = compute_mipo(&extension_field, &extension_field.add(i, extension_field.one()), &poly_ring);
    let x = poly_ring.unknown();
    let expected = poly_ring.add(poly_ring.pow(&x, 2), poly_ring.add(
        poly_ring.mul(x.clone(), poly_ring.from_z(-2)), 
        poly_ring.from_z(2)
    ));
    assert!(poly_ring.is_eq(&expected, &mipo));

    let mipo = compute_mipo(&extension_field, &extension_field.from_z(2), &poly_ring);
    let expected = poly_ring.add_ref(poly_ring.from_z(-2), &x);
    assert!(poly_ring.is_eq(&expected, &mipo));
}

#[test]
fn test_norm() {
    let base_field = r64::RING;
    let extension_field: FiniteExtensionImpl<_> = FiniteExtensionImpl::new(base_field, Vector::new(vec![r64::RING.from_z(-1), r64::RING.from_z(0)]), "x");
    let i = extension_field.generator();
    let norm = extension_field.norm(extension_field.add(i, extension_field.one()));
    assert_eq!(r64::from(2), norm);
}

#[test]
fn test_trace() {
    let base_field = r64::RING;
    let extension_field: FiniteExtensionImpl<_> = FiniteExtensionImpl::new(base_field, Vector::new(vec![r64::RING.from_z(-1), r64::RING.from_z(0)]), "x");
    let i = extension_field.generator();
    let trace = extension_field.trace(extension_field.add(i, extension_field.one()));
    assert_eq!(r64::from(2), trace);
}

#[test]
fn test_trace_finite_field() {
    let base_field = F2;
    let extension_field = F16;
    let a = extension_field.generator();
    let trace = extension_field.trace(extension_field.add(a, extension_field.one()));
    assert_eq!(F2.from_z(0), trace);
    assert!(false);
}