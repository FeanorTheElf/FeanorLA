use super::prelude::*;
use super::integer::bigint_soo::*;

use std::marker::PhantomData;

///
/// Trait for rings S such that there is a canonical embedding (injective
/// ring homomorphism) R -> S. Deciding whether an embedding is canonical is left
/// to the implementor, but we require that all composed canonical embeddings
/// between the same rings are the same, i.e. if there are rings R1 = S1,
/// Rn = Sm and R2, ..., R(n - 1) and S2, ..., S(m - 1) with canonical embeddings
/// R1 -> ... -> Rn and S1 -> ... -> Sm, then the composition of both embedding
/// sequences must give the same value on all inputs.
/// 
/// To get the embedding, use the global function `embedding()`.
/// 
/// # Note on blanket implementations
/// 
/// Sadly, we cannot provide all blanket implementations for references,
/// as the implementations
/// ```text
///     impl CanonicalEmbeddingInfo<&R> for S
///     impl CanonicalEmbeddingInfo<R> for &S
/// ```
/// leave two possibilities for &S: CanonicalEmbeddingInfo<&R> and we do
/// not yet have lattice traits.
///
/// So we use the following rationale:
///  - we definitely want to have `impl CanonicalEmbeddingInfo<&R> for &S`, so
///    we are already quite limited options for further blanket implementations.
///  - we can still make blanket implementations for somewhat concrete types (that
///    do not contain references for sure); up to now, these are `StaticRing` and
///    `BigIntSOORing`, as they are required for the definition of the Integers.
///  - we require `CanonicalIsomorphismInfo<R> for R` for each ring. However, we 
///    do not provide a blanket implementation, so each ring implementation must
///    implement this in the standard way
///  - to implement [`crate::ring::Ring`], it is required to implement
///    `CanonicalIsomorphismInfo<Self>`, so this is always ensured, however without
///    providing a blanket implementation.
///
/// The same holds for `CanonicalIsomorphismInfo`.
/// 
/// # Note on the relationship with [`crate::ring::Ring::from_z()`]
/// 
/// In many cases, the functionality given by [`crate::ring::RingExtension::from()`] or the `from_z`-family
/// [`crate::ring::Ring::from_z()`], [`crate::ring::Ring::from_z_big()`] and [`crate::ring::Ring::from_z_gen()`] is the
/// same. In these cases, you should provide the actual implementation in `CanonicalEmbeddingInfo`, and then delegate
/// from the appropriate functions. In particular, this helps to avoid infinite loops of the form
/// ```text
/// IntegerRing1::from_z_gen() -> IntegerRing2::preimage() -> IntegerRing1::from_z_gen() -> ...
/// ```
/// 
/// I have made this decision mostly on the fact that `CanonicalEmbeddingInfo` is the "general" location
/// for this kind of functionality, and the others are "special cases". Furthermore, the implementation
/// as trait for `CanonicalEmbeddingInfo` is more flexible, e.g. w.r.t. specialization.
/// 
pub trait CanonicalEmbeddingInfo<R>: RingBase
    where R: RingBase
{
    fn has_embedding(&self, from: &R) -> RingPropValue;
    fn embed(&self, from: &R, el: R::El) -> Self::El;
}

///
/// Trait for rings S that are canonically isomorphic.
/// 
/// The definition of this trait is asymmetric, i.e. we can
/// have R: CanonicalIsomorphismInfo<S> without S: CanonicalIsomorphismInfo<R>,
/// to allow implementing this trait for new types.
/// 
/// # Note on blanket implementations
/// 
/// Sadly, we cannot provide blanket implementations for all reference-to-ring
/// combinations. For more details, see the notes on [`CanonicalEmbeddingInfo`].
/// 
pub trait CanonicalIsomorphismInfo<R>: CanonicalEmbeddingInfo<R>
    where R: RingBase
{
    fn has_isomorphism(&self, from: &R) -> RingPropValue;
    fn preimage(&self, from: &R, el: Self::El) -> R::El;
}

impl<R, T> CanonicalEmbeddingInfo<StaticRing<T>> for R
    where T: RingEl, R: RingDecorator, R::DecoratedRing: CanonicalEmbeddingInfo<StaticRing<T>>
{
    fn has_embedding(&self, from: &StaticRing<T>) -> RingPropValue {
        self.decorated_ring().has_embedding(from)
    }

    fn embed(&self, from: &StaticRing<T>, el: T) -> Self::El {
        self.decorated_ring().embed(from, el)
    }
}

impl<R> CanonicalEmbeddingInfo<BigIntSOORing> for R
    where R: RingDecorator, R::DecoratedRing: CanonicalEmbeddingInfo<BigIntSOORing>
{
    fn has_embedding(&self, from: &BigIntSOORing) -> RingPropValue {
        self.decorated_ring().has_embedding(from)
    }

    fn embed(&self, from: &BigIntSOORing, el: BigIntSOO) -> Self::El {
        self.decorated_ring().embed(from, el)
    }
}

impl<R, S> CanonicalEmbeddingInfo<R> for S
    where R: RingDecorator, S: RingDecorator, S::DecoratedRing: CanonicalEmbeddingInfo<R::DecoratedRing>
{
    fn has_embedding(&self, from: &R) -> RingPropValue {
        self.decorated_ring().has_embedding(from.decorated_ring())
    }

    fn embed(&self, from: &R, el: <R as RingBase>::El) -> Self::El {
        self.decorated_ring().embed(from.decorated_ring(), el)
    }
}

impl<R, S> CanonicalIsomorphismInfo<R> for S
    where R: RingDecorator, S: RingDecorator, S::DecoratedRing: CanonicalIsomorphismInfo<R::DecoratedRing>
{
    fn has_isomorphism(&self, from: &R) -> RingPropValue {
        self.decorated_ring().has_isomorphism(from.decorated_ring())
    }

    fn preimage(&self, from: &R, el: Self::El) -> <R as RingBase>::El {
        self.decorated_ring().preimage(from.decorated_ring(), el)
    }
}

impl<R, T> CanonicalIsomorphismInfo<StaticRing<T>> for R
    where T: RingEl, R: RingDecorator, R::DecoratedRing: CanonicalIsomorphismInfo<StaticRing<T>>
{
    fn has_isomorphism(&self, from: &StaticRing<T>) -> RingPropValue {
        self.decorated_ring().has_embedding(from)
    }

    fn preimage(&self, from: &StaticRing<T>, el: Self::El) -> T {
        self.decorated_ring().preimage(from, el)
    }
}

impl<R> CanonicalIsomorphismInfo<BigIntSOORing> for R
    where R: RingDecorator, R::DecoratedRing: CanonicalIsomorphismInfo<BigIntSOORing>
{
    fn has_isomorphism(&self, from: &BigIntSOORing) -> RingPropValue {
        self.decorated_ring().has_isomorphism(from)
    }

    fn preimage(&self, from: &BigIntSOORing, el: Self::El) -> BigIntSOO {
        self.decorated_ring().preimage(from, el)
    }
}

///
/// Creates an embedding between two rings that can be canonically
/// embedded into each other.
/// 
/// # Example
/// 
/// ```
/// # use feanor_la::prelude::*;
/// # use feanor_la::rational::primitive_rational::*;
/// let integers = i64::RING;
/// let rationals = r64::RING;
/// let f = embedding(integers, rationals);
/// assert_eq!(r64::new(3, 1), f(3));
/// ```
/// 
pub fn embedding<R, S>(from: R, to: S) -> StandardEmbedding<R, S> 
    where R: Ring, S: Ring + CanonicalEmbeddingInfo<R>
{
    assert!(to.has_embedding(&from).can_use());
    StandardEmbedding { base: from, ring: to }
}

///
/// Creates an isomorphism and its inverse between two canonically
/// isomorphic rings.
/// 
pub fn isomorphism<R, S>(from: R, to: S) -> (impl Fn(R::El) -> S::El, impl Fn(S::El) -> R::El)
    where R: Ring, S: CanonicalIsomorphismInfo<R>
{
    assert!(to.has_isomorphism(&from).can_use());
    let from_copy = from.clone();
    let to_copy = to.clone();
    (move |x| to_copy.embed(&from_copy, x), move |x| to.preimage(&from, x))
}

///
/// Creates the natural, not necessarily injective homomorphism
/// `Z -> R` that exists for every ring R.
/// 
pub fn z_hom<'a, R>(to: &'a R) -> impl 'a + Copy + Fn(i64) -> R::El 
    where R: Ring
{
    move |x| to.from_z(x)
}

///
/// The canonical embedding between two rings that can be
/// canonically embedded into each other, i.e. implement the
/// appropriate version of [`crate::embedding::CanonicalEmbeddingInfo`].
/// 
/// It is recommended to create instances by using the global
/// function [`crate::embedding::embedding()`].
/// 
#[derive(Clone, Copy)]
pub struct StandardEmbedding<R, S>
    where R: Ring, S: Ring + CanonicalEmbeddingInfo<R>
{
    base: R,
    ring: S
}

impl<R, S> PartialEq for StandardEmbedding<R, S>
    where R: Ring + PartialEq, S: Ring + CanonicalEmbeddingInfo<R> + PartialEq
{
    fn eq(&self, rhs: &Self) -> bool {
        assert_eq!(self.base, rhs.base);
        assert_eq!(self.ring, rhs.ring);
        true
    }
}

impl<R, S> FnOnce<(El<R>,)> for StandardEmbedding<R, S>
    where R: Ring, S: Ring + CanonicalEmbeddingInfo<R>
{
    type Output = El<S>;

    extern "rust-call" fn call_once(
        mut self, 
        (x, ): (R::El, )
    ) -> Self::Output {
        self.call_mut((x, ))
    }
}

impl<R, S> FnMut<(El<R>,)> for StandardEmbedding<R, S>
    where R: Ring, S: Ring + CanonicalEmbeddingInfo<R>
{
    extern "rust-call" fn call_mut(
        &mut self, 
        (x, ): (R::El, )
    ) -> Self::Output {
        self.call((x, ))
    }
}

impl<R, S> Fn<(El<R>,)> for StandardEmbedding<R, S>
    where R: Ring, S: Ring + CanonicalEmbeddingInfo<R>
{
    extern "rust-call" fn call(
        &self, 
        (x, ): (R::El, )
    ) -> Self::Output {
        self.ring.embed(&self.base, x)
    }
}

///
/// An embedding of two rings that is defined as the composition
/// of two given embeddings.
/// 
/// It is recommended to create instances using the global function
/// [`crate::embedding::compose()`].
/// 
#[derive(Clone, Copy)]
pub struct ComposedEmbedding<F, G, R, S, T>
    where R: Ring, S: Ring, T: Ring, F: Fn(R::El) -> S::El, G: Fn(S::El) -> T::El
{
    f: F,
    g: G,
    r: PhantomData<R>,
    s: PhantomData<S>,
    t: PhantomData<T>
}

impl<F, G, R, S, T> PartialEq for ComposedEmbedding<F, G, R, S, T>
    where R: Ring + PartialEq, S: Ring + PartialEq, T: Ring + PartialEq, F: Fn(R::El) -> S::El + PartialEq, G: Fn(S::El) -> T::El + PartialEq
{
    fn eq(&self, rhs: &Self) -> bool {
        assert_eq!(self.r, rhs.r);
        assert_eq!(self.s, rhs.s);
        assert_eq!(self.t, rhs.t);
        self.f == rhs.f && self.g == rhs.g
    }
}

impl<F, G, R, S, T> FnOnce<(R::El, )> for ComposedEmbedding<F, G, R, S, T>
    where R: Ring, S: Ring, T: Ring, F: Fn(R::El) -> S::El, G: Fn(S::El) -> T::El
{
    type Output = T::El;

    extern "rust-call" fn call_once(
        mut self, 
        (x, ): (R::El, )
    ) -> Self::Output {
        self.call_mut((x, ))
    }
}

impl<F, G, R, S, T> FnMut<(R::El, )> for ComposedEmbedding<F, G, R, S, T>
    where R: Ring, S: Ring, T: Ring, F: Fn(R::El) -> S::El, G: Fn(S::El) -> T::El
{
    extern "rust-call" fn call_mut(
        &mut self, 
        (x, ): (R::El, )
    ) -> Self::Output {
        self.call((x, ))
    }
}

impl<F, G, R, S, T> Fn<(R::El, )> for ComposedEmbedding<F, G, R, S, T>
    where R: Ring, S: Ring, T: Ring, F: Fn(R::El) -> S::El, G: Fn(S::El) -> T::El
{
    extern "rust-call" fn call(
        &self, 
        (x, ): (R::El, )
    ) -> Self::Output {
        (self.g)((self.f)(x))
    }
}

///
/// Composes two arbitrary embeddings of rings.
/// 
/// Due to the generality of this function, it is often
/// required to provide type arguments
/// 
/// There is no check that the given maps are injective, or
/// even ring homomorphisms. If wrong input is passed, the
/// resulting map might not be an embedding as well.
/// 
/// # Example
/// 
/// ```
/// # use feanor_la::prelude::*;
/// # use feanor_la::rational::primitive_rational::*;
/// # use feanor_la::poly::uni_var::*;
/// let integers = i64::RING;
/// let rationals = r64::RING;
/// let rational_polys = PolyRingImpl::adjoint(rationals, "X");
/// let f = rationals.embedding();
/// let g = rational_polys.embedding();
/// let h: ComposedEmbedding<_, _, StaticRing<i64>, StaticRing<r64>, PolyRingImpl<StaticRing<r64>>> = compose(g, f);
/// assert_eq!(rational_polys.from_z(3), h(3));
/// ```
/// 
pub fn compose<F, G, R, S, T>(g: G, f: F) -> ComposedEmbedding<F, G, R, S, T>
    where R: Ring, S: Ring, T: Ring, F: Fn(R::El) -> S::El, G: Fn(S::El) -> T::El
{
    ComposedEmbedding {
        f: f, g: g, r: PhantomData, s: PhantomData, t: PhantomData
    }
}
