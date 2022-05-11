use super::ring::*;
use super::primitive::*;
use super::bigint::*;
use super::ring_property::*;

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
///    `BigIntRing`, as they are required for the definition of the Integers.
/// 
/// New ring implentations are strongly encouraged to implement
/// ```text
///     impl CanonicalEmbeddingInfo<R> for R
/// ```
/// However, be careful not to violate the constraints of CanonicalEmbedding.
/// Strictly speaking, rings with nontrivial automorphisms are not canonically
/// isomorphic in the strict sense.
/// However, it might still be possible to choose compatible isomorphisms for
/// all ring implementations in the isomorphism class, and therefore implement
/// `CanonicalEmbedding` resp. `CanonicalIsomorphism` correctly.
///
/// The same holds for `CanonicalIsomorphismInfo`.
/// 
pub trait CanonicalEmbeddingInfo<R>: Ring
    where R: Ring
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
/// combinations. For more details, see the notes on `CanonicalEmbeddingInfo`.
/// 
pub trait CanonicalIsomorphismInfo<R>: CanonicalEmbeddingInfo<R>
    where R: Ring
{
    fn has_isomorphism(&self, from: &R) -> RingPropValue;
    fn preimage(&self, from: &R, el: Self::El) -> R::El;
}

impl<'b, R, T: RingEl> CanonicalEmbeddingInfo<StaticRing<T>> for &'b R
    where R: CanonicalEmbeddingInfo<StaticRing<T>>
{
    fn has_embedding(&self, from: &StaticRing<T>) -> RingPropValue {
        (**self).has_embedding(from)
    }

    fn embed(&self, from: &StaticRing<T>, el: T) -> Self::El {
        (**self).embed(from, el)
    }
}

impl<'b, R> CanonicalEmbeddingInfo<BigIntRing> for &'b R
    where R: CanonicalEmbeddingInfo<BigIntRing>
{
    fn has_embedding(&self, from: &BigIntRing) -> RingPropValue {
        (**self).has_embedding(from)
    }

    fn embed(&self, from: &BigIntRing, el: BigInt) -> Self::El {
        (**self).embed(from, el)
    }
}

impl<'a, 'b, R, S> CanonicalEmbeddingInfo<&'a R> for &'b S
    where R: Ring, S: CanonicalEmbeddingInfo<R>
{
    fn has_embedding(&self, from: &&'a R) -> RingPropValue {
        (**self).has_embedding(&**from)
    }

    fn embed(&self, from: &&'a R, el: R::El) -> Self::El {
        (**self).embed(&**from, el)
    }
}

impl<'a, 'b, R, S> CanonicalIsomorphismInfo<&'a R> for &'b S
    where R: CanonicalIsomorphismInfo<S>, S: CanonicalIsomorphismInfo<R>
{
    fn has_isomorphism(&self, from: &&'a R) -> RingPropValue {
        (**self).has_isomorphism(&**from)
    }

    fn preimage(&self, from: &&'a R, el: Self::El) -> R::El {
        (**self).preimage(&**from, el)
    }
}

impl<'b, R, T: RingEl> CanonicalIsomorphismInfo<StaticRing<T>> for &'b R
    where R: CanonicalIsomorphismInfo<StaticRing<T>>
{
    fn has_isomorphism(&self, from: &StaticRing<T>) -> RingPropValue {
        (**self).has_embedding(from)
    }

    fn preimage(&self, from: &StaticRing<T>, el: Self::El) -> T {
        (**self).preimage(from, el)
    }
}

impl<'b, R> CanonicalIsomorphismInfo<BigIntRing> for &'b R
    where R: CanonicalIsomorphismInfo<BigIntRing>
{
    fn has_isomorphism(&self, from: &BigIntRing) -> RingPropValue {
        (**self).has_isomorphism(from)
    }

    fn preimage(&self, from: &BigIntRing, el: Self::El) -> BigInt {
        (**self).preimage(from, el)
    }
}

pub fn embedding<R, S>(from: R, to: S) -> impl Fn(R::El) -> S::El 
    where R: Ring, S: CanonicalEmbeddingInfo<R>
{
    assert!(to.has_embedding(&from).can_use());
    move |x| to.embed(&from, x)
}

pub fn isomorphism<R, S>(from: R, to: S) -> (impl Fn(R::El) -> S::El, impl Fn(S::El) -> R::El)
    where R: Ring, S: CanonicalIsomorphismInfo<R>
{
    assert!(to.has_isomorphism(&from).can_use());
    let from_copy = from.clone();
    let to_copy = to.clone();
    (move |x| to_copy.embed(&from_copy, x), move |x| to.preimage(&from, x))
}

pub fn z_hom<'a, R>(to: &'a R) -> impl 'a + Fn(i64) -> R::El 
    where R: Ring
{
    move |x| to.from_z(x)
}

pub struct ComposedEmbedding<F, G, R, S, T>
    where R: Ring, S: Ring, T: Ring, F: Fn(R::El) -> S::El, G: Fn(S::El) -> T::El
{
    f: F,
    g: G,
    r: PhantomData<R>,
    s: PhantomData<S>,
    t: PhantomData<T>
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

pub fn compose<F, G, R, S, T>(g: G, f: F) -> ComposedEmbedding<F, G, R, S, T>
    where R: Ring, S: Ring, T: Ring, F: Fn(R::El) -> S::El, G: Fn(S::El) -> T::El
{
    ComposedEmbedding {
        f: f, g: g, r: PhantomData, s: PhantomData, t: PhantomData
    }
}
