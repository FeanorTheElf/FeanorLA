use super::ring::*;

use std::marker::PhantomData;

pub trait CanonicalEmbeddingInfo<R>: Ring
    where R: Ring
{
    fn has_embedding(&self, from: &R) -> RingPropValue;
    fn embed(&self, from: &R, el: R::El) -> Self::El;
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

pub fn embedding<R, S>(from: R, to: S) -> impl Fn(R::El) -> S::El 
    where R: Ring, S: CanonicalEmbeddingInfo<R>
{
    assert!(to.has_embedding(&from).can_use());
    move |x| to.embed(&from, x)
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