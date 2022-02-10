use super::ring::*;
use super::bigint::*;
use super::primitive::*;

use std::marker::PhantomData;

pub trait CanonicalEmbeddingInfo<R>: Ring
    where R: Ring
{
    type Embedding: Fn(R::El) -> Self::El;

    fn has_embedding(&self, from: &R) -> RingPropValue;
    fn embedding(self, from: R) -> Self::Embedding;
}

pub trait CanonicalZEmbeddingInfo<R>: Ring
    where R: Ring
{
    type Embedding: Fn(R::El) -> Self::El;

    fn z_embedding(self, from: R) -> Self::Embedding;
}

pub struct BigIntEmbedding<R>
    where R: Ring
{
    target: R
}

impl<R> FnOnce<(BigInt, )> for BigIntEmbedding<R> 
    where R: Ring
{
    type Output = R::El;

    extern "rust-call" fn call_once(
        mut self, 
        (x, ): (BigInt, )
    ) -> Self::Output {
        self.call_mut((x, ))
    }
}

impl<R> FnMut<(BigInt, )> for BigIntEmbedding<R> 
    where R: Ring
{
    extern "rust-call" fn call_mut(
        &mut self, 
        (x, ): (BigInt, )
    ) -> Self::Output {
        self.call((x, ))
    }
}

impl<R> Fn<(BigInt, )> for BigIntEmbedding<R> 
    where R: Ring
{
    extern "rust-call" fn call(
        &self, 
        (x, ): (BigInt, )
    ) -> Self::Output {
        self.target.from_z_big(x)
    }
}

impl<R> CanonicalZEmbeddingInfo<BigIntRing> for R
    where R: Ring
{
    type Embedding = BigIntEmbedding<R>;

    fn z_embedding(self, _from: BigIntRing) -> Self::Embedding {
        BigIntEmbedding {
            target: self
        }
    }
}

pub struct IntEmbedding<R>
    where R: Ring
{
    target: R
}

impl<R> FnOnce<(i64, )> for IntEmbedding<R> 
    where R: Ring
{
    type Output = R::El;

    extern "rust-call" fn call_once(
        mut self, 
        (x, ): (i64, )
    ) -> Self::Output {
        self.call_mut((x, ))
    }
}

impl<R> FnMut<(i64, )> for IntEmbedding<R> 
    where R: Ring
{
    extern "rust-call" fn call_mut(
        &mut self, 
        (x, ): (i64, )
    ) -> Self::Output {
        self.call((x, ))
    }
}

impl<R> Fn<(i64, )> for IntEmbedding<R> 
    where R: Ring
{
    extern "rust-call" fn call(
        &self, 
        (x, ): (i64, )
    ) -> Self::Output {
        self.target.from_z(x)
    }
}

impl<R> CanonicalZEmbeddingInfo<StaticRing<i64>> for R
    where R: Ring
{
    type Embedding = IntEmbedding<R>;

    fn z_embedding(self, _from: StaticRing<i64>) -> Self::Embedding {
        IntEmbedding {
            target: self
        }
    }
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