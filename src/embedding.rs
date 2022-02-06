use super::ring::*;
use super::bigint::*;
use super::primitive::*;

pub trait CanonicalEmbeddingInfo<R>: Ring
    where R: Ring
{
    type Embedding: Fn(R::El) -> Self::El;

    fn has_embedding(&self, from: &R) -> RingPropValue;
    fn embedding(self, from: R) -> Self::Embedding;
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

impl<R> CanonicalEmbeddingInfo<BigIntRing> for R
    where R: Ring
{
    type Embedding = BigIntEmbedding<R>;

    fn has_embedding(&self, _from: &BigIntRing) -> RingPropValue {
        RingPropValue::True
    }

    fn embedding(self, _from: BigIntRing) -> Self::Embedding {
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

impl<R> CanonicalEmbeddingInfo<StaticRing<i64>> for R
    where R: Ring
{
    type Embedding = IntEmbedding<R>;

    fn has_embedding(&self, _from: &StaticRing<i64>) -> RingPropValue {
        RingPropValue::True
    }

    fn embedding(self, _from: StaticRing<i64>) -> Self::Embedding {
        IntEmbedding {
            target: self
        }
    }
}