use super::super::prelude::*;

#[derive(Clone)]
pub struct ExtensionWrapper<R, S, F>
    where R: Ring, S: Ring, F: Fn(El<R>) -> El<S> + Clone
{
    base: R,
    ring: S,
    embedding: F
}

impl<R, S, F> std::fmt::Debug for ExtensionWrapper<R, S, F> 
    where R: Ring, S: Ring, F: Fn(El<R>) -> El<S> + Clone
{
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "Extension ring {:?} of {:?}", self.ring, self.base)
    }
}

impl<R, S, F> RingBase for ExtensionWrapper<R, S, F> 
    where R: Ring, S: Ring, F: Fn(El<R>) -> El<S> + Clone
{
    type El = El<S>;

    fn add_ref(&self, lhs: Self::El, rhs: &Self::El) -> Self::El { self.ring.add_ref(lhs, rhs) }
    fn mul_ref(&self, lhs: &Self::El, rhs: &Self::El) -> Self::El { self.ring.mul_ref(lhs, rhs) }
    fn add_assign(&self, lhs: &mut Self::El, rhs: Self::El) { self.ring.add_assign(lhs, rhs) }
    fn mul_assign(&self, lhs: &mut Self::El, rhs: Self::El) { self.ring.mul_assign(lhs, rhs) }
    fn neg(&self, val: Self::El) -> Self::El { self.ring.neg(val) }
    fn zero(&self) -> Self::El { self.ring.zero() }
    fn one(&self) -> Self::El { self.ring.one() }
    fn eq(&self, lhs: &Self::El, rhs: &Self::El) -> bool { self.ring.eq(lhs, rhs) }
    fn unspecified_element(&self) -> Self::El { self.ring.unspecified_element() }
    fn sub_ref_fst(&self, lhs: &Self::El, rhs: Self::El) -> Self::El { self.ring.sub_ref_fst(lhs, rhs) }
    fn sub_ref_snd(&self, lhs: Self::El, rhs: &Self::El) -> Self::El { self.ring.sub_ref_snd(lhs, rhs) }
    fn add(&self, lhs: Self::El, rhs: Self::El) -> Self::El { self.ring.add(lhs, rhs) }
    fn mul(&self, lhs: Self::El, rhs: Self::El) -> Self::El { self.ring.mul(lhs, rhs) }
    fn sub(&self, lhs: Self::El, rhs: Self::El) -> Self::El { self.ring.sub(lhs, rhs) }
    fn pow(&self, basis: &Self::El, exp: u32) -> Self::El { self.ring.pow(basis, exp) }
    fn pow_big(&self, basis: &Self::El, exp: &BigInt) -> Self::El { self.ring.pow_big(basis, exp) }
    fn from_z(&self, x: i64) -> Self::El { self.ring.from_z(x) }
    fn from_z_big(&self, x: &BigInt) -> Self::El { self.ring.from_z_big(x) }
    fn is_zero(&self, val: &Self::El) -> bool { self.ring.is_zero(val) }
    fn is_one(&self, val: &Self::El) -> bool { self.ring.is_one(val) }
    fn is_neg_one(&self, val: &Self::El) -> bool { self.ring.is_neg_one(val) }
    fn is_integral(&self) -> RingPropValue { self.ring.is_integral() }
    fn characteristic(&self) -> BigInt { self.ring.characteristic() }
    fn is_field(&self) -> RingPropValue { self.ring.is_field() }
    fn is_noetherian(&self) -> bool { self.ring.is_noetherian() }
    fn div(&self, lhs: Self::El, rhs: &Self::El) -> Self::El { self.ring.div(lhs, rhs) }
    fn format(&self, el: &Self::El, f: &mut std::fmt::Formatter, in_prod: bool) -> std::fmt::Result { self.ring.format(el, f, in_prod) }
    fn format_in_brackets(&self, el: &Self::El, f: &mut std::fmt::Formatter) -> std::fmt::Result { self.ring.format_in_brackets(el, f) }
}

impl<R, S, F> CanonicalEmbeddingInfo<ExtensionWrapper<R, S, F>> for ExtensionWrapper<R, S, F> 
    where R: Ring, S: Ring, F: Fn(El<R>) -> El<S> + Clone
{
    fn has_embedding(&self, from: &ExtensionWrapper<R, S, F>) -> RingPropValue {
        self.ring.has_embedding(&from.ring)
    }

    fn embed(&self, from: &ExtensionWrapper<R, S, F>, el: El<ExtensionWrapper<R, S, F>>) -> El<Self> {
        self.ring.embed(&from.ring, el)
    }
}