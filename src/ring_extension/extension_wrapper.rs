use super::super::prelude::*;
use super::super::wrapper::*;

#[derive(Clone)]
pub struct ExtensionWrapper<R, S, F>
    where R: Ring + PartialEq, S: Ring + PartialEq, F: Fn(El<R>) -> El<S> + Clone + PartialEq
{
    base: R,
    ring: S,
    embedding: F
}

impl<R, S, F> std::fmt::Debug for ExtensionWrapper<R, S, F> 
    where R: Ring + PartialEq, S: Ring + PartialEq, F: Fn(El<R>) -> El<S> + Clone + PartialEq
{
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "Extension ring {:?} of {:?}", self.ring, self.base)
    }
}

impl<R, S, F> RingBase for ExtensionWrapper<R, S, F> 
    where R: Ring + PartialEq, S: Ring + PartialEq, F: Fn(El<R>) -> El<S> + Clone + PartialEq
{
    type El = El<S>;

    fn add_ref(&self, lhs: Self::El, rhs: &Self::El) -> Self::El { self.ring.add_ref(lhs, rhs) }
    fn mul_ref(&self, lhs: &Self::El, rhs: &Self::El) -> Self::El { self.ring.mul_ref(lhs, rhs) }
    fn add_assign(&self, lhs: &mut Self::El, rhs: Self::El) { self.ring.add_assign(lhs, rhs) }
    fn mul_assign(&self, lhs: &mut Self::El, rhs: Self::El) { self.ring.mul_assign(lhs, rhs) }
    fn neg(&self, val: Self::El) -> Self::El { self.ring.neg(val) }
    fn zero(&self) -> Self::El { self.ring.zero() }
    fn one(&self) -> Self::El { self.ring.one() }
    fn is_eq(&self, lhs: &Self::El, rhs: &Self::El) -> bool { self.ring.is_eq(lhs, rhs) }
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
    where R: Ring + PartialEq, S: Ring + PartialEq, F: Fn(El<R>) -> El<S> + Clone + PartialEq
{
    fn has_embedding(&self, from: &ExtensionWrapper<R, S, F>) -> RingPropValue {
        RingPropValue::True & (self.ring == from.ring && self.base == from.base && self.embedding == from.embedding)
    }

    fn embed(&self, from: &ExtensionWrapper<R, S, F>, el: El<ExtensionWrapper<R, S, F>>) -> El<Self> {
        assert!(self.has_embedding(from).can_use());
        self.ring.embed(&from.ring, el)
    }
}

impl<R, S, F> CanonicalIsomorphismInfo<ExtensionWrapper<R, S, F>> for ExtensionWrapper<R, S, F> 
    where R: Ring + PartialEq, S: Ring + PartialEq, F: Fn(El<R>) -> El<S> + Clone + PartialEq
{
    fn has_isomorphism(&self, from: &ExtensionWrapper<R, S, F>) -> RingPropValue {
        RingPropValue::True & (self.ring == from.ring && self.base == from.base && self.embedding == from.embedding)
    }

    fn preimage(&self, from: &ExtensionWrapper<R, S, F>, el: El<ExtensionWrapper<R, S, F>>) -> El<Self> {
        assert!(self.has_isomorphism(from).can_use());
        self.ring.preimage(&from.ring, el)
    }
}

impl<R, S, F> CanonicalEmbeddingInfo<R> for ExtensionWrapper<R, S, F> 
    where R: Ring + PartialEq, S: Ring + PartialEq, F: Fn(El<R>) -> El<S> + Clone + PartialEq
{
    fn has_embedding(&self, from: &R) -> RingPropValue {
        RingPropValue::True & (self.base == *from)
    }

    fn embed(&self, from: &R, el: El<R>) -> El<Self> {
        assert!(self.has_embedding(from).can_use());
        (self.embedding)(el)
    }
}

impl<R, S, F> ExtensionWrapper<R, S, F> 
    where R: Ring + PartialEq, S: Ring + PartialEq, F: Fn(El<R>) -> El<S> + Clone + PartialEq
{
    pub fn new(base: R, ring: S, embedding: F) -> Self {
        ExtensionWrapper {
            base, ring, embedding
        }
    }

    pub fn wrapped_ring(&self) -> &S {
        &self.ring
    }

    pub fn into_wrapped_ring(&self) -> &S {
        &self.ring
    }

    pub fn wrapping_embedding<'a>(&'a self) -> StandardEmbedding<&'a S, &'a S> {
        embedding(&self.ring, &self.ring)
    }
}

impl<R, S, F> RingElWrapper<ExtensionWrapper<R, S, F>>
    where R: Ring + PartialEq, S: Ring + PartialEq, F: Fn(El<R>) -> El<S> + Clone + PartialEq
{
    pub fn unwrap(self) -> RingElWrapper<S> {
        let (el, ring) = self.destruct();
        return ring.into_wrapped_ring().bind_by_value(el);
    }
}

#[cfg(test)]
use super::super::fq::zn_big::*;
#[cfg(test)]
use super::simple_extension::*;

#[allow(non_snake_case)]
#[test]
fn test_no_embedding() {
    let R1 = Zn::new(BigInt::from(5));
    let R2 = Zn::new(BigInt::from(7));
    let S1 = ExtensionWrapper::new(R1.clone(), R1.clone(), embedding(&R1, &R1));
    let S2 = ExtensionWrapper::new(R2.clone(), R2.clone(), embedding(&R1, &R1));
    assert!(S1.has_embedding(&S1).can_use());
    assert!(!S1.has_embedding(&S2).can_use());

    type GaussianIntegers = SimpleRingExtension<StaticRing<i64>, ConstVector2<-1, 0>, VectorArray<i64, 2>>;

    #[derive(Debug, Clone, Copy, PartialEq)]
    enum GaussianIntegersAutomorphism {
        Identity, Conjugation
    }

    impl FnOnce<(El<GaussianIntegers>, )> for GaussianIntegersAutomorphism {

        type Output = El<GaussianIntegers>;

        extern "rust-call" fn call_once(
            mut self, 
            (x, ): (El<GaussianIntegers>, )
        ) -> Self::Output {
            self.call_mut((x, ))
        }
    }

    impl FnMut<(El<GaussianIntegers>, )> for GaussianIntegersAutomorphism {

        extern "rust-call" fn call_mut(
            &mut self, 
            (x, ): (El<GaussianIntegers>, )
        ) -> Self::Output {
            self.call((x, ))
        }
    }

    impl Fn<(El<GaussianIntegers>,)> for GaussianIntegersAutomorphism {

        extern "rust-call" fn call(
            &self, 
            (x, ): (El<GaussianIntegers>, )
        ) -> Self::Output {
            match self {
                GaussianIntegersAutomorphism::Identity => x,
                GaussianIntegersAutomorphism::Conjugation => Vector::from_array([x[0], -x[1]])
            }
        }
    }

    let R = GaussianIntegers::new(i64::RING, Vector::new(ConstVector2::INSTANCE), "i");
    let S = R.clone();
    let S1 = ExtensionWrapper::new(R.clone(), S.clone(), GaussianIntegersAutomorphism::Identity);
    let S2 = ExtensionWrapper::new(R.clone(), S.clone(), GaussianIntegersAutomorphism::Conjugation);
    assert!(S1.has_embedding(&S1).can_use());
    assert!(!S1.has_embedding(&S2).can_use());
}