use super::super::prelude::*;
use super::ops::*;
use super::factoring;
use super::super::la::vec::*;
use super::super::fq::*;
use super::super::fq::zn_big::*;
use super::*;

use vector_map::VecMap;

#[derive(Debug)]
pub struct PolyRingImpl<R>
    where R: Ring
{
    base_ring: R,
    var_name: &'static str
}

impl<R> PolyRingImpl<R>
    where R: CanonicalIsomorphismInfo<R>
{
    pub const fn adjoint(base_ring: R, var_name: &'static str) -> Self {
        PolyRingImpl {
            base_ring, var_name
        }
    }

    pub fn evaluation_hom<'a>(&'a self, x: R::El) -> impl 'a + Fn(El<Self>) -> R::El {
        move |poly| self.evaluate_at(&poly, x.clone(), self.base_ring())
    }

    pub fn lift_hom<F, S>(&self, hom: F) -> impl Fn(El<PolyRingImpl<S>>) -> El<Self>
        where S: Ring, F: Fn(S::El) -> R::El
    {
        move |poly| {
            Vector::new(poly.raw_data().into_iter().map(&hom).collect())
        }
    }

    ///
    /// Performs polynomial division, so computes q and r such that lhs = q * rhs + r
    /// with deg(r) < deg(lhs). q is returned and r is contained in lhs after the function 
    /// returned. 
    /// Note that this function has to compute divisions by the leading coefficient of rhs, 
    /// which must be given via a function object. Special cases, e.g when the base ring is a
    /// field and this is just field division can be accessed via the corresponding function
    /// (in this case, `PolyRing::div`). Errors from this function are forwarded, and this is
    /// the only case in which this function returns Err(()).
    /// 
    pub fn poly_division<F>(&self, lhs: &mut El<Self>, rhs: &El<Self>, div_lc: F) -> Result<El<Self>, ()>
        where F: FnMut(&R::El) -> Result<R::El, ()>
    {
        poly_division(&self.base_ring, lhs.as_mut(), rhs.as_ref(), div_lc)
    }
}

impl<R> Clone for PolyRingImpl<R>
    where R: Ring
{
    fn clone(&self) -> Self {
        PolyRingImpl {
            var_name: self.var_name,
            base_ring: self.base_ring.clone()
        }
    }
}

impl<R> Copy for PolyRingImpl<R>
    where R: Ring + Copy
{}

impl<R> PolyRing for PolyRingImpl<R>
    where R: Ring
{
    fn lc(&self, x: &El<Self>) -> Option<El<Self::BaseRing>> {
        if let Some(d) = self.deg(x) {
            Some(x[d].clone())
        } else {
            None
        }
    }

    fn unknown_name(&self) -> &str {
        self.var_name
    }

    fn deg(&self, x: &El<Self>) -> Option<usize> {
        let mut it = x.iter().enumerate().filter(|(_, x)| !self.base_ring().is_zero(x)).map(|(i, _)| i);
        <_ as Iterator>::last(&mut it)
    }

    fn coefficient_at(&self, x: &El<Self>, i: usize) -> El<Self::BaseRing> {
        x[i].clone()
    }

    fn scale(&self, x: &mut El<Self>, coeff: &El<Self::BaseRing>) {
        for i in 0..x.len() {
            *x.at_mut(i) = self.base_ring().mul_ref(x.at(i), &coeff);
        }
    }

    fn evaluate_at<S: Ring + CanonicalEmbeddingInfo<R>>(&self, f: &El<Self>, x: El<S>, ring: &S) -> El<S> {
        poly_evaluate(&self.base_ring, x, f.as_ref(), ring)
    }

    fn derive(&self, el: El<Self>) -> El<Self> {
        let mut result = el.into_owned();
        poly_formal_derivative(&self.base_ring, result.as_mut());
        return result;
    }
    
    fn unknown(&self) -> El<Self> {
        let mut result = Vec::with_capacity(2);
        result.push(self.base_ring.zero());
        result.push(self.base_ring.one());
        return Vector::new(result);
    }
}

impl<R> CanonicalEmbeddingInfo<R> for PolyRingImpl<R> 
    where R: CanonicalIsomorphismInfo<R>
{
    fn has_embedding(&self, from: &R) -> RingPropValue {
        self.base_ring().has_embedding(from)
    }

    fn embed(&self, from: &R, el: R::El) -> Self::El {
        assert!(self.has_embedding(from).can_use());
        self.from(self.base_ring().embed(from, el))
    }
}

impl<R> CanonicalEmbeddingInfo<R> for PolyRingImpl<&R> 
    where R: CanonicalIsomorphismInfo<R>
{

    fn has_embedding(&self, from: &R) -> RingPropValue {
        self.base_ring().has_embedding(&from)
    }

    fn embed(&self, from: &R, el: R::El) -> Self::El {
        assert!(self.has_embedding(from).can_use());
        self.from(self.base_ring().embed(&from, el))
    }
}

impl<R> CanonicalEmbeddingInfo<&R> for PolyRingImpl<R> 
    where R: CanonicalIsomorphismInfo<R>
{

    fn has_embedding(&self, from: &&R) -> RingPropValue {
        self.base_ring().has_embedding(*from)
    }

    fn embed(&self, from: &&R, el: R::El) -> Self::El {
        assert!(self.has_embedding(from).can_use());
        self.from(self.base_ring().embed(*from, el))
    }
}

impl<R> CanonicalEmbeddingInfo<PolyRingImpl<R>> for PolyRingImpl<R>
    where R: CanonicalIsomorphismInfo<R>
{
    fn has_embedding(&self, from: &PolyRingImpl<R>) -> RingPropValue {
        self.base_ring().has_embedding(from.base_ring())
    }

    fn embed(&self, from: &PolyRingImpl<R>, el: Self::El) -> Self::El {
        assert!(self.has_embedding(from).can_use());
        Vector::new(el.raw_data().into_iter().map(|x| self.base_ring().embed(from.base_ring(), x)).collect())
    }
}

impl<R> CanonicalIsomorphismInfo<PolyRingImpl<R>> for PolyRingImpl<R>
    where R: CanonicalIsomorphismInfo<R>
{
    fn has_isomorphism(&self, from: &PolyRingImpl<R>) -> RingPropValue {
        self.base_ring().has_isomorphism(from.base_ring())
    }

    fn preimage(&self, from: &PolyRingImpl<R>, el: Self::El) -> Self::El {
        assert!(self.has_isomorphism(from).can_use());
        Vector::new(el.raw_data().into_iter().map(|x| self.base_ring().preimage(from.base_ring(), x)).collect())
    }
}

impl<R> RingBase for PolyRingImpl<R>
    where R: CanonicalIsomorphismInfo<R>
{
    type El = Vector<Vec<R::El>, R::El>;

    fn add_ref(&self, lhs: Self::El, rhs: &Self::El) -> Self::El {
        poly_add(&self.base_ring, lhs, rhs.as_ref())
    }

    fn mul_ref(&self, lhs: &Self::El, rhs: &Self::El) -> Self::El {
        let result_len = poly_degree(self.base_ring(), lhs.as_ref()).unwrap_or(0) + poly_degree(self.base_ring(), rhs.as_ref()).unwrap_or(0) + 1;
        let mut result = Vector::new((0..result_len).map(|_| self.base_ring().zero()).collect::<Vec<_>>());
        poly_mul(&self.base_ring, lhs.as_ref(), rhs.as_ref(), &mut result);
        return result;
    }

    fn neg(&self, mut val: Self::El) -> Self::El {
        for i in 0..val.len() {
            let value = std::mem::replace(val.at_mut(i), self.base_ring.unspecified_element());
            *val.at_mut(i) = self.base_ring.neg(value);
        }
        return val;
    }

    fn zero(&self) -> Self::El {
        let result = Vec::new();
        return Vector::new(result);
    }

    fn one(&self) -> Self::El {
        let result = self.from(self.base_ring.one());
        return result;
    }

    fn is_eq(&self, lhs: &Self::El, rhs: &Self::El) -> bool {
        let (shorter, longer) = if lhs.len() <= rhs.len() {
            (lhs, rhs)
        } else {
            (rhs, lhs)
        };
        for i in 0..shorter.len() {
            if !self.base_ring.is_eq(&shorter[i], &longer[i]) {
                return false;
            }
        }
        for i in shorter.len()..longer.len() {
            if !self.base_ring.is_zero(&longer[i]) {
                return false;
            }
        }
        return true;
    }

    fn is_zero(&self, val: &Self::El) -> bool {
        val.iter().all(|x| self.base_ring.is_zero(x))
    }

    fn is_integral(&self) -> RingPropValue {
        self.base_ring.is_integral()
    }

    fn characteristic(&self) -> StdInt {
        self.base_ring().characteristic()
    }
    
    fn is_noetherian(&self) -> bool {
        self.base_ring.is_noetherian()
    }

    fn is_field(&self) -> RingPropValue {
        RingPropValue::False
    }

    fn from_z(&self, x: i64) -> Self::El {
        self.from(self.base_ring().from_z(x))
    }

    fn from_z_big(&self, x: &StdInt) -> Self::El {
        self.from(self.base_ring().from_z_big(x))
    }
    
    fn div(&self, _lhs: Self::El, _rhs: &Self::El) -> Self::El {
        panic!("Not a field!")
    }

    fn format(&self, el: &El<Self>, f: &mut std::fmt::Formatter, in_prod: bool) -> std::fmt::Result {
        if in_prod {
            self.format_in_brackets(el, f)
        } else {
            poly_format(&self.base_ring, el.as_ref(), f, self.var_name)
        }
    }
}

impl<R: Ring> RingExtension for PolyRingImpl<R> {
    
    type BaseRing = R;
    type Embedding = StandardEmbedding<R, Self>;

    fn is_extension(&self) -> RingPropValue {
        RingPropValue::True
    }

    fn base_ring(&self) -> &R {
        &self.base_ring
    }

    fn embedding(&self) -> Self::Embedding {
        embedding(self.base_ring().clone(), self.clone())
    }

    fn from(&self, el: El<R>) -> El<Self> {
        Vector::new(vec![el])
    }
}

impl<R> DivisibilityInfoRing for PolyRingImpl<R> 
    where R: Ring
{
    default fn is_divisibility_computable(&self) -> RingPropValue {
        self.base_ring.is_field()
    }

    default fn quotient(&self, lhs: &Self::El, rhs: &Self::El) -> Option<Self::El> {
        assert!(self.is_divisibility_computable().can_use());
        assert!(!self.is_zero(rhs));
        let lc = self.lc(rhs).unwrap();
        let lc_inv = self.base_ring().div(self.base_ring().one(), &lc);
        let mut p = lhs.clone();
        let result = self.poly_division(
            &mut p, 
            &rhs, 
            |x| Ok(self.base_ring().mul_ref(&x, &lc_inv))
        ).ok()?;
        if self.is_zero(&p) {
            return Some(result);
        } else {
            return None;
        }
    }

    default fn is_unit(&self, el: &Self::El) -> bool {
        assert!(self.is_divisibility_computable().can_use());
        self.deg(el) == Some(0)
    }
}

impl<R> DivisibilityInfoRing for PolyRingImpl<R> 
    where R: DivisibilityInfoRing
{
    fn is_divisibility_computable(&self) -> RingPropValue {
        self.base_ring.is_divisibility_computable()
    }

    fn quotient(&self, lhs: &Self::El, rhs: &Self::El) -> Option<Self::El> {
        assert!(self.is_divisibility_computable().can_use());
        assert!(!self.is_zero(rhs));
        let lc = self.lc(rhs).unwrap();
        let mut p = lhs.clone();
        let result = self.poly_division(
            &mut p, 
            &rhs, 
            |x| self.base_ring.quotient(x, &lc).ok_or(())
        ).ok()?;
        if self.is_zero(&p) {
            return Some(result);
        } else {
            return None;
        }
    }

    fn is_unit(&self, el: &Self::El) -> bool {
        assert!(self.is_divisibility_computable().can_use());
        self.deg(el).map(|d| d == 0).unwrap_or(false) && self.base_ring.is_unit(&el[0])
    }
}

impl<R> EuclideanInfoRing for PolyRingImpl<R> 
    where R: Ring + CanonicalIsomorphismInfo<R>
{
    fn is_euclidean(&self) -> RingPropValue {
        if self.base_ring.is_field().can_use() {
            return RingPropValue::True;
        } else {
            return RingPropValue::Unknown;
        }
    }

    fn euclidean_deg(&self, el: Self::El) -> StdInt {
        assert!(self.is_euclidean().can_use());
        self.deg(&el).map(|x| x + 1)
            .map(|x| StdInt::RING.from_z(x as i64))
            .unwrap_or(StdInt::zero())
    }

    fn euclidean_div_rem(&self, mut lhs: Self::El, rhs: &Self::El) -> (Self::El, Self::El) {
        assert!(self.is_euclidean().can_use());
        assert!(self.base_ring().is_field().can_use());
        assert!(!self.is_zero(&rhs));
        let rhs_lc = self.lc(&rhs).unwrap();
        let rhs_lc_inv = self.base_ring.div(self.base_ring.one(), &rhs_lc);
        let q = self.poly_division(
            &mut lhs, 
            &rhs, 
            |x| Ok(self.base_ring.mul_ref(x, &rhs_lc_inv))
        ).unwrap();
        return (q, lhs);
    }
}

impl<R> PartialEq for PolyRingImpl<R> 
    where R: Ring + PartialEq
{
    fn eq(&self, rhs: &Self) -> bool {
        self.base_ring() == rhs.base_ring()
    }
}

impl<R> UfdInfoRing for PolyRingImpl<R>
    where R: FiniteRing + DivisibilityInfoRing
{
    fn is_ufd(&self) -> RingPropValue {
        if self.base_ring.is_field().can_use() && self.base_ring().characteristic().is_odd() {
            return RingPropValue::True;
        } else {
            return RingPropValue::Unknown;
        }
    }

    default fn is_prime(&self, el: &Self::El) -> bool {
        assert!(self.is_ufd().can_use());
        return factoring::is_prime(self, el);
    }

    default fn calc_factor(&self, el: &Self::El) -> Option<Self::El> {
        let factorization = self.factor(el.clone());
        for (factor, _) in factorization.into_iter() {
            if self.deg(el) == self.deg(factor.val()) {
                return None;
            } else if !self.is_unit(factor.val()) {
                return Some(factor.into_val());
            }
        }
        return None;
    }

    default fn factor<'a>(&'a self, el: Self::El) -> VecMap<RingElWrapper<&'a Self>, usize> {
        assert!(self.is_ufd().can_use());
        assert!(!self.is_zero(&el));
        return factoring::factor_complete(self, el);
    }
}

impl<R> UfdInfoRing for PolyRingImpl<R>
    where R: IntegerQuotientRing + DivisibilityInfoRing
{
    fn is_prime(&self, el: &Self::El) -> bool {
        return factoring::is_prime(self, el);
    }

    ///
    /// If the base ring is a field, this is just the normal factor routine. However,
    /// this also tries to factor the polynomial if the polynomial ring is not a UFD,
    /// but a polynomial ring over `Z/p^eZ`. In this case, there is not a unique factorization
    /// anymore, but in some cases, there is a unique "canonical" factorization based on
    /// Hensel's lemma. However, in such a situation, the function will fail on some inputs.
    /// 
    fn factor<'a>(&'a self, el: Self::El) -> VecMap<RingElWrapper<&'a Self>, usize> {
        assert!(!self.is_zero(&el));
        let base_ring_prime_power = self.base_ring().characteristic().is_prime_power();
        if self.is_ufd().can_use() {
            return factoring::factor_complete(self, el);
        } else if let Some((p, _)) = base_ring_prime_power {
            if p.is_odd() {
                let integer_ring = self.base_ring().lifting_ring();
                let prime_field = Zn::new(&integer_ring, integer_ring.from_z_big(&p));
                return factoring::factor_lifted(self, el, &PolyRingImpl::adjoint(prime_field, "X")).unwrap();
            } else {
                unimplemented!()
            }
        } else {
            unimplemented!()
        }
    }
}

#[cfg(test)]
use super::super::fq::fq_small::*;
#[cfg(test)]
use super::super::integer::bigint::*;
#[cfg(test)]
use test::Bencher;

#[test]
fn test_poly_arithmetic() {
    let ring = WrappingRing::new(PolyRingImpl::adjoint(i32::RING, "X"));
    let x = ring.unknown();

    let x2_2x_1 = (&x * &x) + (&x + &x) + 1;
    let x_1_square = (&x + 1) * (&x + 1);

    assert_eq!(x2_2x_1, x_1_square);
    assert!(x2_2x_1 != x);
    assert_eq!(ring.zero(), x2_2x_1 - x_1_square);
}

#[test]
fn test_format() {
    let ring = WrappingRing::new(PolyRingImpl::adjoint(i32::RING, "X"));
    let x = ring.unknown();

    let poly = &x * &x * &x + &x * &x * 2 - 1;
    assert_eq!("-1 + 2 * X^2 + X^3", format!("{}", poly));
}

#[test]
fn test_poly_div() {
    let ring = WrappingRing::new(PolyRingImpl::adjoint(i32::RING, "X"));
    let x = ring.unknown();

    let mut p = &x * &x * &x + &x * &x + &x + 1;
    let q = &x + 1;
    let expected = &x * &x + 1;
    let result = ring.wrap(ring.wrapped_ring().poly_division(p.val_mut(), q.val(), |x| Ok(*x)).unwrap());
    assert_eq!(ring.zero(), p);
    assert_eq!(expected, result);
}

#[test]
fn test_quotient() {
    let ring = WrappingRing::new(PolyRingImpl::adjoint(i32::RING, "X"));
    let x = ring.unknown();

    let p = ring.one();
    let q = &x + 1;
    assert_eq!(None, ring.quotient(&p, &q));
}

#[test]
fn test_poly_degree() {
    let ring = WrappingRing::new(PolyRingImpl::adjoint(i32::RING, "X"));
    let x = ring.unknown();

    let p = &x * &x * &x + 4;
    assert_eq!(Some(3), p.deg());
}

#[test]
fn test_factor_zn() {
    let coeff_ring = Zn::new(BigInt::RING, BigInt::RING.from_z(3));
    let ring = WrappingRing::new(PolyRingImpl::adjoint(&coeff_ring, "X"));
    let x = ring.unknown();

    let p = x.clone().pow(9) - &x;
    let mut expected = VecMap::new();
    expected.insert(x.clone(), 1);
    expected.insert(&x + 1, 1);
    expected.insert(&x + 2, 1);
    expected.insert(&x * &x + &x + 2, 1);
    expected.insert(&x * &x + &x * 2 + 2, 1);
    expected.insert(&x * &x + 1, 1);
    let factorization = p.factor();
    assert_eq!(expected, factorization);
}

#[test]
fn test_factor_fq() {
    let coeff_ring = F49.clone();
    let ring = PolyRingImpl::adjoint(coeff_ring, "X");
    let ring = WrappingRing::new(&ring);
    let x = ring.unknown();

    let f = x.pow(2) + &x * 6 + 3;
    let factor = <_ as Iterator>::next(&mut ring.factor(f).into_iter()).unwrap().0.into_val();
    let coeff_ring_gen = -factor.coefficient_at(0) / factor.coefficient_at(1);
    let a = ring.embedding()(coeff_ring_gen.clone_ring()); // up to field automorphism, this is the generator picked by sage

    let p = x.clone().pow(14) - &x;
    let mut expected = VecMap::new();
    expected.insert(x.clone(), 1);
    expected.insert(&x + 6, 1);
    expected.insert(x.pow(6) + (&a * 3 + 6) * x.pow(5) + x.pow(4) * 2 + (&a * 3 + 5) * x.pow(3) + x.pow(2) * 2 + (&a * 3 + 6) * &x + 1, 1);
    expected.insert(x.pow(6) + (&a * 4 + 2) * x.pow(5) + x.pow(4) * 2 + (&a * 4 + 1) * x.pow(3) + x.pow(2) * 2 + (&a * 4 + 2) * &x + 1, 1);
    let factorization = (**ring.wrapped_ring()).factor(p.into_val());
    assert_eq!(expected, factorization);
}

#[test]
fn test_is_prime() {
    let coeff_ring = Zn::new(BigInt::RING, BigInt::RING.from_z(3));
    let ring = WrappingRing::new(PolyRingImpl::adjoint(&coeff_ring, "X"));
    let x = ring.unknown();

    let p = &x + 1;
    let q = &x * &x * &x + &x * 2 + 2;
    let a = &x * &x + 2;
    let b = &x * &x + &x * 2 + 1;
    assert_eq!(true, ring.is_prime(&p));
    assert_eq!(true, ring.is_prime(&q));
    assert_eq!(false, ring.is_prime(&a));
    assert_eq!(false, ring.is_prime(&b));
}

#[test]
fn test_evaluate() {
    let poly_ring = WrappingRing::new(PolyRingImpl::adjoint(i32::RING, "X"));
    let base_ring = poly_ring.base_ring();
    let i = base_ring.wrapping_embedding();
    let x = poly_ring.unknown();
    let f = x.pow(4) + x.pow(2) * 3 - x + 7;
    assert_eq!(f(i(2)), 16 + 12 - 2 + 7);

    let poly_ring = WrappingRing::new(PolyRingImpl::adjoint(i32::RING, "X"));
    let x = poly_ring.unknown();
    let f = x.pow(4) + x.pow(2) * 3 - x + 7;
    assert_eq!(f(i(2)), 16 + 12 - 2 + 7);
}

#[cfg(test)]
mod f1369 {
    use super::super::super::fq::fq_small::define_fq::*;

    type F37El = ZnElImpl<37, true>;
    gen_const_vector!(ConstVector2F37; F37El; V0, V1);
    type F1369MipoType = ConstVector2F37<{F37El::project(2)}, {F37El::project(33)}>;
    const F1369_MIPO: Vector<F1369MipoType, F37El> = Vector::new(F1369MipoType::INSTANCE);
    type F1369Type = FiniteExtensionImpl<StaticRing<F37El>, F1369MipoType, VectorArray<F37El, 2>>;
    pub static F1369: F1369Type = F1369Type::new(F37El::RING, F1369_MIPO, "α");
}

#[cfg(test)]
use f1369::F1369;

#[bench]
fn bench_poly_multiplication(b: &mut Bencher) {
    let poly_ring = PolyRingImpl::adjoint(&F1369, "X");
    let poly_ring = WrappingRing::new(&poly_ring);
    let a = poly_ring.embedding()(poly_ring.base_ring().clone_ring().generator());
    let x = poly_ring.unknown();
    let f = (0..=100).map(|i| x.pow(i) * i as i64).sum::<RingElWrapper<&_>>();
    let g = (0..=100).map(|i| x.pow(i) * (100 - i) as i64 * &a).sum::<RingElWrapper<&_>>();
    let h = (0..=100).map(|n| x.pow(n) * &a * ((100 - n) * n * (n + 1) / 2 + n * (n + 1) * (2 * n + 1) / 6) as i64).sum::<RingElWrapper<&_>>()
         + (0..100).map(|n| x.pow(200 - n) * &a * ((100 - n) * n * (n + 1) / 2 + n * (n + 1) * (2 * n + 1) / 6) as i64).sum::<RingElWrapper<&_>>();
    b.iter(|| {
        assert_eq!(h, &f * &g);
    });
}
