use super::super::prelude::*;
use super::super::wrapper::*;
use super::ops::*;
use super::factoring;
use super::super::fq::*;

use vector_map::VecMap;

#[derive(Debug, Clone)]
pub struct PolyRing<R>
    where R: Ring
{
    base_ring: R,
    var_name: &'static str
}

impl<R> PolyRing<R>
    where R: CanonicalIsomorphismInfo<R>
{
    pub fn derive(&self, el: El<Self>) -> El<Self> {
        let mut result = el.into_owned();
        poly_formal_derivative(&self.base_ring, result.as_mut());
        return result;
    }

    pub const fn adjoint(base_ring: R, var_name: &'static str) -> Self {
        PolyRing {
            base_ring, var_name
        }
    }

    pub fn evaluate<S>(&self, poly: &El<Self>, value: S::El, ring: &S) -> S::El 
        where S: Ring + CanonicalEmbeddingInfo<R>
    {
        poly_evaluate(&self.base_ring, value, poly.as_ref(), ring)
    }

    pub fn from(&self, el: R::El) -> El<Self> {
        let mut result = Vec::with_capacity(1);
        result.push(el);
        return Vector::new(result);
    }

    pub fn evaluation_hom<'a>(&'a self, x: R::El) -> impl 'a + FnMut(El<Self>) -> R::El {
        move |poly| self.evaluate(&poly, x.clone(), self.base_ring())
    }

    pub fn lift_hom<F, S>(&self, (ring_ext, mut hom): (S, F)) -> (PolyRing<S>, impl FnMut(El<Self>) -> El<PolyRing<S>>)
        where S: CanonicalIsomorphismInfo<S>, F: FnMut(R::El) -> S::El
    {
        (
            PolyRing::adjoint(ring_ext, self.var_name),
            move |poly| {
                Vector::new(poly.raw_data().into_iter().map(&mut hom).collect())
            }
        )
    }

    pub fn deg(&self, el: &El<Self>) -> Option<usize> {
        poly_degree(&self.base_ring, el.as_ref())
    }

    pub fn lc<'a>(&self, el: &'a El<Self>) -> Option<&'a R::El> {
        self.deg(el).map(|i| &el[i])
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

    pub fn unknown(&self) -> El<Self> {
        let mut result = Vec::with_capacity(2);
        result.push(self.base_ring.zero());
        result.push(self.base_ring.one());
        return Vector::new(result);
    }

    pub fn base_ring(&self) -> &R {
        &self.base_ring
    }

    pub fn normalize(&self, f: El<Self>) -> (El<Self>, R::El) {
        assert!(self.base_ring().is_field().can_use());
        let lc = self.lc(&f).unwrap().clone();
        let lc_inv = self.base_ring().div(self.base_ring().one(), &lc);
        return (self.mul(f, self.from(lc_inv)), lc);
    }

    pub fn scale(&self, mut f: El<Self>, factor: &El<R>) -> El<Self> {
        f.scale(factor, self.base_ring());
        return f;
    }
}

impl<R> CanonicalEmbeddingInfo<R> for PolyRing<R> 
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

impl<R> CanonicalEmbeddingInfo<R> for PolyRing<&R> 
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

impl<R> CanonicalEmbeddingInfo<&R> for PolyRing<R> 
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

impl<R> CanonicalEmbeddingInfo<PolyRing<R>> for PolyRing<R>
    where R: CanonicalIsomorphismInfo<R>
{
    fn has_embedding(&self, from: &PolyRing<R>) -> RingPropValue {
        self.base_ring().has_embedding(from.base_ring())
    }

    fn embed(&self, from: &PolyRing<R>, el: Self::El) -> Self::El {
        assert!(self.has_embedding(from).can_use());
        Vector::new(el.raw_data().into_iter().map(|x| self.base_ring().embed(from.base_ring(), x)).collect())
    }
}

impl<R> CanonicalIsomorphismInfo<PolyRing<R>> for PolyRing<R>
    where R: CanonicalIsomorphismInfo<R>
{
    fn has_isomorphism(&self, from: &PolyRing<R>) -> RingPropValue {
        self.base_ring().has_isomorphism(from.base_ring())
    }

    fn preimage(&self, from: &PolyRing<R>, el: Self::El) -> Self::El {
        assert!(self.has_isomorphism(from).can_use());
        Vector::new(el.raw_data().into_iter().map(|x| self.base_ring().preimage(from.base_ring(), x)).collect())
    }
}

impl<R> RingBase for PolyRing<R>
    where R: CanonicalIsomorphismInfo<R>
{
    type El = Vector<Vec<R::El>, R::El>;

    fn add_ref(&self, lhs: Self::El, rhs: &Self::El) -> Self::El {
        poly_add(&self.base_ring, lhs, rhs.as_ref())
    }

    fn mul_ref(&self, lhs: &Self::El, rhs: &Self::El) -> Self::El {
        poly_mul(&self.base_ring, lhs.as_ref(), rhs.as_ref())
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

    fn characteristic(&self) -> BigInt {
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

    fn from_z_big(&self, x: &BigInt) -> Self::El {
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

impl<R> DivisibilityInfoRing for PolyRing<R> 
    where R: DivisibilityInfoRing + CanonicalIsomorphismInfo<R>
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
            |x| self.base_ring.quotient(x, lc).ok_or(())
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

impl<R> EuclideanInfoRing for PolyRing<R> 
    where R: Ring + CanonicalIsomorphismInfo<R>
{
    fn is_euclidean(&self) -> RingPropValue {
        if self.base_ring.is_field().can_use() {
            return RingPropValue::True;
        } else {
            return RingPropValue::Unknown;
        }
    }

    fn euclidean_deg(&self, el: Self::El) -> BigInt {
        assert!(self.is_euclidean().can_use());
        self.deg(&el).map(|x| (x + 1) as i64).map(BigInt::from).unwrap_or(BigInt::ZERO)
    }

    fn euclidean_div_rem(&self, mut lhs: Self::El, rhs: &Self::El) -> (Self::El, Self::El) {
        assert!(self.is_euclidean().can_use());
        assert!(self.base_ring().is_field().can_use());
        assert!(!self.is_zero(&rhs));
        let rhs_lc = self.lc(&rhs).unwrap();
        let rhs_lc_inv = self.base_ring.div(self.base_ring.one(), rhs_lc);
        let q = self.poly_division(
            &mut lhs, 
            &rhs, 
            |x| Ok(self.base_ring.mul_ref(x, &rhs_lc_inv))
        ).unwrap();
        return (q, lhs);
    }
}

impl<R> PartialEq for PolyRing<R> 
    where R: Ring + PartialEq
{
    fn eq(&self, rhs: &Self) -> bool {
        self.base_ring() == rhs.base_ring()
    }
}

impl<R> RingElWrapper<PolyRing<R>>
    where R: Ring
{
    pub fn lc(&self) -> Option<El<WrappingRing<R>>> {
        self.parent_ring().lc(self.val()).map(|x| self.parent_ring().base_ring().clone().bind_by_value(x.clone()))
    }

    pub fn deg(&self) -> Option<usize> {
        self.parent_ring().deg(self.val())
    }

    pub fn coefficient_at(&self, i: usize) -> RingElWrapper<R> {
        self.parent_ring().base_ring().clone().bind_by_value(self.val()[i].clone())
    }

    pub fn scaled(self, factor: &RingElWrapper<R>) -> Self {
        let (el, ring) = self.destruct();
        ring.bind_by_value(ring.scale(el, factor.val()))
    }
}

trait Evaluatable<S: Ring>: Ring {

    fn evaluate_at(&self, f: &El<Self>, x: El<S>) -> El<S>;
}

impl<R, S> Evaluatable<WrappingRing<S>> for PolyRing<R>
    where S: Ring + CanonicalEmbeddingInfo<R>, R: Ring
{
    fn evaluate_at(&self, f: &El<Self>, x: El<WrappingRing<S>>) -> El<WrappingRing<S>> { 
        let poly_vec = f.as_ref();
        let (el, ring) = x.destruct();
        let coeff_ring = self.base_ring();
        let result: S::El = poly_evaluate(coeff_ring, el, poly_vec, &ring);
        return ring.bind_by_value(result);
    }
}

impl<S, P> Evaluatable<S> for &P
    where S: Ring, P: Evaluatable<S>
{
    fn evaluate_at(&self, f: &El<Self>, x: El<S>) -> El<S> { 
        (**self).evaluate_at(f, x)
    }
}

impl<S, P> FnOnce<(RingElWrapper<S>, )> for RingElWrapper<P>
    where S: Ring, P: Evaluatable<WrappingRing<S>>
{
    type Output = RingElWrapper<S>;

    extern "rust-call" fn call_once(
        mut self, 
        (x, ): (RingElWrapper<S>, )
    ) -> Self::Output {
        self.call_mut((x, ))
    }
}

impl<S, P> FnMut<(RingElWrapper<S>, )> for RingElWrapper<P>
    where S: Ring, P: Evaluatable<WrappingRing<S>>
{
    extern "rust-call" fn call_mut(
        &mut self, 
        (x, ): (RingElWrapper<S>, )
    ) -> Self::Output {
        self.call((x, ))
    }
}

impl<S, P> Fn<(RingElWrapper<S>, )> for RingElWrapper<P>
    where S: Ring, P: Evaluatable<WrappingRing<S>>
{
    extern "rust-call" fn call(
        &self, 
        (x, ): (RingElWrapper<S>, )
    ) -> Self::Output {
        self.parent_ring().evaluate_at(self.val(), x)
    }
}

impl<R> UfdInfoRing for PolyRing<R>
    where R: FiniteRing + DivisibilityInfoRing
{
    fn is_ufd(&self) -> RingPropValue {
        if self.base_ring.is_field().can_use() && self.base_ring().characteristic().is_odd() {
            return RingPropValue::True;
        } else {
            return RingPropValue::Unknown;
        }
    }

    fn is_prime(&self, el: &Self::El) -> bool {
        assert!(self.is_ufd().can_use());
        if self.is_zero(el) {
            return false;
        }
        let d = self.deg(el).unwrap();
        let sqrfree_part = factoring::poly_squarefree_part(self.base_ring(), el.as_ref().into_owned());
        if self.deg(&sqrfree_part) != Some(d) {
            return false;
        }
        let distinct_degree_factorization = factoring::distinct_degree_factorization(
            self.base_ring(), 
            &self.base_ring().size(), 
            sqrfree_part
        );
        if d >= distinct_degree_factorization.len() || self.is_unit(&distinct_degree_factorization[d]) {
            return false;
        }
        return true;
    }

    fn calc_factor(&self, el: &Self::El) -> Option<Self::El> {
        assert!(self.is_ufd().can_use());
        assert!(!self.is_zero(el));
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

    fn factor<'a>(&'a self, mut el: Self::El) -> VecMap<RingElWrapper<&'a Self>, usize> {
        assert!(self.is_ufd().can_use());
        assert!(!self.is_zero(&el));
        let mut result = VecMap::new();
        let mut unit = self.base_ring().one();
        while !self.is_unit(&el) {
            let sqrfree_part = factoring::poly_squarefree_part(self.base_ring(), el.clone());
            for (d, el) in factoring::distinct_degree_factorization(
                self.base_ring(), 
                &self.base_ring().size(), 
                sqrfree_part.clone()
            ).into_iter().enumerate() {
                let mut stack = Vec::new();
                stack.push(el);
                while let Some(el) = stack.pop() {
                    let (el, scaling) = self.normalize(el);
                    unit = self.base_ring().mul(unit, scaling);
                    if self.is_one(&el) {
                        continue;
                    } else if self.deg(&el) == Some(d) {
                        let wrapped_el = self.bind(el);
                        if let Some(power) = result.get_mut(&wrapped_el) {
                            *power += 1;
                        } else {
                            debug_assert!(!self.is_unit(wrapped_el.val()));
                            result.insert(wrapped_el, 1);
                        }
                    } else {
                        let factor = factoring::cantor_zassenhaus(
                            self.base_ring(), 
                            &self.base_ring().size(), 
                            el.clone(), 
                            d
                        );
                        stack.push(self.quotient(&el, &factor).unwrap());
                        stack.push(factor);
                    }
                }
            }
            el = self.quotient(&el, &sqrfree_part).unwrap();
        }
        unit = self.base_ring().mul_ref(&unit, el.at(0));
        debug_assert!(self.base_ring().is_unit(&unit));
        if !self.base_ring().is_one(&unit) {
            result.insert(self.bind(self.from(unit)), 1);
        }
        return result;
    }
}

#[cfg(test)]
use super::super::fq::zn_big::*;
#[cfg(test)]
use super::super::fq::fq_small::*;
#[cfg(test)]
use test::Bencher;

#[test]
fn test_poly_arithmetic() {
    let ring = PolyRing::adjoint(i32::RING, "X");
    let x = ring.bind(ring.unknown());

    let x2_2x_1 = (&x * &x) + (&x + &x) + 1;
    let x_1_square = (&x + 1) * (&x + 1);

    assert_eq!(x2_2x_1, x_1_square);
    assert!(x2_2x_1 != x);
    assert_eq!(ring.bind(ring.zero()), x2_2x_1 - x_1_square);
}

#[test]
fn test_format() {
    let ring = PolyRing::adjoint(i32::RING, "X");
    let x = ring.bind(ring.unknown());

    let poly = &x * &x * &x + &x * &x * 2 - 1;
    assert_eq!("-1 + 2 * X^2 + X^3", format!("{}", poly));
}

#[test]
fn test_poly_div() {
    let ring = PolyRing::adjoint(i32::RING, "X");
    let x = ring.bind(ring.unknown());

    let mut p = &x * &x * &x + &x * &x + &x + 1;
    let q = &x + 1;
    let expected = &x * &x + 1;
    let result = ring.bind(ring.poly_division(p.val_mut(), q.val(), |x| Ok(*x)).unwrap());
    assert_eq!(ring.bind(ring.zero()), p);
    assert_eq!(expected, result);
}

#[test]
fn test_quotient() {
    let ring = PolyRing::adjoint(i32::RING, "X");
    let x = ring.bind(ring.unknown());

    let p = ring.bind(ring.one());
    let q = &x + 1;
    assert_eq!(None, ring.quotient(p.val(), q.val()));
}

#[test]
fn test_poly_degree() {
    let ring = PolyRing::adjoint(i32::RING, "X");
    let x = ring.bind(ring.unknown());

    let p = &x * &x * &x + 4;
    assert_eq!(Some(3), ring.deg(p.val()));
}

#[test]
fn test_factor() {
    let coeff_ring = Zn::new(BigInt::from(3));
    let ring = PolyRing::adjoint(coeff_ring, "X");
    let x = ring.bind(ring.unknown());

    let p = x.clone().pow(9) - &x;
    let mut expected = VecMap::new();
    expected.insert(x.clone(), 1);
    expected.insert(&x + 1, 1);
    expected.insert(&x + 2, 1);
    expected.insert(&x * &x + &x + 2, 1);
    expected.insert(&x * &x + &x * 2 + 2, 1);
    expected.insert(&x * &x + 1, 1);
    let factorization = ring.factor(p.into_val());
    assert_eq!(expected, factorization);
}

#[test]
fn test_factor_fq() {
    let coeff_ring = F49.clone();
    let ring = PolyRing::adjoint(&coeff_ring, "X");
    let x = ring.bind(ring.unknown());

    let f = x.pow(2) + &x * 6 + 3;
    let factor = <_ as Iterator>::next(&mut ring.factor(f.into_val()).into_iter()).unwrap().0.clone_ring();
    let coeff_ring_gen = -factor.coefficient_at(0) / factor.coefficient_at(1);
    let a = ring.bind(ring.from(coeff_ring_gen.into_val())); // up to field automorphism, this is the generator picked by sage

    let p = x.clone().pow(14) - &x;
    let mut expected = VecMap::new();
    expected.insert(x.clone(), 1);
    expected.insert(&x + 6, 1);
    expected.insert(x.pow(6) + (&a * 3 + 6) * x.pow(5) + x.pow(4) * 2 + (&a * 3 + 5) * x.pow(3) + x.pow(2) * 2 + (&a * 3 + 6) * &x + 1, 1);
    expected.insert(x.pow(6) + (&a * 4 + 2) * x.pow(5) + x.pow(4) * 2 + (&a * 4 + 1) * x.pow(3) + x.pow(2) * 2 + (&a * 4 + 2) * &x + 1, 1);
    let factorization = ring.factor(p.into_val());
    assert_eq!(expected, factorization);
}

#[test]
fn test_is_prime() {
    let coeff_ring = Zn::new(BigInt::from(3));
    let ring = PolyRing::adjoint(coeff_ring, "X");
    let x = ring.bind(ring.unknown());

    let p = &x + 1;
    let q = &x * &x * &x + &x * 2 + 2;
    let a = &x * &x + 2;
    let b = &x * &x + &x * 2 + 1;
    assert_eq!(true, ring.is_prime(p.val()));
    assert_eq!(true, ring.is_prime(q.val()));
    assert_eq!(false, ring.is_prime(a.val()));
    assert_eq!(false, ring.is_prime(b.val()));
}

#[test]
fn test_evaluate() {
    let i = embedding(i64::RING, i64::RING.bind_ring_by_value());
    let poly_ring = PolyRing::adjoint(i64::RING, "X").bind_ring_by_value();
    let x = poly_ring.wrapped_ring().clone().bind_by_value(poly_ring.wrapped_ring().unknown());
    let f = x.pow(4) + x.pow(2) * 3 - x + 7;
    assert_eq!(f(i(2)), 16 + 12 - 2 + 7);

    let poly_ring = PolyRing::adjoint(i64::RING, "X");
    let x = poly_ring.bind(poly_ring.unknown());
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
    type F1369Type = SimpleRingExtension<StaticRing<F37El>, F1369MipoType, VectorArray<F37El, 2>>;
    pub static F1369: F1369Type = F1369Type::new(F37El::RING, F1369_MIPO, "α");
}

#[cfg(test)]
use f1369::F1369;

#[bench]
fn bench_poly_multiplication(b: &mut Bencher) {
    let poly_ring = PolyRing::adjoint(F1369.clone(), "x");
    let a = poly_ring.bind(poly_ring.from(poly_ring.base_ring().generator()));
    let x = poly_ring.bind(poly_ring.unknown());
    let f = (0..=100).map(|i| x.pow(i) * i as i64).sum::<RingElWrapper<&_>>();
    let g = (0..=100).map(|i| x.pow(i) * (100 - i) as i64 * &a).sum::<RingElWrapper<&_>>();
    let h = (0..=100).map(|n| x.pow(n) * &a * ((100 - n) * n * (n + 1) / 2 + n * (n + 1) * (2 * n + 1) / 6) as i64).sum::<RingElWrapper<&_>>()
         + (0..100).map(|n| x.pow(200 - n) * &a * ((100 - n) * n * (n + 1) / 2 + n * (n + 1) * (2 * n + 1) / 6) as i64).sum::<RingElWrapper<&_>>();
    b.iter(|| {
        assert_eq!(h, &f * &g);
    });
}
