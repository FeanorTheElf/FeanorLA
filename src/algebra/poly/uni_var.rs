use super::super::super::ring::*;
use super::super::super::la::vec::*;
use super::super::super::bigint::*;
use super::super::super::embedding::*;
use super::ops::*;
use super::factoring;
use super::super::super::wrapper::*;

use vector_map::VecMap;

#[derive(Debug, Clone)]
pub struct PolyRing<R>
    where R: Ring
{
    base_ring: R,
    var_name: &'static str
}

impl<R> PolyRing<R>
    where R: Ring
{
    pub fn derive(&self, el: <Self as Ring>::El) -> <Self as Ring>::El {
        let mut result = el.into_owned();
        poly_formal_derivative(&self.base_ring, result.as_mut());
        return result;
    }

    pub const fn adjoint(base_ring: R, var_name: &'static str) -> Self {
        PolyRing {
            base_ring, var_name
        }
    }

    pub fn evaluate(&self, poly: &<Self as Ring>::El, value: R::El) -> R::El {
        poly_evaluate(&self.base_ring, value, poly.as_ref())
    }

    pub fn from(&self, el: R::El) -> <Self as Ring>::El {
        let mut result = Vec::with_capacity(1);
        result.push(el);
        return Vector::new(result);
    }

    pub fn evaluation_hom<'a>(&'a self, x: R::El) -> impl 'a + FnMut(<Self as Ring>::El) -> R::El {
        move |poly| self.evaluate(&poly, x.clone())
    }

    pub fn lift_hom<F, S>(&self, (ring_ext, mut hom): (S, F)) -> (PolyRing<S>, impl FnMut(<Self as Ring>::El) -> <PolyRing<S> as Ring>::El)
        where S: Ring, F: FnMut(R::El) -> S::El
    {
        (
            PolyRing::adjoint(ring_ext, self.var_name),
            move |poly| {
                Vector::new(poly.raw_data().into_vec().into_iter().map(&mut hom).collect())
            }
        )
    }

    pub fn deg(&self, el: &<Self as Ring>::El) -> Option<usize> {
        poly_degree(&self.base_ring, el.as_ref())
    }

    pub fn lc<'a>(&self, el: &'a <Self as Ring>::El) -> Option<&'a R::El> {
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
    pub fn poly_division<F>(&self, lhs: &mut <Self as Ring>::El, rhs: &<Self as Ring>::El, div_lc: F) -> Result<<Self as Ring>::El, ()>
        where F: FnMut(&R::El) -> Result<R::El, ()>
    {
        poly_division(&self.base_ring, lhs.as_mut(), rhs.as_ref(), div_lc)
    }

    pub fn unknown(&self) -> <Self as Ring>::El {
        let mut result = Vec::with_capacity(2);
        result.push(self.base_ring.zero());
        result.push(self.base_ring.one());
        return Vector::new(result);
    }

    pub fn base_ring(&self) -> &R {
        &self.base_ring
    }

    pub fn normalize(&self, f: <Self as Ring>::El) -> (<Self as Ring>::El, R::El) {
        assert!(self.base_ring().is_field().can_use());
        let lc = self.lc(&f).unwrap().clone();
        let lc_inv = self.base_ring().div(self.base_ring().one(), &lc);
        return (self.mul(f, self.from(lc_inv)), lc);
    }
}

impl<R> CanonicalEmbeddingInfo<R> for PolyRing<R> 
    where R: Ring
{

    fn has_embedding(&self, _from: &R) -> RingPropValue {
        RingPropValue::True
    }

    fn embed(&self, _from: &R, el: R::El) -> Self::El {
        self.from(el)
    }
}

impl<R> CanonicalEmbeddingInfo<PolyRing<R>> for PolyRing<R>
    where R: Ring
{
    fn has_embedding(&self, _from: &PolyRing<R>) -> RingPropValue {
        RingPropValue::True
    }

    fn embed(&self, _from: &PolyRing<R>, el: Self::El) -> Self::El {
        el
    }
}

impl<R> CanonicalIsomorphismInfo<PolyRing<R>> for PolyRing<R>
    where R: Ring
{
    fn has_isomorphism(&self, _from: &PolyRing<R>) -> RingPropValue {
        RingPropValue::True
    }

    fn preimage(&self, _from: &PolyRing<R>, el: Self::El) -> Self::El {
        el
    }
}

impl<R> Ring for PolyRing<R>
    where R: Ring
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
            take_mut::take_or_recover(
                val.at_mut(i), 
                || self.base_ring.unspecified_element(), 
                |v| self.base_ring.neg(v)
            );
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

    fn eq(&self, lhs: &Self::El, rhs: &Self::El) -> bool {
        let (shorter, longer) = if lhs.len() <= rhs.len() {
            (lhs, rhs)
        } else {
            (rhs, lhs)
        };
        for i in 0..shorter.len() {
            if !self.base_ring.eq(&shorter[i], &longer[i]) {
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
        panic!("Not a field")
    }

    fn format(&self, el: &<Self as Ring>::El, f: &mut std::fmt::Formatter, in_prod: bool) -> std::fmt::Result {
        if in_prod {
            self.format_in_brackets(el, f)
        } else {
            poly_format(&self.base_ring, el.as_ref(), f, self.var_name)
        }
    }
}

impl<R> DivisibilityInfoRing for PolyRing<R> 
    where R: DivisibilityInfoRing
{
    fn is_divisibility_computable(&self) -> bool {
        self.base_ring.is_divisibility_computable()
    }

    fn quotient(&self, lhs: &Self::El, rhs: &Self::El) -> Option<Self::El> {
        assert!(self.is_divisibility_computable());
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
        assert!(self.is_divisibility_computable());
        self.deg(el).map(|d| d == 0).unwrap_or(false) && self.base_ring.is_unit(&el[0])
    }
}

impl<R> EuclideanInfoRing for PolyRing<R> 
    where R: Ring
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

impl<R> RingElWrapper<PolyRing<R>>
    where R: Ring
{
    pub fn lc(&self) -> <WrappingRing<&R> as Ring>::El {
        self.base_ring().base_ring().bind(self.base_ring().lc(self.val()).unwrap().clone())
    }
}

use super::super::fq::*;

impl<R> UfdInfoRing for PolyRing<R>
    where R: FiniteRing
{
    fn is_ufd(&self) -> RingPropValue {
        if self.base_ring.is_field().can_use() && self.base_ring().characteristic() != 2 && self.base_ring().size() == self.base_ring().characteristic() {
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
            &self.base_ring().characteristic(), 
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
                &self.base_ring().characteristic(), 
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
                            &self.base_ring().characteristic(), 
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
use super::super::super::primitive::*;
#[cfg(test)]
use super::super::fq::zn_big::*;

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