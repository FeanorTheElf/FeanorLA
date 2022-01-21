use super::super::super::alg::*;
use super::super::super::la::vec::*;
use super::ops::*;
use super::super::primality::*;

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

    pub fn from(&self, el: R::El) -> <Self as Ring>::El {
        let mut result = Vec::with_capacity(1);
        result.push(el);
        return Vector::new(result);
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

    fn is_integral(&self) -> bool {
        self.base_ring.is_integral()
    }

    fn is_euclidean(&self) -> bool {
        self.base_ring.is_field()
    }

    fn is_field(&self) -> bool {
        false
    }
    
    fn euclidean_div_rem(&self, mut lhs: Self::El, rhs: &Self::El) -> (Self::El, Self::El) {
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

    ///
    /// Calculates the polynomial division.
    /// This function panics if the division function panics when called with the leading 
    /// coefficient of rhs as second argument, or if lhs is not divisble by rhs.
    /// 
    /// Therefore, if the base ring division works for all divisible pairs of arguments,
    /// this is also the case for this function.
    /// 
    fn div(&self, mut lhs: Self::El, rhs: &Self::El) -> Self::El {
        assert!(!self.is_zero(&rhs));
        let rhs_lc = self.lc(&rhs).unwrap();
        if self.base_ring.is_field() {
            let rhs_lc_inv = self.base_ring.div(self.base_ring.one(), rhs_lc);
            self.poly_division(
                &mut lhs, 
                &rhs, 
                |x| Ok(self.base_ring.mul_ref(x, &rhs_lc_inv))
            ).unwrap()
        } else {
            self.poly_division(
                &mut lhs, 
                &rhs, 
                |x| Ok(self.base_ring.div(x.clone(), rhs_lc))
            ).unwrap()
        }
    }

    fn format(&self, el: &<Self as Ring>::El, f: &mut std::fmt::Formatter, in_prod: bool) -> std::fmt::Result {
        if in_prod {
            self.format_in_brackets(el, f)
        } else {
            poly_format(&self.base_ring, el.as_ref(), f, self.var_name)
        }
    }
}

impl<R> DivisibilityInformationRing for PolyRing<R> 
    where R: DivisibilityInformationRing
{
    fn is_divisibility_computable(&self) -> bool {
        self.base_ring.is_divisibility_computable()
    }

    fn quotient(&self, lhs: &Self::El, rhs: &Self::El) -> Option<Self::El> {
        assert!(!self.is_zero(rhs));
        let lc = self.lc(rhs).unwrap();
        return poly_division(
            &self.base_ring, 
            lhs.clone().into_owned(), 
            rhs.as_ref(), 
            |x| self.base_ring.quotient(x, lc).ok_or(())
        ).ok();
    }
}

#[cfg(test)]
use super::super::super::alg_env::*;

#[test]
fn test_poly_arithmetic() {
    let ring = PolyRing::adjoint(i32::RING, "X");
    let x = ring.bind::<RingAxiomsIntegralRing>(ring.unknown());

    let x2_2x_1 = (&x * &x) + (&x + &x) + 1;
    let x_1_square = (&x + 1) * (&x + 1);

    assert_eq!(x2_2x_1, x_1_square);
    assert!(x2_2x_1 != x);
    assert_eq!(ring.bind(ring.zero()), x2_2x_1 - x_1_square);
}

#[test]
fn test_format() {
    let ring = PolyRing::adjoint(i32::RING, "X");
    let x = ring.bind::<RingAxiomsIntegralRing>(ring.unknown());

    let poly = &x * &x * &x + &x * &x * 2 - 1;
    assert_eq!("-1 + 2 * X^2 + X^3", format!("{}", poly));
}

#[test]
fn test_poly_div() {
    let ring = PolyRing::adjoint(i32::RING, "X");
    let x = ring.bind::<RingAxiomsIntegralRing>(ring.unknown());

    let mut p = &x * &x * &x + &x * &x + &x + 1;
    let q = &x + 1;
    let expected = &x * &x + 1;
    let result = ring.bind(ring.poly_division(p.val_mut().unwrap(), q.val().unwrap(), |x| Ok(*x)).unwrap());
    assert_eq!(ring.bind(ring.zero()), p);
    assert_eq!(expected, result);
}

#[test]
fn test_poly_degree() {
    let ring = PolyRing::adjoint(i32::RING, "X");
    let x = ring.bind::<RingAxiomsIntegralRing>(ring.unknown());

    let p = &x * &x * &x + 4;
    assert_eq!(Some(3), ring.deg(p.unwrap()));
}