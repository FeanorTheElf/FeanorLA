use super::super::alg::*;
use super::super::alg_env::*;
use super::eea::eea;

use std::ops::{AddAssign, MulAssign, SubAssign, DivAssign, Add, Mul, Sub, Div, Neg};

use std::collections::BTreeMap;
use std::collections::btree_map::Entry;

#[derive(Debug, Clone)]
pub struct MultivariatePolyRing<R>
    where R: Ring
{
    base_ring: R,
    var_names: Vec<&'static str>
}

impl<R> MultivariatePolyRing<R>
    where R: Ring
{
    pub fn from(&self, el: R::El) -> <Self as Ring>::El {
        let mut result = BTreeMap::new();
        result.insert(Vec::new(), el);

        self.assert_valid(&result);
        return result;
    }

    fn assert_valid(&self, el: &<Self as Ring>::El) {
        for key in el.keys() {
            debug_assert!(key.last().is_none() || *key.last().unwrap() != 0);
            debug_assert!(key.len() <= self.var_names.len());
        }
    }

    pub fn new(base_ring: R) -> Self {
        MultivariatePolyRing {
            base_ring: base_ring,
            var_names: Vec::new()
        }
    }

    pub fn adjoint(&mut self, var: &'static str) -> <Self as Ring>::El {
        assert!(!self.var_names.contains(&var));
        let mut result = BTreeMap::new();
        result.insert(
            (0..self.var_names.len()).map(|_| 0).chain(std::iter::once(1)).collect(), 
            self.base_ring.one()
        );
        self.var_names.push(var);

        self.assert_valid(&result);
        return result;
    }

    ///
    /// Evaluates the given polynomial at the given values, i.e. calculates the
    /// value in the base ring of the expression one gets from replacing the
    /// i-th unknown by the i-th given value in the formal sum representation
    /// of this polynomial.
    /// 
    /// # Complexity
    /// 
    /// Complexity is O(V * N) where V is the number of variables and N is the
    /// number of monomials in the given polynomial. The number of multiplications
    /// in the base ring is O(V * P) where P is the highest power of any variable
    /// occuring in the polynomial.
    /// 
    pub fn evaluate_at(&self, mut poly: <Self as Ring>::El, values: &[R::El]) -> R::El {
        self.assert_valid(&poly);

        if poly.len() == 0 {
            return self.base_ring.zero();
        }
        for var in 0..self.var_names.len() {
            let mut accumulator = self.base_ring.one();
            let mut current_power = 0;
            let mut reduced_poly = BTreeMap::new();
            for (mut key, val) in poly.into_iter() {
                if key.len() <= var {
                    reduced_poly.insert(key, val);
                    continue;
                }
                while key[var] > current_power {
                    current_power += 1;
                    accumulator = self.base_ring.mul_ref(&accumulator, &values[var]);
                }
                assert!(key[var] == current_power);
                key[var] = 0;
                let added_value = self.base_ring.mul_ref(&val, &accumulator);
                match reduced_poly.entry(key) {
                    Entry::Occupied(e) => {
                        take_mut::take_or_recover(
                            e.into_mut(), 
                            || self.base_ring.unspecified_element(), 
                            |v| self.base_ring.add(v, added_value)
                        );
                    },
                    Entry::Vacant(e) => {
                        e.insert(added_value);
                    }
                }
            }
            poly = reduced_poly;
        }
        return poly.into_iter()
            .map(|(_key, coeff)| coeff)
            .fold(self.base_ring.zero(), |a, b| self.base_ring.add(a, b));
    }
}

impl<R> Ring for MultivariatePolyRing<R>
    where R: Ring
{
    type El = BTreeMap<Vec<usize>, R::El>;

    fn add_ref(&self, mut lhs: Self::El, rhs: &Self::El) -> Self::El {
        self.assert_valid(&lhs);
        self.assert_valid(&rhs);

        for entry in rhs.iter() {
            if let Some(coeff) = lhs.get_mut(entry.0) {
                take_mut::take_or_recover(
                    coeff, 
                    || self.base_ring.unspecified_element(), 
                    |c| self.base_ring.add_ref(c, entry.1)
                );
            } else {
                lhs.insert(entry.0.clone(), entry.1.clone());
            }
        }

        self.assert_valid(&lhs);
        return lhs;
    }

    fn mul_ref(&self, lhs: &Self::El, rhs: &Self::El) -> Self::El {
        self.assert_valid(&lhs);
        self.assert_valid(&rhs);

        let mut result = self.zero();
        for (lhs_vars, lhs_coeff) in lhs.into_iter() {
            for (rhs_vars, rhs_coeff) in rhs.iter() {
                let mut new_vars = lhs_vars.clone();
                let common_index = rhs_vars.len().min(lhs_vars.len());
                for i in 0..common_index {
                    new_vars[i] += rhs_vars[i];
                }
                new_vars.extend(&rhs_vars[common_index..]);

                let new_coeff = self.base_ring.mul_ref(lhs_coeff, rhs_coeff);

                match result.entry(new_vars) {
                    Entry::Occupied(e) => {
                        take_mut::take_or_recover(
                            e.into_mut(), 
                            || self.base_ring.unspecified_element(), 
                            |c| self.base_ring.add(c, new_coeff)
                        );
                    },
                    Entry::Vacant(e) => {
                        e.insert(new_coeff);
                    }
                };
            }
        }
        
        self.assert_valid(&result);
        return result;
    }

    fn neg(&self, mut val: Self::El) -> Self::El {
        self.assert_valid(&val);

        for (_, coeff) in val.iter_mut() {
            take_mut::take_or_recover(
                coeff, 
                || self.base_ring.unspecified_element(), 
                |c| self.base_ring.neg(c)
            );
        }

        self.assert_valid(&val);
        return val;
    }

    fn zero(&self) -> Self::El {
        let result = BTreeMap::new();

        self.assert_valid(&result);
        return result;
    }

    fn one(&self) -> Self::El {
        let result = self.from(self.base_ring.one());

        self.assert_valid(&result);
        return result;
    }

    fn eq(&self, lhs: &Self::El, rhs: &Self::El) -> bool {
        self.assert_valid(&lhs);
        self.assert_valid(&rhs);

        if lhs.len() != rhs.len() {
            return false;
        }
        for ((lhs_key, lhs_coeff), (rhs_key, rhs_coeff)) in lhs.iter().zip(rhs.iter()) {
            if lhs_key != rhs_key || !self.base_ring.eq(lhs_coeff, rhs_coeff) {
                return false;
            }
        }
        return true;
    }

    fn is_integral(&self) -> bool {
        self.base_ring.is_integral()
    }

    fn is_euclidean(&self) -> bool {
        if self.var_names.len() <= 1 {
            unimplemented!()
        }
        false
    }

    fn is_field(&self) -> bool {
        if self.var_names.len() == 0 {
            unimplemented!()
        }
        false
    }
    
    fn euclidean_div_rem(&self, _lhs: Self::El, _rhs: &Self::El) -> (Self::El, Self::El) {
        panic!("not euclidean")
    }

    fn div(&self, lhs: Self::El, rhs: &Self::El) -> Self::El {
        panic!("not a field")
    }

    fn format(&self, el: &<Self as Ring>::El, f: &mut std::fmt::Formatter, in_prod: bool) -> std::fmt::Result {
        let mut it = el.iter();
        if in_prod {
            return self.format_in_brackets(el, f);
        } else if let Some(first) = it.next() {

            let write_part = |key: &Vec<usize>, coeff, f: &mut std::fmt::Formatter| {
                self.base_ring.format(coeff, f, true)?;
                if key.len() > 0 {
                    write!(f, " *")?;
                    for (i, pow) in key.iter().enumerate() {
                        if *pow > 0 {
                            write!(f, " {}^{}", self.var_names[i], pow)?;
                        }
                    }
                }
                return Ok(());
            };
            
            write_part(first.0, first.1, f)?;
            for el in it {
                write!(f, " + ");
                write_part(el.0, el.1, f)?;
            }
            return Ok(());
        } else {
            return self.base_ring.format(&self.base_ring.zero(), f, in_prod);
        }
    }
}

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
    pub const fn adjoint(base_ring: R, var_name: &'static str) -> Self {
        PolyRing {
            base_ring, var_name
        }
    }

    pub fn from(&self, el: R::El) -> <Self as Ring>::El {
        let mut result = Vec::with_capacity(1);
        result.push(el);
        return result;
    }

    pub fn deg(&self, el: &<Self as Ring>::El) -> Option<usize> {
        el.iter().enumerate().filter(|(_i, x)| !self.base_ring.is_zero(x)).map(|(i, _x)| i).max()
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
    pub fn poly_division<F>(&self, lhs: &mut <Self as Ring>::El, rhs: &<Self as Ring>::El, mut div_lc: F) -> Result<<Self as Ring>::El, ()>
        where F: FnMut(&R::El) -> Result<R::El, ()>
    {
        assert!(!self.is_zero(rhs));
        let rhs_deg = self.deg(rhs).unwrap();
        let mut result = Vec::new();
        while self.deg(lhs) >= self.deg(rhs) {
            let lhs_deg = self.deg(lhs).unwrap();
            let coeff = div_lc(&lhs[lhs_deg])?;
            let pow = lhs_deg - rhs_deg;
            result.push(coeff);
            for i in 0..=rhs_deg {
                take_mut::take_or_recover(
                    &mut lhs[i + pow], 
                    || self.base_ring.unspecified_element(), 
                    |v| self.base_ring.sub_ref_snd(v, &rhs[i])
                );
            }
            assert!(self.base_ring.is_zero(&lhs[lhs_deg]));
        }
        result.reverse();
        return Ok(result);
    }

    pub fn unknown(&self) -> <Self as Ring>::El {
        let mut result = Vec::with_capacity(2);
        result.push(self.base_ring.zero());
        result.push(self.base_ring.one());
        return result;
    }
}

impl<R> Ring for PolyRing<R>
    where R: Ring
{
    type El = Vec<R::El>;

    fn add_ref(&self, mut lhs: Self::El, rhs: &Self::El) -> Self::El {
        for i in 0..rhs.len() {
            if lhs.len() <= i {
                lhs.push(rhs[i].clone());
            } else {
                take_mut::take_or_recover(
                    &mut lhs[i], 
                    || self.base_ring.unspecified_element(), 
                    |v| self.base_ring.add_ref(v, &rhs[i])
                );
            }
        }
        return lhs;
    }

    fn mul_ref(&self, lhs: &Self::El, rhs: &Self::El) -> Self::El {
        let mut result = Vec::with_capacity(lhs.len() + rhs.len());
        for i in 0..(lhs.len() + rhs.len()) {
            let mut val = self.base_ring.zero();
            for j in (rhs.len().max(i + 1) - rhs.len())..lhs.len().min(i + 1) {
                val = self.base_ring.add(val, self.base_ring.mul_ref(&lhs[j], &rhs[i - j]));
            }
            result.push(val);
        }
        return result;
    }

    fn neg(&self, mut val: Self::El) -> Self::El {
        for i in 0..val.len() {
            take_mut::take_or_recover(
                &mut val[i], 
                || self.base_ring.unspecified_element(), 
                |v| self.base_ring.neg(v)
            );
        }
        return val;
    }

    fn zero(&self) -> Self::El {
        let result = Vec::new();
        return result;
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
    
    fn euclidean_div_rem(&self, _lhs: Self::El, _rhs: &Self::El) -> (Self::El, Self::El) {
        panic!("not euclidean")
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
            return self.format_in_brackets(el, f);
        } else if self.is_zero(el) {
            return self.base_ring.format(&self.base_ring.zero(), f, in_prod);
        } else {
            let mut it = el.iter().enumerate().filter(|(_i, x)| !self.base_ring.is_zero(x)).rev();

            let print_part = |pow, coeff, formatter: &mut std::fmt::Formatter| {
                self.base_ring.format(coeff, formatter, true)?;
                if pow > 0 {
                    write!(formatter, " * {}^{}", self.var_name, pow)?;
                }
                return Ok(());
            };

            let (fst_pow, fst_coeff) = it.next().unwrap();
            print_part(fst_pow, fst_coeff, f)?;
            for (pow, coeff) in it {
                write!(f, " + ")?;
                print_part(pow, coeff, f)?;
            }
            return Ok(());
        }
    }
}

#[test]
fn test_binomial_formula() {
    let mut ring = MultivariatePolyRing::new(StaticRing::<i32>::RING);
    let x = ring.adjoint("X");
    let y = ring.adjoint("Y");
    let two = ring.from(2);
    
    let x2_2xy_y2 = fixed_ring_env!{ &ring; x, y, two; {
        (x * x) + (two * x * y) + (y * y)
    }};
    let x_y_square = fixed_ring_env!{ &ring; x, y; {
        (x + y) * (x + y)
    }};

    assert!(ring.eq(&x2_2xy_y2, &x_y_square));
    assert!(!ring.eq(&x2_2xy_y2, &x));
}

#[test]
fn test_eq() {
    let mut ring = MultivariatePolyRing::new(StaticRing::<i32>::RING);
    let one = ring.from(1);
    let _x = ring.adjoint("x");
    let one_prime = ring.from(1);
    assert!(ring.eq(&one, &one_prime));
}

#[test]
fn test_evaluate_at() {
    let mut ring = MultivariatePolyRing::new(StaticRing::<i32>::RING);
    let x = ring.adjoint("X");
    let y = ring.adjoint("Y");
    let thirteen = ring.from(13);

    // the polynomial x^2 y + 2 x y^2 + x + 13
    let poly = fixed_ring_env!{ &ring; x, y, thirteen; {
        (x * x * y) + (x * y * y + x * y * y) + x + thirteen
    }};

    assert_eq!(14, ring.evaluate_at(poly.clone(), &[1, 0]));
    assert_eq!(12 + 36 + 2 + 13, ring.evaluate_at(poly, &[2, 3]));
}

#[test]
fn test_assumption_option_ord() {
    // in PolyRing::poly_div we rely on this behavior, 
    // and it is useful in general for degree comparison
    assert!(None < Some(0));
    assert!(Some(0) < Some(1));
}

#[test]
fn test_poly_arithmetic() {
    let ring = PolyRing::adjoint(StaticRing::<i32>::RING, "X");
    let x = ring.unknown();
    let one = ring.one();

    let x2_2x_1 = fixed_ring_env!{ &ring; x, one; {
        (x * x) + (x + x) + one
    }};
    let x_one_square = fixed_ring_env!{ &ring; x, one; {
        (x + one) * (x + one)
    }};

    assert!(ring.eq(&x2_2x_1, &x_one_square));
    assert!(!ring.eq(&x2_2x_1, &x));
    assert!(ring.is_zero(&ring.sub(x2_2x_1, x_one_square)));
}

#[test]
fn test_format() {
    let ring = PolyRing::adjoint(StaticRing::<i32>::RING, "X");
    let x = ring.unknown();
    let one = ring.one();

    let poly = fixed_ring_env!{ &ring; x, one; {
        x * x * x + (one + one) * x * x - one
    }};
    assert_eq!("1 * X^3 + 2 * X^2 + -1", format!("{}", ring.display(&poly)));
}

#[test]
fn test_format_multivar_poly_ring() {
    let mut ring = MultivariatePolyRing::new(StaticRing::<i32>::RING);
    let x = ring.adjoint("X");
    let y = ring.adjoint("Y");
    let one = ring.one();

    let poly = fixed_ring_env!{ &ring; x, y, one; {
        x * x * x - y + (one + one) * y * x - one
    }};
    assert_eq!("-1 + -1 * Y^1 + 2 * X^1 Y^1 + 1 * X^3", format!("{}", ring.display(&poly)));
}