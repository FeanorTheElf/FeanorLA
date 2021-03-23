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
                    accumulator = self.base_ring.mul_ref(accumulator, &values[var]);
                }
                assert!(key[var] == current_power);
                key[var] = 0;
                let added_value = self.base_ring.mul_ref(val, &accumulator);
                match reduced_poly.entry(key) {
                    Entry::Occupied(e) => {
                        take_mut::take_or_recover(e.into_mut(), || self.base_ring.unspecified_element(), |v| self.base_ring.add(v, added_value));
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
                take_mut::take_or_recover(coeff, || self.base_ring.unspecified_element(), |c| self.base_ring.add_ref(c, entry.1));
            } else {
                lhs.insert(entry.0.clone(), entry.1.clone());
            }
        }

        self.assert_valid(&lhs);
        return lhs;
    }

    fn mul_ref(&self, lhs: Self::El, rhs: &Self::El) -> Self::El {
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

                let new_coeff = self.base_ring.mul_ref(lhs_coeff.clone(), rhs_coeff);

                match result.entry(new_vars) {
                    Entry::Occupied(e) => {
                        take_mut::take_or_recover(e.into_mut(), || self.base_ring.unspecified_element(), |c| self.base_ring.add(c, new_coeff));
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
            take_mut::take_or_recover(coeff, || self.base_ring.unspecified_element(), |c| self.base_ring.neg(c));
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
    
    fn euclidean_div_rem(&self, _lhs: Self::El, _rhs: Self::El) -> (Self::El, Self::El) {
        panic!("not euclidean")
    }

    fn div(&self, _lhs: Self::El, _rhs: Self::El) -> Self::El {
        panic!("not a field")
    }

    fn format(&self, el: &<Self as Ring>::El, f: &mut std::fmt::Formatter) -> std::fmt::Result { 
        let mut it = el.iter();
        if let Some(first) = it.next() {
            self.base_ring.format(first.1, f)?;
            for el in it {
                write!(f, " + ")?;
                self.base_ring.format(el.1, f)?;
                write!(f, " *")?;
                for (i, pow) in el.0.iter().enumerate() {
                    write!(f, " {}^{}", self.var_names[i], pow)?;
                }
            }
        } else {
            self.base_ring.format(&self.base_ring.zero(), f)?;
        }
        return Ok(());
    }
}

#[test]
fn test_binomial_formula() {
    let mut ring = MultivariatePolyRing::new(StaticRing::<i32>::RING);
    let x = ring.adjoint("X");
    let y = ring.adjoint("Y");
    
    let x2_2xy_y2 = fixed_ring_env!{ &ring; x, y; {
        (x * x) + (x * y + x * y) + (y * y)
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

    // the polynomial x^2 y + 2 x y^2 + x + 13
    let poly = ring.add(fixed_ring_env!{ &ring; x, y; {
        (x * x * y) + (x * y * y + x * y * y) + x
    }}, ring.from(13));

    assert_eq!(14, ring.evaluate_at(poly.clone(), &[1, 0]));
    assert_eq!(12 + 36 + 2 + 13, ring.evaluate_at(poly, &[2, 3]));
}