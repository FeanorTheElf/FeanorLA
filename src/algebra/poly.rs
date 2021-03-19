use super::super::alg::*;
use super::super::alg_env::*;
use super::eea::eea;

use std::ops::{AddAssign, MulAssign, SubAssign, DivAssign, Add, Mul, Sub, Div, Neg};

use std::collections::HashMap;
use std::collections::hash_map::Entry;

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
        let mut result = HashMap::new();
        result.insert(Vec::new(), el);
        return result;
    }

    pub fn new(base_ring: R) -> Self {
        MultivariatePolyRing {
            base_ring: base_ring,
            var_names: Vec::new()
        }
    }

    pub fn adjoint(&mut self, var: &'static str) -> <Self as Ring>::El {
        assert!(!self.var_names.contains(&var));
        self.var_names.push(var);
        let mut result = HashMap::new();
        result.insert(vec![self.var_names.len() - 1], self.base_ring.one());
        return result;
    }
}

impl<R> Ring for MultivariatePolyRing<R>
    where R: Ring
{
    type El = HashMap<Vec<usize>, R::El>;

    fn add_ref(&self, mut lhs: Self::El, rhs: &Self::El) -> Self::El {
        for entry in rhs.iter() {
            if let Some(coeff) = lhs.get_mut(entry.0) {
                take_mut::take(coeff, |c| self.base_ring.add_ref(c, entry.1));
            } else {
                lhs.insert(entry.0.clone(), entry.1.clone());
            }
        }
        return lhs;
    }

    fn mul_ref(&self, lhs: Self::El, rhs: &Self::El) -> Self::El {
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
                        take_mut::take(e.into_mut(), |c| self.base_ring.add(c, new_coeff));
                    },
                    Entry::Vacant(e) => {
                        e.insert(new_coeff);
                    }
                };
            }
        }
        return result;
    }

    fn neg(&self, mut val: Self::El) -> Self::El {
        for (_, coeff) in val.iter_mut() {
            take_mut::take(coeff, |c| self.base_ring.neg(c));
        }
        return val;
    }

    fn zero(&self) -> Self::El {
        HashMap::new()
    }

    fn one(&self) -> Self::El {
        self.from(self.base_ring.one())
    }

    fn eq(&self, lhs: &Self::El, rhs: &Self::El) -> bool {
        if lhs.len() != rhs.len() {
            return false;
        }
        for (key, lhs_coeff) in lhs.iter() {
            if let Some(rhs_coeff) = rhs.get(key) {
                if !self.base_ring.eq(lhs_coeff, rhs_coeff) {
                    return false;
                }
            } else {
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
    
    ///
    /// May panic if the ring is not euclidean. The first result is the quotient
    /// and the second result the remainder. The inequality
    ///  lhs = quo * rhs + rem
    /// must always hold.
    /// 
    fn euclidean_div_rem(&self, _lhs: Self::El, _rhs: Self::El) -> (Self::El, Self::El) {
        panic!("not euclidean")
    }

    ///
    /// May panic if the ring is not a field. If it does not panic, the result
    /// must be valid. For a non-field ring, it therefore must panic if rhs does not
    /// divide lhs, and if it does, it may either compute the correct quotient but may
    /// also panic nevertheless.
    /// 
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