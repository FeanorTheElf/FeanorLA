use super::super::prelude::*;
use super::uni_var::*;
use super::*;

use std::ops::Index;
use std::collections::BTreeMap;
use std::collections::btree_map::Entry;

#[derive(Debug, Clone)]
pub struct MultivariatePolyRing<R>
    where R: Ring
{
    base_ring: R,
    var_names: Vec<&'static str>
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Var(usize);

impl<R> MultivariatePolyRing<R>
    where R: Ring
{
    fn truncate_element_zeros(el: &mut Vec<usize>) {
        el.truncate(
            el.len() - el.iter().rev().take_while(|x| **x == 0).count()
        );
    }

    pub fn get_var(&self, name: &'static str) -> Var {
        Var(self.var_names.iter().enumerate().filter(|(_, x)| **x == name).next().unwrap().0)
    }

    pub fn get_indeterminate(&self, Var(index): Var) -> El<Self> {
        let mut result = BTreeMap::new();
        result.insert(
            (0..index).map(|_| 0).chain(std::iter::once(1)).collect(), 
            self.base_ring.one()
        );
        return result;
    }

    pub fn derive(&self, el: &El<Self>, variable: Var) -> El<Self> {
        let var = variable.0;
        let mut result = BTreeMap::new();
        for (key, coeff) in el {
            if var < key.len() && key[var] > 0 {
                let mut new_key = key.clone();
                let new_coeff = self.base_ring.mul_ref(&self.base_ring.from_z(key[var] as i64), coeff);
                new_key[var] -= 1;
                Self::truncate_element_zeros(&mut new_key);
                result.insert(new_key, new_coeff);
            }
        }
        return result;
    }

    pub fn gradient(&self, el: &El<Self>) -> Vector<VectorOwned<El<Self>>, El<Self>> {
        Vector::from_fn(self.var_names.len(), |i| self.derive(el, Var(i)))
    }

    pub fn elevate_var<'a>(&'a self, var: Var) -> (
        PolyRingImpl<&'a MultivariatePolyRing<R>>, 
        impl Fn(El<Self>) -> El<PolyRingImpl<&'a MultivariatePolyRing<R>>>,
        impl Fn(El<PolyRingImpl<&'a MultivariatePolyRing<R>>>) -> El<Self>
    ) {
        (
            self.elevate_var_ring(var),
            move |x| self.elevate_var_element(var, x),
            move |x| self.de_elevate_var(var, x)
        )
    }

    fn elevate_var_ring(&self, var: Var) -> PolyRingImpl<&MultivariatePolyRing<R>> {
        PolyRingImpl::adjoint(self, self.var_names[var.0])
    }

    fn elevate_var_element(&self, variable: Var, x: El<Self>) -> El<PolyRingImpl<&MultivariatePolyRing<R>>> {
        self.assert_valid(&x);

        let var = variable.0;
        let mut result = Vec::new();
        for (mut key, coeff) in x.into_iter() {
            let pow = *key.get(var).unwrap_or(&0);
            if var < key.len() {
                key[var] = 0;
                Self::truncate_element_zeros(&mut key);
            }
            result.resize_with(result.len().max(pow + 1), || self.zero());
            let v = &mut result[pow];
            debug_assert!(v.get(&key).is_none());
            v.insert(key, coeff);
        }
        for coeff in &result {
            self.assert_valid(coeff);
        }
        return Vector::new(result);
    }

    fn de_elevate_var(&self, variable: Var, x: El<PolyRingImpl<&MultivariatePolyRing<R>>>) -> El<Self> {
        let var = variable.0;
        let mut result = BTreeMap::new();
        for (pow, coeff) in x.raw_data().into_iter().enumerate() {
            result.extend(coeff.into_iter().map(|(mut key, c)| {
                if pow > 0 {
                    key.resize(key.len().max(var + 1), 0);
                    key[var] = pow;
                }
                return (key, c);
            }));
        }
        
        self.assert_valid(&result);
        return result;
    }

    fn assert_valid(&self, el: &El<Self>) {
        //
        // We require that the power vector has no trailing zeros, as otherwise this
        // would allow extremely blown-up polynomials (like n * X = 1 * X + ... + 1 * X) 
        // and makes code simpler.
        // Note that we do not require that the coefficients are non-zero, so the following
        // are three perfectly valid representations of zero: {} and {[]: 0} and 
        // {[]: 0, [0, 1]: 0}
        //
        for key in el.keys() {
            debug_assert!(key.last().is_none() || *key.last().unwrap() != 0);
            debug_assert!(key.len() <= self.var_names.len());
        }
    }

    fn nonzero_monomials<'a>(&'a self, el: &'a El<Self>) -> impl 'a + Iterator<Item = (&'a Vec<usize>, &'a El<R>)> {
        self.assert_valid(el);
        el.iter().filter(move |(_, coeff)| !self.base_ring.is_zero(coeff))
    }

    pub fn as_constant<'a>(&self, el: &El<Self>) -> Option<R::El> {
        self.assert_valid(el);
        if self.is_zero(el) {
            return Some(self.base_ring.zero());
        } else {
            if self.nonzero_monomials(el).any(|(key, _)| *key != Vec::new()) {
                return None;
            } else {
                return  Some(el.get(&Vec::new()).cloned().unwrap_or_else(|| self.base_ring.zero()));
            }
        }
    }

    pub fn new(base_ring: R) -> Self {
        MultivariatePolyRing {
            base_ring: base_ring,
            var_names: Vec::new()
        }
    }

    pub fn adjoint(&mut self, var: &'static str) -> El<Self> {
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
    /// of this polynomial. This algorithm is optimized for base rings in which
    /// multiplication is expensive.
    /// 
    /// # Complexity
    /// 
    /// Complexity is O(V * N + V * P * T) where V is the number of variables, N is
    /// the number of monomials in the given polynomial, P is the highest power of a
    /// variable and T is the complexity of a multiplication in the base ring.
    /// 
    pub fn evaluate_at<V>(&self, mut poly: El<Self>, values: &V) -> R::El 
        where V: Index<usize, Output = R::El>
    {
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
                        self.base_ring.add_assign(e.into_mut(), added_value);
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

    pub fn evaluate_matrix_at<M, V>(&self, m: Matrix<M, El<Self>>, vars: &V) -> Matrix<MatrixOwned<R::El>, R::El> 
        where M: MatrixView<El<Self>>, V: Index<usize, Output = R::El>
    {
        Matrix::from_fn(m.row_count(), m.col_count(), |i, j| self.evaluate_at(m.at(i, j).clone(), vars))
    }

    pub fn evaluate_vector_at<M, V>(&self, m: Vector<M, El<Self>>, vars: &V) -> Vector<VectorOwned<R::El>, R::El> 
        where M: VectorView<El<Self>>, V: Index<usize, Output = R::El>
    {
        Vector::from_fn(m.len(), |i| self.evaluate_at(m.at(i).clone(), vars))
    }
}

impl<R> RingBase for MultivariatePolyRing<R>
    where R: Ring
{
    ///
    /// This contains all monomials as the power vector of the variables and
    /// the corresonding coefficients.
    /// For invariants, see `assert_valid()`.
    /// 
    type El = BTreeMap<Vec<usize>, R::El>;

    fn add_ref(&self, mut lhs: Self::El, rhs: &Self::El) -> Self::El {
        self.assert_valid(&lhs);
        self.assert_valid(&rhs);

        for entry in rhs.iter() {
            if let Some(coeff) = lhs.get_mut(entry.0) {
                self.base_ring.add_assign(coeff, entry.1.clone());
            } else {
                lhs.insert(entry.0.clone(), entry.1.clone());
            }
        }

        self.assert_valid(&lhs);
        return lhs;
    }

    ///
    /// # Complexity
    /// 
    /// The complexity is O(N^2 * V + N^2 * T) where N is the number of monomials
    /// in lhs resp. rhs, V is the count of variables and T is the complexity of
    /// a multiplication in the base ring of two coefficients.
    /// 
    fn mul_ref(&self, lhs: &Self::El, rhs: &Self::El) -> Self::El {
        self.assert_valid(&lhs);
        self.assert_valid(&rhs);

        let mut result = self.zero();
        for (lhs_vars, lhs_coeff) in lhs.iter() {
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
                        self.base_ring.add_assign(e.into_mut(), new_coeff);
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
            let value = std::mem::replace(coeff, self.base_ring.unspecified_element());
            *coeff = self.base_ring.neg(value);
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

    fn is_eq(&self, lhs: &Self::El, rhs: &Self::El) -> bool {
        self.assert_valid(&lhs);
        self.assert_valid(&rhs);

        let cmp = |a: &Self::El, b: &Self::El| {
            for (key, lhs_coeff) in a.iter() {
                if let Some(rhs_coeff) = b.get(key) {
                    if !self.base_ring.is_eq(lhs_coeff, rhs_coeff) {
                        return false;
                    }
                } else {
                    if !self.base_ring.is_zero(lhs_coeff) {
                        return false;
                    }       
                }
            }
            return true;
        };

        return cmp(lhs, rhs) && cmp(rhs, lhs);
    }

    fn from_z(&self, x: i64) -> Self::El {
        self.from(self.base_ring.from_z(x))
    }

    fn from_z_big(&self, x: &BigInt) -> Self::El {
        self.from(self.base_ring.from_z_big(x))
    }

    fn characteristic(&self) -> BigInt {
        self.base_ring.characteristic()
    }

    fn is_integral(&self) -> RingPropValue {
        self.base_ring.is_integral()
    }

    fn is_noetherian(&self) -> bool {
        self.base_ring.is_noetherian()
    }

    fn is_field(&self) -> RingPropValue {
        if self.var_names.len() == 0 {
            self.base_ring.is_field()
        } else {
            RingPropValue::False
        }
    }
    
    fn div(&self, _lhs: Self::El, _rhs: &Self::El) -> Self::El {
        panic!("Not a field!")
    }

    fn format(&self, el: &El<Self>, f: &mut std::fmt::Formatter, in_prod: bool) -> std::fmt::Result {
        self.assert_valid(el);
        if self.is_zero(el) {
            return self.base_ring.format(&self.base_ring.zero(), f, in_prod);
        } else if let Some(c) = self.as_constant(el) {
            return self.base_ring.format(&c, f, in_prod);
        } else if in_prod && self.nonzero_monomials(el).next().is_some() && !self.base_ring.is_one(&self.nonzero_monomials(el).next().unwrap().1) {
            return self.format_in_brackets(el, f);
        } else if in_prod && self.nonzero_monomials(el).skip(1).next().is_some() {
            return self.format_in_brackets(el, f);
        } else {
            let mut monomial_it = self.nonzero_monomials(el);
            let first = monomial_it.next().unwrap();

            let write_part = |key: &Vec<usize>, coeff, f: &mut std::fmt::Formatter| {
                if !self.base_ring.is_one(coeff) {
                    self.base_ring.format(coeff, f, true)?;
                    if key.len() > 0 {
                        write!(f, " * ")?;
                    }
                }
                if key.len() > 0 {
                    for (i, pow) in key.iter().enumerate() {
                        if *pow == 1 {
                            write!(f, "{}", self.var_names[i])?;
                        } else if *pow > 1 {
                            write!(f, "{}^{}", self.var_names[i], pow)?;
                        }
                    }
                }
                return Ok(());
            };
            
            write_part(first.0, first.1, f)?;
            for el in monomial_it {
                write!(f, " + ")?;
                write_part(el.0, el.1, f)?;
            }
            return Ok(());
        }
    }
}

impl<R: Ring> RingExtension for MultivariatePolyRing<R> {
    
    type BaseRing = R;
    type Embedding = StandardEmbedding<R, Self>;

    fn base_ring(&self) -> &R {
        &self.base_ring
    }

    fn embedding(&self) -> Self::Embedding {
        embedding(self.base_ring().clone(), self.clone())
    }

    fn from(&self, el: El<R>) -> El<Self> {
        let mut result = BTreeMap::new();
        result.insert(Vec::new(), el);
        self.assert_valid(&result);
        return result;
    }
}

impl<R> DivisibilityInfoRing for MultivariatePolyRing<R> 
    where R: DivisibilityInfoRing
{
    fn is_divisibility_computable(&self) -> RingPropValue {
        self.base_ring.is_divisibility_computable()
    }

    fn quotient(&self, lhs: &Self::El, rhs: &Self::El) -> Option<Self::El> {
        assert!(!self.is_zero(rhs));
        if let Some(division_var) = self.nonzero_monomials(rhs)
            .filter_map(|(key, _coeff)| 
                key.iter().enumerate().filter(|(_i, pow)| **pow != 0).map(|(i, _pow)| i).next()
            ).min() 
        {
            let ring = self.elevate_var_ring(Var(division_var));
            let lhs_new = self.elevate_var_element(Var(division_var), lhs.clone());
            let rhs_new = self.elevate_var_element(Var(division_var), rhs.clone());
            let result = ring.quotient(&lhs_new, &rhs_new)?;
            return Some(self.de_elevate_var(Var(division_var), result));
        } else {
            let mut result = lhs.clone();
            // rhs is only a scalar
            let (key, scalar) = self.nonzero_monomials(rhs).next().unwrap();
            debug_assert_eq!(Vec::<usize>::new(), *key);
            for coeff in result.values_mut() {
                let new_coeff = self.base_ring.quotient(&coeff, &scalar)?;
                *coeff = new_coeff;
            }
            return Some(result);
        }
    }
}

impl<R> CanonicalEmbeddingInfo<R> for MultivariatePolyRing<R>
    where R: Ring
{
    fn has_embedding(&self, from: &R) -> RingPropValue {
        self.base_ring.has_embedding(from)
    }

    fn embed(&self, from: &R, el: R::El) -> Self::El {
        assert!(self.has_embedding(from).can_use());
        self.from(self.base_ring.embed(from, el))
    }
}

impl<R> CanonicalEmbeddingInfo<R> for MultivariatePolyRing<&R> 
    where R: Ring
{

    fn has_embedding(&self, from: &R) -> RingPropValue {
        self.base_ring.has_embedding(&from)
    }

    fn embed(&self, from: &R, el: R::El) -> Self::El {
        assert!(self.has_embedding(from).can_use());
        self.from(self.base_ring.embed(&from, el))
    }
}

impl<R> CanonicalEmbeddingInfo<&R> for MultivariatePolyRing<R> 
    where R: Ring
{

    fn has_embedding(&self, from: &&R) -> RingPropValue {
        self.base_ring.has_embedding(*from)
    }

    fn embed(&self, from: &&R, el: R::El) -> Self::El {
        assert!(self.has_embedding(from).can_use());
        self.from(self.base_ring.embed(*from, el))
    }
}

impl<R> CanonicalEmbeddingInfo<PolyRingImpl<R>> for MultivariatePolyRing<R>
    where R: CanonicalIsomorphismInfo<R>
{
    fn has_embedding(&self, _from: &PolyRingImpl<R>) -> RingPropValue {
        // It might seem natural to have the canonical embedding f(X) -> f(X),
        // but this would contradict the general commutativity of canonical embeddings;
        // The problem here is that we have defined k[X] and k[Y] to be canonically
        // isomorphic as univariate polynomial rings, but not as multivariate polynomial
        // rings.
        RingPropValue::False
    }

    fn embed(&self, _from: &PolyRingImpl<R>, _el: El<PolyRingImpl<R>>) -> Self::El {
        panic!("No embedding defined!")
    }
}

impl<R> CanonicalEmbeddingInfo<MultivariatePolyRing<R>> for MultivariatePolyRing<R>
    where R: Ring
{
    fn has_embedding(&self, from: &MultivariatePolyRing<R>) -> RingPropValue {
        self.base_ring.has_embedding(&from.base_ring) & if from.var_names.iter().all(|name| self.var_names.contains(name)) {
            RingPropValue::True
        } else {
            RingPropValue::False
        }
    }

    fn embed(&self, from: &MultivariatePolyRing<R>, el: El<Self>) -> El<Self> {
        assert!(self.has_embedding(from).can_use());
        let var_mapping = from.var_names.iter().map(|v| 
            self.var_names.iter().enumerate().find(|(_, v2)| v == *v2).unwrap().0
        ).collect::<Vec<_>>();
        let mut key_recycling = Vec::new();
        let result = el.into_iter()
            .map(|(key, val)| (key, self.base_ring.embed(&from.base_ring, val)))
            .map(|(mut key, val)| {
                if key.len() == 0 {
                    return (Vec::new(), val);
                }
                key_recycling.clear();
                key_recycling.resize((0..key.len()).map(|i| var_mapping[i] + 1).max().unwrap(), 0);
                for (i, power) in key.iter().enumerate() {
                    key_recycling[var_mapping[i]] = *power;
                }
                std::mem::swap(&mut key_recycling, &mut key);
                return (key, val);
            }).collect();
        self.assert_valid(&result);
        return result;
    }
}

impl<R> CanonicalIsomorphismInfo<MultivariatePolyRing<R>> for MultivariatePolyRing<R>
    where R: Ring
{
    fn has_isomorphism(&self, from: &MultivariatePolyRing<R>) -> RingPropValue {
        self.has_embedding(from) & (self.var_names.len() == from.var_names.len())
    }

    fn preimage(&self, from: &MultivariatePolyRing<R>, el: Self::El) -> Self::El {
        from.embed(self, el)
    }
}

#[cfg(test)]
use super::super::wrapper::*;

#[test]
fn test_binomial_formula() {
    let mut ring = MultivariatePolyRing::new(i32::RING);
    let x = ring.adjoint("X");
    let y = ring.adjoint("Y");
    let x = ring.bind(x);
    let y = ring.bind(y);
    let two = ring.bind(ring.from(2));
    
    let x2_2xy_y2 = (&x * &x) + (two * &x * &y) + (&y * &y);
    let x_y_square = (&x + &y) * (&x + &y);

    assert_eq!(x2_2xy_y2, x_y_square);
    assert!(x2_2xy_y2 != x);
}

#[test]
fn test_eq() {
    let mut ring = MultivariatePolyRing::new(i32::RING);
    let one = ring.from(1);
    let _x = ring.adjoint("x");
    let one_prime = ring.from(1);
    assert!(ring.is_eq(&one, &one_prime));
}

#[test]
fn test_evaluate_at() {
    let mut ring = MultivariatePolyRing::new(i32::RING);
    let x = ring.adjoint("X");
    let y = ring.adjoint("Y");
    let x = ring.bind(x);
    let y = ring.bind(y);
    let thirteen = y.ring().from_z(13);

    // the polynomial x^2 y + 2 x y^2 + x + 13
    let poly = (&x * &x * &y) + (&x * &y * &y + &x * &y * &y) + &x + thirteen;

    assert_eq!(14, ring.evaluate_at(poly.val().clone(), &vec![1, 0]));
    assert_eq!(12 + 36 + 2 + 13, ring.evaluate_at(poly.val().clone(), &vec![2, 3]));
}

#[test]
fn test_assumption_option_ord() {
    // in PolyRing::poly_div we rely on this behavior, 
    // and it is useful in general for degree comparison
    assert!(None < Some(0));
    assert!(Some(0) < Some(1));
}

#[test]
fn test_format_multivar_poly_ring() {
    let mut ring = MultivariatePolyRing::new(i32::RING);
    let x = ring.adjoint("X");
    let y = ring.adjoint("Y");
    let x = ring.bind(x);
    let y = ring.bind(y);
    let one = y.ring().from_z(1);

    let poly = &x * &x * &x - &y + (&one + &one) * &y * &x - &one;
    assert_eq!("-1 + -1 * Y + 2 * XY + X^3", format!("{}", poly));
}

#[test]
fn test_div_multivar_poly_ring() {
    let mut ring = MultivariatePolyRing::new(i32::RING);
    let x = ring.adjoint("X");
    let y = ring.adjoint("Y");
    let x = ring.bind(x);
    let y = ring.bind(y);

    let a = &x + &y;
    let b = &x * &x * 2;
    let c = &y + 1;
    let d = (&x + &y + 1) * (&x - 1);
    let p = &a * &b * &c * &d;
    let q = &a * &c;
    let expected = &b * &d;
    let result = ring.quotient(p.val(), q.val()).unwrap();
    assert_eq!(expected, ring.bind(result));
}

#[test]
fn test_elevate_var() {
    let mut ring = MultivariatePolyRing::new(i32::RING);
    let x = ring.adjoint("X");
    let y = ring.adjoint("Y");
    let x = ring.bind(x);
    let y = ring.bind(y);

    let p = &x * &y + &y * &y * &x * 2 + 1 + &x;

    let uni_ring = ring.elevate_var_ring(Var(1));

    let uni_y = uni_ring.bind(uni_ring.unknown());
    let uni_x = uni_ring.bind(uni_ring.from(x.val().clone()));
    let uni_one = uni_ring.bind(uni_ring.one());

    let expected = &uni_x * &uni_y + &uni_y * &uni_y * &uni_x * (&uni_one + &uni_one) + &uni_one + &uni_x;

    let actual = uni_ring.bind(ring.elevate_var_element(Var(1), p.val().clone()));
    assert_eq!(expected, actual);

    let original = ring.bind(ring.de_elevate_var(Var(1), actual.val().clone()));
    assert_eq!(p, original);
}

#[test]
fn test_gradient() {
    let mut ring = MultivariatePolyRing::new(i32::RING);
    let x = ring.adjoint("X");
    let y = ring.adjoint("Y");
    let x = ring.bind(x);
    let y = ring.bind(y);

    let p = &x + &x * &x * &y + &x * &y + &x * &x * &x * 2;
    let dx = (&y + 1) + &x * &y * 2 + &x * &x * 6;
    let dy = &x * &x + &x;

    assert_eq!(dx, ring.bind(ring.derive(p.val(), ring.get_var("X"))));
    assert_eq!(dy, ring.bind(ring.derive(p.val(), ring.get_var("Y"))));
}

#[test]
fn test_quotient() {
    let mut ring = MultivariatePolyRing::new(i32::RING);
    let x = ring.adjoint("X");
    ring.adjoint("Y");
    let x = ring.bind(x);

    assert_eq!(None, ring.quotient(&ring.one(), (&x + 1).val()));
}

#[test]
fn test_is_unit() {
    let mut ring = MultivariatePolyRing::new(i32::RING);
    let x = ring.adjoint("X");
    ring.adjoint("Y");
    let x = ring.bind(x);

    assert_eq!(false, ring.is_unit(&(x + 1).val()));
    assert_eq!(false, ring.is_unit(&ring.from_z(2)));
    assert_eq!(true, ring.is_unit(&ring.from_z(-1)));
}

#[test]
fn test_poly_ring_embedding() {
    let mut ring1 = MultivariatePolyRing::new(i32::RING);
    let x1 = ring1.adjoint("X");
    let y1 = ring1.adjoint("Y");

    let mut ring2 = MultivariatePolyRing::new(i32::RING);
    let y2 = ring2.adjoint("Y");
    let x2 = ring2.adjoint("X");

    let (f, fi) = isomorphism(&ring1, &ring2);
    assert!(!ring1.is_eq(&x1, &x2));
    assert!(ring2.is_eq(&x2, &f(x1)));
    assert!(ring1.is_eq(&y1, &fi(y2)));
}