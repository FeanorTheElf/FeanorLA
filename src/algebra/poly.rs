use super::super::alg::*;
use super::super::la::mat::*;

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

    pub fn derive(&self, el: &<Self as Ring>::El, variable: Var) -> <Self as Ring>::El {
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

    pub fn gradient(&self, el: &<Self as Ring>::El) -> Vector<VectorOwned<<Self as Ring>::El>, <Self as Ring>::El> {
        Vector::from_fn(self.var_names.len(), |i| self.derive(el, Var(i)))
    }

    fn elevate_var_ring(&self, var: Var) -> PolyRing<&MultivariatePolyRing<R>> {
        PolyRing::adjoint(self, self.var_names[var.0])
    }

    fn elevate_var(&self, variable: Var, x: <Self as Ring>::El) -> <PolyRing<&MultivariatePolyRing<R>> as Ring>::El
    {
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
        return result;
    }

    fn de_elevate_var(&self, variable: Var, x: <PolyRing<&MultivariatePolyRing<R>> as Ring>::El) -> <Self as Ring>::El {
        let var = variable.0;
        let mut result = BTreeMap::new();
        for (pow, coeff) in x.into_iter().enumerate() {
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

    pub fn from(&self, el: R::El) -> <Self as Ring>::El {
        let mut result = BTreeMap::new();
        result.insert(Vec::new(), el);

        self.assert_valid(&result);
        return result;
    }

    fn assert_valid(&self, el: &<Self as Ring>::El) {
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
    /// of this polynomial. This algorithm is optimized for base rings in which
    /// multiplication is expensive.
    /// 
    /// # Complexity
    /// 
    /// Complexity is O(V * N + V * P * T) where V is the number of variables, N is
    /// the number of monomials in the given polynomial, P is the highest power of a
    /// variable and T is the complexity of a multiplication in the base ring.
    /// 
    pub fn evaluate_at<V>(&self, mut poly: <Self as Ring>::El, values: &V) -> R::El 
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

    pub fn evaluate_matrix_at<M, V>(&self, m: Matrix<M, <Self as Ring>::El>, vars: &V) -> Matrix<MatrixOwned<R::El>, R::El> 
        where M: MatrixView<<Self as Ring>::El>, V: Index<usize, Output = R::El>
    {
        Matrix::from_fn(m.row_count(), m.col_count(), |i, j| self.evaluate_at(m.at(i, j).clone(), vars))
    }

    pub fn evaluate_vector_at<M, V>(&self, m: Vector<M, <Self as Ring>::El>, vars: &V) -> Vector<VectorOwned<R::El>, R::El> 
        where M: VectorView<<Self as Ring>::El>, V: Index<usize, Output = R::El>
    {
        Vector::from_fn(m.len(), |i| self.evaluate_at(m.at(i).clone(), vars))
    }
}

impl<R> Ring for MultivariatePolyRing<R>
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

        let cmp = |a: &Self::El, b: &Self::El| {
            for (key, lhs_coeff) in a.iter() {
                if let Some(rhs_coeff) = b.get(key) {
                    if !self.base_ring.eq(lhs_coeff, rhs_coeff) {
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
        panic!("Not euclidean!")
    }

    fn div(&self, mut lhs: Self::El, rhs: &Self::El) -> Self::El {
        assert!(!self.is_zero(rhs));
        if let Some(division_var) = rhs.iter()
            .filter_map(|(key, _coeff)| 
                key.iter().enumerate().filter(|(_i, pow)| **pow != 0).map(|(i, _pow)| i).next()
            ).min() 
        {
            let ring = self.elevate_var_ring(Var(division_var));
            let lhs_new = self.elevate_var(Var(division_var), lhs);
            let rhs_new = self.elevate_var(Var(division_var), rhs.clone());
            let result = ring.div(lhs_new, &rhs_new);
            return self.de_elevate_var(Var(division_var), result);
        } else {
            // rhs is only a scalar
            debug_assert!(rhs.len() == 1);
            let (key, scalar) = rhs.iter().next().unwrap();
            debug_assert_eq!(Vec::<usize>::new(), *key);
            for coeff in lhs.values_mut() {
                take_mut::take_or_recover(
                    coeff, 
                    || self.base_ring.unspecified_element(), 
                    |v| self.base_ring.div(v, &scalar)
                );
            }
            return lhs;
        }
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
                write!(f, " + ")?;
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
    pub fn derive(&self, el: &<Self as Ring>::El) -> <Self as Ring>::El {
        let mut result = Vec::with_capacity(el.len());
        for i in 1..el.len() {
            result.push(self.base_ring.mul_ref(&self.base_ring.from_z(i as i64), &el[i]));
        }
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
        result.resize_with(self.deg(lhs).unwrap_or(0), || self.base_ring.zero());
        while self.deg(lhs) >= self.deg(rhs) {
            let lhs_deg = self.deg(lhs).unwrap();
            let coeff = div_lc(&lhs[lhs_deg])?;
            let pow = lhs_deg - rhs_deg;
            result[pow] = coeff;
            for i in 0..=rhs_deg {
                take_mut::take_or_recover(
                    &mut lhs[i + pow], 
                    || self.base_ring.unspecified_element(), 
                    |v| self.base_ring.sub(v, self.base_ring.mul_ref(&result[pow], &rhs[i]))
                );
            }
            if !self.base_ring.is_zero(&lhs[lhs_deg]) {
                panic!("Passed division function yielded the wrong result!");
            }
        }
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
        panic!("Not a euclidean domain!")
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

#[cfg(test)]
use super::super::alg_env::*;

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

    assert_eq!(14, ring.evaluate_at(poly.clone(), &vec![1, 0]));
    assert_eq!(12 + 36 + 2 + 13, ring.evaluate_at(poly, &vec![2, 3]));
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
    assert_eq!("1 * X^3 + 2 * X^2 + -1", format!("{}", display_ring_el(&ring, &poly)));
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
    assert_eq!("-1 + -1 * Y^1 + 2 * X^1 Y^1 + 1 * X^3", format!("{}", display_ring_el(&ring, &poly)));
}

#[test]
fn test_poly_div() {
    let ring = PolyRing::adjoint(StaticRing::<i32>::RING, "X");
    let x = ring.unknown();
    let one = ring.one();

    let mut p = fixed_ring_env!{ &ring; x, one; {
        x * x * x + x * x + x + one
    }};
    let q = fixed_ring_env!{ &ring; x, one; {
        x + one
    }};
    let expected = fixed_ring_env!{ &ring; x, one; {
        x * x + one
    }};
    let result = ring.poly_division(&mut p, &q, |x| Ok(*x)).unwrap();
    println!("{} ?= 0", display_ring_el(&ring, &p));
    assert!(ring.is_zero(&p));
    println!("{} ?= {}", display_ring_el(&ring, &expected), display_ring_el(&ring, &result));
    assert!(ring.eq(&expected, &result));
}

#[test]
fn test_div_multivar_poly_ring() {
    let mut ring = MultivariatePolyRing::new(StaticRing::<i32>::RING);
    let x = ring.adjoint("X");
    let y = ring.adjoint("Y");
    let one = ring.one();

    let a = fixed_ring_env!{ &ring; x, y; {
        x + y
    }};
    let b = fixed_ring_env!{ &ring; x, one; {
        (one + one) * x * x
    }};
    let c = fixed_ring_env!{ &ring; y, one; {
        y + one
    }};
    let d = fixed_ring_env!{ &ring; x, y, one; {
        (x + y + one) * (x - one)
    }};
    let p = fixed_ring_env!{ &ring; a, b, c, d; {
        a * b * c * d
    }};
    let q = fixed_ring_env!{ &ring; a, c; {
        a * c
    }};
    let expected = fixed_ring_env!{ &ring; b, d; {
        b * d
    }};
    let result = ring.div(p, &q);
    println!("{} ?= {}", display_ring_el(&ring, &expected), display_ring_el(&ring, &result));
    assert!(ring.eq(&expected, &result));
}

#[test]
fn test_elevate_var() {
    let mut ring = MultivariatePolyRing::new(StaticRing::<i32>::RING);
    let x = ring.adjoint("X");
    let y = ring.adjoint("Y");
    let one = ring.one();

    let p: BTreeMap<Vec<usize>, <StaticRing::<i32> as Ring>::El> = fixed_ring_env!{ &ring; x, y, one; {
        x * y + y * y * x * (one + one) + one + x
    }};

    let uni_ring = ring.elevate_var_ring(Var(1));

    let uni_y = uni_ring.unknown();
    let uni_x = uni_ring.from(x);
    let uni_one = uni_ring.one();

    let expected = fixed_ring_env!{ &uni_ring; uni_y, uni_x, uni_one; {
        uni_x * uni_y + uni_y * uni_y * uni_x * (uni_one + uni_one) + uni_one + uni_x
    }};

    let actual = ring.elevate_var(Var(1), p.clone());
    println!("{} ?= {}", display_ring_el(&uni_ring, &expected), display_ring_el(&uni_ring, &actual));
    assert!(uni_ring.eq(&expected, &actual));

    let original = ring.de_elevate_var(Var(1), actual);
    println!("{} ?= {}", display_ring_el(&ring, &p), display_ring_el(&ring, &original));
    assert!(ring.eq(&p, &original));
}

#[test]
fn test_gradient() {
    let mut ring = MultivariatePolyRing::new(StaticRing::<i32>::RING);
    let x = ring.adjoint("X");
    let y = ring.adjoint("Y");
    let one = ring.one();

    let p = fixed_ring_env!{ &ring; x, y, one; {
        x + x * x * y + x * y + (one + one) * x * x * x
    }};
    let dx = fixed_ring_env!{ &ring; x, y, one; {
        (one + y) + (one + one) * x * y + (one + one) * (one + one + one) * x * x
    }};
    let dy = fixed_ring_env!{ &ring; x; {
        x * x + x
    }};
    assert_eq!(dx, ring.derive(&p, ring.get_var("X")));
    assert_eq!(dy, ring.derive(&p, ring.get_var("Y")));
}