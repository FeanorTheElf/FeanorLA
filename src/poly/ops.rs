use super::super::multiplication::karatsuba;
use super::super::prelude::*;
use super::super::la::vec::*;

pub fn poly_degree<V, R>(coeff_ring: &R, poly_coeffs: Vector<V, R::El>) -> Option<usize>
    where R: Ring, V: VectorView<R::El>
{
    poly_coeffs.iter().enumerate().filter(|(_, x)| !coeff_ring.is_zero(x)).map(|(i, _)| i).max()
}

pub fn poly_formal_derivative<R, V>(coeff_ring: &R, mut el: Vector<V, R::El>)
    where R: Ring, V: VectorViewMut<R::El>
{
    let mut tmp: R::El = coeff_ring.zero();
    if let Some(deg) = poly_degree(coeff_ring, el.as_ref()) {
        for i in (0..=deg).rev() {
            let old_coeff = std::mem::replace(el.at_mut(i), tmp);
            tmp = coeff_ring.mul(old_coeff, coeff_ring.from_z(i as i64));
        }
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
pub fn poly_division<V, W, R, F>(coeff_ring: &R, mut lhs: Vector<V, R::El>, rhs: Vector<W, R::El>, mut div_lc: F) -> Result<Vector<VectorOwned<R::El>, R::El>, ()>
    where R: Ring, V: VectorViewMut<R::El>, W: VectorView<R::El>, F: FnMut(&R::El) -> Result<R::El, ()>
{
    assert!(poly_degree(coeff_ring, rhs.as_ref()).is_some());

    let rhs_deg = poly_degree(coeff_ring, rhs.as_ref()).unwrap();
    let mut result = Vec::new();
    result.resize_with(poly_degree(coeff_ring, lhs.as_ref()).unwrap_or(0) + 1, || coeff_ring.zero());
    while poly_degree(coeff_ring, lhs.as_ref()) >= poly_degree(coeff_ring, rhs.as_ref()) {
        let lhs_deg = poly_degree(coeff_ring, lhs.as_ref()).unwrap();
        let coeff = div_lc(lhs.at(lhs_deg))?;
        let pow = lhs_deg - rhs_deg;
        result[pow] = coeff;
        for i in 0..=rhs_deg {
            coeff_ring.add_assign(lhs.at_mut(i + pow), coeff_ring.neg(coeff_ring.mul_ref(&result[pow], &rhs[i])));
        }
        if !coeff_ring.is_zero(&lhs[lhs_deg]) {
            panic!("Passed division function yielded the wrong result!");
        }
    }
    return Ok(Vector::new(result));
}

pub fn poly_add<R, V, W>(coeff_ring: &R, lhs: Vector<V, R::El>, rhs: Vector<W, R::El>) -> Vector<VectorOwned<R::El>, R::El>
    where R: Ring, V: VectorView<R::El>, W: VectorView<R::El>
{
    fn poly_add_assign<R, V, W>(coeff_ring: &R, mut base: Vector<W, R::El>, add: Vector<V, R::El>) -> Vector<W, R::El>
        where R: Ring, V: VectorView<R::El>, W: VectorViewMut<R::El>
    {
        for i in 0..add.len() {
            coeff_ring.add_assign(base.at_mut(i), add.at(i).clone());
        }
        return base;
    }
    if lhs.len() > rhs.len() {
        poly_add_assign(coeff_ring, lhs.into_owned(), rhs)
    } else {
        poly_add_assign(coeff_ring, rhs.into_owned(), lhs)
    }
}

pub fn poly_mul<R, U, V, W>(coeff_ring: &R, lhs: Vector<V, R::El>, rhs: Vector<W, R::El>, out: &mut Vector<U, R::El>) 
    where R: Ring, U: VectorViewMut<R::El>, V: VectorView<R::El>, W: VectorView<R::El>
{
    // do not just use lhs.len() and rhs.len(), as in some situations with huge amount of
    // multiplications, this can blow up the resulting vector terribly (e.g. when calculating powers) 
    let lhs_deg = poly_degree(coeff_ring, lhs.as_ref());
    let rhs_deg = poly_degree(coeff_ring, rhs.as_ref());
    if lhs_deg.is_none() || rhs_deg.is_none() {
        return;
    }
    let lhs_deg = lhs_deg.unwrap();
    let rhs_deg = rhs_deg.unwrap();
    // this is only the degree of result if the product of the leading coeffs is nonzero
    // e.g. if coeff_ring is integral
    // the code works anyway
    let result_deg = lhs_deg + rhs_deg;
    assert!(out.len() >= result_deg);
    unimplemented!()
}

pub fn poly_format<R, V>(coeff_ring: &R, el: Vector<V, R::El>, f: &mut std::fmt::Formatter, var_name: &str) -> std::fmt::Result
    where R: Ring, V: VectorView<R::El>
{
    if poly_degree(coeff_ring, el.as_ref()).is_none() {
        return coeff_ring.format(&coeff_ring.zero(), f, false);
    } else {
        let mut monomial_it = el.iter().enumerate().filter(|(_i, x)| !coeff_ring.is_zero(x));

        let print_monomial = |pow, coeff, formatter: &mut std::fmt::Formatter| {
            if !coeff_ring.is_one(coeff) || pow == 0 {
                coeff_ring.format(coeff, formatter, true)?;
                if pow > 0 {
                    write!(formatter, " * ")?;
                }
            }
            if pow == 1 {
                write!(formatter, "{}", var_name)?;
            } else if pow > 0 {
                write!(formatter, "{}^{}", var_name, pow)?;
            }
            return Ok(());
        };

        let (fst_pow, fst_coeff) = monomial_it.next().unwrap();
        print_monomial(fst_pow, fst_coeff, f)?;
        for (pow, coeff) in monomial_it {
            write!(f, " + ")?;
            print_monomial(pow, coeff, f)?;
        }
        return Ok(());
    }
}

pub fn poly_evaluate<R, S, V>(coeff_ring: &R, value: S::El, poly: Vector<V, R::El>, ring: &S) -> S::El
    where S: Ring + CanonicalEmbeddingInfo<R>, R: Ring, V: VectorView<R::El>
{
    let incl = embedding(coeff_ring, ring);
    let mut result = incl(poly.at(0).clone());
    let mut current_power = value.clone();
    for i in 1..poly.len() {
        result = ring.add(result, ring.mul(incl(poly.at(i).clone()), current_power.clone()));
        current_power = ring.mul(current_power, value.clone());
    }
    return result;
}

#[test]
fn test_poly_mul() {
    let mut a = Vector::from_array([1]);
    poly_mul(&i64::RING, Vector::from_array([1]).as_ref(), Vector::from_array([1]).as_ref(), &mut a);
    assert_eq!(Vector::from_array([2]), a);
}