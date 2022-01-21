use super::super::super::alg::*;
use super::super::super::la::vec::*;
use super::super::super::la::vector_view::*;

pub fn poly_degree<V, R>(coeff_ring: &R, poly_coeffs: Vector<V, R::El>) -> Option<usize>
    where R: Ring, V: VectorView<R::El>
{
    poly_coeffs.iter().enumerate().filter(|(_, x)| !coeff_ring.is_zero(x)).map(|(i, _)| i).next()
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
            take_mut::take_or_recover(
                lhs.at_mut(i + pow), 
                || coeff_ring.unspecified_element(), 
                |v| coeff_ring.sub(v, coeff_ring.mul_ref(&result[pow], &rhs[i]))
            );
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
    fn poly_add_assign<R, V>(coeff_ring: &R, mut base: Vector<VectorOwned<R::El>, R::El>, add: Vector<V, R::El>) -> Vector<VectorOwned<R::El>, R::El>
        where R: Ring, V: VectorView<R::El>
    {
        for i in 0..add.len() {
            take_mut::take_or_recover(
                base.at_mut(i), 
                || coeff_ring.unspecified_element(), 
                |v| coeff_ring.add_ref(v, add.at(i))
            );
        }
        return base;
    }
    if lhs.len() > rhs.len() {
        poly_add_assign(coeff_ring, lhs.into_owned(), rhs)
    } else {
        poly_add_assign(coeff_ring, rhs.into_owned(), lhs)
    }
}

pub fn poly_mul<R, V, W>(coeff_ring: &R, lhs: Vector<V, R::El>, rhs: Vector<W, R::El>) -> Vector<VectorOwned<R::El>, R::El>
    where R: Ring, V: VectorView<R::El>, W: VectorView<R::El>
{
    let mut result = Vec::with_capacity(lhs.len() + rhs.len());
    for i in 0..(lhs.len() + rhs.len()) {
        let mut val = coeff_ring.zero();
        for j in (rhs.len().max(i + 1) - rhs.len())..lhs.len().min(i + 1) {
            val = coeff_ring.add(val, coeff_ring.mul_ref(&lhs[j], &rhs[i - j]));
        }
        result.push(val);
    }
    return Vector::new(result);
}