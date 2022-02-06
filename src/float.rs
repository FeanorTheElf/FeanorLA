use super::primitive::FieldEl;
use super::la::mat::*;

pub trait ApproxEq<Rhs = Self, D = Self> {
    ///
    /// Checks whether self is approximately equal to rhs, with
    /// a relative or absolute error of at most delta.
    /// 
    /// # Details
    /// 
    /// For values far away from 0, the relative error should be 
    /// bounded by delta, i.e. 2 |self - rhs| / (|self| + |rhs|).
    /// For values close to 0, this does not work as the quotient
    /// might become very small in absolute value. In this case, the
    /// absolute error must be bounded instead.
    /// 
    /// Note that therefore `a.approx_eq(b, d)` yields a much more
    /// precise result as `(a - b).approx_eq(0, d)`, as in the first
    /// case, the relative error may be used which respects the size
    /// of a and b in the comparision, while the second must yield the
    /// same values for `2.approx_eq(1, d)` and for 
    /// `100001.approx_eq(100000)`. Note that because of this, comparing
    /// greater values does in principle not require a bigger delta.
    /// However, if those greater values are the result of longer
    /// computations, then clearly, this is the case.
    /// 
    fn approx_eq(&self, rhs: &Rhs, delta: &D) -> bool;
}

impl ApproxEq for f64 {

    fn approx_eq(&self, rhs: &Self, delta: &f64) -> bool {
        if (self.abs() + rhs.abs()) < 100. * *delta {
            (self - rhs).abs() < *delta
        } else {
            (self - rhs).abs() / (self.abs() + rhs.abs()) < *delta
        }
    }
}

impl ApproxEq for f32 {

    fn approx_eq(&self, rhs: &Self, delta: &f32) -> bool {
        if (self.abs() + rhs.abs()) < 100. * *delta {
            (self - rhs).abs() < *delta
        } else {
            (self - rhs).abs() / (self.abs() + rhs.abs()) < *delta
        }
    }
}

impl<M, N, T, D> ApproxEq<Matrix<N, T>, D> for Matrix<M, T>
    where M: MatrixView<T>, N: MatrixView<T>, T: ApproxEq<T, D>
{
    fn approx_eq(&self, rhs: &Matrix<N, T>, delta: &D) -> bool {
        for row in 0..self.row_count() {
            for col in 0..self.col_count() {
                if !self.at(row, col).approx_eq(rhs.at(row, col), delta) {
                    return false;
                }
            }
        }
        return true;
    }
}

impl<V, W, T, D> ApproxEq<Vector<W, T>, D> for Vector<V, T>
    where V: VectorView<T>, W: VectorView<T>, T: ApproxEq<T, D>
{
    fn approx_eq(&self, rhs: &Vector<W, T>, delta: &D) -> bool {
        for i in 0..self.len() {
            if !self.at(i).approx_eq(rhs.at(i), delta) {
                return false;
            }
        }
        return true;
    }
}

///
/// For types that may have rounding errors
/// 
pub trait Float: Clone + ApproxEq + PartialEq + PartialOrd + FieldEl {

    fn sqrt(&self) -> Self;
}

impl Float for f32 {

    fn sqrt(&self) -> Self {
        (*self).sqrt()
    }
}

impl Float for f64 {

    fn sqrt(&self) -> Self {
        (*self).sqrt()
    }
}
