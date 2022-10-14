use super::matrix_view::*;
use super::super::ring::*;

pub trait MatrixScale<M>: Ring 
    where M: MatrixViewMut<Self::El>
{
    fn scale_matrix(&self, l: &Self::El, a: &mut M);
    fn negate_matrix(&self, a: &mut M);
}

impl<R: Ring, M: MatrixViewMut<Self::El>> MatrixScale<M> for R {
    
    default fn scale_matrix(&self, l: &Self::El, a: &mut M)
    {
        for i in 0..a.row_count() {
            for j in 0..a.col_count() {
                *a.at_mut(i, j) = self.mul_ref(a.at(i, j), l);
            }
        }
    }

    default fn negate_matrix(&self, a: &mut M) {
        for i in 0..a.row_count() {
            for j in 0..a.col_count() {
                let value = std::mem::replace(a.at_mut(i, j), self.unspecified_element());
                *a.at_mut(i, j) = self.neg(value);
            }
        }
    }
}

pub trait MatrixEq<M, N>: Ring 
    where M: MatrixView<Self::El>, N: MatrixView<Self::El>
{
    fn eq_matrix(&self, a: M, b: N) -> bool;
}

impl<R: Ring, M: MatrixView<Self::El>, N: MatrixView<Self::El>> MatrixEq<M, N> for R {

    default fn eq_matrix(&self, a: M, b: N) -> bool {
        assert_eq!(a.row_count(), b.row_count());
        assert_eq!(a.col_count(), b.col_count());
        for row in 0..a.row_count() {
            for col in 0..a.col_count() {
                if !self.is_eq(a.at(row, col), b.at(row, col)) {
                    return false;
                }
            }
        }
        return true;
    }
}

pub trait MatrixFrobenius<M>: Ring
    where M: MatrixView<Self::El>
{
    fn calc_matrix_frobenius_norm_square(&self, a: M) -> Self::El;
}

impl<R, M> MatrixFrobenius<M> for R 
    where R: Ring, M: MatrixView<R::El>
{
    default fn calc_matrix_frobenius_norm_square(&self, a: M) -> Self::El
    {
        let mut result = self.zero();
        for i in 0..a.row_count() {
            for j in 0..a.col_count() {
                result = self.add(result, self.mul_ref(a.at(i, j), a.at(i, j)));
            }
        }
        return result;
    }
}
