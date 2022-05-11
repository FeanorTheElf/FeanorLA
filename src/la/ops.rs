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

pub trait MatrixMul<M, N>: Ring 
    where M: MatrixView<Self::El>, N: MatrixView<Self::El>
{
    fn mul_matrix(&self, a: M, b: N) -> MatrixOwned<Self::El>;
}

impl<R: Ring, M: MatrixView<Self::El>, N: MatrixView<Self::El>> MatrixMul<M, N> for R {

    default fn mul_matrix(&self, a: M, b: N) -> MatrixOwned<Self::El> {
        assert_eq!(a.col_count(), b.row_count());
        debug_assert!(a.col_count() > 0);
        MatrixOwned::from_fn(a.row_count(), b.col_count(), |i, j| {
            let mut it = (0..a.col_count()).map(|k| 
                self.mul_ref(a.at(i, k), b.at(k, j))
            );
            let initial = it.next().unwrap();
            it.fold(initial, |a, b| self.add(a, b))
        })
    }
}

pub trait MatrixAddAssign<M, N>: Ring 
    where M: MatrixViewMut<Self::El>, N: MatrixView<Self::El>
{
    fn add_assign_matrix(&self, a: &mut M, b: N);
    fn sub_assign_matrix(&self, a: &mut M, b: N);
}

impl<R: Ring, M: MatrixViewMut<Self::El>, N: MatrixView<Self::El>> MatrixAddAssign<M, N> for R {
    
    default fn add_assign_matrix(&self, a: &mut M, b: N) {
        assert_eq!(a.row_count(), b.row_count());
        assert_eq!(a.col_count(), b.col_count());
        for row in 0..a.row_count() {
            for col in 0..a.col_count() {
                self.add_assign(a.at_mut(row, col), b.at(row, col).clone());
            }
        }
    }

    default fn sub_assign_matrix(&self, a: &mut M, b: N) {
        assert_eq!(a.row_count(), b.row_count());
        assert_eq!(a.col_count(), b.col_count());
        for row in 0..a.row_count() {
            for col in 0..a.col_count() {
                self.add_assign(a.at_mut(row, col), self.neg(b.at(row, col).clone()));
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
                if !self.eq(a.at(row, col), b.at(row, col)) {
                    return false;
                }
            }
        }
        return true;
    }
}

pub trait MatrixAssign<M, N>: Clone 
    where M: MatrixViewMut<Self>, N: MatrixView<Self>
{
    fn assign_matrix(a: &mut M, b: N);
}

impl<T: Clone, M: MatrixViewMut<T>, N: MatrixView<T>> MatrixAssign<M, N> for T {
    
    default fn assign_matrix(a: &mut M, b: N) {
        assert_eq!(a.row_count(), b.row_count());
        assert_eq!(a.col_count(), b.col_count());
        for row in 0..a.row_count() {
            for col in 0..a.col_count() {
                *a.at_mut(row, col) = b.at(row, col).clone();
            }
        }
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
