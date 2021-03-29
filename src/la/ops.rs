use super::matrix_view::*;
use super::matrix_owned::*;
use super::super::alg::*;

pub trait MatrixScale<M>: Ring 
    where M: MatrixViewMut<Self::El>
{
    fn scale_matrix(&self, l: &Self::El, a: &mut M);
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
                take_mut::take_or_recover(
                    a.at_mut(row, col), 
                    || self.unspecified_element(), 
                    |v| self.add_ref(v, b.at(row, col))
                );
            }
        }
    }

    default fn sub_assign_matrix(&self, a: &mut M, b: N) {
        assert_eq!(a.row_count(), b.row_count());
        assert_eq!(a.col_count(), b.col_count());
        for row in 0..a.row_count() {
            for col in 0..a.col_count() {
                take_mut::take_or_recover(
                    a.at_mut(row, col), 
                    || self.unspecified_element(), 
                    |v| self.sub_ref_snd(v, b.at(row, col))
                );
            }
        }
    }
}

pub trait MatrixAddAssignProduct<M, N, P>: Ring 
    where M: MatrixView<Self::El>, N: MatrixView<Self::El>, P: MatrixViewMut<Self::El>
{
    fn add_assign_product(&self, a: M, b: N, alpha: &Self::El, c: &mut P, gamma: &Self::El);
}

impl<R: Ring, M: MatrixView<Self::El>, N: MatrixView<Self::El>, P: MatrixViewMut<Self::El>> MatrixAddAssignProduct<M, N, P> for R {
    
    default fn add_assign_product(&self, a: M, b: N, alpha: &Self::El, c: &mut P, gamma: &Self::El) {
        assert_eq!(a.row_count(), c.row_count());
        assert_eq!(b.col_count(), c.col_count());
        assert_eq!(a.col_count(), b.row_count());
        for row in 0..a.row_count() {
            for col in 0..b.col_count() {
                let c_part = self.mul_ref(gamma, c.at(row, col));
                let mut a_part = self.zero();
                for k in 0..a.col_count() {
                    a_part = self.add(a_part, self.mul_ref(a.at(row, k), b. at(k, col)));
                }
                *c.at_mut(row, col) = self.add(self.mul_ref(&a_part, alpha), c_part);
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