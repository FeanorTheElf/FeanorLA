#![allow(non_snake_case)]
use super::super::la::mat::*;
use super::super::float::*;

pub fn cholesky<M, T>(A: &mut Matrix<M, T>)
    where M: MatrixViewMut<T>, T: Float
{
    assert_eq!(A.row_count(), A.col_count());
    let n = A.row_count();
    for j in 0..n {
        // compute the diagonal element L_jj
        let mut row_square_sum = A.at(j, j).clone();
        for k in 0..j {
            row_square_sum -= A.at(j, k).clone() * A.at(j, k).clone();
        }
        *A.at_mut(j, j) = row_square_sum.sqrt();
        // compute the column below jj
        for i in (j + 1)..n {
            let mut rows_inner_prod = A.at(i, j).clone();
            for k in 0..j {
                rows_inner_prod -= A.at(i, k).clone() * A.at(j, k).clone();
            }
            *A.at_mut(i, j) = rows_inner_prod / A.at(j, j).clone();
        }
        // clear upper half of A
        for i in 0..j {
            *A.at_mut(i, j) = T::zero();
        }
    }
}

#[test]
fn test_cholesky() {
    #[rustfmt::skip]
    let mut A = Matrix::from_array([[1.,  2.,  4.], 
                                    [2., 13., 23.],
                                    [4., 23., 77.]]);

    #[rustfmt::skip]
    let B = Matrix::from_array([[1., 0., 0.], 
                                [2., 3., 0.],
                                [4., 5., 6.]]);

    cholesky(&mut A);
    assert_eq!(A, B);
}