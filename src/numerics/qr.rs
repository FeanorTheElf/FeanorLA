#![allow(non_snake_case)]
use super::super::la::mat::*;
use super::super::la::algorithms::MatrixFrobenius;
use super::super::alg::*;
use super::super::float::*;
use std::ops::AddAssign;

fn two<T>() -> T
    where T: AddAssign + One
{
    let mut two = T::one();
    two += T::one();
    two
}

fn householder_left<V, M, T>(y: Vector<V, T>, mut A: Matrix<M, T>)
    where V: VectorView<T>, M: MatrixViewMut<T>, T: Float + Clone
{
    assert_eq!(y.len(), A.row_count());
    let yTA = Matrix::row_vec(y.as_ref()) * A.as_ref();
    for i in 0..A.row_count() {
        for j in 0..A.col_count() {
            *A.at_mut(i, j) -= two::<T>() * y.at(i).clone() * yTA.at(0, j).clone();
        }
    }
}

fn householder_right<V, M, T>(y: Vector<V, T>, mut A: Matrix<M, T>)
    where V: VectorView<T>, M: MatrixViewMut<T>, T: Float + Clone
{
    assert_eq!(y.len(), A.col_count());
    let Ay = A.as_ref() * Matrix::col_vec(y.as_ref());
    for i in 0..A.row_count() {
        for j in 0..A.col_count() {
            *A.at_mut(i, j) -= two::<T>() * Ay.at(i, 0).clone() * y.at(j).clone();
        }
    }
}

fn abs<T>(val: T) -> T
    where T: Float
{
    if !(val < T::zero()) {
        val
    } else {
        -val
    }
}

fn sgn<T>(val: &T) -> T
    where T: Float
{
    if !(*val < T::zero()) {
        T::one()
    } else {
        -T::one()
    }
}

///
/// Given a matrix A, calculates an orthonormal matirx Q and an 
/// upper triangle matrix R such that A = QR. The matrix Q is returned
/// and R is assigned to the input matrix.
/// 
pub fn qr_decompose<M, T>(A: &mut Matrix<M, T>) -> Matrix<MatrixOwned<T>, T>
    where M: MatrixViewMut<T>, T: Float + Root + std::fmt::Display
{
    let mut Q = Matrix::identity(A.row_count(), A.row_count());
    let mut y_base = Vector::zero(A.row_count()).to_owned();
    let half = T::one() / two::<T>();
    let two = two::<T>();

    for k in 0..A.col_count().min(A.row_count() - 1) {
        let mut y = y_base.subvector_mut(k..);
        let x = A.submatrix(k.., k..=k);
        let gamma = two.clone() * StaticRing::<T>::RING.calc_matrix_frobenius_norm_square(x.as_ref()).sqrt();

        // by choosing this correctly, the addition for y1 involves two positive
        // numbers, preventing catastrophic cancellation
        let sgn = sgn(x.at(0, 0));
        let y1 = (half.clone() + abs(x.at(0, 0).clone()) / gamma.clone()).sqrt();
        *y.at_mut(0) = y1.clone();
        let factor = sgn.clone() / (gamma * y1);
        for i in 1..y.len() {
            *y.at_mut(i) = factor.clone() * x.at(i, 0).clone();
        }
        householder_left(y.as_ref(), A.submatrix_mut(k.., k..));
        householder_right(y.as_ref(), Q.submatrix_mut(.., k..));
    }
    return Q;
}

#[cfg(test)]
use super::macros::ApproxEq;

#[test]
fn test_qr() {
    #[rustfmt::skip]
    let mut m = Matrix::from_array([[1., 0., 1.], 
                                    [1., 1., 1.], 
                                    [0., 2., 2.]]);
    
    #[rustfmt::skip]
    let q = Matrix::from_array([[-0.707,  0.236,  0.667], 
                                [-0.707, -0.236, -0.667], 
                                [0.,     -0.943,  0.333]]);

    #[rustfmt::skip]
    let r = Matrix::from_array([[-1.414, -0.707, -1.414], 
                                [0.,     -2.121, -1.886], 
                                [0.,      0.,     0.667]]);
                                
    let actual_q = qr_decompose(&mut m);
    assert_approx_eq!(q, &actual_q, 0.005);
    assert_approx_eq!(r, &m, 0.005);
}

#[test]
fn test_qr_rows_ge_cols() {
    #[rustfmt::skip]
    let mut m = Matrix::from_array([[1., 0.], 
                                    [1., 1.], 
                                    [0., 2.]]);

    #[rustfmt::skip]
    let q = Matrix::from_array([[-0.707,  0.236,  0.666], 
                                [-0.707, -0.236, -0.6666], 
                                [0.,     -0.943,  0.333]]);
    
    #[rustfmt::skip]
    let r = Matrix::from_array([[-1.414, -0.707], 
                                [0.,     -2.121], 
                                [0.,      0.]]);

    let actual_q = qr_decompose(&mut m);
    assert_approx_eq!(r, &m, 0.005);
    assert_approx_eq!(q, &actual_q, 0.005);
}

#[test]
fn test_qr_cols_ge_rows() {
    #[rustfmt::skip]
    let mut m = Matrix::from_array([[1., 0., 3.], 
                                    [1., 1., 2.]]);

    #[rustfmt::skip]
    let q = Matrix::from_array([[-0.707, -0.707],
                                [-0.707,  0.707]]);
    
    #[rustfmt::skip]
    let r = Matrix::from_array([[-1.414, -0.707, -3.535],
                                [ 0.   ,  0.707, -0.707]]);

    let actual_q = qr_decompose(&mut m);
    assert_approx_eq!(r, &m, 0.005);
    assert_approx_eq!(q, &actual_q, 0.005);
}