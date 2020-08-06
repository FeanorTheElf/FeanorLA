use super::prelude::*;
use std::ops::{AddAssign, Neg};

fn two<T>() -> T
    where T: AddAssign + One
{
    let mut two = T::one();
    two += T::one();
    two
}

fn householder_left<T>(y: VectorRef<T>, mut A: MatrixRefMut<T>)
    where T: Field + Clone + std::fmt::Debug
{
    assert_eq!(y.len(), A.rows());
    let mut yyTA = Matrix::from_func(A.rows(), A.cols(), |row, col| {
        let mut yTA = y[0].clone() * A[0][col].clone();
        for k in 1..y.len() {
            yTA += y[k].clone() * A[k][col].clone();
        }
        return yTA * y[row].clone();
    });
    yyTA.scal(two());
    A -= yyTA;
}

fn householder_right<T>(y: VectorRef<T>, mut A: MatrixRefMut<T>)
    where T: Field + Clone
{
    assert_eq!(y.len(), A.cols());
    let mut AyyT = Matrix::from_func(A.rows(), A.cols(), |row, col| {
        let mut yTA = y[0].clone() * A[row][0].clone();
        for k in 1..y.len() {
            yTA += y[k].clone() * A[row][k].clone();
        }
        return yTA * y[col].clone();
    });
    AyyT.scal(two());
    A -= AyyT;
}

fn abs<T>(val: T) -> T
    where T: Neg<Output = T> + Zero + PartialOrd
{
    if !(val < T::zero()) {
        val
    } else {
        -val
    }
}

fn sgn<T>(val: &T) -> T
    where T: Neg<Output = T> + Zero + PartialOrd + One
{
    if !(*val < T::zero()) {
        T::one()
    } else {
        -T::one()
    }
}

pub fn qr_decompose<T>(mut A: MatrixRefMut<T>) -> Matrix<T>
    where T: Field + Root + Clone + PartialOrd + std::fmt::Debug
{
    let mut Q = Matrix::identity(A.rows());
    let mut y_base: Vector<T> = Vector::zero(A.rows());
    let two = T::one() + T::one();
    let half = T::one() / two.clone();
    for k in 0..(A.cols().min(A.rows()) - 1) {
        let mut y = y_base.get_mut(k..);
        let x = A.get((k.., k..=k));
        let gamma = two.clone() * x.frobenius_square().sqrt();
        let sgn = sgn(&x[0][0]);
        let y1 = (half.clone() + abs(x[0][0].clone()) / gamma.clone()).sqrt();
        y[0] = y1.clone();
        for i in 1..y.len() {
            y[i] = sgn.clone() * x[i][0].clone() / (gamma.clone() * y1.clone());
        }
        householder_left(y.as_const(), A.get_mut((k.., k..)));
        householder_right(y.as_const(), Q.get_mut((.., k..)));
    }
    return Q;
}

#[cfg(test)]
use super::macros::ApproxEq;

#[test]
fn test_qr() {
    let mut m = matlab![1., 0., 1.; 1., 1., 1.; 0., 2., 2.];
    let q = matlab![-0.707, 0.236, 0.667; -0.707, -0.236, -0.667; 0., -0.943, 0.333];
    let r = matlab![-1.414, -0.707, -1.414; 0., -2.121, -1.886; 0., 0., 0.667];
    let actual_q = qr_decompose(m.as_mut());
    assert_approx_eq!(q, &actual_q, 0.005);
    assert_approx_eq!(r, &m, 0.005);
}