use super::super::ring::*;
use super::super::la::mat::*;
use super::super::algebra::fractions::*;
use super::super::embedding::*;

pub trait MatrixDeterminant<M>: Ring
    where M: MatrixView<Self::El>
{
    fn matrix_determinant(&self, matrix: Matrix<M, Self::El>) -> Self::El;
}

fn compute_det<F>(field: F, mut work_matrix: Matrix<MatrixOwned<F::El>, F::El>) -> F::El
    where F: Ring
{
    assert!(field.is_field().can_use());
    let mut det_factor_inv = field.one();
    let mut negated = false;
    let result = work_matrix.gaussion_elimination_half(
        |_, a, ()| { field.mul_assign(&mut det_factor_inv, a) },
        |_, _, ()| { negated = !negated; },
        |_, _, _, ()| {},
        &mut (),
        &field
    );
    if let Ok(()) = result {
        let value = work_matrix.into_nonmain_diag(0).into_owned().raw_data().into_vec().into_iter().fold(field.one(), |a, b| field.mul(a, b));
        let scaled = field.div(value, &det_factor_inv);
        if negated {
            return field.neg(scaled);
        } else {
            return scaled;
        }
    } else {
        return field.zero();
    }
}

impl<R, M> MatrixDeterminant<M> for R
    where R: DivisibilityInfoRing + CanonicalIsomorphismInfo<R>, M: MatrixView<R::El>
{
    default fn matrix_determinant(&self, matrix: Matrix<M, R::El>) -> R::El {
        if self.is_field().can_use() {
            compute_det(self, matrix.into_owned())
        } else if self.is_integral().can_use() {
            assert!(self.is_divisibility_computable());
            let field = FieldOfFractions::new(self);
            let incl = embedding(self, field);
            let work_matrix = Matrix::from_fn(matrix.row_count(), matrix.col_count(), |i, j| incl(matrix.at(i, j).clone()));
            let result = compute_det(&field, work_matrix);
            field.in_base_ring(&result).unwrap()
        } else {
            unimplemented!()
        }
    }
}

impl<M, T> Matrix<M, T>
    where M: MatrixView<T>
{
    pub fn det<R>(&self, ring: &R) -> R::El
        where R: Ring<El = T> + DivisibilityInfoRing + CanonicalIsomorphismInfo<R>
    {
        ring.matrix_determinant(self.as_ref())
    }
}

#[cfg(test)]
use super::super::primitive::*;

#[test]
fn test_det() {
    let a = Matrix::from_array([[1, 2], [3, 4]]);
    assert_eq!(-2, a.det(&i32::RING));
}