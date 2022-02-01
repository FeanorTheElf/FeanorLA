use super::super::alg::*;
use super::mat::*;
use super::super::algebra::primality::*;
use super::super::algebra::fractions::*;

pub trait MatrixDeterminant<M>: Ring
    where M: MatrixView<Self::El>
{
    fn matrix_determinant(&self, matrix: Matrix<M, Self::El>) -> Self::El;
}

fn compute_det<F>(field: F, mut work_matrix: Matrix<MatrixOwned<F::El>, F::El>) -> F::El
    where F: Ring
{
    assert!(field.is_field());
    let mut det_factor = field.one();
    let mut negated = false;
    let result = work_matrix.gaussion_elimination_half(
        |_, a, ()| {
            take_mut::take_or_recover(&mut det_factor, || field.unspecified_element(), |x| field.div(x, &a))
        },
        |_, _, ()| { negated = !negated; },
        |_, _, _, ()| {},
        &mut (),
        &field
    );
    if let Ok(()) = result {
        let value = work_matrix.into_nonmain_diag(0).into_owned().raw_data().into_vec().into_iter().fold(field.one(), |a, b| field.mul(a, b));
        let scaled = field.mul(value, det_factor);
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
    where R: DivisibilityInformationRing, M: MatrixView<R::El>
{
    default fn matrix_determinant(&self, matrix: Matrix<M, R::El>) -> R::El {
        if self.is_field() {
            compute_det(self, matrix.into_owned())
        } else if self.is_integral() {
            assert!(self.is_divisibility_computable());
            let (field, mut incl) = self.field_of_fractions();
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
        where R: Ring<El = T> + DivisibilityInformationRing
    {
        ring.matrix_determinant(self.as_ref())
    }
}

#[test]
fn test_det() {
    let a = Matrix::from_array([[1, 2], [3, 4]]);
    assert_eq!(-2, i32::RING.matrix_determinant(a));
}