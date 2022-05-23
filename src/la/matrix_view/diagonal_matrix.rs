use super::*;
use super::constant_value_matrix::*;
use super::super::vector_view::*;

fn get_matrix_index(diagonal_index: i64, i: usize) -> (usize, usize) {
    if diagonal_index >= 0 {
        (i as usize, i + diagonal_index as usize)
    } else {
        (i + ((-diagonal_index) as usize), i as usize)
    }
}

#[derive(Debug)]
pub struct DiagonalMatrix<V, T>
    where V: VectorView<T>
{
    diagonal: V,
    diagonal_index: i64,
    zero: T,
}

impl<V, T> DiagonalMatrix<V, T>
where V: VectorView<T>
{
    pub fn new(diagonal: V, diagonal_index: i64, zero: T) -> DiagonalMatrix<V, T> {
        let diag_len = diagonal.len();
        let result = DiagonalMatrix {
            diagonal, diagonal_index, zero
        };
        assert_eq!((result.row_count(), result.col_count()), get_matrix_index(diagonal_index, diag_len));
        return result;
    }
}

impl<V, T> MatrixView<T> for DiagonalMatrix<V, T>
where V: VectorView<T>
{
    fn row_count(&self) -> usize {
        if self.diagonal_index >= 0 {
            self.diagonal.len()
        } else {
            self.diagonal.len() + (-self.diagonal_index) as usize
        }
    }

    fn col_count(&self) -> usize {
        if self.diagonal_index >= 0 {
            self.diagonal.len() + self.diagonal_index as usize
        } else {
            self.diagonal.len()
        }
    }

    fn at(&self, i: usize, j: usize) -> &T {
        self.assert_row_in_range(i);
        self.assert_col_in_range(j);
        if self.diagonal_index >= 0 && i + self.diagonal_index as usize == j {
            debug_assert_eq!((i, j), get_matrix_index(self.diagonal_index, i));
            return self.diagonal.at(i);
        } else if self.diagonal_index < 0 && i == j + (-self.diagonal_index) as usize {
            debug_assert_eq!((i, j), get_matrix_index(self.diagonal_index, j));
            return self.diagonal.at(j);
        } else {
            return &self.zero;
        }
    }

    fn into_owned(self) -> MatrixOwned<T> 
        where T: Clone
    {
        let mut result = MatrixConstant::new(self.row_count(), self.col_count(), self.zero).into_owned();
        for index in 0..self.diagonal.len() {
            let (i, j) = get_matrix_index(self.diagonal_index, index);
            *result.at_mut(i, j) = self.diagonal.at(index).clone();
        }
        return result;
    }
}
