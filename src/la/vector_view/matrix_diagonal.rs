use super::super::matrix_view::*;
use super::*;

use std::marker::PhantomData;

fn get_matrix_index(diagonal_index: i64, i: usize) -> (usize, usize) {
    if diagonal_index >= 0 {
        (i as usize, i + diagonal_index as usize)
    } else {
        (i + ((-diagonal_index) as usize), i as usize)
    }
}

#[derive(Debug)]
pub struct MatrixDiagonal<M, T>
    where M: MatrixView<T>
{
    matrix: M,
    diagonal: i64,
    element: PhantomData<T>
}

impl<M, T> MatrixDiagonal<M, T>
where M: MatrixView<T>
{
    pub fn new(matrix: M, diagonal: i64) -> MatrixDiagonal<M, T> {
        assert!(diagonal <= matrix.col_count() as i64);
        assert!(-diagonal <= matrix.row_count() as i64);
        MatrixDiagonal {
            matrix: matrix,
            diagonal: diagonal,
            element: PhantomData
        }
    }
}

impl<M, T> VectorView<T> for MatrixDiagonal<M, T>
where M: MatrixView<T>
{
    type Subvector = Subvector<Self, T>;

    fn len(&self) -> usize {
        if self.diagonal >= 0 {
            usize::min(self.matrix.row_count(), self.matrix.col_count() - self.diagonal as usize)
        } else {
            usize::min(self.matrix.row_count() - ((-self.diagonal) as usize), self.matrix.col_count())
        }
    }

    fn at(&self, index: usize) -> &T {
        self.assert_in_range(index);
        let (i, j) = get_matrix_index(self.diagonal, index);
        self.matrix.at(i, j)
    }

    fn create_subvector(self, from: usize, to: usize) -> Self::Subvector {
        Subvector::new(from, to, self)
    }
}

impl<M, T> VectorViewMut<T> for MatrixDiagonal<M, T>
where M: MatrixViewMut<T>
{
    type SubvectorMut = Subvector<Self, T>;

    fn at_mut(&mut self, index: usize) -> &mut T {
        self.assert_in_range(index);
        let (i, j) = get_matrix_index(self.diagonal, index);
        self.matrix.at_mut(i, j)
    }

    fn swap(&mut self, fst: usize, snd: usize) {
        self.assert_in_range(fst);
        self.assert_in_range(snd);

        if fst != snd {
            let (fst_i, fst_j) = get_matrix_index(self.diagonal, fst);
            let (snd_i, snd_j) = get_matrix_index(self.diagonal, snd);
            self.matrix.swap((fst_i, fst_j), (snd_i, snd_j));
        }
    }

    fn cast_subvector(subvector: Self::Subvector) -> Self::SubvectorMut {
        subvector
    }
}