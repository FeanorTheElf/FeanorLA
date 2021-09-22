use super::matrix_view::*;

use std::marker::PhantomData;

#[derive(Debug)]
pub struct MatrixTranspose<T, M>
    where M: MatrixView<T> 
{
    matrix: M,
    element: PhantomData<T>
}

impl<T, M> MatrixTranspose<T, M>
    where M: MatrixView<T> 
{
    pub fn new(matrix: M) -> Self {
        MatrixTranspose { 
            matrix: matrix,
            element: PhantomData
        }
    }
}

impl<T, M> Copy for MatrixTranspose<T, M> 
    where M: MatrixView<T> + Copy {}

impl<T, M> Clone for MatrixTranspose<T, M>
    where M: MatrixView<T> + Clone
{
    fn clone(&self) -> Self {
        MatrixTranspose::new(self.matrix.clone())
    }
}

impl<T, M> MatrixView<T> for MatrixTranspose<T, M> 
    where M: MatrixView<T>
{
    fn row_count(&self) -> usize {
        self.matrix.col_count()
    }

    fn col_count(&self) -> usize {
        self.matrix.row_count()
    }
    
    fn at(&self, i: usize, j: usize) -> &T {
        self.matrix.at(i, j)
    }
}

impl<T, M> MatrixViewMut<T> for MatrixTranspose<T, M> 
    where M: MatrixViewMut<T>
{
    fn at_mut(&mut self, i: usize, j: usize) -> &mut T {
        self.assert_row_in_range(i);
        self.assert_col_in_range(j);
        self.matrix.at_mut(j, i)
    }

    fn swap(&mut self, (fst_i, fst_j): (usize, usize), (snd_i, snd_j): (usize, usize)) {
        self.assert_row_in_range(fst_i);
        self.assert_col_in_range(fst_j);
        self.assert_row_in_range(snd_i);
        self.assert_col_in_range(snd_j);
        self.matrix.swap((fst_j, fst_i), (snd_j, snd_i));
    }
}
