use super::*;
use std::marker::PhantomData;

#[derive(Debug)]
pub struct Submatrix<M, T>
    where M: MatrixView<T>
{
    rows_begin: usize,
    rows_end: usize,
    cols_begin: usize,
    cols_end: usize,
    matrix: M,
    element: PhantomData<T>
}
impl<M, T> Submatrix<M, T>
    where M: MatrixView<T>
{
    pub fn new(rows_begin: usize, rows_end: usize, cols_begin: usize, cols_end: usize, matrix: M) -> Self {
        assert!(rows_begin <= rows_end);
        assert!(cols_begin <= cols_end);
        assert!(rows_end <= matrix.row_count());
        assert!(cols_end <= matrix.col_count());
        Submatrix {
            rows_begin, rows_end, cols_end, cols_begin, matrix, element: PhantomData
        }
    }
}

impl<M, T> MatrixView<T> for Submatrix<M, T>
    where M: MatrixView<T>
{
    fn row_count(&self) -> usize {
        self.rows_end - self.rows_begin
    }

    fn col_count(&self) -> usize {
        self.cols_end - self.cols_begin
    }

    fn at(&self, row: usize, col: usize) -> &T {
        self.matrix.at(row + self.rows_begin, col + self.cols_begin)
    }
}

impl<M, T> MatrixViewMut<T> for Submatrix<M, T> 
    where M: MatrixViewMut<T>
{
    fn at_mut(&mut self, row: usize, col: usize) -> &mut T {
        self.matrix.at_mut(row + self.rows_begin, col + self.cols_begin)
    }
    
    fn swap(&mut self, fst: (usize, usize), snd: (usize, usize)) {
        assert!(fst != snd);
        self.assert_row_in_range(fst.0);
        self.assert_row_in_range(snd.0);
        self.assert_col_in_range(fst.1);
        self.assert_col_in_range(snd.1);
        self.matrix.swap((fst.0 + self.rows_begin, fst.1 + self.cols_begin), (snd.0 + self.rows_begin, snd.1 + self.cols_begin))
    }
}

impl<M, T> Clone for Submatrix<M, T> 
    where M: MatrixView<T> + Clone
{
    fn clone(&self) -> Self {
        Submatrix {
            rows_begin: self.rows_begin,
            rows_end: self.rows_end,
            cols_begin: self.cols_begin,
            cols_end: self.cols_end,
            matrix: self.matrix.clone(),
            element: PhantomData
        }
    }
}

impl<M, T> Copy for Submatrix<M, T> 
    where M: MatrixView<T> + Copy {}
