use super::matrix_view::*;

#[derive(Debug)]
pub struct MatrixConstant<T>
{
    rows: usize,
    cols: usize,
    constant: T
}

impl<T> Copy for MatrixConstant<T>
    where T: Copy
{}

impl<T> Clone for MatrixConstant<T>
    where T: Clone
{
    fn clone(&self) -> Self {
        Self::new(self.rows, self.cols, self.constant.clone())
    }
}

impl<T> MatrixConstant<T>
{
    pub fn new(rows: usize, cols: usize, constant: T) -> Self {
        MatrixConstant {
            rows, cols, constant
        }
    }
}

impl<T> MatrixView<T> for MatrixConstant<T> {

    fn row_count(&self) -> usize {
        self.rows
    }

    fn col_count(&self) -> usize {
        self.cols
    }

    fn at(&self, row: usize, col: usize) -> &T {
        self.assert_row_in_range(row);
        self.assert_col_in_range(col);
        &self.constant
    }
}
