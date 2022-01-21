use super::matrix_view::*;
use super::vector_view::*;

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

#[derive(Debug)]
pub struct VectorConstant<T>
{
    len: usize,
    constant: T
}

impl<T> Copy for VectorConstant<T>
    where T: Copy
{}

impl<T> Clone for VectorConstant<T>
    where T: Clone
{
    fn clone(&self) -> Self {
        Self::new(self.len, self.constant.clone())
    }
}

impl<T> VectorConstant<T>
{
    pub fn new(len: usize, constant: T) -> Self {
        VectorConstant {
            len, constant
        }
    }
}

impl<T> VectorView<T> for VectorConstant<T> {
    
    fn len(&self) -> usize {
        self.len
    }

    fn at(&self, i: usize) -> &T {
        self.assert_in_range(i);
        &self.constant
    }
}

#[derive(Debug)]
pub struct VectorUnit<T>
{
    len: usize,
    index: usize,
    zero: T,
    one: T
}

impl<T> Copy for VectorUnit<T>
    where T: Copy
{}

impl<T> Clone for VectorUnit<T>
    where T: Clone
{
    fn clone(&self) -> Self {
        Self::new(self.len, self.index, self.zero.clone(), self.one.clone())
    }
}

impl<T> VectorUnit<T>
{
    pub fn new(len: usize, index: usize, zero: T, one: T) -> Self {
        assert!(index < len);
        VectorUnit {
            len, index, zero, one
        }
    }
}

impl<T> VectorView<T> for VectorUnit<T> {
    
    fn len(&self) -> usize {
        self.len
    }

    fn at(&self, i: usize) -> &T {
        self.assert_in_range(i);
        if i == self.index {
            return &self.one;
        } else {
            return &self.zero;
        }
    }
}