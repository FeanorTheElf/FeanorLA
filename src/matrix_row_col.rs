use super::matrix_view::*;
use super::vector_view::*;
use super::vec::*;
use std::marker::PhantomData;

pub struct MatrixRow<'a, T, M>
    where M: MatrixView<T> 
{
    row: usize,
    matrix: &'a M,
    item: PhantomData<T>
}

impl<'a, T, M> MatrixRow<'a, T, M>
    where M: MatrixView<T> 
{
    pub fn new(matrix: &'a M, row: usize) -> Self {
        MatrixRow {
            matrix: matrix,
            row: row,
            item: PhantomData
        }
    }
}

impl<'a, T, M> Copy for MatrixRow<'a, T, M> 
    where M: MatrixView<T> {}

impl<'a, T, M> Clone for MatrixRow<'a, T, M>
    where M: MatrixView<T> 
{
    fn clone(&self) -> Self {
        *self
    }
}

impl<'a, T, M> VectorView<T> for MatrixRow<'a, T, M> 
    where M: MatrixView<T>
{
    fn len(&self) -> usize {
        self.matrix.col_count()
    }

    fn at(&self, i: usize) -> &T {
        self.assert_in_range(i);
        self.matrix.at(self.row, i)
    }
}

pub struct MatrixRowIter<'a, T, M>
    where M: MatrixView<T> 
{
    current: MatrixRow<'a, T, M>
}

impl<'a, T, M> MatrixRowIter<'a, T, M>
    where M: MatrixView<T> 
{
    pub fn new(matrix: &'a M) -> Self {
        MatrixRowIter {
            current: MatrixRow::new(matrix, 0)
        }
    }
}

impl<'a, T, M> Iterator for MatrixRowIter<'a, T, M> 
    where M: MatrixView<T>
{
    type Item = Vector<MatrixRow<'a, T, M>, T>;

    fn next(&mut self) -> Option<Self::Item> {
        let result = self.current;
        if result.row < result.matrix.row_count() {
            self.current.row += 1;
            return Some(Vector::new(result));
        } else {
            return None;
        }
    }
}

impl<'a, T, M> Copy for MatrixRowIter<'a, T, M> 
    where M: MatrixView<T> {}

impl<'a, T, M> Clone for MatrixRowIter<'a, T, M>
    where M: MatrixView<T> 
{
    fn clone(&self) -> Self {
        *self
    }
}

pub struct MatrixCol<'a, T, M>
    where M: MatrixView<T> 
{
    col: usize,
    matrix: &'a M,
    item: PhantomData<T>
}

impl<'a, T, M> MatrixCol<'a, T, M>
    where M: MatrixView<T> 
{
    pub fn new(matrix: &'a M, col: usize) -> Self {
        MatrixCol {
            col: col,
            matrix: matrix,
            item: PhantomData
        }
    }
}

impl<'a, T, M> Copy for MatrixCol<'a, T, M> 
    where M: MatrixView<T> {}

impl<'a, T, M> Clone for MatrixCol<'a, T, M>
    where M: MatrixView<T> 
{
    fn clone(&self) -> Self {
        *self
    }
}

impl<'a, T, M> VectorView<T> for MatrixCol<'a, T, M> 
    where M: MatrixView<T>
{
    fn len(&self) -> usize {
        self.matrix.row_count()
    }

    fn at(&self, i: usize) -> &T {
        self.assert_in_range(i);
        self.matrix.at(i, self.col)
    }
}

pub struct MatrixColIter<'a, T, M>
    where M: MatrixView<T> 
{
    current: MatrixCol<'a, T, M>
}

impl<'a, T, M> MatrixColIter<'a, T, M>
    where M: MatrixView<T> 
{
    pub fn new(matrix: &'a M) -> Self {
        MatrixColIter {
            current: MatrixCol::new(matrix, 0)
        }
    }
}

impl<'a, T, M> Iterator for MatrixColIter<'a, T, M> 
    where M: MatrixView<T>
{
    type Item = Vector<MatrixCol<'a, T, M>, T>;

    fn next(&mut self) -> Option<Self::Item> {
        let result = self.current;
        if result.col < result.matrix.col_count() {
            self.current.col += 1;
            return Some(Vector::new(result));
        } else {
            return None;
        }
    }
}

impl<'a, T, M> Copy for MatrixColIter<'a, T, M> 
    where M: MatrixView<T> {}

impl<'a, T, M> Clone for MatrixColIter<'a, T, M>
    where M: MatrixView<T> 
{
    fn clone(&self) -> Self {
        *self
    }
}