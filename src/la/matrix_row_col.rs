use super::matrix_view::*;
use super::vector_view::*;
use super::vec::*;
use std::marker::PhantomData;

#[derive(Debug)]
pub struct MatrixRow<T, M>
    where M: MatrixView<T> 
{
    row: usize,
    matrix: M,
    item: PhantomData<T>
}

impl<T, M> MatrixRow<T, M>
    where M: MatrixView<T> 
{
    pub fn new(matrix: M, row: usize) -> Self {
        MatrixRow {
            matrix: matrix,
            row: row,
            item: PhantomData
        }
    }
}

impl<'a, T, M> Copy for MatrixRow<T, M> 
    where M: MatrixView<T> + Copy {}

impl<'a, T, M> Clone for MatrixRow<T, M>
    where M: MatrixView<T> + Clone
{
    fn clone(&self) -> Self {
        MatrixRow::new(self.matrix.clone(), self.row)
    }
}

impl<'a, T, M> VectorView<T> for MatrixRow<T, M> 
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

pub struct MatrixRowIter<T, M>
    where M: MatrixView<T> + Copy
{
    current: MatrixRow<T, M>
}

impl<T, M> MatrixRowIter<T, M>
    where M: MatrixView<T> + Copy
{
    pub fn new(matrix: M) -> Self {
        MatrixRowIter {
            current: MatrixRow::new(matrix, 0)
        }
    }
}

impl<T, M> Iterator for MatrixRowIter<T, M> 
    where M: MatrixView<T> + Copy
{
    type Item = Vector<MatrixRow<T, M>, T>;

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

impl<T, M> Copy for MatrixRowIter<T, M> 
    where M: MatrixView<T> + Copy {}

impl<T, M> Clone for MatrixRowIter<T, M>
    where M: MatrixView<T> + Copy
{
    fn clone(&self) -> Self {
        *self
    }
}

#[derive(Debug)]
pub struct MatrixCol<T, M>
    where M: MatrixView<T> 
{
    col: usize,
    matrix: M,
    item: PhantomData<T>
}

impl<T, M> MatrixCol<T, M>
    where M: MatrixView<T> 
{
    pub fn new(matrix: M, col: usize) -> Self {
        MatrixCol {
            col: col,
            matrix: matrix,
            item: PhantomData
        }
    }
}

impl<T, M> Copy for MatrixCol<T, M> 
    where M: MatrixView<T> + Copy {}

impl<T, M> Clone for MatrixCol<T, M>
    where M: MatrixView<T> + Clone
{
    fn clone(&self) -> Self {
        MatrixCol::new(self.matrix.clone(), self.col)
    }
}

impl<T, M> VectorView<T> for MatrixCol<T, M> 
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

pub struct MatrixColIter<T, M>
    where M: MatrixView<T> + Copy
{
    current: MatrixCol<T, M>
}

impl<T, M> MatrixColIter<T, M>
    where M: MatrixView<T> + Copy
{
    pub fn new(matrix: M) -> Self {
        MatrixColIter {
            current: MatrixCol::new(matrix, 0)
        }
    }
}

impl<'a, T, M> Iterator for MatrixColIter<T, M> 
    where M: MatrixView<T> + Copy
{
    type Item = Vector<MatrixCol<T, M>, T>;

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

impl<'a, T, M> Copy for MatrixColIter<T, M> 
    where M: MatrixView<T> + Copy {}

impl<'a, T, M> Clone for MatrixColIter<T, M>
    where M: MatrixView<T> + Copy
{
    fn clone(&self) -> Self {
        *self
    }
}