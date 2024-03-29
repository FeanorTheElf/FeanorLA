use super::super::matrix_view::*;
use super::*;
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

impl<T, M> Copy for MatrixRow<T, M> 
    where M: MatrixView<T> + Copy {}

impl<T, M> Clone for MatrixRow<T, M>
    where M: MatrixView<T> + Clone
{
    fn clone(&self) -> Self {
        MatrixRow::new(self.matrix.clone(), self.row)
    }
}

impl<T, M> VectorView<T> for MatrixRow<T, M> 
    where M: MatrixView<T>
{
    type Subvector = Subvector<Self, T>;

    fn len(&self) -> usize {
        self.matrix.col_count()
    }

    fn at(&self, i: usize) -> &T {
        self.assert_in_range(i);
        self.matrix.at(self.row, i)
    }

    fn create_subvector(self, from: usize, to: usize) -> Self::Subvector {
        Subvector::new(from, to, self)
    }
}

impl<T, M> VectorViewMut<T> for MatrixRow<T, M> 
    where M: MatrixViewMut<T>
{
    type SubvectorMut = Subvector<Self, T>;

    fn at_mut(&mut self, i: usize) -> &mut T {
        self.assert_in_range(i);
        self.matrix.at_mut(self.row, i)
    }

    fn swap(&mut self, fst: usize, snd: usize) {
        self.assert_in_range(fst);
        self.assert_in_range(snd);
        self.matrix.swap((self.row, fst), (self.row, snd));
    }
    
    fn cast_subvector(subvector: Self::Subvector) -> Self::SubvectorMut {
        subvector
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
    type Item = MatrixRow<T, M>;

    fn next(&mut self) -> Option<Self::Item> {
        let result = self.current;
        if result.row < result.matrix.row_count() {
            self.current.row += 1;
            return Some(result);
        } else {
            return None;
        }
    }
}

impl<T, M> std::iter::FusedIterator for MatrixRowIter<T, M> 
    where M: MatrixView<T> + Copy {}

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
    type Subvector = Subvector<Self, T>;

    fn len(&self) -> usize {
        self.matrix.row_count()
    }

    fn at(&self, i: usize) -> &T {
        self.assert_in_range(i);
        self.matrix.at(i, self.col)
    }

    fn create_subvector(self, from: usize, to: usize) -> Self::Subvector {
        Subvector::new(from, to, self)
    }
}

impl<T, M> VectorViewMut<T> for MatrixCol<T, M> 
    where M: MatrixViewMut<T>
{
    type SubvectorMut = Subvector<Self, T>;

    fn at_mut(&mut self, i: usize) -> &mut T {
        self.assert_in_range(i);
        self.matrix.at_mut(i, self.col)
    }

    fn swap(&mut self, fst: usize, snd: usize) {
        self.assert_in_range(fst);
        self.assert_in_range(snd);
        self.matrix.swap((fst, self.col), (snd, self.col));
    }

    fn cast_subvector(subvector: Self::Subvector) -> Self::SubvectorMut {
        subvector
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

impl<T, M> Iterator for MatrixColIter<T, M> 
    where M: MatrixView<T> + Copy
{
    type Item = MatrixCol<T, M>;

    fn next(&mut self) -> Option<Self::Item> {
        let result = self.current;
        if result.col < result.matrix.col_count() {
            self.current.col += 1;
            return Some(result);
        } else {
            return None;
        }
    }
}

impl<T, M> Copy for MatrixColIter<T, M> 
    where M: MatrixView<T> + Copy {}

impl<T, M> Clone for MatrixColIter<T, M>
    where M: MatrixView<T> + Copy
{
    fn clone(&self) -> Self {
        *self
    }
}

impl<T, M> std::iter::FusedIterator for MatrixColIter<T, M> 
    where M: MatrixView<T> + Copy {}
