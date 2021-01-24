use super::matrix_view::*;
use super::vector_view::*;
use super::vector::*;
use std::marker::PhantomData;

#[derive(Debug)]
pub struct MatrixRef<'a, M, T>
    where M: MatrixView<T>
{
    rows_begin: usize,
    rows_end: usize,
    cols_begin: usize,
    cols_end: usize,
    matrix: &'a M,
    element: PhantomData<T>
}

#[derive(Debug)]
pub struct MatrixRefMut<'a, M, T>
    where M: MatrixViewMut<T>
{
    rows_begin: usize,
    rows_end: usize,
    cols_begin: usize,
    cols_end: usize,
    matrix: &'a mut M,
    element: PhantomData<T>
}

impl<'a, M, T> MatrixView<T> for MatrixRef<'a, M, T>
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

impl<'a, M, T> MatrixView<T> for MatrixRefMut<'a, M, T> 
    where M: MatrixViewMut<T>
{
    fn row_count(&self) -> usize {
        self.rows_end - self.rows_begin
    }

    fn col_count(&self) -> usize {
        self.cols_end - self.cols_begin
    }

    fn at(&self, row: usize, col: usize) -> &T {
        self.matrix.at(row + self.rows_begin, col + self.rows_end)
    }
}

impl<'a, M, T> MatrixViewMut<T> for MatrixRefMut<'a, M, T> 
    where M: MatrixViewMut<T>
{
    fn at_mut(&mut self, row: usize, col: usize) -> &mut T {
        self.matrix.at_mut(row + self.rows_begin, col + self.rows_end)
    }
}

impl<'a, M, T> Clone for MatrixRef<'a, M, T> 
    where M: MatrixView<T>
{
    fn clone(&self) -> Self {
        *self
    }
}

impl<'a, M, T> Copy for MatrixRef<'a, M, T> 
    where M: MatrixView<T> {}

pub struct MatrixRefMutRowIter<'a, M, T>
    where M: LifetimeMatrixMutRowIter<'a, T>
{
    current: M::RowIter,
    from_col: usize,
    to_col: usize,
    to_yield: usize
}

impl<'a, M, T> Iterator for MatrixRefMutRowIter<'a, M, T>
    where M: LifetimeMatrixMutRowIter<'a, T>
{
    type Item = VectorRestriction<M::RowRef, T>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.to_yield == 0 {
            return None;
        } else {
            self.to_yield -= 1;
            return self.current.next().map(|r| VectorRestriction::restrict(r, self.from_col, self.to_col));
        }
    }
}

impl<'a, 'b, M: 'a, T: 'a> LifetimeMatrixMutRowIter<'a, T> for MatrixRefMut<'b, M, T> 
    where M: LifetimeMatrixMutRowIter<'a, T>
{
    type RowRef = VectorRestriction<M::RowRef, T>;
    type RowIter = MatrixRefMutRowIter<'a, M, T>;

    fn rows_mut(&'a mut self) -> Self::RowIter {
        let mut it = self.matrix.rows_mut();
        for _i in 0..self.rows_begin {
            it.next();
        }
        MatrixRefMutRowIter {
            current: it,
            from_col: self.cols_begin,
            to_col: self.cols_end,
            to_yield: self.rows_end - self.rows_begin
        }
    }

    fn get_row_mut(&'a mut self, i: usize) -> Self::RowRef {
        self.assert_row_in_range(i);
        VectorRestriction::restrict(self.matrix.get_row_mut(i + self.rows_begin), self.cols_begin, self.cols_end)
    }
}

#[cfg(test)]
use super::matrix_owned::MatrixOwned;

#[test]
fn test_row_iter() {
    let a = MatrixOwned::from_fn(4, 4, |i, j| i + 4 * j);
    let b = MatrixRef {
        rows_begin: 1,
        rows_end: 3,
        cols_begin: 2,
        cols_end: 4,
        element: PhantomData,
        matrix: &a
    };
    let mut it = b.rows();
    let r1 = it.next().unwrap();
    let r2 = it.next().unwrap();
    assert!(it.next().is_none());

    assert_eq!(2, r1.len());
    assert_eq!(2, r2.len());
    assert_eq!(9, *r1.at(0));
    assert_eq!(13, *r1.at(1));
    assert_eq!(10, *r2.at(0));
    assert_eq!(14, *r2.at(1));

    assert_eq!(3, *a.get_row(3).at(0));
}

#[test]
fn test_row_iter_mut() {
    let mut a = MatrixOwned::from_fn(4, 4, |i, j| i + 4 * j);
    let mut b = MatrixRefMut {
        rows_begin: 1,
        rows_end: 3,
        cols_begin: 2,
        cols_end: 4,
        element: PhantomData,
        matrix: &mut a
    };
    let mut it = b.rows_mut();
    let mut r1 = it.next().unwrap();
    let mut r2 = it.next().unwrap();
    assert!(it.next().is_none());

    assert_eq!(2, r1.len());
    assert_eq!(9, *r1.at(0));
    assert_eq!(13, *r1.at(1));
    assert_eq!(14, *r2.at(1));
    assert_eq!(9, *r1.at_mut(0));
    assert_eq!(14, *r2.at_mut(1));

    *r1.at_mut(1) = 20;
    *a.get_row_mut(3).at_mut(2) = 21;
    assert_eq!(20, *a.at(1, 3));
    assert_eq!(21, *a.at(3, 2));
}