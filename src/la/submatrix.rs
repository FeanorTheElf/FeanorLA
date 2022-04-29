use super::matrix_view::*;
use super::subvector::*;
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
    type Item = Subvector<M::RowRef, T>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.to_yield == 0 {
            return None;
        } else {
            self.to_yield -= 1;
            return self.current.next().map(|r| Subvector::new(self.from_col, self.to_col, r));
        }
    }
}

impl<'a, M, T: 'a> LifetimeMatrixMutRowIter<'a, T> for Submatrix<M, T> 
    where M: LifetimeMatrixMutRowIter<'a, T>
{
    type RowRef = Subvector<M::RowRef, T>;
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
        Subvector::new(self.cols_begin, self.cols_end, self.matrix.get_row_mut(i + self.rows_begin))
    }
}

impl<M, T: 'static> MatrixMutRowIter<T> for Submatrix<M, T> 
    where M: MatrixMutRowIter<T>
{}

#[cfg(test)]
use super::vector_view::*;

#[test]
fn test_row_iter_mut() {
    let mut a = MatrixOwned::from_fn(4, 4, |i, j| i + 4 * j);
    let mut b = Submatrix {
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