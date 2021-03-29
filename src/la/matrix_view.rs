use super::vector_view::*;

pub trait MatrixView<T>: Sized {

    fn row_count(&self) -> usize;
    fn col_count(&self) -> usize;
    fn at(&self, row: usize, col: usize) -> &T;

    fn assert_row_in_range(&self, row: usize) {
        assert!(row < self.row_count(), "Row index {} out of range 0..{}", row, self.row_count());
    }

    fn assert_col_in_range(&self, col: usize) {
        assert!(col < self.col_count(), "Column index {} out of range 0..{}", col, self.col_count());
    }

    fn to_owned(self) -> MatrixOwned<T> 
        where T: Clone
    {
        MatrixOwned::from_fn(self.row_count(), self.col_count(), |i, j| self.at(i, j).clone())
    }
}

pub trait FromFnCreateable<T> : MatrixView<T> {

    fn from_fn<F>(rows: usize, cols: usize, f: F) -> Self
        where F: FnMut(usize, usize) -> T;
}

pub trait MatrixViewMut<T>: MatrixView<T> {

    fn at_mut(&mut self, row: usize, col: usize) -> &mut T;
    fn swap(&mut self, fst: (usize, usize), snd: (usize, usize));
}

pub trait LifetimeMatrixMutRowIter<'a, T>: MatrixViewMut<T> {

    type RowRef: VectorViewMut<T>;
    type RowIter: Iterator<Item = Self::RowRef>;

    fn rows_mut(&'a mut self) -> Self::RowIter;

    fn get_row_mut(&'a mut self, index: usize) -> Self::RowRef {
        self.assert_row_in_range(index);
        self.rows_mut().nth(index).unwrap()
    }
}

pub trait MatrixMutRowIter<T>: for<'a> LifetimeMatrixMutRowIter<'a, T> {}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MatrixOwned<T> {
    cols: usize,
    data: Box<[T]>,
}

impl<T> MatrixView<T> for MatrixOwned<T> {
    
    fn row_count(&self) -> usize {
        self.data.len() / self.cols
    }

    fn col_count(&self) -> usize {
        self.cols
    }

    fn at(&self, row: usize, col: usize) -> &T {
        self.assert_col_in_range(col);
        self.assert_row_in_range(row);
        &self.data[row * self.cols + col]
    }

    fn to_owned(self) -> MatrixOwned<T> {
        self
    }
}

impl<T> MatrixViewMut<T> for MatrixOwned<T> {
    fn at_mut(&mut self, row: usize, col: usize) -> &mut T {
        self.assert_col_in_range(col);
        self.assert_row_in_range(row);
        &mut self.data[row * self.cols + col]
    }

    fn swap(&mut self, fst: (usize, usize), snd: (usize, usize)) {
        self.assert_row_in_range(fst.0);
        self.assert_row_in_range(snd.0);
        self.assert_col_in_range(fst.1);
        self.assert_col_in_range(snd.1);
        if fst == snd {
            return;
        }
        self.data.swap(fst.1 + fst.0 * self.cols, snd.1 + snd.0 * self.cols);
    }
}


impl<T> MatrixOwned<T> {

    pub fn from_data(data: Box<[T]>, cols: usize) -> MatrixOwned<T> {
        assert!(data.len() > 0, "Cannot create matrix with zero elements");
        assert!(
            data.len() % cols == 0,
            "Data length must be a multiple of column count, but got {} and {}",
            data.len(),
            cols
        );
        MatrixOwned {
            data: data,
            cols: cols
        }
    }

    pub fn from_array<const R: usize, const C: usize>(
        array: [[T; C]; R]
    ) -> MatrixOwned<T> 
    {
        let data = std::array::IntoIter::new(array).flat_map(|row| 
            std::array::IntoIter::new(row)
        ).collect::<Vec<T>>().into_boxed_slice();
        Self::from_data(data, C)
    }
}

impl<T> FromFnCreateable<T> for MatrixOwned<T> {

    fn from_fn<F>(rows: usize, cols: usize, mut f: F) -> Self
        where F: FnMut(usize, usize) -> T
    {
        let mut data = Vec::with_capacity(rows * cols);
        for row in 0..rows {
            data.extend((0..cols).map(|col| f(row, col)));
        }
        return Self::from_data(data.into_boxed_slice(), cols);
    }
}

#[derive(Debug)]
pub struct OwnedMatrixRowMutRef<'b, T> {
    data: &'b mut [T]
}

impl<'b, T> VectorView<T> for OwnedMatrixRowMutRef<'b, T> {
    
    fn len(&self) -> usize {
        self.data.len()
    }

    fn at(&self, index: usize) -> &T {
        &self.data[index]
    }
}

impl<'b, T> VectorViewMut<T> for OwnedMatrixRowMutRef<'b, T> {
    
    fn at_mut(&mut self, index: usize) -> &mut T {
        &mut self.data[index]
    }

    fn swap(&mut self, i: usize, j: usize) {
        self.assert_in_range(i);
        self.assert_in_range(j);
        self.data.swap(i, j);
    }
}

pub struct OwnedMatrixRowMutIter<'b, T> {
    rows: std::slice::ChunksExactMut<'b, T>
}

impl<'b, T> Iterator for OwnedMatrixRowMutIter<'b, T> {
    type Item = OwnedMatrixRowMutRef<'b, T>;

    fn next(&mut self) -> Option<Self::Item> {
        self.rows.next().map(|r| OwnedMatrixRowMutRef {
            data: r
        })
    }
}

impl<'b, T: 'b> LifetimeMatrixMutRowIter<'b, T> for MatrixOwned<T> {

    type RowRef = OwnedMatrixRowMutRef<'b, T>;
    type RowIter = OwnedMatrixRowMutIter<'b, T>;
    
    fn rows_mut(&'b mut self) -> OwnedMatrixRowMutIter<'b, T> {
        OwnedMatrixRowMutIter {
            rows: self.data.chunks_exact_mut(self.col_count())
        }
    }
}

impl<T: 'static> MatrixMutRowIter<T> for MatrixOwned<T> {}

#[test]
fn test_from_fn() {
    let a = MatrixOwned::from_array([[0, 1], [1, 2]]);
    let b = MatrixOwned::from_fn(2, 2, |i, j| i + j);
    assert_eq!(a, b);
}

#[test]
fn test_mut_row_iter() {
    let mut a = MatrixOwned::from_array([[0, 1], [3, 4]]);
    let mut it = a.rows_mut();
    let mut r1 = it.next().unwrap();
    let mut r2 = it.next().unwrap();
    assert!(it.next().is_none());
    assert_eq!(0, *r1.at_mut(0));
    assert_eq!(4, *r2.at_mut(1));
    *r1.at_mut(1) = 10;
    assert_eq!(10, *a.at(0, 1));
}