use super::mat::*;
use super::super::prelude::*;

pub trait MatFn<T>: Sized {

    fn row_count(&self) -> usize;
    fn col_count(&self) -> usize;
    fn at(&self, i: usize, j: usize) -> T;

    fn add_to<R, M>(self, dst: &mut Matrix<M, T>, ring: &R)
        where M: MatrixViewMut<T>, R: Ring<El = T>
    {
        assert_eq!(self.row_count(), dst.row_count());
        assert_eq!(self.col_count(), dst.col_count());
        for row in 0..self.row_count() {
            for col in 0..self.col_count() {
                ring.add_assign(dst.at_mut(row, col), self.at(row, col));
            }
        }
    }

    fn assign_to<M>(self, dst: &mut Matrix<M, T>)
        where M: MatrixViewMut<T>
    {
        assert_eq!(self.row_count(), dst.row_count());
        assert_eq!(self.col_count(), dst.col_count());
        for row in 0..self.row_count() {
            for col in 0..self.col_count() {
                *dst.at_mut(row, col) = self.at(row, col);
            }
        }
    }

    fn compute(self) -> Matrix<MatrixOwned<T>, T> {
        Matrix::from_fn(self.row_count(), self.col_count(), |i, j| self.at(i, j))
    }

    fn scaled<R>(self, coeff: El<R>, ring: R) -> MatrixScaled<R, Self> 
        where R: Ring<El = T>
    {
        MatrixScaled::new(self, coeff, ring)
    }
}

impl<M, T> MatFn<T> for Matrix<M, T>
    where T: Clone, M: MatrixView<T>
{
    fn row_count(&self) -> usize {
        Matrix::<M, T>::row_count(self)
    }

    fn col_count(&self) -> usize {
        Matrix::<M, T>::col_count(self)
    }

    fn at(&self, i: usize, j: usize) -> T {
        Matrix::<M, T>::at(self, i, j).clone()
    }

    default fn add_to<R, N>(self, dst: &mut Matrix<N, T>, ring: &R)
        where N: MatrixViewMut<T>, R: Ring<El = T>
    {
        assert_eq!(self.row_count(), dst.row_count());
        assert_eq!(self.col_count(), dst.col_count());
        for row in 0..self.row_count() {
            for col in 0..self.col_count() {
                ring.add_assign_ref(dst.at_mut(row, col), self.at(row, col));
            }
        }
    }

    fn compute(self) -> Matrix<MatrixOwned<T>, T> {
        self.into_owned()
    }
}

pub struct MatrixSum<R, M1, M2>
    where R: Ring, M1: MatFn<El<R>>, M2: MatFn<El<R>>
{
    ring: R,
    lhs: M1,
    rhs: M2
}

impl<R, M1, M2> MatrixSum<R, M1, M2>
    where R: Ring, M1: MatFn<El<R>>, M2: MatFn<El<R>>
{
    pub fn new(lhs: M1, rhs: M2, ring: R) -> Self {
        assert_eq!(lhs.col_count(), rhs.row_count());
        MatrixSum { ring, lhs, rhs }
    }
}

impl<R, M1, M2> MatFn<El<R>> for MatrixSum<R, M1, M2>
    where R: Ring, M1: MatFn<El<R>>, M2: MatFn<El<R>>
{
    fn row_count(&self) -> usize {
        self.lhs.row_count()
    }

    fn col_count(&self) -> usize {
        self.lhs.col_count()
    }

    fn at(&self, i: usize, j: usize) -> El<R> {
        self.ring.add(self.lhs.at(i, j), self.rhs.at(i, j))
    }
    
    fn add_to<R2, M>(self, dst: &mut Matrix<M, El<R>>, ring: &R2)
        where M: MatrixViewMut<El<R>>, R2: Ring<El = El<R>>
    {
        assert_eq!(self.row_count(), dst.row_count());
        assert_eq!(self.col_count(), dst.col_count());
        self.lhs.add_to(dst, ring);
        self.rhs.add_to(dst, ring);
    }

    fn assign_to<M>(self, dst: &mut Matrix<M, El<R>>)
        where M: MatrixViewMut<El<R>>
    {
        assert_eq!(self.row_count(), dst.row_count());
        assert_eq!(self.col_count(), dst.col_count());
        self.lhs.assign_to(dst);
        self.rhs.add_to(dst, &self.ring);
    }
}

pub struct MatrixNeg<R, M> 
    where R: Ring, M: MatFn<El<R>>
{
    ring: R,
    val: M
}

impl<R, M> MatrixNeg<R, M>
    where R: Ring, M: MatFn<El<R>>
{
    pub fn new(val: M, ring: R) -> Self {
        MatrixNeg { ring, val }
    }
}

impl<R, M> MatFn<El<R>> for MatrixNeg<R, M>
    where R: Ring, M: MatFn<El<R>>
{
    fn row_count(&self) -> usize {
        self.val.row_count()
    }

    fn col_count(&self) -> usize {
        self.val.col_count()
    }

    fn at(&self, i: usize, j: usize) -> El<R> {
        self.ring.neg(self.val.at(i, j))
    }
}

pub struct MatrixProd<R, M1, M2>
    where R: Ring, M1: MatrixView<El<R>>, M2: MatrixView<El<R>>
{
    ring: R,
    lhs: Matrix<M1, El<R>>,
    rhs: Matrix<M2, El<R>>
}

impl<R, M1, M2> MatrixProd<R, M1, M2>
    where R: Ring, M1: MatrixView<El<R>>, M2: MatrixView<El<R>>
{
    pub fn new(lhs: Matrix<M1, El<R>>, rhs: Matrix<M2, El<R>>, ring: R) -> Self {
        assert_eq!(lhs.col_count(), rhs.row_count());
        MatrixProd { ring, lhs, rhs }
    }
}

impl<R, M1, M2> MatFn<El<R>> for MatrixProd<R, M1, M2>
    where R: Ring, M1: MatrixView<El<R>>, M2: MatrixView<El<R>>
{
    fn row_count(&self) -> usize {
        self.lhs.row_count()
    }

    fn col_count(&self) -> usize {
        self.rhs.col_count()
    }

    fn at(&self, i: usize, j: usize) -> El<R> {
        self.ring.sum(
            (0..self.lhs.col_count()).map(|k| self.ring.mul_ref(self.lhs.at(i, k), self.rhs.at(k, j)))
        )
    }
}

pub struct MatrixScaled<R, M> 
    where R: Ring, M: MatFn<El<R>>
{
    ring: R,
    coeff: El<R>,
    val: M
}

impl<R, M> MatrixScaled<R, M>
    where R: Ring, M: MatFn<El<R>>
{
    pub fn new(val: M, coeff: El<R>, ring: R) -> Self {
        MatrixScaled { ring, coeff, val }
    }
}

impl<R, M> MatFn<El<R>> for MatrixScaled<R, M>
    where R: Ring, M: MatFn<El<R>>
{
    fn row_count(&self) -> usize {
        self.val.row_count()
    }

    fn col_count(&self) -> usize {
        self.val.col_count()
    }

    fn at(&self, i: usize, j: usize) -> El<R> {
        self.ring.mul_ref(&self.val.at(i, j), &self.coeff)
    }
}
