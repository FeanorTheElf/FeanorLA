use super::super::alg::*;

pub trait VectorView<T>: Sized {
    
    fn len(&self) -> usize;
    fn at(&self, index: usize) -> &T;

    fn assert_in_range(&self, index: usize) {
        assert!(index < self.len(), "Vector index {} out of range 0..{}", index, self.len());
    }

    fn to_owned(self) -> VectorOwned<T>
        where T: Clone
    {
        VectorOwned::from_fn(self.len(), |i| self.at(i).clone())
    }
}

pub trait VectorViewMut<T>: VectorView<T> {
    fn at_mut(&mut self, index: usize) -> &mut T;
    fn swap(&mut self, i: usize, j: usize);
}

#[derive(Debug, Clone)]
pub struct VectorOwned<T> {
    data: Box<[T]>,
}

impl<T> VectorView<T> for VectorOwned<T> {

    fn len(&self) -> usize {
        self.data.len()
    }

    fn at(&self, index: usize) -> &T {
        &self.data[index]
    }
}

impl<T> VectorViewMut<T> for VectorOwned<T> {

    fn at_mut(&mut self, index: usize) -> &mut T {
        &mut self.data[index]
    }

    fn swap(&mut self, i: usize, j: usize) {
        self.assert_in_range(i);
        self.assert_in_range(j);
        if i == j {
            return;
        }
        self.data.swap(i, j);
    }
}

impl<T> VectorOwned<T> {

    pub fn new(data: Box<[T]>) -> VectorOwned<T> {
        VectorOwned { data }
    }

    pub fn from_fn<F>(len: usize, mut f: F) -> Self
        where F: FnMut(usize) -> T
    {
        Self::new((0..len).map(|i| f(i)).collect::<Vec<_>>().into_boxed_slice())
    }

    pub fn from_array<const L: usize>(data: [T; L]) -> Self {
        VectorOwned::new(std::array::IntoIter::new(data).collect::<Vec<_>>().into_boxed_slice())
    }
}

impl<T> VectorOwned<T>
    where T: Zero 
{
    pub fn zero(len: usize) -> VectorOwned<T> {
        VectorOwned::new((0..len).map(|_| T::zero()).collect::<Vec<T>>().into_boxed_slice())
    }
}

impl<T> VectorOwned<T>
    where T: Zero + One 
{
    pub fn unit_vector(index: usize, len: usize) -> VectorOwned<T> {
        assert!(index < len, "Expected index of 1-entry in unit vector to be within 0 and {}, got {}", len, index);
        VectorOwned::new((0..len).map(|i| if i == index { T::one() } else { T::zero() }).collect::<Vec<T>>().into_boxed_slice())
    }
}