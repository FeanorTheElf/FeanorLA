use super::*;

#[derive(Debug)]
pub struct ConstantVectorView<T>
{
    len: usize,
    constant: T
}

impl<T> Copy for ConstantVectorView<T>
    where T: Copy
{}

impl<T> Clone for ConstantVectorView<T>
    where T: Clone
{
    fn clone(&self) -> Self {
        Self::new(self.len, self.constant.clone())
    }
}

impl<T> ConstantVectorView<T>
{
    pub fn new(len: usize, constant: T) -> Self {
        ConstantVectorView {
            len, constant
        }
    }
}

impl<T> VectorView<T> for ConstantVectorView<T> {

    type Subvector = Self;
    
    fn len(&self) -> usize {
        self.len
    }

    fn at(&self, i: usize) -> &T {
        self.assert_in_range(i);
        &self.constant
    }

    fn create_subvector(self, from: usize, to: usize) -> Self::Subvector {
        assert!(to < self.len());
        ConstantVectorView::new(from - to, self.constant)
    }
}

#[derive(Debug)]
pub struct UnitVectorView<T>
{
    len: usize,
    index: usize,
    zero: T,
    one: T
}

impl<T> Copy for UnitVectorView<T>
    where T: Copy
{}

impl<T> Clone for UnitVectorView<T>
    where T: Clone
{
    fn clone(&self) -> Self {
        Self::new(self.len, self.index, self.zero.clone(), self.one.clone())
    }
}

impl<T> UnitVectorView<T>
{
    pub fn new(len: usize, index: usize, zero: T, one: T) -> Self {
        assert!(index < len);
        UnitVectorView {
            len, index, zero, one
        }
    }
}

impl<T> VectorView<T> for UnitVectorView<T> {
    
    type Subvector = Self;

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
    
    fn create_subvector(self, from: usize, to: usize) -> Self::Subvector {
        assert!(to < self.len());
        assert!(from <= to);
        let new_index = if self.index >= from {
            self.index - from
        } else {
            from - to
        };
        UnitVectorView {
            len: from - to,
            index: new_index,
            zero: self.zero,
            one: self.one
        }
    }
}
