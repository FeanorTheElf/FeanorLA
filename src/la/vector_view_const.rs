pub use super::vector_view::*;
pub use super::super::{discard, gen_const_vector};

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

    fn subvector(self, from: usize, to: usize) -> Self::Subvector {
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
    
    fn subvector(self, from: usize, to: usize) -> Self::Subvector {
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

#[macro_export]
macro_rules! discard {
    ($value:expr) => {
        0
    };
}

#[macro_export]
macro_rules! gen_const_vector {
    ($name:ident; $value_type:ty; $($value:ident),*) => {
        #[derive(Debug, Clone, Copy)]
        pub struct $name<$(const $value: $value_type),*> {}

        impl<$(const $value: $value_type),*> $name<$($value),*> {

            pub const LEN: usize = [$(discard!($value)),*].len();
            // using LEN in the below expression results in a warning that might become a hard error later
            // I guess this is a compiler bug, here we directly inline the expression 
            // `[$(discard!($value)),*].len()` as a workaround
            pub const DATA: [$value_type; [$(discard!($value)),*].len()] = [$($value),*];

            pub const INSTANCE: Self = Self {};
        }

        impl<$(const $value: $value_type),*> VectorView<$value_type> for $name<$($value),*> {

            type Subvector = Subvector<Self, $value_type>;

            fn len(&self) -> usize {
                Self::LEN
            }

            fn at(&self, i: usize) -> &$value_type {
                &Self::DATA[i]
            }

            fn subvector(self, from: usize, to: usize) -> Self::Subvector {
                Subvector::new(from, to, self)
            }
        }
    };
}

gen_const_vector!(ConstVector1; i64; V0);
gen_const_vector!(ConstVector2; i64; V0, V1);
gen_const_vector!(ConstVector3; i64; V0, V1, V2);
gen_const_vector!(ConstVector4; i64; V0, V1, V2, V3);
gen_const_vector!(ConstVector5; i64; V0, V1, V2, V3, V5);
gen_const_vector!(ConstVector6; i64; V0, V1, V2, V3, V4, V5);
gen_const_vector!(ConstVector7; i64; V0, V1, V2, V3, V4, V5, V6);
gen_const_vector!(ConstVector8; i64; V0, V1, V2, V3, V4, V5, V6, V7);

#[test]
fn test_create_const_vector() {
    type V = ConstVector4<0, 1, 2, 3>;
    assert_eq!(1, V::DATA[1]);
}