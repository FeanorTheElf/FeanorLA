use super::*;
pub use super::super::super::{discard, gen_const_vector};

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
            
            fn create_subvector(self, from: usize, to: usize) -> Self::Subvector {
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