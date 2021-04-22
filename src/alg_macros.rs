pub use super::alg::*;

///
/// Implements the basic ring operations for some type. The implementation
/// of these ring operations taken from a ring object, that is easily obtainable
/// from each function. Usually, this will be a constant.
/// 
/// It is not recommended to use this for copyable ring elements, as it also
/// implements all value/ref-combinations for arithmetic operations.
/// 
/// In some way, this is the converse to `StaticRing`. This creates a `RingEl`-type
/// (i.e. a type for ring elements in a ring that is completely known at compile-time)
/// from a given ring object.
/// `StaticRing` on the other hand provides the ring object that represents the
/// ring of such `RingEl`-types.
/// 
macro_rules! impl_ring_el {
    ($t:ty; $ring_constant:expr; $axioms:ty) => {
        impl std::ops::Add for $t {
            type Output = $t;

            fn add(self, rhs: $t) -> Self::Output {
                ($ring_constant).add(self, rhs)
            }
        }

        impl std::ops::Add<&$t> for $t {
            type Output = $t;

            fn add(self, rhs: &$t) -> Self::Output {
                ($ring_constant).add_ref(self, rhs)
            }
        }

        impl std::ops::Add<$t> for &$t {
            type Output = $t;

            fn add(self, rhs: $t) -> Self::Output {
                ($ring_constant).add_ref(rhs, self)
            }
        }

        impl std::ops::Add<&$t> for &$t {
            type Output = $t;

            fn add(self, rhs: &$t) -> Self::Output {
                ($ring_constant).add_ref(self.clone(), rhs)
            }
        }

        impl std::ops::AddAssign for $t {
            fn add_assign(&mut self, rhs: $t) {
                take_mut::take_or_recover(self, || ($ring_constant).unspecified_element(), |v| ($ring_constant).add(v, rhs));
            }
        }

        impl std::ops::AddAssign<&$t> for $t {
            fn add_assign(&mut self, rhs: &$t) {
                take_mut::take_or_recover(self, || ($ring_constant).unspecified_element(), |v| ($ring_constant).add_ref(v, rhs));
            }
        }

        impl std::ops::Mul for $t {
            type Output = $t;

            fn mul(self, rhs: $t) -> Self::Output {
                ($ring_constant).mul(self, rhs)
            }
        }

        impl std::ops::Mul<&$t> for $t {
            type Output = $t;

            fn mul(self, rhs: &$t) -> Self::Output {
                ($ring_constant).mul_ref(&self, rhs)
            }
        }

        impl std::ops::Mul<$t> for &$t {
            type Output = $t;

            fn mul(self, rhs: $t) -> Self::Output {
                ($ring_constant).mul_ref(&rhs, self)
            }
        }

        impl std::ops::Mul<&$t> for &$t {
            type Output = $t;

            fn mul(self, rhs: &$t) -> Self::Output {
                ($ring_constant).mul_ref(self, rhs)
            }
        }

        impl std::ops::MulAssign for $t {

            fn mul_assign(&mut self, rhs: $t) {
                take_mut::take_or_recover(self, || ($ring_constant).unspecified_element(), |v| ($ring_constant).mul(v, rhs));
            }
        }

        impl std::ops::MulAssign<&$t> for $t {

            fn mul_assign(&mut self, rhs: &$t) {
                *self = ($ring_constant).mul_ref(self, rhs);
            }
        }

        impl std::ops::Sub for $t {
            type Output = $t;

            fn sub(self, rhs: $t) -> Self::Output {
                ($ring_constant).sub(self, rhs)
            }
        }

        impl std::ops::Sub<&$t> for $t {
            type Output = $t;

            fn sub(self, rhs: &$t) -> Self::Output {
                ($ring_constant).sub_ref_snd(self, rhs)
            }
        }

        impl std::ops::Sub<$t> for &$t {
            type Output = $t;

            fn sub(self, rhs: $t) -> Self::Output {
                ($ring_constant).sub_ref_fst(self, rhs)
            }
        }

        impl std::ops::Sub<&$t> for &$t {
            type Output = $t;

            fn sub(self, rhs: &$t) -> Self::Output {
                ($ring_constant).sub_ref_snd(self.clone(), rhs)
            }
        }

        impl std::ops::SubAssign for $t {

            fn sub_assign(&mut self, rhs: $t) {
                take_mut::take_or_recover(self, || ($ring_constant).unspecified_element(), |v| ($ring_constant).sub(v, rhs));
            }
        }

        impl std::ops::SubAssign<&$t> for $t {

            fn sub_assign(&mut self, rhs: &$t) {
                take_mut::take_or_recover(self, || ($ring_constant).unspecified_element(), |v| ($ring_constant).sub_ref_snd(v, rhs));
            }
        }
        impl PartialEq for $t {
            fn eq(&self, rhs: &$t) -> bool {
                ($ring_constant).eq(self, rhs)
            }
        }

        impl std::ops::Neg for $t {
            type Output = $t;

            fn neg(self) -> Self::Output {
                ($ring_constant).neg(self)
            }
        }

        impl Zero for $t {
            fn zero() -> $t {
                ($ring_constant).zero()
            }
        }

        impl One for $t {
            fn one() -> $t {
                ($ring_constant).one()
            }
        }

        impl std::fmt::Display for $t {
            fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
                ($ring_constant).format(self, f, false)
            }
        }

        impl RingEl for $t {
            type Axioms = $axioms;
        }
    };
    ($t:ty; $ring_constant:expr; $axioms:ty; $($t_gen_constraits:tt)*) => {

        impl <$($t_gen_constraits)*> std::ops::Add for $t {
            type Output = $t;

            fn add(self, rhs: $t) -> Self::Output {
                ($ring_constant).add(self, rhs)
            }
        }

        impl <$($t_gen_constraits)*> std::ops::Add<&$t> for $t {
            type Output = $t;

            fn add(self, rhs: &$t) -> Self::Output {
                ($ring_constant).add_ref(self, rhs)
            }
        }

        impl <$($t_gen_constraits)*> std::ops::Add<$t> for &$t {
            type Output = $t;

            fn add(self, rhs: $t) -> Self::Output {
                ($ring_constant).add_ref(rhs, self)
            }
        }

        impl <$($t_gen_constraits)*> std::ops::Add<&$t> for &$t {
            type Output = $t;

            fn add(self, rhs: &$t) -> Self::Output {
                ($ring_constant).add_ref(self.clone(), rhs)
            }
        }

        impl <$($t_gen_constraits)*> std::ops::AddAssign for $t {
            fn add_assign(&mut self, rhs: $t) {
                take_mut::take_or_recover(self, || ($ring_constant).unspecified_element(), |v| ($ring_constant).add(v, rhs));
            }
        }

        impl <$($t_gen_constraits)*> std::ops::AddAssign<&$t> for $t {
            fn add_assign(&mut self, rhs: &$t) {
                take_mut::take_or_recover(self, || ($ring_constant).unspecified_element(), |v| ($ring_constant).add_ref(v, rhs));
            }
        }

        impl <$($t_gen_constraits)*> std::ops::Mul for $t {
            type Output = $t;

            fn mul(self, rhs: $t) -> Self::Output {
                ($ring_constant).mul(self, rhs)
            }
        }

        impl <$($t_gen_constraits)*> std::ops::Mul<&$t> for $t {
            type Output = $t;

            fn mul(self, rhs: &$t) -> Self::Output {
                ($ring_constant).mul_ref(&self, rhs)
            }
        }

        impl <$($t_gen_constraits)*> std::ops::Mul<$t> for &$t {
            type Output = $t;

            fn mul(self, rhs: $t) -> Self::Output {
                ($ring_constant).mul_ref(&rhs, self)
            }
        }

        impl <$($t_gen_constraits)*> std::ops::Mul<&$t> for &$t {
            type Output = $t;

            fn mul(self, rhs: &$t) -> Self::Output {
                ($ring_constant).mul_ref(self, rhs)
            }
        }

        impl <$($t_gen_constraits)*> std::ops::MulAssign for $t {

            fn mul_assign(&mut self, rhs: $t) {
                take_mut::take_or_recover(self, || ($ring_constant).unspecified_element(), |v| ($ring_constant).mul(v, rhs));
            }
        }

        impl <$($t_gen_constraits)*> std::ops::MulAssign<&$t> for $t {

            fn mul_assign(&mut self, rhs: &$t) {
                *self = ($ring_constant).mul_ref(self, rhs);
            }
        }

        impl <$($t_gen_constraits)*> std::ops::Sub for $t {
            type Output = $t;

            fn sub(self, rhs: $t) -> Self::Output {
                ($ring_constant).sub(self, rhs)
            }
        }

        impl <$($t_gen_constraits)*> std::ops::Sub<&$t> for $t {
            type Output = $t;

            fn sub(self, rhs: &$t) -> Self::Output {
                ($ring_constant).sub_ref_snd(self, rhs)
            }
        }

        impl <$($t_gen_constraits)*> std::ops::Sub<$t> for &$t {
            type Output = $t;

            fn sub(self, rhs: $t) -> Self::Output {
                ($ring_constant).sub_ref_fst(self, rhs)
            }
        }

        impl <$($t_gen_constraits)*> std::ops::Sub<&$t> for &$t {
            type Output = $t;

            fn sub(self, rhs: &$t) -> Self::Output {
                ($ring_constant).sub_ref_snd(self.clone(), rhs)
            }
        }

        impl <$($t_gen_constraits)*> std::ops::SubAssign for $t {

            fn sub_assign(&mut self, rhs: $t) {
                take_mut::take_or_recover(self, || ($ring_constant).unspecified_element(), |v| ($ring_constant).sub(v, rhs));
            }
        }

        impl <$($t_gen_constraits)*> std::ops::SubAssign<&$t> for $t {

            fn sub_assign(&mut self, rhs: &$t) {
                take_mut::take_or_recover(self, || ($ring_constant).unspecified_element(), |v| ($ring_constant).sub_ref_snd(v, rhs));
            }
        }

        impl <$($t_gen_constraits)*> PartialEq for $t {
            fn eq(&self, rhs: &$t) -> bool {
                ($ring_constant).eq(self, rhs)
            }
        }

        impl <$($t_gen_constraits)*> std::ops::Neg for $t {
            type Output = $t;

            fn neg(self) -> Self::Output {
                ($ring_constant).neg(self)
            }
        }

        impl <$($t_gen_constraits)*> Zero for $t {
            fn zero() -> $t {
                ($ring_constant).zero()
            }
        }

        impl <$($t_gen_constraits)*> One for $t {
            fn one() -> $t {
                ($ring_constant).one()
            }
        }

        impl <$($t_gen_constraits)*> std::fmt::Display for $t {
            fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
                ($ring_constant).format(self, f, false)
            }
        }

        impl <$($t_gen_constraits)*> RingEl for $t {
            type Axioms = $axioms;
        }
    };
}

///
/// See `impl_ring_el!`. This should be used additionally the former macro.
/// 
macro_rules! impl_euclidean_ring_el {
    ($t:ty; $ring_constant:expr) => {
        
        impl std::ops::Div for $t {
            type Output = $t;

            fn div(self, rhs: $t) -> Self::Output {
                ($ring_constant).euclidean_div(self, &rhs)
            }
        }

        impl std::ops::Rem for $t {
            type Output = $t;

            fn rem(self, rhs: $t) -> Self::Output {
                ($ring_constant).euclidean_rem(self, &rhs)
            }
        }

        impl std::ops::RemAssign for $t {
            fn rem_assign(&mut self, rhs: $t) {
                take_mut::take_or_recover(self, || ($ring_constant).unspecified_element(), |v| ($ring_constant).euclidean_rem(v, &rhs));
            }
        }

        impl std::ops::DivAssign for $t {
            fn div_assign(&mut self, rhs: $t) {
                take_mut::take_or_recover(self, || ($ring_constant).unspecified_element(), |v| ($ring_constant).euclidean_div(v, &rhs));
            }
        }

        impl EuclideanRingEl for $t {

            fn div_rem(&mut self, rhs: $t) -> $t {
                let mut result: Option<$t> = None;
                take_mut::take_or_recover(self, || ($ring_constant).unspecified_element(), |v| {
                    let (quo, rem) = ($ring_constant).euclidean_div_rem(v, &rhs);
                    result = Some(quo);
                    return rem;
                });
                return result.unwrap();
            }
        }
    };
    ($t:ty; $ring_constant:expr; $($t_gen_constraits:tt)*) => {
        
        impl <$($t_gen_constraits)*> std::ops::Div for $t {
            type Output = $t;

            fn div(self, rhs: $t) -> Self::Output {
                ($ring_constant).euclidean_div(self, &rhs)
            }
        }

        impl <$($t_gen_constraits)*> std::ops::Rem for $t {
            type Output = $t;

            fn rem(self, rhs: $t) -> Self::Output {
                ($ring_constant).euclidean_rem(self, &rhs)
            }
        }

        impl <$($t_gen_constraits)*> std::ops::RemAssign for $t {
            fn rem_assign(&mut self, rhs: $t) {
                take_mut::take_or_recover(self, || ($ring_constant).unspecified_element(), |v| ($ring_constant).euclidean_rem(v, &rhs));
            }
        }

        impl <$($t_gen_constraits)*> std::ops::DivAssign for $t {
            fn div_assign(&mut self, rhs: $t) {
                take_mut::take_or_recover(self, || ($ring_constant).unspecified_element(), |v| ($ring_constant).euclidean_div(v, &rhs));
            }
        }

        impl <$($t_gen_constraits)*> EuclideanRingEl for $t {

            fn div_rem(&mut self, rhs: $t) -> $t {
                let mut result: Option<$t> = None;
                take_mut::take_or_recover(self, || ($ring_constant).unspecified_element(), |v| {
                    let (quo, rem) = ($ring_constant).euclidean_div_rem(v, &rhs);
                    result = Some(quo);
                    return rem;
                });
                return result.unwrap();
            }
        }
    };
}

///
/// See `impl_ring_el!`. This should be used additionally the former macro.
/// 
macro_rules! impl_field_ring_el {
    ($t:ty; $ring_constant:expr) => {
        
        impl std::ops::Div for $t {
            type Output = $t;

            fn div(self, rhs: $t) -> Self::Output {
                ($ring_constant).euclidean_div(self, &rhs)
            }
        }

        impl std::ops::DivAssign for $t {
            fn div_assign(&mut self, rhs: $t) {
                take_mut::take_or_recover(self, || ($ring_constant).unspecified_element(), |v| ($ring_constant).euclidean_div(v, &rhs));
            }
        }

        impl FieldEl for $t {}
    };
    ($t:ty; $ring_constant:expr; $($t_gen_constraits:tt)*) => {
        
        impl <$($t_gen_constraits)*> std::ops::Div for $t {
            type Output = $t;

            fn div(self, rhs: $t) -> Self::Output {
                ($ring_constant).euclidean_div(self, &rhs)
            }
        }

        impl <$($t_gen_constraits)*> std::ops::DivAssign for $t {
            fn div_assign(&mut self, rhs: $t) {
                take_mut::take_or_recover(self, || ($ring_constant).unspecified_element(), |v| ($ring_constant).euclidean_div(v, &rhs));
            }
        }

        impl <$($t_gen_constraits)*> FieldEl for $t {}
    };
}