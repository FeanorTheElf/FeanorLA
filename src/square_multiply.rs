use super::prelude::*;
use super::integer::*;

pub fn abs_square_and_multiply<T, F, G, I>(base: &T, power: &El<I>, int_ring: I, mut multiply: F, mut multiply_ref: G, identity: T) -> T
    where I: IntegerRing, F: FnMut(T, T) -> T, G: FnMut(&T, &T) -> T, T: Clone
{
    if int_ring.is_zero(&power) {
        return identity;
    } else if int_ring.is_one(&power) {
        return base.clone();
    }

    let mut result = identity;
    for i in (0..(int_ring.abs_log2_floor(power) + 1)).rev() {
        if int_ring.abs_is_bit_set(power, i) {
            result = multiply(multiply_ref(&result, base), result);
        } else {
            result = multiply_ref(&result, &result);
        }
    }
    return result;
}