use super::super::alg::*;
use super::bigint::*;
use super::zn::*;
use super::factoring_algorithms;

impl<T> DivisibilityInformationRing for StaticRingImpl<RingAxiomsField, T> 
    where T: FieldEl
{
    fn is_divisibility_computable(&self) -> bool { 
        true
    }

    fn is_divisible_by(&self, lhs: &Self::El, rhs: &Self::El) -> bool {
        self.is_zero(lhs) || !self.is_zero(rhs)
    }

    fn quotient(&self, lhs: &Self::El, rhs: &Self::El) -> Option<Self::El> {
        if self.is_divisible_by(lhs, rhs) {
            Some(self.div(lhs.clone(), rhs))
        } else {
            None
        }
    }

    fn is_unit(&self, el: &Self::El) -> bool {
        !self.is_zero(el)
    }
}

impl DivisibilityInformationRing for StaticRing<i32> {
    
    fn is_divisibility_computable(&self) -> bool { true }

    fn quotient(&self, lhs: &Self::El, rhs: &Self::El) -> Option<Self::El> { 
        if lhs % rhs != 0 {
            None
        } else {
            Some(lhs / rhs)
        }
    }

    fn is_unit(&self, el: &Self::El) -> bool {
        *el == 1 || *el == -1
    }
}

impl DivisibilityInformationRing for StaticRing<i64> {
    
    fn is_divisibility_computable(&self) -> bool { true }

    fn quotient(&self, lhs: &Self::El, rhs: &Self::El) -> Option<Self::El> { 
        if lhs % rhs != 0 {
            None
        } else {
            Some(lhs / rhs)
        }
    }

    fn is_unit(&self, el: &Self::El) -> bool {
        *el == 1 || *el == -1
    }
}

use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};

use oorandom;

///
/// Miller-Rabin primality test.
/// 
/// If n is a prime, this returns true.
/// If n is not a prime, this returns false with probability greater or 
/// equal than 1 - 4^(-k).
/// 
/// Complexity O(k log(n)^3)
/// 
/// # Randomness
/// 
/// Note that the randomness used for this function is derived only from
/// the input, hence it will always yield the same output on the same input.
/// Technically, it follows that the probability of a wrong output is greater
/// than 4^(-k) on some outputs (as it is either 0 or 1), but of course
/// this is not helpful. To be completely precise: If the seed of the used
/// PRNG would be random, then the probability of a wrong output is at 
/// most 4^(-k).
/// 
#[allow(non_snake_case)]
pub fn miller_rabin(n: &BigInt, k: usize) -> bool {
    if *n <= 2 {
        return *n == 2;
    }

    let mut hasher = DefaultHasher::new();
    n.hash(&mut hasher);
    let mut rng = oorandom::Rand32::new(hasher.finish());
    let n_minus_one = n.clone() - 1;
    let s = n_minus_one.highest_dividing_power_of_two();
    let d = n_minus_one.clone() >> s;
    let Z_nZ = FactorRingZ::new(n.clone());

    // Admitted, there is no calculation behind this choice
    const STATISTICAL_DISTANCE_ERROR_BOUND: usize = 5;

    for _i in 0..k {
        let a = Z_nZ.project(BigInt::get_uniformly_random_oorandom(
            &mut rng, &n_minus_one, STATISTICAL_DISTANCE_ERROR_BOUND
        ) + 1);
        let mut current = Z_nZ.pow_big(&a, &d);
        let mut miller_rabin_condition = Z_nZ.is_one(&current);
        for _r in 0..s {
            miller_rabin_condition |= Z_nZ.is_neg_one(&current);
            if miller_rabin_condition {
                break;
            }
            current = Z_nZ.mul(current.clone(), current);
        }
        if Z_nZ.is_zero(&current) || !miller_rabin_condition {
            return false;
        }
    }
    return true;
}

impl DivisibilityInformationRing for BigIntRing {

    fn quotient(&self, lhs: &Self::El, rhs: &Self::El) -> Option<BigInt> {
        let (quo, rem) = self.euclidean_div_rem(lhs.clone(), rhs);
        if rem == 0 {
            return Some(quo);
        } else {
            return None;
        }
    }

    fn is_divisibility_computable(&self) -> bool {
        true
    }

    fn is_unit(&self, el: &Self::El) -> bool {
        self.is_one(el) || self.is_neg_one(el)
    }
}

impl FactoringInformationRing for BigIntRing {

    fn is_ufd(&self) -> bool {
        true
    }

    fn is_prime(&self, el: &Self::El) -> bool {
        miller_rabin(el, 10)
    }

    fn calc_factor(&self, el: &mut Self::El) -> Option<Self::El> {
        // honestly, this is much too small, especially given my very slow implementation of the QS
        // however, now that it exists, I also want to use it :)
        // and all in all, the whole function is really slow
        const QUADRATIC_SIEVE_BOUND: i64 = 1000000000000;
        const IS_PRIME_ERROR_BOUND: usize = 8;

        let n = el.clone().abs();
        
        if n < QUADRATIC_SIEVE_BOUND {
            let n_int = n.to_int().unwrap();
            let potential_divisors = factoring_algorithms::gen_primes((n_int as f64).sqrt() as i64 + 1);
            for p in potential_divisors {
                if n_int % p == 0 {
                    return Some(BigInt::from(p))
                }
            }
            return None;
        } else {
            if miller_rabin(&n, IS_PRIME_ERROR_BOUND) {
                return None;
            } else {
                for i in 2..n.log2_floor() {
                    if n.clone().root_floor(i).pow(i as u32) == n {
                        let root = n.root_floor(i);
                        return Some(root);
                    }
                }
                return Some(factoring_algorithms::sieve::quadratic_sieve(&n));
            }
        }
    }
}