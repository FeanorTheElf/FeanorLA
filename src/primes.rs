use super::sieve::*;
use super::alg::*;
use super::bigint::*;
use super::quotient::*;

use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};

use oorandom;

///
/// Miller-Rabin primality test.
/// 
/// If n is a prime, this returns true.
/// If n is not a prime, this returns false with probability greater or equal than 1 - 4^k
/// 
pub fn is_prime(n: &BigInt, k: usize) -> bool {
    let mut hasher = DefaultHasher::new();
    n.hash(&mut hasher);
    let mut rng = oorandom::Rand32::new(hasher.finish());
    let n_minus_one = n.clone() - 1;
    let s = n_minus_one.highest_dividing_power_of_two();
    let d = n_minus_one >> s;
    let ring = QuotientRingZ::new(n.clone());

    // Admitted, there is no calculation behind this choice
    const STATISTICAL_DISTANCE_BOUND: usize = 5;

    for _i in 0..k {
        let a = ring.project(BigInt::get_uniformly_random(|| ((rng.rand_u32() as u64) << 32) | (rng.rand_u32() as u64), &n, STATISTICAL_DISTANCE_BOUND));
        let mut current = ring.pow_big(a.clone(), d.clone());
        let mut miller_rabin_condition = ring.is_one(&current);
        for _r in 0..s {
            miller_rabin_condition |= ring.is_neg_one(&current);
            if miller_rabin_condition {
                break;
            }
            current = ring.mul(current.clone(), current);
        }
        if !miller_rabin_condition {
            return false;
        }
    }
    return true;
}

pub fn factor(mut n: BigInt) -> Vec<BigInt> {
    const QUADRATIC_SIEVE_BOUND: i64 = 1000000000000;
    const IS_PRIME_ERROR_BOUND: usize = 8;

    n = n.abs();
    
    if n < QUADRATIC_SIEVE_BOUND as u64 {
        let mut n_int = n.to_int().unwrap();
        let potential_divisors = gen_primes((n_int as f64).sqrt() as i64 + 1);
        let mut result = Vec::new();
        for p in potential_divisors {
            while n_int % p == 0 {
                n_int /= p;
                result.push(BigInt::from(p));
            }
        }
        if n != 1u64 {
            result.push(BigInt::from(n_int));
        }
        return result;
    } else {
        if is_prime(&n, IS_PRIME_ERROR_BOUND) {
            return vec![n];
        } else {
            let first_factor = quadratic_sieve(&n);
            let other_factor = n.div_rem(&first_factor);
            debug_assert!(first_factor != 1u64 && other_factor != 1u64);
            let mut result = factor(first_factor);
            result.extend(factor(other_factor).into_iter());
            return result;
        }
    }
}

#[test]
fn test_is_prime() {
    assert_eq!(true, is_prime(&BigInt::from(97u64), 5));
    assert_eq!(true, is_prime(&BigInt::from(65537u64), 5));
    assert_eq!(false, is_prime(&BigInt::from(91u64), 5));
    assert_eq!(false, is_prime(&BigInt::from(3589u64), 5));
    assert_eq!(true, is_prime(&"125347695556857218151067929142222014224821".parse::<BigInt>().unwrap(), 5));
    assert_eq!(false, is_prime(&"3724981".parse::<BigInt>().unwrap(), 5));
    assert_eq!(false, is_prime(&"3724981832962".parse::<BigInt>().unwrap(), 5));
    let f6 = BigInt::power_of_two(64) + 1;
    assert_eq!(false, is_prime(&f6, 5));
}

#[test]
fn test_factor() {
    assert_eq!(vec![BigInt::from(2u64), BigInt::from(2u64), BigInt::from(3u64), BigInt::from(37u64)], factor(BigInt::from(12 * 37 as u64)));
    assert_eq!(vec![BigInt::from(641u64), BigInt::from(6700417u64)], factor(BigInt::from(4294967297u64)));
    let f6 = BigInt::power_of_two(64) + 1;
    assert_eq!(vec![BigInt::from(67280421310721u64), BigInt::from(274177u64)], factor(f6));
    assert_eq!(vec![BigInt::from(237689u64), BigInt::from(717653u64)], factor(BigInt::from(237689 * 717653 as u64)));
}