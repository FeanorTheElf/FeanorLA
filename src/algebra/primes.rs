use super::super::alg::*;
use super::factoring;
use super::bigint::*;
use super::zn::*;

use std::collections::HashMap;
use std::collections::hash_map::DefaultHasher;
use std::hash::{Hash, Hasher};

use oorandom;

///
/// Miller-Rabin primality test.
/// 
/// If n is a prime, this returns true.
/// If n is not a prime, this returns false with probability greater or 
/// equal than 1 - 4^k
/// 
/// Complexity O(k log(n)^3)
/// 
pub fn is_prime(n: &BigInt, k: usize) -> bool {
    let mut hasher = DefaultHasher::new();
    n.hash(&mut hasher);
    let mut rng = oorandom::Rand32::new(hasher.finish());
    let n_minus_one = n.clone() - 1;
    let s = n_minus_one.highest_dividing_power_of_two();
    let d = n_minus_one >> s;
    let ring = FactorRingZ::new(n.clone());

    // Admitted, there is no calculation behind this choice
    const STATISTICAL_DISTANCE_ERROR_BOUND: usize = 5;

    for _i in 0..k {
        let a = ring.project(BigInt::get_uniformly_random(
            || ((rng.rand_u32() as u64) << 32) | (rng.rand_u32() as u64), &n, 
            STATISTICAL_DISTANCE_ERROR_BOUND
        ));
        let mut current = ring.pow_big(&a, &d);
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
    // honestly, this is much too small, especially given my very slow implementation of the QS
    // however, now that it exists, I also want to use it :)
    // and all in all, the whole function is really slow
    const QUADRATIC_SIEVE_BOUND: i64 = 1000000000000;
    const IS_PRIME_ERROR_BOUND: usize = 8;

    n = n.abs();
    
    if n < QUADRATIC_SIEVE_BOUND {
        let mut n_int = n.to_int().unwrap();
        let potential_divisors = factoring::gen_primes((n_int as f64).sqrt() as i64 + 1);
        let mut result = Vec::new();
        for p in potential_divisors {
            while n_int % p == 0 {
                n_int /= p;
                result.push(BigInt::from(p));
            }
        }
        if n_int != 1 {
            result.push(BigInt::from(n_int));
        }
        return result;
    } else {
        if is_prime(&n, IS_PRIME_ERROR_BOUND) {
            return vec![n];
        } else {
            for i in 2..n.log2_floor() {
                if n.clone().root_floor(i).pow(i as u32) == n {
                    let root = n.root_floor(i);
                    let factors = factor(root);
                    return factors.iter().cycle().cloned().take(factors.len() * i).collect();
                }
            }
            let first_factor = factoring::sieve::quadratic_sieve(&n);
            let other_factor = n.euclidean_div_rem(&first_factor);
            debug_assert!(first_factor != 1 && other_factor != 1);
            let mut result = factor(first_factor);
            result.extend(factor(other_factor).into_iter());
            return result;
        }
    }
}

pub fn factor_grouped(n: BigInt) -> HashMap<BigInt, usize> {
    let factors = factor(n);
    let mut result = HashMap::new();
    for p in factors {
        *result.entry(p).or_insert(0) += 1;
    }
    return result;
}

#[test]
fn test_is_prime() {
    assert_eq!(true, is_prime(&BigInt::from(97), 5));
    assert_eq!(true, is_prime(&BigInt::from(65537), 5));
    assert_eq!(false, is_prime(&BigInt::from(91), 5));
    assert_eq!(false, is_prime(&BigInt::from(3589), 5));
    assert_eq!(true, is_prime(&"125347695556857218151067929142222014224821".parse::<BigInt>().unwrap(), 5));
    assert_eq!(false, is_prime(&"3724981".parse::<BigInt>().unwrap(), 5));
    assert_eq!(false, is_prime(&"3724981832962".parse::<BigInt>().unwrap(), 5));
    let f6 = BigInt::power_of_two(64) + 1;
    assert_eq!(false, is_prime(&f6, 5));
}

#[test]
fn test_factor() {
    assert_eq!(vec![
        BigInt::from(2), 
        BigInt::from(2), 
        BigInt::from(3), 
        BigInt::from(37)
    ], factor(BigInt::from(12 * 37)));

    assert_eq!(vec![
        BigInt::from(641), 
        BigInt::from(6700417)
    ], factor("4294967297".parse::<BigInt>().unwrap()));

    assert_eq!(vec![
        BigInt::from(237689), 
        BigInt::from(717653)
    ], factor(BigInt::from(237689) * BigInt::from(717653)));
}

#[test]
fn test_factor_perfect_power() {
    let n = BigInt::from(7681).pow(32);
    let factors = factor(n);
    assert_eq!(std::iter::repeat(7681).map(BigInt::from).take(32).collect::<Vec<_>>(), factors);
}

#[cfg(release)]
#[test]
fn test_factor_release() {
    let f6 = BigInt::power_of_two(64) + 1;
    assert_eq!(vec![
        BigInt::from(67280421310721), 
        BigInt::from(274177)
    ], factor(f6));
}