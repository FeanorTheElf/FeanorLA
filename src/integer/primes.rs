use super::super::fq::zn_big::*;
use super::super::ring::*;
use super::super::combinatorics::iters::*;
use super::super::integer::*;
use super::super::square_multiply::*;
use super::super::factoring_algorithms;

use std::collections::hash_map::DefaultHasher;
use std::hash::Hasher;

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
pub fn miller_rabin<I>(ring: I, n: &El<I>, k: usize) -> bool 
    where I: IntegerRing + HashableElRing
{
    if ring.is_leq(n, &ring.from_z(2)) {
        return ring.is_eq(n, &ring.from_z(2));
    }

    let mut hasher = DefaultHasher::new();
    ring.hash(&mut hasher, n);
    let mut rng = oorandom::Rand32::new(hasher.finish());
    let n_minus_one = ring.sub_ref_fst(n, ring.one());
    let s = ring.highest_dividing_power_of_two(&n_minus_one);
    let d = ring.euclidean_div_pow_2(n_minus_one.clone(), s as u64);
    let Z_nZ = Zn::new(&ring, n.clone());

    for _i in 0..k {
        let a = Z_nZ.project(ring.add(ring.get_uniformly_random_oorandom(&mut rng, &n_minus_one), ring.one()));
        let mut current = abs_square_and_multiply(&a, &d, &ring, |a, b| Z_nZ.mul(a, b), |a, b| Z_nZ.mul_ref(a, b), Z_nZ.one());
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

pub fn calc_factor(el: &BigInt) -> Option<BigInt> {
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
        if miller_rabin(&BigInt::RING, &n, IS_PRIME_ERROR_BOUND) {
            return None;
        } else {
            for i in 2..BigInt::RING.abs_log2_floor(&n) {
                if BigInt::RING.root_floor(&n, i).pow(i as u32) == n {
                    let root = BigInt::RING.root_floor(&n, i);
                    return Some(root);
                }
            }
            return Some(factoring_algorithms::sieve::quadratic_sieve(&n));
        }
    }
}

pub fn gen_primes() -> impl Iterator<Item = i64> {
    let mut found_primes = vec![2];
    std::iter::once(2).chain(
        condense((3..).step_by(2), move |x| {
            for p in &found_primes {
                if x % p == 0 {
                    return None;
                } else if p * p >= x {
                    found_primes.push(x);
                    return Some(x);
                }
            }
            unreachable!()
        })
    )
}

#[test]
fn test_calc_factor() {
    assert_eq!(Some(BigInt::from(3)), calc_factor(&BigInt::from(81)));
}