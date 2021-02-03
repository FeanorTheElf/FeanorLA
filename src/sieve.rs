use super::bigint::*;
use super::alg::*;

use std::collections::BTreeMap;

///
/// Generates all primes <= bound
/// 
fn gen_primes(bound: i64) -> Vec<i64> {
    assert!(bound > 13);
    let mut numbers = (0..bound).map(|_| true).collect::<Vec<_>>();
    let mut result = Vec::with_capacity(2 * ((bound as f64) / (bound as f64).ln()) as usize);
    for i in 2..(bound as usize) {
        if !numbers[i] {
            continue;
        }
        result.push(i as i64);
        for k in 2..(bound as usize / i) {
            numbers[k * i] = false;
        }
    }
    return result;
}

fn next_around_zero(x: i64) -> i64 {
    if x < 0 {
        return -x;
    } else {
        return -x - 1;
    }
}

fn around_zero_iter() -> impl Iterator<Item = i64> {
    let mut x = 0;
    return std::iter::from_fn(move || {
        x = next_around_zero(x);
        return Some(x);
    });
}

fn check_smooth(mut k: BigInt, factor_base: &Vec<i64>) -> Option<BTreeMap<i64, u32>> {
    let mut result = BTreeMap::new();
    let mut tmp = BigInt::zero();
    for p in factor_base {
        let mut dividing_power = 0;
        tmp.assign(&k);
        loop {
            let (quo, rem) = tmp.div_rem_small(*p as i64);
            tmp = quo;
            if rem != 0 {
                break;
            } else {
                dividing_power += 1;
                k.assign(&tmp);
            }
        }
        if dividing_power > 0 {
            result.insert(*p, dividing_power);
        }
    }
    if k == BigInt::one() {
        return Some(result);
    } else {
        return None;
    }
}

fn quadratic_sieve(n: BigInt) {
    let n_float = n.to_float_approx();
    let smoothness_bound_float = (0.5 * n_float.ln().sqrt() * n_float.ln().ln().sqrt()).exp();
    assert!(smoothness_bound_float < i64::MAX as f64);
    let smoothness_bound = smoothness_bound_float as i64;
    let factor_base = gen_primes(smoothness_bound);
    println!("factor_base: {:?}", factor_base);
    let m = BigInt::from_float_approx(n_float.sqrt());
    let mut relations: Vec<BTreeMap<i64, u32>> = Vec::with_capacity(factor_base.len() + 1);

    for d in around_zero_iter() {
        let mut k = m.clone();
        k += d;
        k = k.clone() * k;
        k -= &n;
        if let Some(rel) = check_smooth(k.clone(), &factor_base) {
            println!("found relation {} ~ {:?}", k, rel);
            relations.push(rel);
            if relations.len() > smoothness_bound as usize {
                break;
            }
        }
    }
}

#[test]
fn experiment() {
    let f5 = BigInt::power_of_two(32) + BigInt::one();
    println!("{}", f5);
    quadratic_sieve(f5);
    assert!(false);
}