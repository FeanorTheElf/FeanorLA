use super::bigint::*;

///
/// Generates at least n primes {2, 3, 5, ...}
/// 
fn gen_primes(n: usize) -> Vec<u64> {
    assert!(n > 13);
    // by the prime number theorem, this is the approximate bound up to which we have to try numbers
    let prime_number_bound = |x: f64| x * x.ln() + x * x.ln().ln();
    let bound = prime_number_bound(n as f64).ceil() as usize;
    let mut numbers = (0..bound).map(|_| true).collect::<Vec<_>>();
    let mut result = Vec::with_capacity(2 * n);
    for i in 2..bound {
        if !numbers[i] {
            continue;
        }
        result.push(i as u64);
        for k in 2..(bound / i) {
            numbers[k * i] = false;
        }
    }
    return result;
}

fn quadratic_sieve(n: BigInt) {
    let n_float = n.to_float_approx();
    let factor_base_size = (0.5 * n_float.ln().sqrt() * n_float.ln().ln().sqrt()).exp() as usize;
    let factor_base = gen_primes(factor_base_size);
    let m = BigInt::from_float_approx(n_float.sqrt());
}