pub mod sieve;
pub mod fermat;

///
/// Generates all primes <= bound
/// 
pub fn gen_primes(bound: i64) -> Vec<i64> {
    assert!(bound >= 0);
    // the prime density formulas used later to estimate the vector 
    // size do not work for very small values 
    if bound <= 13 {
        return [2, 3, 5, 7, 11].iter()
            .map(|p| *p)
            .filter(|p| *p < bound)
            .collect::<Vec<_>>();
    }
    let mut numbers = (0..bound).map(|_| true).collect::<Vec<_>>();
    let mut result = Vec::with_capacity(
        2 * ((bound as f64) / (bound as f64).ln()) as usize
    );
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