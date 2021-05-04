use super::super::alg::*;
use super::super::la::mat::*;
use super::super::la::algorithms::MatrixKernelBase;
use super::bigint::*;
use super::zn::*;
use super::eea::*;

///
/// Generates all primes <= bound
/// 
pub fn gen_primes(bound: i64) -> Vec<i64> {
    assert!(bound >= 0);
    // the prime density formulas used later to estimate the vector 
    // size do not work for very small values 
    if bound <= 13 {
        return [2, 5, 7, 11].iter()
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

type RelVec = Vector<VectorOwned<i32>, i32>;

fn check_smooth(mut k: BigInt, factor_base: &Vec<i64>) -> Option<RelVec> {
    let mut result = Vector::zero(factor_base.len()).to_owned();
    assert!(factor_base[0] == -1);
    if k < 0 {
        *result.at_mut(0) = 1;
        k = -k;
    }
    let mut tmp = BigInt::zero();
    for i in 1..factor_base.len() {
        let mut dividing_power = 0;
        tmp.assign(&k);
        loop {
            let (quo, rem) = tmp.div_rem_small(factor_base[i] as i64);
            tmp = quo;
            if rem != 0 {
                break;
            } else {
                dividing_power += 1;
                k.assign(&tmp);
            }
        }
        *result.at_mut(i) = dividing_power;
    }
    if k == BigInt::one() {
        return Some(result);
    } else {
        return None;
    }
}

fn collect_relations<I>(
    m: &BigInt, 
    n: &BigInt, 
    factor_base: &Vec<i64>, 
    relations: &mut Vec<(BigInt, RelVec)>, 
    count: usize, 
    delta_it: &mut I
)
    where I: Iterator<Item = i64>
{
    let mut k = m.clone();
    loop {
        // here consider (m + d)^2 - n
        // as d is small in absolute value, have that
        // (m + d)^2 - n = m^2 - n + 2md + d^2 is small in absolute value
        // as m is the rounded square root of n, so also m^2 - n is small
        // in absolute value.
        // This maximizes the chances to get a B-smooth square
        let d = delta_it.next().unwrap();
        k.assign(&m);
        k += d;
        let mut square = k.clone();
        square *= &k;
        square -= n;
        if let Some(rel) = check_smooth(square, &factor_base) {
            println!("found relation: {}, relation count: {}", &k, relations.len());
            relations.push((k.clone(), rel));
            // we need enough relations to get a matrix with nontrivial kernel
            // modulo 2; this is the case for sure if the it has more columns
            // than rows
            if relations.len() == count {
                return;
            }
        }
    }
}

type F2 = ZnEl<2>;

///
/// Checks if the congruent square given by choosing exactly 
/// the relations given by sol is a real congruent square that 
/// yields a factor. If it does, the factor is returned
/// 
fn check_congruent_square<V>(
    n: &BigInt, factor_base: &Vec<i64>, 
    relations: &Vec<(BigInt, RelVec)>, 
    sol: Vector<V, F2>
) -> Result<BigInt, ()>
    where V: VectorView<ZnEl<2>>
{
    let mut x = BigInt::one();
    let mut y_powers = Vector::zero(factor_base.len()).to_owned();
    for (i, rel) in relations.iter().enumerate() {
        if *sol.at(i) == F2::ONE {
            x *= &rel.0;
            y_powers += rel.1.as_ref();
        }
    }

    let mut y = BigInt::one();
    for i in 0..factor_base.len() {
        let power = *y_powers.at(i);
        debug_assert!(power % 2 == 0);
        y *= BigInt::from(factor_base[i]).pow(power as u32 / 2);
    }

    let factor = gcd(&BigInt::RING, n.clone(), x.clone() - y.clone());
    if factor != BigInt::one() && factor != *n {
        return Ok(factor);
    } else {
        return Err(());
    }
}

///
/// Uses the quadratic sieve algorithm to find a nontrivial factor of n. 
/// Use only for composite numbers, as it will not terminate for primes.
/// 
pub fn quadratic_sieve(n: &BigInt) -> BigInt {
    assert!(*n >= 2);
    let n_float = n.to_float_approx();
    
    // we choose a factor base that consists of all primes <= B, where this
    // is B. The concrete value L_n(1/2, 1/2) is asymptotically optimal for 
    // performance.
    //
    // this factor base is used to describe multiplicative relations in Z/nZ
    // by the use of linear algebra. Concretely, if some integer factors over
    // the factor base modulo n, then we can describe it via the vector in which
    // each entry corresponds to the power of the prime factor dividing k.
    //
    // However, checking if k factors modulo n is very hard (in particular, as Z/nZ
    // is almost a field, all k do so). Therefore, we only consider integers that
    // factor over the factor base as integers, which are then called B-smooth.
    // In particular, we can learn about the concrete structure of Z/nZ by considering
    // different factorizations over the factor base modulo n (e.g. if k and k * n
    // are both B-smooth). Having enough of these relations allows us to solve
    // the factorization problem by linear algebra.
    //
    // Somewhat more performant is the quadratic sieve used here:
    // Look for integers k such that k is a square modulo n and factors over
    // the factor base. Using linear algebra over F2, we can then find a so
    // called congruent square (i.e. x^2 = y^2 but x != +/- y) which yields
    // a factor of n
    let smoothness_bound_float = (
        0.5 * n_float.ln().sqrt() * n_float.ln().ln().sqrt()
    ).exp();
    assert!(smoothness_bound_float < i64::MAX as f64);
    let smoothness_bound = smoothness_bound_float as i64;
    let factor_base = {
        let mut result = gen_primes(smoothness_bound);
        // also add -1 to the factor base for factoring negative numbers
        result.insert(0, -1);
        result
    };
    dbg!("factor_base: {:?}", &factor_base);

    // we search for B-smooth squares modulo n, and the probability that
    // a given square is B-smooth is better the smaller it is.
    // We get the smallest squares modulo n by considering integers k near
    // sqrt(n) and then looking at k^2 - n
    let m = BigInt::from_float_approx(n_float.sqrt());
    let mut relations: Vec<(BigInt, RelVec)> = Vec::with_capacity(
        factor_base.len() + 1
    );
    let mut delta_it = around_zero_iter();
    let mut count = factor_base.len() + 1;

    loop {

        collect_relations(&m, n, &factor_base, &mut relations, count, &mut delta_it);

        // if we find a vector in the kernel of this matrix, this vector
        // leads to a congruent square, as then the relation matrix (not taken
        // modulo 2) times the vector has even coefficients for each prime,
        // so is a square. The other square we get as product of the 
        // corresponding relation squares modulo n, since we chose them
        // exactly so that they are squares
        let matrix = Matrix::from_fn(factor_base.len(), relations.len(), |r, c| 
            F2::project(*relations[c].1.at(r) as i64)
        );
        let solutions = StaticRing::<F2>::RING.calc_matrix_kernel_space(matrix).unwrap();

        for i in 0..solutions.col_count() {

            if let Ok(factor) = check_congruent_square(
                n, &factor_base, &relations, solutions.col(i)
            ) {
                return factor;
            }
        }

        // increase count and hope we find suitable factors with more relations
        count += (count as f64).ln() as usize;
        count += 1;
    }
}

#[test]
fn test_quadratic_sieve() {
    let f5 = BigInt::power_of_two(32) + 1;
    let factor = quadratic_sieve(&f5);
    assert!(factor != f5);
    assert!(factor != 1);
    assert_eq!(f5 % factor, 0);
}
