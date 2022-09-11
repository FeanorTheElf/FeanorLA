use super::super::ring::*;
use super::super::la::mat::*;
use super::super::la::inversion::*;
use super::super::integer::*;
use super::super::primitive::*;
use super::super::fq::zn_small::*;
use super::super::eea::*;
use super::super::wrapper::*;
use super::gen_primes;

type Int = RingElWrapper<BigIntRing>;
const INT_RING: WrappingRing<BigIntRing> = BigInt::WRAPPED_RING; 

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

///
/// Checks whether the given number is smooth w.r.t. the given factor base,
/// and if this is the case, the corresponding factorization is returned.
/// This is done using basic trial division, so the complexity is
/// O(b) divisions of integers of similar size as k, where b is the length
/// of the factor base.
/// 
fn check_smooth(mut k: Int, factor_base: &Vec<i64>) -> Option<RelVec> {
    let mut result = Vector::zero(factor_base.len()).into_owned();
    assert!(factor_base[0] == -1);
    if k < 0 {
        *result.at_mut(0) = 1;
        k = -k;
    }
    let mut tmp = INT_RING.zero();
    for i in 1..factor_base.len() {
        let mut dividing_power = 0;
        tmp.val_mut().assign(k.val());
        loop {
            let quo = INT_RING.quotient(&tmp, &INT_RING.from_z(factor_base[i] as i64));
            if let Some(quo) = quo {
                dividing_power += 1;
                k.val_mut().assign(quo.val());
                tmp = quo;
            } else {
                break;
            }
        }
        *result.at_mut(i) = dividing_power;
    }
    if k.is_one() {
        return Some(result);
    } else {
        return None;
    }
}

///
/// Collects relations, i.e. squares modulo n that factor over the given
/// factor base. Given enough of those relations, one can then compute the
/// factorization of n.
/// Found relations are inserted into the given vector, and the function
/// terminates after count relations have been found.
/// 
fn collect_relations<J>(
    m: &Int, 
    n: &Int, 
    factor_base: &Vec<i64>, 
    relations: &mut Vec<(Int, RelVec)>, 
    count: usize, 
    delta_it: &mut J
)
    where J: Iterator<Item = i64>
{
    assert!(relations.len() < count);
    let mut k = m.clone();
    loop {
        // here consider (m + d)^2 - n
        // as d is small in absolute value, have that
        // (m + d)^2 - n = m^2 - n + 2md + d^2 is small in absolute value
        // as m is the rounded square root of n, so also m^2 - n is small
        // in absolute value.
        // This maximizes the chances to get a B-smooth square
        let d = delta_it.next().unwrap();
        k.val_mut().assign(m.val());
        k += d;
        let mut square = k.clone();
        square = square * k.clone();
        square = square - n.clone();

        // the probability that the found square is B-smooth is given by
        // the Canfield-Erdös-Pomerance theorem, concretely it is u^(-u(1 + o(1)))
        // for u = log(M)/log(B), where M is the size of the found square
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
/// yields a factor. If it does, the factor is returned.
/// 
fn check_congruent_square<V>(
    n: &Int, 
    factor_base: &Vec<i64>, 
    relations: &Vec<(Int, RelVec)>, 
    sol: Vector<V, F2>
) -> Result<Int, ()>
    where V: VectorView<ZnEl<2>>
{
    let mut x = INT_RING.one();
    let mut y_powers = Vector::zero(factor_base.len()).into_owned();
    for (i, rel) in relations.iter().enumerate() {
        if *sol.at(i) == F2::ONE {
            x = x * rel.0.clone();
            y_powers += rel.1.as_ref();
        }
    }

    let mut y = INT_RING.one();
    for i in 0..factor_base.len() {
        let power = *y_powers.at(i);
        debug_assert!(power % 2 == 0);
        y = y * INT_RING.from_z(factor_base[i]).pow(power as u32 / 2);
    }

    let factor = gcd(&INT_RING, n.clone(), x.clone() - y.clone());
    if factor != INT_RING.one() && factor != *n {
        return Ok(factor);
    } else {
        return Err(());
    }
}

///
/// Uses the quadratic sieve algorithm to find a nontrivial factor of n. 
/// Do not use this for perfect powers, as it will not terminate for prime powers.
/// The running time complexity is L_n(1/2, 1).
/// 
/// # Note
/// 
/// This is only a proof-of-concept implementation, and very slow in practice.
/// 
pub fn quadratic_sieve(n: &Int) -> Int {
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
    let m = n.ring().from_float_approx(n_float.sqrt()).unwrap();
    let mut relations: Vec<(Int, RelVec)> = Vec::with_capacity(
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
        let solutions = F2::RING.calc_matrix_kernel_space(matrix).unwrap();

        for i in 0..solutions.col_count() {

            if let Ok(factor) = check_congruent_square(
                n, &factor_base, &relations, solutions.col(i)
            ) {
                return factor;
            }
        }

        // increase count and hope we find suitable factors with more relations
        count += (count as f64).ln() as usize + 1;
    }
}

#[test]
fn test_quadratic_sieve() {
    let f5 = INT_RING.from(BigInt::RING.mul_pow_2(BigInt::RING.one(), 32)) + 1;
    let factor = quadratic_sieve(&f5);
    assert!(factor != f5);
    assert!(factor != 1);
    assert!(factor.divides(&f5));
}
