use super::super::la::vec::*;

use std::cmp::Ordering;

type BlockInt = u64;
type DoubleBlockInt = u128;
const BLOCK_BITS: u32 = u64::BITS;

fn expand(x: &mut Vector<Vec<BlockInt>, BlockInt>, len: usize) {
	if len > x.len() {
		let mut base_vec = std::mem::replace(x, Vector::new(Vec::new())).raw_data();
		base_vec.resize(len, 0);
		*x = Vector::new(base_vec);
	}
}

#[cfg(test)]
fn truncate_zeros(x: Vector<Vec<BlockInt>, BlockInt>) -> Vector<Vec<BlockInt>, BlockInt> {
	let mut data = x.raw_data();
	data.truncate(data.len() - data.iter().rev().take_while(|a| **a == 0).count());
	return Vector::new(data);
}

fn assign<V>(x: &mut Vector<Vec<BlockInt>, BlockInt>, rhs: Vector<V, BlockInt>) 
	where V: VectorView<BlockInt>
{
	x.assign(Vector::zero(x.len()));
	expand(x, rhs.len());
	x.subvector_mut(..rhs.len()).assign(rhs);
}

pub fn bigint_add<W>(lhs: &mut Vector<Vec<BlockInt>, BlockInt>, rhs: Vector<W, BlockInt>, block_offset: usize)
	where W: VectorView<BlockInt>
{
	let mut buffer: bool = false;
	let mut i = 0;
	while i < rhs.len() || buffer {
		let rhs_val = *rhs.get(i).unwrap_or(&0);
		let j = i + block_offset;
		expand(lhs, j + 1);
		let (sum, overflow) = lhs[j].overflowing_add(rhs_val);
		if buffer {
			let (carry_sum, carry_overflow) = sum.overflowing_add(1);
			*lhs.at_mut(j) = carry_sum;
			buffer = overflow || carry_overflow;
		} else {
			*lhs.at_mut(j) = sum;
			buffer = overflow;
		}
		i += 1;
	}
}

pub fn highest_set_block<V>(x: Vector<V, BlockInt>) -> Option<usize> 
	where V: VectorView<BlockInt>
{
	for i in (0..x.len()).rev() {
		if x[i] != 0 {
			return Some(i);
		}
	}
	return None;
}

pub fn bigint_cmp<V, W>(lhs: Vector<V, BlockInt>, rhs: Vector<W, BlockInt>) -> Ordering 
	where V: VectorView<BlockInt>, W: VectorView<BlockInt>
{
	match (highest_set_block(lhs.as_ref()), highest_set_block(rhs.as_ref())) {
		(None, None) => return Ordering::Equal,
		(Some(_), None) => return Ordering::Greater,
		(None, Some(_)) => return Ordering::Less,
		(Some(x), Some(y)) => match x.cmp(&y) {
			Ordering::Less => return Ordering::Less,
			Ordering::Greater => return Ordering::Greater,
			Ordering::Equal => {
				for i in (0..=x).rev() {
					match lhs[i].cmp(&rhs[i]) {
						Ordering::Less => return Ordering::Less,
						Ordering::Greater => return Ordering::Greater,
						_ => {}
					}
				}
				return Ordering::Equal;
			}
		}
	};
}

pub fn bigint_cmp_small<V>(lhs: Vector<V, BlockInt>, rhs: DoubleBlockInt) -> Ordering 
	where V: VectorView<BlockInt>
{
	match highest_set_block(lhs.as_ref()) {
	   None => 0.cmp(&rhs),
	   Some(0) => (lhs[0] as DoubleBlockInt).cmp(&rhs),
	   Some(1) => (((lhs[1] as DoubleBlockInt) << BLOCK_BITS) | (lhs[0] as DoubleBlockInt)).cmp(&rhs),
	   Some(_) => Ordering::Greater,
	}
}

///
/// Calculate self -= rhs * (1 << BLOCK_BITS)^block_offset
/// 
/// This will panic if the subtraction would result in a negative number
/// 
pub fn bigint_sub<V, W>(lhs: &mut Vector<V, BlockInt>, rhs: Vector<W, BlockInt>, block_offset: usize) 
	where V: VectorViewMut<BlockInt>, W: VectorView<BlockInt>
{
	assert!(bigint_cmp(lhs.as_ref(), rhs.as_ref()) != Ordering::Less);

	if let Some(rhs_high) = highest_set_block(rhs.as_ref()) {
		let mut buffer: bool = false;
		let mut i = 0;
		while i <= rhs_high || buffer {
			let rhs_val = *rhs.get(i).unwrap_or(&0);
			let j = i + block_offset;
			debug_assert!(j < lhs.len());
			let (difference, overflow) = lhs[j].overflowing_sub(rhs_val);
			if buffer {
				let (carry_difference, carry_overflow) = difference.overflowing_sub(1);
				*lhs.at_mut(j) = carry_difference;
				buffer = overflow || carry_overflow;
			} else {
				*lhs.at_mut(j) = difference;
				buffer = overflow;
			}
			i += 1;
		}
	}
}

pub fn bigint_mul<V, W>(lhs: Vector<V, BlockInt>, rhs: Vector<W, BlockInt>) -> Vector<Vec<BlockInt>, BlockInt> 
	where V: VectorView<BlockInt>, W: VectorView<BlockInt>
{
	let mut result = Vector::new(Vec::with_capacity(
		highest_set_block(lhs.as_ref()).unwrap_or(0) + 
		highest_set_block(rhs.as_ref()).unwrap_or(0) + 2
	));
	if let Some(d) = highest_set_block(rhs.as_ref()) {
		let mut val = Vector::new(Vec::new());
		for i in 0..=d {
			assign(&mut val, lhs.as_ref());
			bigint_mul_small(&mut val, rhs[i]);
			bigint_add(&mut result, val.as_ref(), i);
		}
	}
	return result;
}

///
/// Complexity O(log(n))
/// 
pub fn bigint_mul_small(lhs: &mut Vector<Vec<BlockInt>, BlockInt>, factor: BlockInt) {
	if let Some(d) = highest_set_block(lhs.as_ref()) {
		let mut buffer: u64 = 0;
		for i in 0..=d {
			let prod = lhs[i] as u128 * factor as u128 + buffer as u128;
			*lhs.at_mut(i) = (prod & ((1u128 << BLOCK_BITS) - 1)) as u64;
			buffer = (prod >> BLOCK_BITS) as u64;
		}
		expand(lhs, d + 2);
		*lhs.at_mut(d + 1) = buffer;
	}
}

pub fn bigint_add_small(lhs: &mut Vector<Vec<BlockInt>, BlockInt>, rhs: BlockInt) {
	if lhs.len() > 0 {
		let (sum, mut buffer) = lhs[0].overflowing_add(rhs);
		*lhs.at_mut(0) = sum;
		let mut i = 1;
		while buffer {
			expand(lhs, i + 1);
			let (sum, overflow) = lhs[i].overflowing_add(1);
			buffer = overflow;
			*lhs.at_mut(i) = sum;
			i += 1;
		}
	} else {
		expand(lhs, 1);
		*lhs.at_mut(0) = rhs;
	}
}

///
/// Same as division_step, but for self_high == rhs_high == d
/// 
fn division_step_last<V, W>(lhs: &mut Vector<V, BlockInt>, rhs: Vector<W, BlockInt>, d: usize, tmp: &mut Vector<Vec<BlockInt>, BlockInt>) -> u64 
	where V: VectorViewMut<BlockInt>, W: VectorView<BlockInt>
{
	assert!(lhs[d] != 0);
	assert!(lhs[d] != 0);

	let self_high_blocks: u128 = ((lhs[d] as u128) << BLOCK_BITS) | (lhs[d - 1] as u128);
	let rhs_high_blocks: u128 = ((rhs[d] as u128) << BLOCK_BITS) | (rhs[d - 1] as u128);

	if rhs_high_blocks == u128::MAX {
		if bigint_cmp(lhs.as_ref(), rhs.as_ref()) != Ordering::Less {
			bigint_sub(lhs, rhs, 0);
			return 1;
		} else {
			return 0;
		}
	} else {
		let mut quotient = (self_high_blocks / (rhs_high_blocks + 1)) as u64;
		assign(tmp, rhs.as_ref());
		bigint_mul_small(tmp, quotient);
		bigint_sub(lhs, tmp.as_ref(), 0);

		if bigint_cmp(lhs.as_ref(), rhs.as_ref()) != Ordering::Less {
			bigint_sub(lhs, rhs.as_ref(), 0);
			quotient += 1;
		}
		if bigint_cmp(lhs.as_ref(), rhs.as_ref()) != Ordering::Less {
			bigint_sub(lhs, rhs.as_ref(), 0);
			quotient += 1;
		}
		
		debug_assert!(bigint_cmp(lhs.as_ref(), rhs) == Ordering::Less);
		return quotient;
	}
}

///
/// Finds some integer d such that subtracting d * rhs from self clears the top
/// block of self. self will be assigned the value after the subtraction and d
/// will be returned as d = (u * 2 ^ block_bits + l) * 2 ^ (k * block_bits) 
/// where the return value is (u, l, k)
/// 
/// Complexity O(log(n))
/// 
fn division_step<V, W>(
	lhs: &mut Vector<V, BlockInt>, 
	rhs: Vector<W, BlockInt>, 
	lhs_high: usize, 
	rhs_high: usize, 
	tmp: &mut Vector<Vec<BlockInt>, BlockInt>
) -> (u64, u64, usize) 
	where V: VectorViewMut<BlockInt>, W: VectorView<BlockInt>
{
	assert!(lhs_high > rhs_high);
	assert!(lhs[lhs_high] != 0);
	assert!(rhs[rhs_high] != 0);

	// the basic idea is as follows:
	// we know that for a and b, have a - a//(b+1) * b <= b + a/(b+1)
	// Hence, we perform two steps:
	//  - by choosing a and b as the top two blocks of lhs resp. lhs, achieve that a - a//(b+1) * b <= b + 2^k
	//    (where k = BLOCK_BITS); hence, lhs - a//(b+1) * rhs <= rhs + 2^k, and so possibly subtracting rhs
	//    achieves new_lhs <= rhs
	//  - by choosing a as the top two blocks and b as only the top block of lhs resp. lhs (now b < 2^k), achieve that
	//    lhs - a//(b+1) * rhs < 2^k + 2^k = 2 * 2^k, and so after possibly subtracting rhs we find
	//    that the top block of lhs is cleared

	let mut result_upper = 0;
	let mut result_lower = 0;

	// first step
	{
		let lhs_high_blocks = ((lhs[lhs_high] as DoubleBlockInt) << BLOCK_BITS) | (lhs[lhs_high - 1] as DoubleBlockInt);
		let rhs_high_blocks = ((rhs[rhs_high] as DoubleBlockInt) << BLOCK_BITS) | (rhs[rhs_high - 1] as DoubleBlockInt);

		if rhs_high_blocks != DoubleBlockInt::MAX && lhs_high_blocks >= (rhs_high_blocks + 1) {
			let mut quotient = (lhs_high_blocks / (rhs_high_blocks + 1)) as u64;
			debug_assert!(quotient != 0);
			assign(tmp, rhs.as_ref());
			bigint_mul_small(tmp, quotient);
			bigint_sub(lhs, tmp.as_ref(), lhs_high - rhs_high);

			let lhs_high_blocks = ((lhs[lhs_high] as DoubleBlockInt) << BLOCK_BITS) | (lhs[lhs_high - 1] as DoubleBlockInt);

			if lhs_high_blocks > rhs_high_blocks {
				bigint_sub(lhs, rhs.as_ref(), lhs_high - rhs_high);
				quotient += 1;
			}
			result_upper = quotient;
		}

		// this is what we wanted to achieve in the first step
		let lhs_high_blocks = ((lhs[lhs_high] as DoubleBlockInt) << BLOCK_BITS) | (lhs[lhs_high - 1] as DoubleBlockInt);
		debug_assert!(lhs_high_blocks <= rhs_high_blocks);
	}

	// second step
	{
		let lhs_high_blocks = ((lhs[lhs_high] as DoubleBlockInt) << BLOCK_BITS) | (lhs[lhs_high - 1] as DoubleBlockInt);
		let rhs_high_block = rhs[rhs_high] as DoubleBlockInt;

		if lhs[lhs_high] != 0 {
			let mut quotient = (lhs_high_blocks / (rhs_high_block + 1)) as BlockInt;
			assign(tmp, rhs.as_ref());
			bigint_mul_small(tmp, quotient);
			bigint_sub(lhs, tmp.as_ref(), lhs_high - rhs_high - 1);

			if lhs[lhs_high] != 0 {
				bigint_sub(lhs, rhs.as_ref(), lhs_high - rhs_high - 1);
				quotient += 1;
			}
			result_lower = quotient;
		}

		debug_assert!(lhs[lhs_high] == 0);
	}
	return (result_upper, result_lower, lhs_high - rhs_high - 1);
}

///
// Calculates abs(self) = abs(self) % abs(rhs) and returns the quotient
/// of the division abs(self) / abs(rhs). The sign bit of self is ignored
/// and left unchanged.
/// 
/// Complexity O(log(n)^2)
/// 
pub fn bigint_div<V, W>(lhs: &mut Vector<V, BlockInt>, rhs: Vector<W, BlockInt>) -> Vector<Vec<BlockInt>, BlockInt> 
	where V: VectorViewMut<BlockInt>, W: VectorView<BlockInt>
{
	assert!(highest_set_block(rhs.as_ref()).is_some());
	
	if let Some(mut d) = highest_set_block(lhs.as_ref()) {
		let mut tmp = Vector::new(Vec::new());
		let k = highest_set_block(rhs.as_ref()).expect("Division by zero");
		if d < k {
			return Vector::new(Vec::new());
		} else if k == 0 {
			let rem = bigint_div_small(lhs, rhs[0]);
			assign(&mut tmp, lhs.as_ref());
			lhs.assign(Vector::zero(lhs.len()));
			*lhs.at_mut(0) = rem;
			return tmp;
		} else {
			let mut result = Vector::new(Vec::new());
			expand(&mut result, d - k + 1);
			while d > k {
				if lhs[d] != 0 {
					let (quo_upper, quo_lower, quo_power) = division_step(lhs, rhs.as_ref(), d, k, &mut tmp);
					*result.at_mut(quo_power) = quo_lower;
					bigint_add(&mut result, Vector::from_array([quo_upper]), quo_power + 1);
					debug_assert!(lhs[d] == 0);
				}
				d -= 1;
			}
			let quo = division_step_last(lhs, rhs, d, &mut tmp);
			bigint_add(&mut result, Vector::from_array([quo]), 0);
			return result;
		}
	} else {
		return Vector::new(Vec::new());
	}
}

///
/// Calculates self /= divisor and returns the remainder of the division.
/// 
pub fn bigint_div_small<V>(lhs: &mut Vector<V, BlockInt>, rhs: BlockInt) -> BlockInt 
	where V: VectorViewMut<BlockInt>
{
	assert!(rhs != 0);
	let highest_block_opt = highest_set_block(lhs.as_ref());
	if highest_block_opt == Some(0) {
		let (quo, rem) = (lhs[0] / rhs, lhs[0] % rhs);
		*lhs.at_mut(0) = quo;
		return rem;
	} else if let Some(highest_block) = highest_block_opt {
		let (quo, rem) = (lhs[highest_block] / rhs, lhs[highest_block] % rhs);
		let mut buffer = rem as DoubleBlockInt;
		*lhs.at_mut(highest_block) = quo;
		for i in (0..highest_block).rev() {
			buffer = (buffer << BLOCK_BITS) | (lhs[i] as DoubleBlockInt);
			let (quo, rem) = (buffer / rhs as DoubleBlockInt, buffer % rhs as DoubleBlockInt);
			debug_assert!(quo <= BlockInt::MAX as DoubleBlockInt);
			*lhs.at_mut(i) = quo as BlockInt;
			buffer = rem;
		}
		return buffer as BlockInt;
	} else {
		return 0;
	}
}

#[cfg(test)]
use super::bigint::BigInt;

#[test]
fn test_sub() {
    let mut x = "923645871236598172365987287530543".parse::<BigInt>().unwrap().base_u64_repr();
    let y = "58430657823473456743684735863478".parse::<BigInt>().unwrap().base_u64_repr();
    let z = "865215213413124715622302551667065".parse::<BigInt>().unwrap().base_u64_repr();
    bigint_sub(&mut x, y, 0);
    assert_eq!(truncate_zeros(z), truncate_zeros(x));
}

#[test]
fn test_sub_with_carry() {
    let mut x = BigInt::from_str_radix("1000000000000000000", 16).unwrap().base_u64_repr();
    let y = BigInt::from_str_radix("FFFFFFFFFFFFFFFF00", 16).unwrap().base_u64_repr();
    bigint_sub(&mut x, y, 0);
    assert_eq!(Vector::from_array([256]), truncate_zeros(x));
}

#[test]
fn test_add() {
    let mut x = "923645871236598172365987287530543".parse::<BigInt>().unwrap().base_u64_repr();
    let y = "58430657823473456743684735863478".parse::<BigInt>().unwrap().base_u64_repr();
    let z = "982076529060071629109672023394021".parse::<BigInt>().unwrap().base_u64_repr();
    bigint_add(&mut x, y, 0);
    assert_eq!(truncate_zeros(z), truncate_zeros(x));
}

#[test]
fn test_add_with_carry() {
    let mut x = BigInt::from_str_radix("1BC00000000000000BC", 16).unwrap().base_u64_repr();
    let y =  BigInt::from_str_radix("FFFFFFFFFFFFFFFF0000000000000000BC", 16).unwrap().base_u64_repr();
    let z = BigInt::from_str_radix("10000000000000000BC0000000000000178", 16).unwrap().base_u64_repr();
    bigint_add(&mut x, y, 0);
    assert_eq!(truncate_zeros(z), truncate_zeros(x));
}

#[test]
fn test_mul() {
    let x = BigInt::from_str_radix("57873674586797895671345345", 10).unwrap().base_u64_repr();
    let y = BigInt::from_str_radix("21308561789045691782534873921650342768903561413264128756389247568729346542359871235465", 10).unwrap().base_u64_repr();
    let z = BigInt::from_str_radix("1233204770891906354921751949503652431220138020953161094405729272872607166072371117664593787957056214903826660425", 10).unwrap().base_u64_repr();
    assert_eq!(truncate_zeros(z), truncate_zeros(bigint_mul(x, y)));
}

#[test]
fn test_div_no_remainder() {
    let mut x = BigInt::from_str_radix("578435387FF0582367863200000000000000000000", 16).unwrap().base_u64_repr();
    let y = BigInt::from_str_radix("200000000000000000000", 16).unwrap().base_u64_repr();
    let z = BigInt::from_str_radix("2BC21A9C3FF82C11B3C319", 16).unwrap().base_u64_repr();
    let quotient = bigint_div(&mut x, y);
    assert_eq!(Vector::<_, u64>::zero(x.len()), x);
    assert_eq!(truncate_zeros(z), truncate_zeros(quotient));
}

#[test]
fn test_div_with_remainder() {
    let mut x = BigInt::from_str_radix("578435387FF0582367863200000000007651437856", 16).unwrap().base_u64_repr();
    let y = BigInt::from_str_radix("200000000000000000000", 16).unwrap().base_u64_repr();
    let z = BigInt::from_str_radix("2BC21A9C3FF82C11B3C319", 16).unwrap().base_u64_repr();
    let r = BigInt::from_str_radix("7651437856", 16).unwrap().base_u64_repr();
    let quotient = bigint_div(&mut x, y);
    assert_eq!(truncate_zeros(r), truncate_zeros(x));
    assert_eq!(truncate_zeros(z), truncate_zeros(quotient));
}

#[test]
fn test_div_big() {
    let mut x = BigInt::from_str_radix("581239456149785691238569872349872348569871269871234657986123987237865847935698734296434575367565723846982523852347", 10).unwrap().base_u64_repr();
    let y = BigInt::from_str_radix("903852718907268716125180964783634518356783568793426834569872365791233387356325", 10).unwrap().base_u64_repr();
    let q = BigInt::from_str_radix("643068769934649368349591185247155725", 10).unwrap().base_u64_repr();
    let r = BigInt::from_str_radix("265234469040774335115597728873888165088018116561138613092906563355599185141722", 10).unwrap().base_u64_repr();
    let quotient = bigint_div(&mut x, y);
    assert_eq!(truncate_zeros(r), truncate_zeros(x));
    assert_eq!(truncate_zeros(q), truncate_zeros(quotient));
}

#[test]
fn test_div_last_block_overflow() {
    let mut x = BigInt::from_str_radix("3227812347608635737069898965003764842912132241036529391038324195675809527521051493287056691600172289294878964965934366720", 10).unwrap().base_u64_repr();
    let y = BigInt::from_str_radix("302231454903657293676544", 10).unwrap().base_u64_repr();
    let q = BigInt::from_str_radix("10679935179604550411975108530847760573013522611783263849735208039111098628903202750114810434682880", 10).unwrap().base_u64_repr();
    let quotient = bigint_div(&mut x, y);
    assert_eq!(truncate_zeros(q), truncate_zeros(quotient));
    assert_eq!(Vector::<_, u64>::zero(x.len()), x);
}

#[test]
fn test_div_small() {
    let mut x = BigInt::from_str_radix("891023591340178345678931246518793456983745682137459364598623489512389745698237456890239238476873429872346579", 10).unwrap().base_u64_repr();
    let q = BigInt::from_str_radix("255380794307875708133829534685810678413226048190730686328066348384175908769916152734376393945793473738133", 10).unwrap().base_u64_repr();
    bigint_div_small(&mut x, 3489);
    assert_eq!(truncate_zeros(q), truncate_zeros(x));
}