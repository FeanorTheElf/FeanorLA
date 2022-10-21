use super::integer::*;
use super::embedding::*;
use super::primitive::*;
use super::ring_property::*;
use super::square_multiply::abs_square_and_multiply;
use super::wrapper::*;

use vector_map::VecMap;

pub struct RingElDisplay<'a, R: ?Sized> 
    where R: RingBase
{
    ring: &'a R,
    el: &'a R::El
}

impl<'a, R> std::fmt::Display for RingElDisplay<'a, R>
    where R: RingBase
{
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        self.ring.format(self.el, f, false)
    }
}

///
/// Trait to represent a unital, commutative ring as a collection of operations 
/// on the elements. More abstract functionality (global properties, like ideals) 
/// is not provided, this is mainly an interface that can be used for algorithms 
/// that deal with ring elements.
/// 
/// # Further properties
/// 
/// Further properties a ring can have are represented by subtraits, e.g.
/// `EuclideanInfoRing` or `DivisibilityInfoRing`. Note that implementing such a
/// trait means that objects of that type can have the property, not that they
/// necessarily have to. Each of those subtraits provides a function to check
/// whether a ring object has the property at runtime. This is necessary, as
/// the check might require information only available at runtime.
/// 
/// # Value vs Reference design rationale
/// 
/// Types that are cheap to copy can be used with the dd/sub()-functions, 
/// for types that are expensive to copy, one can use the
/// add_ref/sub_ref()-functions. However, as each add/sub() / add_ref/sub_ref() 
/// function returns the result by value, for big types, is is usually 
/// most efficient to clone at least one parameter and potentially reuse 
/// the memory.
///
/// If an operation can improve efficience by consuming both parameters, one should
/// explicitly implement the default-implemented add/sub()-functions.
///
/// For multiplication, the situation is usually different: Mostly, multiplication
/// is given by some kind of operation on a cartesian product of the
/// components of lhs and rhs, with a following reduce. In this case, one usually
/// cannot profit from getting the parameters by value. If this is not true, then
/// one should implement the default-implemented mul()-function
/// 
/// # A note on Clone
/// 
/// Sometimes, the question arises if rings should generally be required to implement
/// Clone. In the end, the decision is to add it, since some very common use cases
/// require Clone. In particular
///  - In order to get `StandardEmbedding` objects (and other kinds of embeddings) that
///    live for static lifetime (required e.g. for `RingExtension`), we must be able to
///    clone source & destination ring
///  - Many functionality of `WrappingRing` requires cloning, even to make it a ring (as
///    elements of a ring are required to be clonable).
/// 
/// I am completely aware that cloning in these situations might lead to a non-obvious
/// problems, e.g. if the ring contains a thread or memory pool. In these cases, we 
/// recommend to implement `Ring` only for objects that contain references or `Rc`'s to
/// the corresponding resource.
/// 
pub trait RingBase : std::fmt::Debug + std::clone::Clone {
    
    type El: Sized + Clone + std::fmt::Debug;

    fn add_ref(&self, lhs: Self::El, rhs: &Self::El) -> Self::El;
    fn mul_ref(&self, lhs: &Self::El, rhs: &Self::El) -> Self::El;
    fn neg(&self, val: Self::El) -> Self::El;
    fn zero(&self) -> Self::El;
    fn one(&self) -> Self::El;
    fn is_eq(&self, lhs: &Self::El, rhs: &Self::El) -> bool;

    ///
    /// Returns any element of the ring.
    /// 
    /// This ring element might be used as drop-in if one needs some 
    /// unspecified element. This is used only in very exceptional cases.
    /// Best do not touch this - chances are, you will never encounter it.
    /// Because of this, it is sensible to choose an element that is very
    /// cheap to construct, probably that is zero.
    /// 
    /// # Example
    /// 
    /// A default implementation of add_assign might move out the value from the
    /// mutable reference, then call add_ref and fill the result in again. However,
    /// if the underlying add-call panics, there is no value to fill in again.
    /// Usually, this is irrelevant, as the variable on which add_assign is called
    /// goes out of scope by the panic, if however the panic is caught (might require
    /// unsafe code due to UnwindSafe ???), this is not the case and the value can
    /// be accessed later. In this case, it will be filled with invalid().
    /// 
    fn unspecified_element(&self) -> Self::El {
        self.zero()
    }

    fn sub_ref_fst(&self, lhs: &Self::El, rhs: Self::El) -> Self::El {
        self.add_ref(self.neg(rhs), lhs)
    }

    fn sub_ref_snd(&self, lhs: Self::El, rhs: &Self::El) -> Self::El {
        self.neg(self.add_ref(self.neg(lhs), rhs))
    }

    fn add(&self, lhs: Self::El, rhs: Self::El) -> Self::El { self.add_ref(lhs, &rhs) }

    fn add_assign(&self, lhs: &mut Self::El, rhs: Self::El) { 
        let value = std::mem::replace(lhs, self.unspecified_element());
        *lhs = self.add(value, rhs);
    }

    fn add_assign_ref(&self, lhs: &mut Self::El, rhs: &Self::El) { 
        let value = std::mem::replace(lhs, self.unspecified_element());
        *lhs = self.add_ref(value, rhs);
    }

    fn add_assign_int(&self, lhs: &mut Self::El, rhs: i64) {
        self.add_assign(lhs, self.from_z(rhs));
    }

    fn sum<I>(&self, data: I) -> Self::El
        where I: Iterator<Item = Self::El>
    {
        data.fold(self.zero(), |a, b| self.add(a, b))
    }

    fn mul(&self, lhs: Self::El, rhs: Self::El) -> Self::El { self.mul_ref(&lhs, &rhs) }

    fn mul_assign(&self, lhs: &mut Self::El, rhs: &Self::El) { 
        *lhs = self.mul_ref(lhs, rhs);
    }

    fn mul_assign_int(&self, lhs: &mut Self::El, rhs: i64) {
        self.mul_assign(lhs, &self.from_z(rhs));
    }

    fn product<I>(&self, data: I) -> Self::El 
        where I: Iterator<Item = Self::El>
    {
        data.fold(self.one(), |a, b| self.mul(a, b))
    }

    fn sub(&self, lhs: Self::El, rhs: Self::El) -> Self::El {
        self.add(lhs, self.neg(rhs))
    }

    fn pow(&self, basis: &Self::El, exp: u32) -> Self::El 
        where Self::El: Clone
    {
        if exp == 0 {
            return self.one();
        // cases 1 and 2 are not necessary for correctness, but will speed up the whole thing
        } else if exp == 1 {
            return basis.clone();
        } else if exp == 2 {
            return self.mul_ref(basis, basis);
        }
        return abs_square_and_multiply(basis, &(exp as i64), i64::RING, |x, y| self.mul(x, y), |x, y| self.mul_ref(x, y), self.one());
    }

    fn pow_big(&self, basis: &Self::El, exp: &StdInt) -> Self::El 
        where Self::El: Clone
    {
        assert!(*exp >= 0);
        return abs_square_and_multiply(basis, exp, StdInt::RING, |x, y| self.mul(x, y), |x, y| self.mul_ref(x, y), self.one());
    }

    fn is_zero(&self, val: &Self::El) -> bool { self.is_eq(val, &self.zero()) }
    fn is_one(&self, val: &Self::El) -> bool { self.is_eq(val, &self.one()) }
    fn is_neg_one(&self, val: &Self::El) -> bool { self.is_eq(val, &self.neg(self.one())) }

    fn characteristic(&self) -> StdInt;
    ///
    /// Returns whether the ring is integral, so if there for all nonzero a, b it holds
    /// that ab != 0.
    /// 
    fn is_integral(&self) -> RingPropValue;
    ///
    /// Returns whether the ring is a field, so whether each nonzero element has a unique 
    /// inverse.
    /// 
    fn is_field(&self) -> RingPropValue;
    fn is_noetherian(&self) -> bool;

    ///
    /// May panic if the ring is not a field (meaning `!is_field().can_use()`). 
    /// If it does not panic, the result must be valid. For a non-field ring, it therefore 
    /// must panic if rhs does not divide lhs, and if it does, it may either compute the 
    /// correct quotient but may also panic nevertheless.
    /// 
    fn div(&self, lhs: Self::El, rhs: &Self::El) -> Self::El;

    fn from_z_big(&self, x: &StdInt) -> Self::El {
        let result = abs_square_and_multiply(&self.one(), x, StdInt::RING, |x, y| self.add(x, y), |x, y| self.add_ref(x.clone(), y), self.zero());
        if *x < 0 {
            return self.neg(result);
        } else {
            return result;
        }
    }

    fn from_z(&self, x: i64) -> Self::El {
        let result = abs_square_and_multiply(&self.one(), &x, i64::RING, |x, y| self.add(x, y), |x, y| self.add_ref(x.clone(), y), self.zero());
        if x < 0 {
            return self.neg(result);
        } else {
            return result;
        }
    }

    fn from_z_gen<I>(&self, x: El<I>, ring: &I) -> Self::El 
        where I: IntegerRing
    {
        let result = abs_square_and_multiply(&self.one(), &x, ring, |x, y| self.add(x, y), |x, y| self.add_ref(x.clone(), y), self.zero());
        if ring.is_neg(&x) {
            return self.neg(result);
        } else {
            return result;
        }
    }

    ///
    /// Writes a textual representation of the element to the formatter.
    /// If in_prod is set, then the representation will be chosen such that
    /// any internal operations bind stronger (or equally strong) than 
    /// potential surrounding multiplications (usually brackets will be placed).
    /// If in_prod is not set, then the representation will be chosen such
    /// that any internal operations bind stronger (or equally strong) than
    /// potential surrounding additions.
    /// 
    /// If an even stronger binding is required (i.e. for surrounding powering),
    /// use format_in_brackets which will unconditionally place brackets.
    /// 
    fn format(&self, el: &Self::El, f: &mut std::fmt::Formatter, _in_prod: bool) -> std::fmt::Result {
        write!(f, "{:?}", el)
    }

    fn format_in_brackets(&self, el: &Self::El, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        write!(f, "(")?;
        self.format(el, f, false)?;
        write!(f, ")")?;
        return Ok(());
    }
    
    fn display<'a>(&'a self, el: &'a Self::El) -> RingElDisplay<'a, Self> {
        RingElDisplay {
            ring: self,
            el: el
        }
    }
}

pub trait Ring: RingBase + CanonicalIsomorphismInfo<Self> {}

impl<R: ?Sized> Ring for R
    where R: RingBase + CanonicalIsomorphismInfo<R>
{}

#[allow(type_alias_bounds)]
pub type El<R: Ring> = <R as RingBase>::El;

///
/// Trait for rings that might have a euclidean division.
/// 
pub trait EuclideanInfoRing: DivisibilityInfoRing {
    
    ///
    /// Returns whether the ring is euclidean, so whether the euclidean division and remainder
    /// functions are implemented and behave correctly.
    /// 
    fn is_euclidean(&self) -> RingPropValue;

    ///
    /// May panic if the ring is not euclidean (meaning `!is_euclidean().can_use()`). 
    /// The first result is the quotient and the second result the remainder. 
    /// The equality
    ///  `lhs == quo * rhs + rem`
    /// must always hold, and furthermore, we must have
    ///  `euclidean_deg(rem) < euclidean_deg(rhs)`
    /// 
    fn euclidean_div_rem(&self, lhs: Self::El, rhs: &Self::El) -> (Self::El, Self::El);

    ///
    /// May panic if the ring is not euclidean (meaning `!is_euclidean().can_use()`).
    /// 
    fn euclidean_rem(&self, lhs: Self::El, rhs: &Self::El) -> Self::El { 
        self.euclidean_div_rem(lhs, rhs).1 
    }

    ///
    /// May panic if the ring is not euclidean (meaning `!is_euclidean().can_use()`).
    /// 
    fn euclidean_div(&self, lhs: Self::El, rhs: &Self::El) -> Self::El {
        self.euclidean_div_rem(lhs, rhs).0
    }

    fn euclidean_deg(&self, el: Self::El) -> StdInt;
}

///
/// Trait for rings that might be a ring extension of another ring, i.e.
/// have a canonical, injective ring homomorphism R' -> R.
/// 
pub trait RingExtension: Ring {

    type BaseRing: Ring;
    type Embedding: Fn(El<Self::BaseRing>) -> El<Self>;

    ///
    /// Checks whether this ring is in fact an extension of the base ring.
    /// 
    fn is_extension(&self) -> RingPropValue;

    fn base_ring(&self) -> &Self::BaseRing;

    fn embedding(&self) -> Self::Embedding;

    fn from(&self, el: El<Self::BaseRing>) -> El<Self>;
}

///
/// Trait for rings in which one might be able to can provide further information about
/// a quotient a/b of ring elements.
/// 
pub trait DivisibilityInfoRing : Ring {

    ///
    /// Returns whether this ring supports computing divisibility information.
    /// 
    fn is_divisibility_computable(&self) -> RingPropValue;

    ///
    /// Checks whether one element divides another.
    /// This may panic if `is_divisibility_computable()` returns false.
    /// 
    fn is_divisible_by(&self, lhs: &Self::El, rhs: &Self::El) -> bool {
        self.quotient(lhs, rhs).is_some()
    }

    ///
    /// Computes a quotient of two elements, if one divides the other.
    /// If this is not the case, None is returned.
    /// This may panic if `is_divisibility_computable()` returns false.
    /// 
    /// # Uniqueness
    /// 
    /// If rhs is not a zero divisor, this is unique. If this is not the
    /// case and at least one x with x * rhs == lhs exists, any one of those
    /// is returned. In particular, taking the quotient of 0/0 may return
    /// any element of the ring.
    /// 
    fn quotient(&self, lhs: &Self::El, rhs: &Self::El) -> Option<Self::El>;

    ///
    /// Checks whether an element in this ring is a unit.
    /// This may panic if `is_ufd()` returns false.
    /// 
    fn is_unit(&self, el: &Self::El) -> bool {
        self.quotient(&self.one(), el).is_some()
    }
}

///
/// Trait for rings that might be able to provide additional information about prime
/// factorizations. These are assumed to be unique for implementors of this
/// trait (provided that `is_ufd().can_use()` is true).
/// 
pub trait UfdInfoRing : DivisibilityInfoRing {

    ///
    /// Determines whether elements in this ring have a unique prime factorization.
    /// 
    fn is_ufd(&self) -> RingPropValue;

    ///
    /// Checks whether an element in this ring is prime.
    /// This may panic if `is_ufd().can_use()` returns false.
    /// 
    fn is_prime(&self, el: &Self::El) -> bool {
        !self.is_unit(el) && self.calc_factor(el).is_none()
    }

    ///
    /// Returns a nontrivial factor of the given element, or None if the element is prime or a unit.
    /// This may panic if `is_ufd().can_use()` returns false.
    /// 
    fn calc_factor(&self, el: &Self::El) -> Option<Self::El>;

    ///
    /// Factors the given element into all prime factors. The returned list contains pairwise
    /// non-associated primes (p, q are associated if there is a unit e with ep = q) together
    /// with the corresponding power dividing the number, and possibly one furtyher unit.
    /// 
    /// Note that the prime factorization is only unique up to association, i.e. multiplication
    /// by units.
    /// 
    /// This may panic if `is_ufd().can_use()` returns false.
    /// 
    fn factor<'a>(&'a self, el: Self::El) -> VecMap<RingElWrapper<&'a Self>, usize> {
        let mut result = VecMap::new();
        let mut stack = Vec::new();
        stack.push(el);
        let wrapping_ring = WrappingRing::new(self);
        while let Some(el) = stack.pop() {
            let wrapped_el = wrapping_ring.from(el);
            if let Some(factor) = self.calc_factor(wrapped_el.val()) {
                stack.push(self.quotient(wrapped_el.val(), &factor).unwrap());
                stack.push(factor);
            } else if let Some(power) = result.get_mut(&wrapped_el) {
                *power += 1;
            } else {
                result.insert(wrapped_el, 1);
            }
        }
        return result;
    }
}

///
/// Trait for rings that are already completely determined by their type, i.e.
/// contain no further runtime data.
/// 
pub trait SingletonRing: Ring {
    fn singleton() -> Self;
}

///
/// Trait for rings whose elements are hashable. Note that the hashing contract
/// of this function is as usual, but relative to the equality notion given by
/// `ring.is_eq(el1, el2)`.
/// 
pub trait HashableElRing: Ring {
    fn hash<H: std::hash::Hasher>(&self, h: &mut H, el: &Self::El);
}

///
/// Trait for ordered rings, i.e. rings whose elements have a total order that is
/// compatible with + and * in the following sense:
/// 
/// If a >= b and c >= d, then a + c >= b + d.
/// 
/// If a >= 0 and b >= 0, then ab >= 0.
/// 
pub trait OrderedRing: Ring {
    fn cmp(&self, lhs: &Self::El, rhs: &Self::El) -> std::cmp::Ordering;

    fn is_leq(&self, lhs: &Self::El, rhs: &Self::El) -> bool {
        self.cmp(lhs, rhs) != std::cmp::Ordering::Greater
    }

    fn is_geq(&self, lhs: &Self::El, rhs: &Self::El) -> bool {
        self.cmp(lhs, rhs) != std::cmp::Ordering::Less
    }

    fn is_lt(&self, lhs: &Self::El, rhs: &Self::El) -> bool {
        self.cmp(lhs, rhs) == std::cmp::Ordering::Less
    }

    fn is_gt(&self, lhs: &Self::El, rhs: &Self::El) -> bool {
        self.cmp(lhs, rhs) == std::cmp::Ordering::Greater
    }

    fn is_pos(&self, x: &Self::El) -> bool {
        self.is_gt(x, &self.zero())
    }

    fn is_neg(&self, x: &Self::El) -> bool {
        self.is_lt(x, &self.zero())
    }
}