use super::integer::*;
use super::embedding::*;
use super::primitive::*;
use super::ring_property::*;
use super::square_multiply::abs_square_and_multiply;
use super::wrapper::*;

use vector_map::VecMap;

///
/// Utility struct providing a nice way to display ring elements.
/// 
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
/// # `RingBase` vs [`Ring`]
/// 
/// See the corresponding remark in [`Ring`].
/// 
/// # Value vs Reference design rationale
/// 
/// Types that are cheap to copy can be used with the add/sub()-functions, 
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
///  - In order to get [`crate::embedding::StandardEmbedding`] objects (and other kinds of embeddings) that
///    live for static lifetime (required e.g. for [`crate::ring::RingExtension`]), we must be able to
///    clone source & destination ring
///  - Much functionality of [`crate::wrapper::WrappingRing`] requires cloning, even to make it a ring (as
///    elements of a ring are required to be clonable).
/// 
/// I am completely aware that cloning in these situations might lead to a non-obvious
/// problems, e.g. if the ring contains a thread or memory pool. In these cases, we 
/// recommend to implement `Ring` only for objects that contain references or [`std::rc::Rc`]'s to
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

    ///
    /// Raises the given element to the given power. The result should be 
    /// equivalent to chaining `self.mul(basis, ...)` the appropriate number
    /// of times.
    /// 
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

    ///
    /// Raises the given element to the given power. This is equivalent to [`pow()`],
    /// but takes an arbitrarily large integer.
    /// 
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

    ///
    /// Returns whether the ring is noetherian, so whether each ascending chain of ideals is
    /// stationary. 
    /// 
    /// Even though computations with ideals are currently not available, this
    /// can have visible consequences. For example, in a noetherian ring it is impossible that
    /// some x is infinitely often divisible by some y, such that x/y^i != x/y^j for i != j.
    /// 
    fn is_noetherian(&self) -> bool;

    ///
    /// May panic if the ring is not a field (meaning `!is_field().can_use()`). 
    /// If it does not panic, the result must be valid. For a non-field ring, it therefore 
    /// must panic if rhs does not divide lhs, and if it does, it may either compute the 
    /// correct quotient but may also panic nevertheless.
    /// 
    fn div(&self, lhs: Self::El, rhs: &Self::El) -> Self::El;

    ///
    /// Maps an integer into the ring. This is analogeous to [`from_z()`], but
    /// supports arbitrarily large integers.
    /// 
    fn from_z_big(&self, x: &StdInt) -> Self::El {
        let result = abs_square_and_multiply(&self.one(), x, StdInt::RING, |x, y| self.add(x, y), |x, y| self.add_ref(x.clone(), y), self.zero());
        if *x < 0 {
            return self.neg(result);
        } else {
            return result;
        }
    }

    ///
    /// Maps an integer into the ring. The result is equivalent to summing up `self.one()`
    /// an appropriate number of times.
    /// 
    fn from_z(&self, x: i64) -> Self::El {
        let result = abs_square_and_multiply(&self.one(), &x, i64::RING, |x, y| self.add(x, y), |x, y| self.add_ref(x.clone(), y), self.zero());
        if x < 0 {
            return self.neg(result);
        } else {
            return result;
        }
    }

    ///
    /// Maps an integer into the ring. This is analogeous to [`from_z()`], but
    /// supports custom implementations of the integer ring.
    /// 
    /// # Relationship with [`CanonicalEmbeddingInfo<I>`]
    /// 
    /// There is a potential overlap with [`CanonicalEmbeddingInfo<I>`]. Consider
    /// this function to be more high-level, i.e. there is no problem with implementing this by
    /// `self.embed(ring, x)` or even `ring.preimage(self, x)` when appropriate.
    /// 
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

    ///
    /// See [`format()`].
    /// 
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

///
/// Trait to represent a unital, commutative ring. More abstract functionality 
/// (global properties, like ideals) is not provided, this is mainly an interface 
/// that can be used for algorithms that deal with ring elements. 
/// 
/// This trait is the central interface for the crate, as everything else are either
/// algorithms working on rings, or ring combinators that provide standard rings and
/// ring constructions.
/// 
/// # Notes
/// 
/// This trait is mainly to be used in algorithms that work in rings. In particular,
/// if you want to create a new ring type, look below.
/// 
/// If you just want to use existing rings, you might want to have a look at [`crate::wrapper::RingElWrapper`] 
/// and [`crate::wrapper::WrappingRing`], which provide a much more convenient interface
/// to computing with them. In particular, it is then possible to use the
/// natural `+`, `*`, `=` and other operations.
/// 
/// ## `Ring` vs [`RingBase`]
/// 
/// As a rule of thumb, use `Ring` everywhere, except when implementing the traits for
/// your own type (see below). The main reason is that `RingBase` defines the basic functionality,
/// and `Ring` adds some (slightly) higher-level functionality on top of that (in particular
/// embeddings & isomorphisms between rings, see [`crate::embedding::CanonicalEmbeddingInfo`]).
/// 
/// ## Implementation
/// 
/// To implement this trait, implement [`RingBase`] and [`crate::embedding::CanonicalIsomorphismInfo<Self>`].
/// For more information, see [`RingBase`].
/// 
/// # Specific rings & properties
/// 
/// Rings can have many properties that are not provided by the basic `Ring` trait (it only contains
/// the query functions [`is_field()`], [`is_noetherian()`] and [`is_integral()`]). For further properties,
/// there exist traits (usually using the naming convention `...Info`) that mark ring type that might
/// have the corresponding property. However, since this property may depend on runtime information, you
/// always have to use the corresponding query function of the specific trait to check.
/// Examples are [`crate::ring::EuclideanInfoRing`] or [`crate::ring::DivisibilityInfoRing`].
/// 
/// We also consider the property of being "a specific ring", e.g. a polynomial ring as a property, and hence
/// there exists a corresponding trait for it (e.g. [`crate::poly::PolyRing`]). This trait is then implemented
/// by all rings that represent the ring (as a purely mathematical object) in question. Usually, this trait
/// comes together with a default implementation (e.g. [`crate::poly::uni_var::PolyRingImpl`]), and the trait
/// also requires all implementations to be canonically isomorphic to this default implementation.
/// 
/// # Examples
/// 
/// ## Integers
/// 
/// The most classical ring is probably the integers Z. There are many different implementations
/// of this ring, with different performance characteristics.
/// ```
/// # use feanor_la::prelude::*;
/// type Ztype = StaticRing<i64>;
/// let Z: Ztype = i64::RING; // StaticRing<_> implements Ring
/// let a: El<Ztype> = 7i64;
/// let b: El<Ztype> = 6i64;
/// assert!(Z.is_eq(&13, &Z.add(a, b)));
/// ```
/// Of course, this example is equivalent to the standard, more convenient notation
/// ```
/// # use feanor_la::prelude::*;
/// let a: i64 = 7;
/// let b: i64 = 6;
/// assert!(13 == a + b);
/// ```
/// To work with arbitrarily large integers, the most important ring is [`crate::integer::bigint_soo::BigIntSOORing`],
/// usually used through [`crate::integer::StdInt`] (note that this uses [`crate::wrapper::RingElWrapper`] internally,
/// for improved convenience).
/// ```
/// # use feanor_la::prelude::*;
/// type Ztype = <StdInt as RingEl>::RingType;
/// let Z: Ztype = StdInt::RING; // StaticRing<_> implements Ring
/// let a: El<Ztype> = StdInt::from(7);
/// let b: El<Ztype> = StdInt::from(6);
/// assert!(Z.is_eq(&StdInt::from(13), &Z.add(a, b)));
/// ```
/// Since `StdInt` is wrapped, we can also use a nice notation
/// ```
/// # use feanor_la::prelude::*;
/// let a: StdInt = StdInt::from(7);
/// let b: StdInt = StdInt::from(6);
/// assert!(StdInt::from(13) == a + b);
/// ```
/// However, if directly use the underlying ring [`crate::integer::bigint_soo::BigIntSOORing`],
/// we cannot do this anymore. Instead, we would have to write
/// ```
/// # use feanor_la::prelude::*;
/// # use feanor_la::integer::bigint_soo::*;
/// type Ztype = BigIntSOORing;
/// let Z = BigIntSOORing::singleton();
/// let a: El<Ztype> = Z.from_z(7);
/// let b: El<Ztype> = Z.from_z(6);
/// assert!(Z.is_eq(&BigIntSOO::RING.from_z(13), &Z.add(a, b)));
/// ```
/// However, the wrapper around `BigIntSOORing` is trivial, so we do not suffer any performance
/// loss when using the above, nicer variant.
/// 
/// ## Other rings
/// For example polynomial rings are provided by [`crate::poly::uni_var::PolyRing`].
/// ```
/// # use feanor_la::prelude::*;
/// # use feanor_la::poly::uni_var::*;
/// # use feanor_la::poly::*;
/// let P = PolyRingImpl::adjoint(i64::RING, "X");
/// let x = P.unknown();
/// let f = P.add(x.clone(), P.one()); // x + 1
/// let g = P.add(P.pow(&x, 2), P.add(P.mul(x, P.from_z(2)), P.one())); // x^2 + 2x + 1
/// assert!(P.is_eq(&g, &P.mul_ref(&f, &f)));
/// ```
/// As this shows, this ring-based interface can easily get very cumbersome. To
/// avoid this, look at [`crate::wrapper::RingElWrapper`].
/// 
pub trait Ring: RingBase + CanonicalIsomorphismInfo<Self> {}

impl<R: ?Sized> Ring for R
    where R: RingBase + CanonicalIsomorphismInfo<R>
{}

#[allow(type_alias_bounds)]
pub type El<R: Ring> = <R as RingBase>::El;

///
/// Trait for rings that might have a euclidean division, that is they provide
/// a function `d: R -> N` and a division procedure `x, y -> x = qy + r` such that
/// `d(r) < d(y)` (if y is nonzero).
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
    /// # Uniqueness
    ///
    /// The euclidean division does not have to be unique. In the case of integers for
    /// example, we have that `-5 = 2 (-3) + 1` and `-5 = 1(-3) + (-2)` are both valid
    /// euclidean division identities. An implementation is allowed to return any of
    /// those.
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

    ///
    /// Returns the canonical ring homomorphism `Self::BaseRing -> Self`.
    /// 
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
    /// This may panic if `is_divisibility_computable().can_use()` returns false.
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
    /// with the corresponding power dividing the number, and possibly one further unit.
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
            let wrapped_el = wrapping_ring.wrap(el);
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
/// For a similar functionality, see [`crate::primitive::StaticRing`].
/// 
/// # Difference of `SingletonRing` and [`crate::primitive::StaticRing`]
/// 
/// Both `SingletonRing` and `StaticRing` can be used for "constant-ring" rings, i.e.
/// ring types that have only one sensible instance. Hence, from an abstract perspective,
/// they are equivalent. However, due to the different kind of implementation, they
/// have significantly different properties.
/// 
/// In particular, `StaticRing` is always used with elements of type implementing
/// [`crate::primitive::RingEl`], and so the elements can be multiplied without having
/// any ring instance at all. Furthermore, a ring instance must be a compile-time constant,
/// thus very small and trivial to copy.
/// 
/// On the other hand, there is no reason why an implementation of `SingletonRing` should be
/// small, it could for example contain a memory pool to reuse memory for elements. In this case,
/// we really require access to the ring object for operations, as creating one on-the-fly every
/// time would be extremely expensive.
/// 
/// As a summary, `StaticRing` is more narrow than `SingletonRing`. The struct `StaticRing` does
/// implement `SingletonRing`, but not every ring implementing `SingletonRing` can be sensibly
/// made into `StaticRing<E>` for some `RingEl`-type `E`.
/// 
/// Note that one could bridge this gap using lazy statics, but that seems like an unnecessary
/// effort at the moment.
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