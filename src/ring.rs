use super::bigint::*;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RingPropValue {
    True, False, NotComputable, Unknown
}

impl RingPropValue {

    pub fn can_use(&self) -> bool {
        *self == RingPropValue::True
    }
}

impl std::fmt::Display for RingPropValue {

    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            RingPropValue::True => write!(f, "true"),
            RingPropValue::False => write!(f, "false"),
            RingPropValue::NotComputable => write!(f, "not computable"),
            RingPropValue::Unknown => write!(f, "unknown")
        }
    }
}

pub struct RingElDisplay<'a, R: ?Sized> 
    where R: Ring
{
    ring: &'a R,
    el: &'a R::El
}

impl<'a, R> std::fmt::Display for RingElDisplay<'a, R>
where R: Ring
{
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        self.ring.format(self.el, f, false)
    }
}

///
/// Trait to represent a commutative ring with one as a collection of operations 
/// on the elements. More abstract functionality (global properties, like ideals) 
/// is not provided, this is mainly an interface that can be used for algorithms 
/// that deal with ring elements.
/// 
pub trait Ring : std::fmt::Debug + std::clone::Clone {
    type El: Sized + Clone + std::fmt::Debug;

    //
    // Design rationale: Types that are cheap to copy can be used with the 
    // add/sub()-functions, for types that are expensive to copy, one can use the 
    // add_ref/sub_ref()-functions. However, as each add/sub() / add_ref/sub_ref() 
    // function returns the result by value, for big types, is is usually 
    // most efficient to clone at least one parameter and potentially reuse 
    // the memory.
    //
    // If an operation can improve efficience by consuming both parameters, one should
    // explicitly implement the default-implemented add/sub()-functions.
    //
    // For multiplication, the situation is usually different: Mostly, multiplication
    // is given by some kind of operation on a cartesian product of the
    // components of lhs and rhs, with a following reduce. In this case, one usually
    // cannot profit from getting the parameters by value. If this is not true, then
    // one should implement the default-implemented mul()-function
    //

    fn add_ref(&self, lhs: Self::El, rhs: &Self::El) -> Self::El;

    ///
    /// Calculates the product of lhs and rhs. Note that multiplication is assumed to
    /// be commutative.
    /// 
    fn mul_ref(&self, lhs: &Self::El, rhs: &Self::El) -> Self::El;
    fn neg(&self, val: Self::El) -> Self::El;
    fn zero(&self) -> Self::El;
    fn one(&self) -> Self::El;
    fn eq(&self, lhs: &Self::El, rhs: &Self::El) -> bool;

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
    /// goes out of scope by the panic, if however the panic is caught (might requrie
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
    fn sum(&self, data: &[Self::El]) -> Self::El {
        data.iter().fold(self.zero(), |a, b| self.add_ref(a, b))
    }
    fn mul(&self, lhs: Self::El, rhs: Self::El) -> Self::El { self.mul_ref(&lhs, &rhs) }
    fn prod(&self, data: &[Self::El]) -> Self::El {
        data.iter().fold(self.one(), |a, b| self.mul_ref(&a, b))
    }

    fn sub(&self, lhs: Self::El, rhs: Self::El) -> Self::El {
        self.add(lhs, self.neg(rhs))
    }

    fn pow(&self, basis: &Self::El, exp: u32) -> Self::El 
        where Self::El: Clone
    {
        self.pow_big(basis, &BigInt::from(exp as i64))
    }

    fn pow_big(&self, basis: &Self::El, exp: &BigInt) -> Self::El 
        where Self::El: Clone
    {
        assert!(*exp >= 0);
        if exp.is_zero() {
            return self.one();
        }

        let mut result = self.one();
        for i in (0..(exp.log2_floor() + 1)).rev() {
            if exp.is_bit_set(i) {
                result = self.mul(self.mul_ref(&result, &basis), result);
            } else {
                result = self.mul_ref(&result, &result);
            }
        }
        return result;
    }

    fn is_zero(&self, val: &Self::El) -> bool { self.eq(val, &self.zero()) }
    fn is_one(&self, val: &Self::El) -> bool { self.eq(val, &self.one()) }
    fn is_neg_one(&self, val: &Self::El) -> bool { self.eq(val, &self.neg(self.one())) }

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

    fn from_z_big(&self, x: BigInt) -> Self::El {
        if x == 0 {
            return self.zero();
        }
        let mut result = self.zero();
        for i in (0..(x.log2_floor() + 1)).rev() {
            if x.is_bit_set(i) {
                result = self.add(self.add_ref(self.one(), &result), result);
            } else {
                result = self.add_ref(result.clone(), &result);
            }
        }
        if x < 0 {
            return self.neg(result);
        } else {
            return result;
        }
    }

    fn from_z(&self, x: i64) -> Self::El {
        if x == 0 {
            return self.zero();
        }
        let mut result = self.zero();
        for i in (0..63).rev() {
            if (x.abs() >> i) & 1 == 1 {
                result = self.add(self.add_ref(self.one(), &result), result);
            } else {
                result = self.add_ref(result.clone(), &result);
            }
        }
        if x < 0 {
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

impl<'a, R: Ring> Ring for &'a R {
    type El = R::El;

    fn add_ref(&self, lhs: Self::El, rhs: &Self::El) -> Self::El { (**self).add_ref(lhs, rhs) }
    fn mul_ref(&self, lhs: &Self::El, rhs: &Self::El) -> Self::El { (**self).mul_ref(lhs, rhs) }
    fn neg(&self, val: Self::El) -> Self::El { (**self).neg(val) }
    fn zero(&self) -> Self::El { (**self).zero() }
    fn one(&self) -> Self::El { (**self).one() }
    fn eq(&self, lhs: &Self::El, rhs: &Self::El) -> bool { (**self).eq(lhs, rhs) }
    fn unspecified_element(&self) -> Self::El { (**self).unspecified_element() }
    fn sub_ref_fst(&self, lhs: &Self::El, rhs: Self::El) -> Self::El { (**self).sub_ref_fst(lhs, rhs) }
    fn sub_ref_snd(&self, lhs: Self::El, rhs: &Self::El) -> Self::El { (**self).sub_ref_snd(lhs, rhs) }
    fn add(&self, lhs: Self::El, rhs: Self::El) -> Self::El { (**self).add(lhs, rhs) }
    fn mul(&self, lhs: Self::El, rhs: Self::El) -> Self::El { (**self).mul(lhs, rhs) }
    fn sub(&self, lhs: Self::El, rhs: Self::El) -> Self::El { (**self).sub(lhs, rhs) }
    fn pow(&self, basis: &Self::El, exp: u32) -> Self::El 
        where Self::El: Clone
    {
        (**self).pow(basis, exp)
    }
    fn pow_big(&self, basis: &Self::El, exp: &BigInt) -> Self::El 
        where Self::El: Clone
    {
        (**self).pow_big(basis, exp)
    }
    fn is_zero(&self, val: &Self::El) -> bool { (**self).is_zero(val) }
    fn is_one(&self, val: &Self::El) -> bool { (**self).is_one(val) }
    fn is_neg_one(&self, val: &Self::El) -> bool { (**self).is_neg_one(val) }
    fn is_integral(&self) -> RingPropValue { (**self).is_integral() }
    fn is_field(&self) -> RingPropValue { (**self).is_field() }
    fn is_noetherian(&self) -> bool { (**self).is_noetherian() }
    fn div(&self, lhs: Self::El, rhs: &Self::El) -> Self::El { (**self).div(lhs, rhs) }
    fn format(&self, el: &Self::El, f: &mut std::fmt::Formatter, in_prod: bool) -> std::fmt::Result { (**self).format(el, f, in_prod) }
    fn format_in_brackets(&self, el: &Self::El, f: &mut std::fmt::Formatter) -> std::fmt::Result { (**self).format_in_brackets(el, f) }
}

pub trait EuclideanInfoRing: Ring {
    
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
    fn euclidean_deg(&self, el: Self::El) -> BigInt;
}

impl<'a, R: EuclideanInfoRing> EuclideanInfoRing for &'a R {
    
    fn is_euclidean(&self) -> RingPropValue { (**self).is_euclidean() }
    fn euclidean_div_rem(&self, lhs: Self::El, rhs: &Self::El) -> (Self::El, Self::El) { (**self).euclidean_div_rem(lhs, rhs) }
    fn euclidean_rem(&self, lhs: Self::El, rhs: &Self::El) -> Self::El { (**self).euclidean_rem(lhs, rhs) }
    fn euclidean_div(&self, lhs: Self::El, rhs: &Self::El) -> Self::El { (**self).euclidean_div(lhs, rhs) }
    fn euclidean_deg(&self, el: Self::El) -> BigInt { (**self).euclidean_deg(el) }
}

pub trait DivisibilityInfoRing : Ring {

    ///
    /// Returns whether this ring supports computing divisibility information.
    /// 
    fn is_divisibility_computable(&self) -> bool;

    ///
    /// Checks whether one element divides another.
    /// This may panic if `is_divisibility_computable()` returns false.
    /// 
    fn is_divisible_by(&self, lhs: &Self::El, rhs: &Self::El) -> bool {
        self.quotient(lhs, rhs).is_some()
    }

    ///
    /// Computes the quotient of two elements, if one divides the other.
    /// If this is not the case, None is returned.
    /// This may panic if `is_divisibility_computable()` returns false.
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

impl<'a, R> DivisibilityInfoRing for &'a R
    where R: DivisibilityInfoRing
{
    fn is_divisibility_computable(&self) -> bool { (**self).is_divisibility_computable() }
    fn is_divisible_by(&self, lhs: &Self::El, rhs: &Self::El) -> bool { (**self).is_divisible_by(lhs, rhs) }
    fn quotient(&self, lhs: &Self::El, rhs: &Self::El) -> Option<Self::El> { (**self).quotient(lhs, rhs) }
    fn is_unit(&self, el: &Self::El) -> bool { (**self).is_unit(el) }
}

pub trait FactoringInfoRing : DivisibilityInfoRing {

    ///
    /// Determines whether elements in this ring have a unique prime factorization.
    /// This may return false even for UFDs, if it is intractable to compute this
    /// information or intractable to compute prime factorizations. 
    /// 
    fn is_ufd(&self) -> RingPropValue;
    ///
    /// Checks whether an element in this ring is prime.
    /// This may panic if `is_ufd().can_use()` returns false.
    /// 
    fn is_prime(&self, el: &Self::El) -> bool;
    ///
    /// Returns a nontrivial factor of the given element, or None if the element is a unit.
    /// This may panic if `is_ufd().can_use()` returns false.
    /// 
    fn calc_factor(&self, el: &mut Self::El) -> Option<Self::El>;
}

impl<'a, R> FactoringInfoRing for &'a R
    where R: FactoringInfoRing
{
    fn is_ufd(&self) -> RingPropValue { (**self).is_ufd() }
    fn is_prime(&self, el: &Self::El) -> bool { (**self).is_prime(el) }
    fn calc_factor(&self, el: &mut Self::El) -> Option<Self::El> { (**self).calc_factor(el) }
}

pub trait SingletonRing: Ring {

    fn singleton() -> Self;
}