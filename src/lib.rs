#![feature(const_generics)]
#![feature(specialization)]
#![feature(test)]
#![feature(array_value_iter)]
#![feature(const_fn)]
#![feature(const_panic)]
#![feature(unboxed_closures)]

#[cfg(test)]
extern crate test;

#[macro_use]
pub mod alg;
pub mod la;
pub mod numerics;
pub mod algebra;
pub mod combinatorics;
