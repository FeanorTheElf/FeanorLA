#![feature(const_generics)]
#![feature(specialization)]
#![feature(test)]
#![feature(const_fn)]
#![feature(const_panic)]
#![feature(unboxed_closures)]
#![feature(trace_macros)]

#[cfg(test)]
extern crate test;

#[macro_use]
pub mod alg;
pub mod float;
#[macro_use]
pub mod alg_env;

pub mod la;
pub mod numerics;
pub mod algebra;
pub mod combinatorics;
