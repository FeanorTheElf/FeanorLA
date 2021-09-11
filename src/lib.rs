#![feature(const_generics)]
#![feature(specialization)]
#![feature(core_intrinsics)]
#![feature(test)]
#![feature(const_panic)]
#![feature(unboxed_closures)]
#![feature(trace_macros)]
#![feature(bench_black_box)]
#![feature(const_fn_trait_bound)]

#[cfg(test)]
extern crate test;
extern crate take_mut;
extern crate oorandom;
extern crate array_init;

pub mod alg;
#[macro_use]
pub mod alg_macros;
pub mod float;
pub mod alg_env;

pub mod la;
pub mod numerics;
pub mod algebra;
pub mod combinatorics;
