#![feature(const_generics)]
#![feature(specialization)]
#![feature(fn_traits)]
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

pub mod bigint;
pub mod ring;
pub mod primitive;
pub mod embedding;
pub mod wrapper;
pub mod float;

#[macro_use]
pub mod la;
pub mod algebra;

pub mod numerics;
pub mod combinatorics;
