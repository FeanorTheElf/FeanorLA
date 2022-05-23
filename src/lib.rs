#![feature(adt_const_params)]
#![feature(generic_const_exprs)]
#![feature(specialization)]
#![feature(generic_associated_types)]
#![feature(fn_traits)]
#![feature(test)]
#![feature(unboxed_closures)]
#![feature(trace_macros)]
#![feature(bench_black_box)]

#[cfg(test)]
extern crate test;
extern crate oorandom;
extern crate vector_map;

pub mod ring_property;
pub mod bigint;
pub mod ring;
pub mod primitive;
pub mod embedding;
pub mod wrapper;
pub mod float;
pub mod prelude;

#[macro_use]
pub mod la;
pub mod algebra;

pub mod numerics;
pub mod combinatorics;
