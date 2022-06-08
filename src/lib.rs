#![feature(adt_const_params)]
#![feature(generic_const_exprs)]
#![feature(specialization)]
#![feature(step_trait)]
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
pub mod ring;
pub mod square_multiply;
pub mod primitive;
pub mod embedding;
pub mod wrapper;
pub mod float;
pub mod prelude;
pub mod karatsuba;
pub mod eea;

#[macro_use]
pub mod la;
pub mod fq;
pub mod integer;
pub mod rational;
pub mod factoring_algorithms;
pub mod p_adic;
pub mod fraction_field;
pub mod poly;
pub mod ring_extension;
pub mod primes;
pub mod diophantine;
pub mod elliptic_curve;
pub mod finite_field_sqrt;
pub mod discrete_log;

pub mod numerics;
pub mod combinatorics;
