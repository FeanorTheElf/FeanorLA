#![feature(adt_const_params)]
#![feature(generic_const_exprs)]
#![feature(specialization)]
#![feature(associated_type_bounds)]
#![feature(step_trait)]
#![feature(result_option_inspect)]
#![feature(fn_traits)]
#![feature(test)]
#![feature(unboxed_closures)]
#![feature(trace_macros)]
#![feature(bench_black_box)]

#[cfg(test)]
extern crate test;
extern crate oorandom;
extern crate vector_map;
extern crate libc;

pub mod ring_property;
pub mod ring;
pub mod primitive;
pub mod embedding;
pub mod wrapper;
pub mod static_ring;
pub mod float;
pub mod ring_decorator;

pub mod square_multiply;
pub mod multiplication;
pub mod eea;
pub mod lattice;

pub mod prelude;

pub mod la;
pub mod integer;
pub mod fq;
pub mod rational;
pub mod factoring_algorithms;
pub mod p_adic;
pub mod fraction_field;
pub mod poly;
pub mod finite_extension;
pub mod elliptic_curve;
pub mod finite_field_sqrt;
pub mod discrete_log;
pub mod number_field;
pub mod extension_wrapper;

pub mod numerics;
pub mod combinatorics;

#[cfg(feature = "mpir")]
pub mod mpir;