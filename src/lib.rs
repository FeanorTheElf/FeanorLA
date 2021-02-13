#![feature(test)]
#![feature(array_value_iter)]
#![feature(const_generics)]
#![feature(const_fn)]
#![feature(const_panic)]
#![feature(unboxed_closures)]
#![allow(non_snake_case)]

#[cfg(test)]
extern crate test;

#[cfg(test)]
#[macro_use]
mod macros;

pub mod alg;
pub mod eea;

pub mod vector_view;
pub mod vector;
pub mod matrix_view;
pub mod matrix_row_col;
pub mod submatrix;
pub mod matrix_owned;
pub mod matrix_vector;

pub mod vec;
pub mod mat;

pub mod prelude;
pub mod rat;
pub mod bigint;
pub mod simplex;
pub mod diophantine;
pub mod qr;
pub mod zn;
pub mod sieve;