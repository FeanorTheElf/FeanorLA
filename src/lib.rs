#![feature(test)]
#![allow(non_snake_case)]
#![feature(array_value_iter)]
#![feature(const_generics)]

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
pub mod submatrix;
pub mod matrix_owned;
pub mod matrix_vector;

pub mod mat;

// pub mod prelude;
// pub mod rat;
// pub mod simplex;
// pub mod diophantine;
// pub mod qr;
