#![feature(test)]
#![allow(non_snake_case)]
#![feature(array_value_iter)]
#![feature(const_generics)]
#![feature(const_generic_impls_guard)]

#[cfg(test)]
extern crate test;

#[cfg(test)]
#[macro_use]
mod macros;

pub mod alg;
pub mod indexed;
pub mod vector;
#[macro_use]
pub mod matrix;
pub mod prelude;
pub mod rat;
pub mod simplex;
pub mod diophantine;
pub mod qr;
