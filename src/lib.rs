#![feature(test)]
#![allow(non_snake_case)]

#[cfg(test)]
extern crate test;

#[cfg(test)]
#[macro_use]
mod macros;

pub mod arith;
pub mod diophantine;
pub mod indexed;
pub mod matrix;
pub mod prelude;
pub mod rat;
pub mod simplex;
