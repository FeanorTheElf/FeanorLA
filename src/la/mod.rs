pub mod vector_view;
pub mod vector;
pub mod matrix_view;
pub mod matrix_owned;
pub mod matrix_row_col;
pub mod matrix_vector;
pub mod submatrix;
pub mod vec;
pub mod mat;

#[cfg(test)]
#[macro_use]
pub mod macros;

pub mod qr;
pub mod simplex;
pub mod prelude;