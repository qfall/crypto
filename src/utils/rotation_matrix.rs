// Copyright © 2023 Marvin Beckmann
//
// This file is part of qFALL-crypto.
//
// qFALL-crypto is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

//! This module includes a specialized implementation called rotation matrices,
//! which find application in ring-based implementations for the special rings of the form
//! `Z[X]/(X^n + 1)`.

use qfall_math::{
    integer::MatZ,
    traits::{Concatenate, GetEntry, GetNumColumns, GetNumRows, SetEntry},
};

/// Takes in a vector and computes and the rotation matrix as follows.
/// For a vector `[[a_1],[a_2],...,[a_n]]` it computes the rotation matrix as
/// `[[a_1, -a_n, ..., -a_2],[a_2, a_1, ..., -a_3],...,[a_n, a_{n-1}, ..., a_1]]`
///
/// It takes in both a column vector or a row vector, but the format of the matrix
/// remains as described above.
///
/// Parameters:
/// - `vec`: The vector for which the rotation matrix will be computed.
///
/// Returns the rotation matrix of `vec` as a [`MatZ`].
///
/// # Examples
/// ```
/// use qfall_crypto::utils::rotation_matrix::rot_minus;
/// use qfall_math::integer::MatZ;
/// use std::str::FromStr;
///
/// let vec = MatZ::from_str("[[1],[5],[-1],[9]]").unwrap();
/// let row_vec = MatZ::from_str("[[1,5,-1,9]]").unwrap();
///
/// let rot_col = rot_minus(&vec);
/// let rot_row = rot_minus(&row_vec);
/// ```
///
/// # Panics ...
/// - if the provided matrix is not a vector.
pub fn rot_minus(vec: &MatZ) -> MatZ {
    let vec = if vec.is_column_vector() {
        vec.clone()
    } else if vec.is_row_vector() {
        vec.transpose()
    } else {
        panic!("The input must be a vector.")
    };

    let mut out = MatZ::new(vec.get_num_rows(), vec.get_num_rows());

    for i in 0..out.get_num_rows() {
        let entry = vec.get_entry(i, 0).unwrap();
        for j in 0..out.get_num_columns() {
            let (row, sign) = match i + j {
                k if k >= out.get_num_rows() => ((k) % out.get_num_rows(), -1),
                k => (k, 1),
            };
            out.set_entry(row, j, sign * &entry).unwrap();
        }
    }
    out
}

/// Takes in a matrix, splits the matrix into separate columns and concatenates the
/// rotation matrices as:
/// `[rot^-(a_1) | rot^-(a_2) | ... | rot^-(a_m)]`
///
/// It takes in both a column vector or a row vector, but the format of the matrix
/// remains as described above.
///
/// Parameters:
/// - `vec`: The vector for which the rotation matrix will be computed.
///
/// # Examples
/// ```
/// use qfall_crypto::utils::rotation_matrix::rot_minus_matrix;
/// use qfall_math::integer::MatZ;
/// use std::str::FromStr;
///
/// let mat = MatZ::from_str("[[1,5,-1,9],[2,3,4,5]]").unwrap();
///
/// let rot_mat = rot_minus_matrix(&mat);
/// ```
pub fn rot_minus_matrix(matrix: &MatZ) -> MatZ {
    let mut vec = Vec::new();
    for i in 0..matrix.get_num_columns() {
        vec.push(rot_minus(&matrix.get_column(i).unwrap()));
    }

    let mut out = vec.first().unwrap().clone();
    for mat in vec.iter().skip(1) {
        out = out.concat_horizontal(mat).unwrap();
    }
    out
}

#[cfg(test)]
mod test_rot_minus {
    use crate::utils::rotation_matrix::{rot_minus, rot_minus_matrix};
    use qfall_math::integer::MatZ;
    use std::str::FromStr;

    /// Ensures that the rotation minus matrix works correctly for a vector.
    #[test]
    fn correct_rotation_matrix_vec() {
        let vec = MatZ::from_str("[[1],[5],[-1],[9]]").unwrap();
        let row_vec = MatZ::from_str("[[1,5,-1,9]]").unwrap();

        let rot_col = rot_minus(&vec);
        let rot_row = rot_minus(&row_vec);

        let cmp_rot =
            MatZ::from_str("[[1, -9, 1, -5],[5, 1, -9, 1],[-1, 5, 1, -9],[9, -1, 5, 1]]").unwrap();
        assert_eq!(rot_col, rot_row);
        assert_eq!(cmp_rot, rot_col)
    }

    /// Ensures that the rotation minus matrix works correctly for matrices.
    #[test]
    fn correct_rotation_matrix_mat() {
        let mat = MatZ::from_str(&format!("[[1,5,-1,9],[{}, 1, 2, 3]]", u64::MAX)).unwrap();

        let rot_mat = rot_minus_matrix(&mat);

        let cmp_rot = MatZ::from_str(&format!(
            "[[1, -{}, 5, -1, -1, -2, 9, -3],\
            [{}, 1, 1, 5, 2, -1, 3, 9]]",
            u64::MAX,
            u64::MAX
        ))
        .unwrap();
        assert_eq!(cmp_rot, rot_mat);
    }

    /// Ensures that the minus rotation matrix for vectors panics if a matrix is provided.
    #[test]
    #[should_panic]
    fn not_vector() {
        let mat = MatZ::from_str("[[1,5,-1,9],[1,2,3,4]]").unwrap();

        let _ = rot_minus(&mat);
    }
}
