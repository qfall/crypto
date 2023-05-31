// Copyright © 2023 Marvin Beckmann
//
// This file is part of qFALL-crypto.
//
// qFALL-crypto is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

//! This module contains an implementation to generate a gadget trapdoor in a
//! classical setting.

use super::gadget_parameters::GadgetParameters;
use qfall_math::{
    error::MathError,
    integer::{MatZ, Z},
    integer_mod_q::MatZq,
    traits::{Concatenate, Pow, SetEntry, Tensor},
};
use std::fmt::Display;

/// Generates a trapdoor according to Algorithm 1 in [\[1\]](<../index.html#:~:text=[1]>).
/// Namely:
/// - Generates the gadget matrix: `G`
/// - Samples the trapdoor `R` from the specified distribution in `params`
/// - Outputs `([a_bar | tag * g - a_bar * r], r)` as a tuple of `(A,R)`,
/// where `R` defines a trapdoor for `A`.
///
/// Parameters:
/// - `params`: all gadget parameters which are required to generate the trapdoor
/// - `a_bar`: the matrix defining the first part of the G-Trapdoor
/// - `tag`: the tag which is hidden within the matrix ´A`
///
/// Returns a a parity-check matrix `a` derived from `a_bar` and its gadget-trapdoor `r`
/// under a give tag `h`.
///
/// # Examples
/// ```
/// use qfall_crypto::sample::g_trapdoor::{gadget_parameters::GadgetParameters, gadget_classical::gen_trapdoor};
/// use qfall_math::integer::Z;
/// use qfall_math::integer_mod_q::{Modulus, MatZq};
///
/// let params = GadgetParameters::init_default(42, &Modulus::try_from(&Z::from(42)).unwrap());
/// let a_bar = MatZq::sample_uniform(42, &params.m_bar, &params.q).unwrap();
/// let tag = MatZq::identity(42, 42, &params.q).unwrap();
///
/// let (a,r) = gen_trapdoor(&params, &a_bar, &tag).unwrap();
/// ```
/// # Errors and Failures
/// - Returns a [`MathError`] of type [`InvalidMatrix`](MathError::InvalidMatrix)
/// or of type [`OutOfBounds`](MathError::OutOfBounds), if `params.k`
/// or `params.n` is either `0`,
/// it is negative or it does not fit into an [`i64`].
pub fn gen_trapdoor(
    params: &GadgetParameters,
    a_bar: &MatZq,
    tag: &MatZq,
) -> Result<(MatZq, MatZ), MathError> {
    let g = gen_gadget_mat(&params.n, &params.k, &params.base)?;
    let r = params
        .distribution
        .sample(&params.m_bar, &(&params.n * &params.k));
    // set A = [\bar A | HG - \bar A R]
    let a = a_bar.concat_horizontal(&(tag * g - a_bar * &r))?;
    Ok((a, r))
}

/// Generates a gadget matrix based on its definition in [\[1\]](<../index.html#:~:text=[1]>).
/// This corresponds to `I_n \oplus g^t` where `g` is a gadget vector for the `base`.
///
/// Parameters:
/// - `n`: the size of the identity matrix, with which the tensor product is defined
/// - `k`: the size of the gadget vector
/// - `base`: the base with which the entries in the gadget vector are defined
///
/// Returns a gadget matrix of size `n*nk` with `base` as the base for the gadget vector.
///
/// # Examples
/// ```
/// use qfall_crypto::sample::g_trapdoor::gadget_classical::gen_gadget_mat;
/// use qfall_math::integer::Z;
///
/// let g = gen_gadget_mat(&Z::from(3), &Z::from(4), &Z::from(2));
/// ```
/// # Errors and Failures
/// - Returns a [`MathError`] of type [`InvalidMatrix`](MathError::InvalidMatrix)
/// or of type [`OutOfBounds`](MathError::OutOfBounds), if `k` or `n` is either `0`,
/// it is negative or it does not fit into an [`i64`].
pub fn gen_gadget_mat(
    n: impl TryInto<i64> + Display + Clone,
    k: impl TryInto<i64> + Display,
    base: &Z,
) -> Result<MatZ, MathError> {
    let gadget_vec = gen_gadget_vec(k, base)?;
    let identity = MatZ::identity(n.clone(), n)?;
    Ok(identity.tensor_product(&gadget_vec.transpose()))
}

/// Generates a gadget vector based on its definition in [\[1\]](<../index.html#:~:text=[1]>).
/// This corresponds to a vector `(base ^0, base^1, ..., base^{k-1})`
///
/// Parameters:
/// - `k`: the size of the gadget vector
/// - `base`: the base with which the entries in the gadget vector are defined
///
/// Returns a gadget vector of length `k` with `base` as its base.
///
/// # Examples
/// ```
/// use qfall_crypto::sample::g_trapdoor::gadget_classical::gen_gadget_vec;
/// use qfall_math::integer::Z;
///
/// let g = gen_gadget_vec(4, &Z::from(2));
/// ```
///
/// # Errors and Failures
/// - Returns a [`MathError`] of type [`InvalidMatrix`](MathError::InvalidMatrix)
/// or of type [`OutOfBounds`](MathError::OutOfBounds), if `k` is either `0`,
/// it is negative or it does not fit into an [`i64`].
pub fn gen_gadget_vec(k: impl TryInto<i64> + Display, base: &Z) -> Result<MatZ, MathError> {
    let mut out = MatZ::new(k, 1).unwrap();
    let mut i: i64 = 0;
    while out.set_entry(i, 0, &base.pow(i)?).is_ok() {
        i += 1;
    }
    Ok(out)
}

#[cfg(test)]
mod test_gen_gadget_vec {
    use crate::sample::g_trapdoor::gadget_classical::gen_gadget_vec;
    use qfall_math::integer::{MatZ, Z};
    use std::str::FromStr;

    /// assure that the gadget vector with base `2` and length `5` works correctly
    #[test]
    fn correctness_base_2() {
        let gadget_vec = gen_gadget_vec(5, &Z::from(2)).unwrap();

        let vec = MatZ::from_str("[[1],[2],[4],[8],[16]]").unwrap();
        assert_eq!(vec, gadget_vec);
    }

    /// assure that the gadget vector with base `5` and length `4` works correctly
    #[test]
    fn correctness_base_5() {
        let gadget_vec = gen_gadget_vec(4, &Z::from(5)).unwrap();

        let vec = MatZ::from_str("[[1],[5],[25],[125]]").unwrap();
        assert_eq!(vec, gadget_vec);
    }
}

#[cfg(test)]
mod test_gen_gadget_mat {
    use super::gen_gadget_mat;
    use qfall_math::integer::{MatZ, Z};
    use std::str::FromStr;

    /// assure that the gadget matrix with gadget vector `[1, 2, 4]^t`(base 3) and
    /// `I_3` works correctly
    #[test]
    fn correctness_base_2_3x3() {
        let gadget_mat = gen_gadget_mat(3, 3, &Z::from(2)).unwrap();

        let mat_str = "[[1, 2, 4, 0, 0, 0, 0, 0, 0],\
                            [0, 0, 0, 1, 2, 4, 0, 0, 0],\
                            [0, 0, 0, 0, 0, 0, 1, 2, 4]]";

        let mat = MatZ::from_str(mat_str).unwrap();
        assert_eq!(mat, gadget_mat);
    }

    /// assure that the gadget matrix with gadget vector `[1, 3, 9, 27, 81]^t`(base 3) and
    /// `I_2` works correctly
    #[test]
    fn correctness_base_3_2x5() {
        let gadget_mat = gen_gadget_mat(2, 5, &Z::from(3)).unwrap();

        let mat_str = "[[1, 3, 9, 27, 81, 0, 0, 0, 0, 0],\
                            [ 0, 0, 0, 0, 0, 1, 3, 9, 27, 81]]";

        let mat = MatZ::from_str(mat_str).unwrap();
        assert_eq!(mat, gadget_mat);
    }
}

#[cfg(test)]
mod test_gen_trapdoor {
    use super::gen_trapdoor;
    use crate::sample::g_trapdoor::{
        gadget_classical::gen_gadget_mat, gadget_parameters::GadgetParameters,
    };
    use qfall_math::{
        integer::{MatZ, Z},
        integer_mod_q::{MatZq, Modulus},
        traits::{Concatenate, GetNumColumns, GetNumRows, SetEntry},
    };

    /// assure that the trapdoor `r` returned from [`gen_trapdoor`] is actually a
    /// trapdoor for `a`
    #[test]
    fn is_trapdoor_without_tag() {
        let modulus = Modulus::try_from(&Z::from(32)).unwrap();
        let params = GadgetParameters::init_default(42, &modulus);
        let a_bar = MatZq::sample_uniform(42, &params.m_bar, &params.q).unwrap();
        let tag = MatZq::identity(42, 42, &params.q).unwrap();

        // call gen_trapdoor to get matrix a and its 'trapdoor' r
        let (a, r) = gen_trapdoor(&params, &a_bar, &tag).unwrap();

        // generate the trapdoor for a from r as trapdoor = [[r],[I]]
        let trapdoor = r
            .concat_vertical(
                &MatZ::identity(a.get_num_columns() - r.get_num_rows(), r.get_num_columns())
                    .unwrap(),
            )
            .unwrap();

        // ensure G = A*trapdoor (definition of a trapdoor)
        let gadget_mat = gen_gadget_mat(&params.n, &params.k, &Z::from(2)).unwrap();
        assert_eq!(
            MatZq::from((&gadget_mat, &modulus)),
            a * MatZq::from((&trapdoor, &modulus))
        );
    }

    /// assure that the trapdoor `r` returned from [`gen_trapdoor`] is actually a
    /// trapdoor for `a`
    #[test]
    fn is_trapdoor_with_tag() {
        let modulus = Modulus::try_from(&Z::from(32)).unwrap();
        let params = GadgetParameters::init_default(42, &modulus);
        let a_bar = MatZq::sample_uniform(42, &params.m_bar, &params.q).unwrap();
        // calculate an invertible tag in Z_q^{n \times n}
        let tag = calculate_invertible_tag(42, &modulus);

        // call gen_trapdoor to get matrix a and its 'trapdoor' r
        let (a, r) = gen_trapdoor(&params, &a_bar, &tag).unwrap();

        // generate the trapdoor for a from r as trapdoor = [[r],[I]]
        let trapdoor = r
            .concat_vertical(
                &MatZ::identity(a.get_num_columns() - r.get_num_rows(), r.get_num_columns())
                    .unwrap(),
            )
            .unwrap();

        // ensure tag*G = A*trapdoor (definition of a trapdoor)
        let gadget_mat = gen_gadget_mat(&params.n, &params.k, &Z::from(2)).unwrap();
        assert_eq!(
            tag * MatZq::from((&gadget_mat, &modulus)),
            a * MatZq::from((&trapdoor, &modulus))
        );
    }

    /// Generates an invertible tag matrix (generates a diagonal matrix)
    fn calculate_invertible_tag(size: i64, modulus: &Modulus) -> MatZq {
        let max_value = Z::from(modulus);
        let mut out = MatZq::identity(size, size, modulus).unwrap();
        // create a diagonal matrix with random values (because it is a diagonal matrix
        // with `1` on the diagonal, it is always invertible)
        for row in 0..size {
            for column in 0..size {
                if row < column {
                    out.set_entry(row, column, Z::sample_uniform(&0, &max_value).unwrap())
                        .unwrap();
                }
            }
        }
        out
    }
}
