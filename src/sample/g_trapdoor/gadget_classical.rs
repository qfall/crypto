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
    integer_mod_q::{MatZq, Zq},
    traits::{Concatenate, GetEntry, GetNumColumns, GetNumRows, Pow, SetEntry, Tensor},
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
/// - `tag`: the tag which is hidden within the matrix `A`
///
/// Returns a a parity-check matrix `a` derived from `a_bar` and its gadget-trapdoor `r`
/// under a give tag `h`.
///
/// # Examples
/// ```
/// use qfall_crypto::sample::g_trapdoor::{gadget_parameters::GadgetParameters, gadget_classical::gen_trapdoor};
/// use qfall_math::integer_mod_q::MatZq;
///
/// let params = GadgetParameters::init_default(42, 42);
/// let a_bar = MatZq::sample_uniform(42, &params.m_bar, &params.q);
/// let tag = MatZq::identity(42, 42, &params.q);
///
/// let (a,r) = gen_trapdoor(&params, &a_bar, &tag).unwrap();
/// ```
///
/// # Errors and Failures
/// - Returns a [`MathError`] of type [`MismatchingMatrixDimension`](MathError::MismatchingMatrixDimension)
/// if the matrices can not be concatenated due to mismatching dimensions.
/// - Returns a [`MathError`] of type [`MismatchingModulus`](MathError::MismatchingModulus)
/// if the matrices can not be concatenated due to mismatching moduli.
///
/// # Panics ...
/// - if `params.k < 1` or it does not fit into an [`i64`].
/// - if `params.n < 1`.
pub fn gen_trapdoor(
    params: &GadgetParameters,
    a_bar: &MatZq,
    tag: &MatZq,
) -> Result<(MatZq, MatZ), MathError> {
    let g = gen_gadget_mat(&params.n, &params.k, &params.base);
    let r = params
        .distribution
        .sample(&params.m_bar, &(&params.n * &params.k));
    // set A = [A_bar | HG - A_bar R]
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
/// let g = gen_gadget_mat(4, 4, &Z::from(2));
/// ```
///
/// # Panics ...
/// - if `k < 1` or it does not fit into an [`i64`].
/// - if `n < 1`.
pub fn gen_gadget_mat(
    n: impl TryInto<i64> + Display + Clone,
    k: impl TryInto<i64> + Display,
    base: &Z,
) -> MatZ {
    let gadget_vec = gen_gadget_vec(k, base);
    let identity = MatZ::identity(n.clone(), n);
    identity.tensor_product(&gadget_vec.transpose())
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
/// # Panics ...
/// - if `k < 1` or it does not fit into an [`i64`].
pub fn gen_gadget_vec(k: impl TryInto<i64> + Display, base: &Z) -> MatZ {
    let mut out = MatZ::new(k, 1);
    for i in 0..out.get_num_rows() {
        out.set_entry(i, 0, &base.pow(i).unwrap()).unwrap();
    }
    out
}

/// Computes an arbitrary solution for `g^t x = value mod q`.
///
/// Parameters:
/// - `value`: the matrix for which a solution has to be computed
/// - `k`: the length of a gadget vector
/// - `base`: the base with which the gadget vector is defined
///
/// Returns an arbitrary solution for `g^tx = value mod q`
///
/// # Examples
/// ```
/// use qfall_math::integer::Z;
/// use qfall_math::integer_mod_q::Zq;
/// use qfall_math::integer::MatZ;
/// use qfall_crypto::sample::g_trapdoor::gadget_classical::find_solution_gadget_vec;
/// use qfall_crypto::sample::g_trapdoor::gadget_classical::gen_gadget_vec;
/// use qfall_math::traits::GetEntry;
/// use std::str::FromStr;
///
/// let k = Z::from(5);
/// let base = Z::from(3);
/// let value = Zq::from((29,125));
///
/// let sol = find_solution_gadget_vec(&value, &k, &base);
///
/// assert_eq!(
///     value.get_value(),
///     (gen_gadget_vec(&k, &base).transpose() * sol)
///         .get_entry(0, 0)
///         .unwrap()
/// )
/// ```
///
/// # Panics ...
/// - if the modulus of the value is greater than `base^k`.
pub fn find_solution_gadget_vec(value: &Zq, k: &Z, base: &Z) -> MatZ {
    if base.pow(k).unwrap() < Z::from(&value.get_mod()) {
        panic!("The modulus is too large, the value is potentially not representable.");
    }

    let mut value = value.get_value();
    let mut out = MatZ::new(k, 1);
    for i in 0..out.get_num_rows() {
        let val_i = value.modulo(base);
        out.set_entry(i, 0, &val_i).unwrap();
        value = (value - val_i).div_exact(base).unwrap();
    }
    out
}

/// Computes an arbitrary solution for `GX = value mod q`.
///
/// Computes a entrywise solution using the structure of the gadget matrix to its
/// advantage and utilizing `find_solution_gadget_vec`.
///
/// Parameters:
/// - `value`: the matrix for which a solution has to be computed
/// - `k`: the length of a gadget vector
/// - `base`: the base with which the gadget vector is defined
///
/// Returns an arbitrary solution for `GX = value mod q`.
///
/// # Examples
/// ```
/// use qfall_math::integer::Z;
/// use qfall_math::integer::MatZ;
/// use qfall_math::integer_mod_q::MatZq;
/// use qfall_crypto::sample::g_trapdoor::gadget_classical::find_solution_gadget_mat;
/// use qfall_crypto::sample::g_trapdoor::gadget_classical::gen_gadget_mat;
/// use std::str::FromStr;
///
/// let k = Z::from(5);
/// let base = Z::from(3);
/// let value = MatZq::from_str("[[1, 42],[2, 30],[3, 12]] mod 125").unwrap();
///
/// let sol = find_solution_gadget_mat(&value, &k, &base);
///
/// assert_eq!(
///     MatZ::from(&value),
///     gen_gadget_mat(3, &k, &base) * sol
/// )
/// ```
///
/// # Panics ...
/// - if the modulus of the value is greater than `base^k`.
pub fn find_solution_gadget_mat(value: &MatZq, k: &Z, base: &Z) -> MatZ {
    let mut out = MatZ::new(k * value.get_num_rows(), value.get_num_columns());
    for i in 0..value.get_num_columns() as usize {
        let mut _out: MatZ = find_solution_gadget_vec(&value.get_entry(0, i).unwrap(), k, base);
        for j in 1..value.get_num_rows() as usize {
            let sol_j = find_solution_gadget_vec(&value.get_entry(j, i).unwrap(), k, base);
            _out = _out.concat_vertical(&sol_j).unwrap();
        }
        out.set_column(i, &_out, 0).unwrap();
    }
    out
}

#[cfg(test)]
mod test_gen_gadget_vec {
    use crate::sample::g_trapdoor::gadget_classical::gen_gadget_vec;
    use qfall_math::integer::{MatZ, Z};
    use std::str::FromStr;

    /// Assure that the gadget vector with base `2` and length `5` works correctly.
    #[test]
    fn correctness_base_2() {
        let gadget_vec = gen_gadget_vec(5, &Z::from(2));

        let vec = MatZ::from_str("[[1],[2],[4],[8],[16]]").unwrap();
        assert_eq!(vec, gadget_vec);
    }

    /// Assure that the gadget vector with base `5` and length `4` works correctly.
    #[test]
    fn correctness_base_5() {
        let gadget_vec = gen_gadget_vec(4, &Z::from(5));

        let vec = MatZ::from_str("[[1],[5],[25],[125]]").unwrap();
        assert_eq!(vec, gadget_vec);
    }
}

#[cfg(test)]
mod test_gen_gadget_mat {
    use super::gen_gadget_mat;
    use qfall_math::integer::{MatZ, Z};
    use std::str::FromStr;

    /// Assure that the gadget matrix with gadget vector `[1, 2, 4]^t`(base 3) and
    /// `I_3` works correctly.
    #[test]
    fn correctness_base_2_3x3() {
        let gadget_mat = gen_gadget_mat(3, 3, &Z::from(2));

        let mat_str = "[[1, 2, 4, 0, 0, 0, 0, 0, 0],\
                            [0, 0, 0, 1, 2, 4, 0, 0, 0],\
                            [0, 0, 0, 0, 0, 0, 1, 2, 4]]";

        let mat = MatZ::from_str(mat_str).unwrap();
        assert_eq!(mat, gadget_mat);
    }

    /// Assure that the gadget matrix with gadget vector `[1, 3, 9, 27, 81]^t`(base 3) and
    /// `I_2` works correctly.
    #[test]
    fn correctness_base_3_2x5() {
        let gadget_mat = gen_gadget_mat(2, 5, &Z::from(3));

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

    /// Assure that the trapdoor `r` returned from [`gen_trapdoor`] is actually a
    /// trapdoor for `a`.
    #[test]
    fn is_trapdoor_without_tag() {
        let params = GadgetParameters::init_default(42, 32);
        let a_bar = MatZq::sample_uniform(42, &params.m_bar, &params.q);
        let tag = MatZq::identity(42, 42, &params.q);

        // call gen_trapdoor to get matrix a and its 'trapdoor' r
        let (a, r) = gen_trapdoor(&params, &a_bar, &tag).unwrap();

        // generate the trapdoor for a from r as trapdoor = [[r],[I]]
        let trapdoor = r
            .concat_vertical(&MatZ::identity(
                a.get_num_columns() - r.get_num_rows(),
                r.get_num_columns(),
            ))
            .unwrap();

        // ensure G = A*trapdoor (definition of a trapdoor)
        let gadget_mat = gen_gadget_mat(&params.n, &params.k, &Z::from(2));
        assert_eq!(
            MatZq::from((&gadget_mat, &params.q)),
            a * MatZq::from((&trapdoor, &params.q))
        );
    }

    /// Assure that the trapdoor `r` returned from [`gen_trapdoor`] is actually a
    /// trapdoor for `a`.
    #[test]
    fn is_trapdoor_with_tag() {
        let modulus = Modulus::from(32);
        let params = GadgetParameters::init_default(42, &modulus);
        let a_bar = MatZq::sample_uniform(42, &params.m_bar, &params.q);
        // calculate an invertible tag in Z_q^{n × n}
        let tag = calculate_invertible_tag(42, &modulus);

        // call gen_trapdoor to get matrix a and its 'trapdoor' r
        let (a, r) = gen_trapdoor(&params, &a_bar, &tag).unwrap();

        // generate the trapdoor for a from r as trapdoor = [[r],[I]]
        let trapdoor = r
            .concat_vertical(&MatZ::identity(
                a.get_num_columns() - r.get_num_rows(),
                r.get_num_columns(),
            ))
            .unwrap();

        // ensure tag*G = A*trapdoor (definition of a trapdoor)
        let gadget_mat = gen_gadget_mat(&params.n, &params.k, &Z::from(2));
        assert_eq!(
            tag * MatZq::from((&gadget_mat, &modulus)),
            a * MatZq::from((&trapdoor, &modulus))
        );
    }

    /// Generates an invertible tag matrix (generates a diagonal matrix)
    fn calculate_invertible_tag(size: i64, modulus: &Modulus) -> MatZq {
        let max_value = Z::from(modulus);
        let mut out = MatZq::identity(size, size, modulus);
        // create a diagonal matrix with random values (because it is a diagonal matrix
        // with `1` on the diagonal, it is always invertible)
        for row in 0..size {
            for column in 0..size {
                if row < column {
                    out.set_entry(row, column, Z::sample_uniform(0, &max_value).unwrap())
                        .unwrap();
                }
            }
        }
        out
    }
}

#[cfg(test)]
mod test_find_solution_gadget {
    use super::find_solution_gadget_vec;
    use crate::sample::g_trapdoor::gadget_classical::{
        find_solution_gadget_mat, gen_gadget_mat, gen_gadget_vec,
    };
    use qfall_math::{
        integer::{MatZ, Z},
        integer_mod_q::{MatZq, Zq},
        traits::GetEntry,
    };
    use std::str::FromStr;

    /// Ensure that the found solution is actually correct.
    #[test]
    fn returns_correct_solution_vec() {
        let k = Z::from(5);
        let base = Z::from(3);
        for i in 0..124 {
            let value = Zq::from((i, 125));

            let sol = find_solution_gadget_vec(&value, &k, &base);

            assert_eq!(
                value.get_value(),
                (gen_gadget_vec(&k, &base).transpose() * sol)
                    .get_entry(0, 0)
                    .unwrap()
            )
        }
    }

    /// Ensure that the found solution is actually correct.
    #[test]
    fn returns_correct_solution_mat() {
        let k = Z::from(5);
        let base = Z::from(3);
        let value = MatZq::from_str("[[1, 42],[2, 40],[3, 90]] mod 125").unwrap();

        let sol = find_solution_gadget_mat(&value, &k, &base);

        println!("{sol}");
        println!("{}", gen_gadget_mat(3, &k, &base));

        assert_eq!(MatZ::from(&value), gen_gadget_mat(3, &k, &base) * sol)
    }
}
