// Copyright © 2023 Marvin Beckmann
//
// This file is part of qFALL-crypto.
//
// qFALL-crypto is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

//! This module contains an implementation to generate a short basis from a G-Trapdoor
//! and its parity check matrix.

use super::{gadget_classical::find_solution_gadget_mat, gadget_parameters::GadgetParameters};
use qfall_math::{
    integer::{MatZ, Z},
    integer_mod_q::{MatZq, Zq},
    traits::{Concatenate, GetNumColumns, GetNumRows, Pow, SetEntry, Tensor},
};

/// Generates a short basis according to [\[1\]](<../index.html#:~:text=[1]>).
/// Also refer to Lemma 5.3 in the eprint version <https://eprint.iacr.org/2011/501.pdf>.
///
/// *Note*: At the moment this function does not support tags other than the identity.
/// The function will panic right now, if anything else was provided.
///
/// The matrix is generated as `[ I | R, 0 | I ] * [ I | 0, W | S ]`
/// where `W` is a solution of `GW = -H^{-1}A [ I | 0 ] mod q`
///
/// Parameters:
/// - `params`: the gadget parameters with which the trapdoor was generated
/// - `tag`: the corresponding tag
/// - `a`: the parity check matrix
/// - `r`: the trapdoor for `a`
///
/// Returns a short basis for the lattice `\Lambda^\perp(a)` using the trapdoor `r`
///
/// # Examples
/// ```
/// use qfall_crypto::sample::g_trapdoor::{gadget_parameters::GadgetParameters,
/// gen_trapdoor_default};
/// use qfall_crypto::sample::g_trapdoor::short_basis_classical::gen_short_basis_for_trapdoor;
/// use qfall_math::{
///     integer::Z,
///     integer_mod_q::{MatZq, Modulus},
///     traits::{GetNumColumns, GetNumRows},
/// };
///
/// let modulus = Modulus::from(127);
/// let params = GadgetParameters::init_default(10, &modulus);
/// let (a, r) = gen_trapdoor_default(&params.n, &modulus);
///
/// let tag = MatZq::identity(&params.n, &params.n, &modulus);
///
/// let short_basis = gen_short_basis_for_trapdoor(&params, &tag, &a, &r);
/// ```
pub fn gen_short_basis_for_trapdoor(
    params: &GadgetParameters,
    tag: &MatZq,
    a: &MatZq,
    r: &MatZ,
) -> MatZ {
    let sa_l = gen_sa_l(r);
    let sa_r = gen_sa_r(params, tag, a);
    sa_l * sa_r
}

/// Computes [ I | R, 0 | I ]
fn gen_sa_l(r: &MatZ) -> MatZ {
    let left = MatZ::identity(r.get_num_rows() + r.get_num_columns(), r.get_num_rows());
    let identity_right_lower = MatZ::identity(r.get_num_columns(), r.get_num_columns());
    let right = r.concat_vertical(&identity_right_lower).unwrap();
    left.concat_horizontal(&right).unwrap()
}

/// Computes `[ I | 0, W | S ]`
fn gen_sa_r(params: &GadgetParameters, tag: &MatZq, a: &MatZq) -> MatZ {
    let s = compute_s(params);
    let w = compute_w(params, tag, a);

    let identity_upper = MatZ::identity(
        w.get_num_columns(),
        w.get_num_columns() + s.get_num_columns(),
    );

    let ws = w.concat_horizontal(&s).unwrap();
    identity_upper.concat_vertical(&ws).unwrap()
}

/// Compute S for `[ I | 0, W | S ]`
fn compute_s(params: &GadgetParameters) -> MatZ {
    let id_k = MatZ::identity(&params.k, &params.k);
    let mut sk = &params.base * id_k;
    for i in 0..(sk.get_num_rows() - 1) {
        sk.set_entry(i + 1, i, Z::MINUS_ONE).unwrap();
    }
    sk = if params.base.pow(&params.k).unwrap() == Z::from(&params.q) {
        // compute s in the special case where the modulus is a power of base
        // i.e. the last column can remain as it is
        sk
    } else {
        // compute s for any arbitrary modulus
        // represent modulus in `base` and set last row accordingly
        let mut q = Z::from(&params.q);
        for i in 0..(sk.get_num_rows()) {
            let q_i = Zq::from((&q, &params.base)).get_value();
            sk.set_entry(i, sk.get_num_columns() - 1, &q_i).unwrap();
            q = q - q_i;
            q = q.div_exact(&params.base).unwrap();
        }
        sk
    };
    MatZ::identity(&params.n, &params.n).tensor_product(&sk)
}

/// Computes `W` with `GW = -H^{-1}A [ I | 0 ] mod q`
fn compute_w(params: &GadgetParameters, tag: &MatZq, a: &MatZq) -> MatZ {
    // TODO invert tag, remove this, once we can invert tags.
    let identity = MatZq::identity(tag.get_num_rows(), tag.get_num_columns(), tag.get_mod());
    assert_eq!(&identity, tag);
    let tag_inv = tag;

    let rhs = Z::MINUS_ONE * tag_inv * (a * MatZ::identity(a.get_num_columns(), &params.m_bar));
    find_solution_gadget_mat(&rhs, &params.k, &params.base)
}

#[cfg(test)]
mod test_gen_short_basis_for_trapdoor {
    use super::gen_short_basis_for_trapdoor;
    use crate::sample::g_trapdoor::{gadget_parameters::GadgetParameters, gen_trapdoor_default};
    use qfall_math::{
        integer_mod_q::{MatZq, Modulus},
        traits::{GetNumColumns, GetNumRows},
    };

    /// ensure that every vector within the returned basis is in `\Lambda^\perp(A)`
    #[test]
    fn is_basis_not_power_tag_identity() {
        for n in [1, 5, 10, 12] {
            let modulus = Modulus::from(127 + 3 * n);
            let params = GadgetParameters::init_default(n, &modulus);
            let (a, r) = gen_trapdoor_default(&params.n, &modulus);

            let tag = MatZq::identity(&params.n, &params.n, &modulus);

            let short_basis = gen_short_basis_for_trapdoor(&params, &tag, &a, &r);

            let zero_vec = MatZq::new(a.get_num_rows(), 1, &modulus);

            for i in 0..short_basis.get_num_columns() {
                assert_eq!(zero_vec, &a * short_basis.get_column(i).unwrap())
            }
        }
    }
}

#[cfg(test)]
mod test_gen_sa {
    use super::gen_sa_l;
    use crate::sample::g_trapdoor::{
        gadget_parameters::GadgetParameters, short_basis_classical::gen_sa_r,
    };
    use qfall_math::{
        integer::MatZ,
        integer_mod_q::{MatZq, Modulus},
    };
    use std::str::FromStr;

    /// Returns a fixed trapdoor and a matrix a for a fixed parameter set
    fn get_fixed_trapdoor_for_tag_identity() -> (GadgetParameters, MatZq, MatZ) {
        let params = GadgetParameters::init_default(2, &Modulus::from(8));

        let a = MatZq::from_str(
            "[\
            [2, 6, 2, 5, 3, 0, 1, 1, 1, 6, 5, 0, 6],\
            [6, 0, 3, 1, 5, 6, 2, 7, 0, 3, 7, 7, 0]] mod 8",
        )
        .unwrap();

        let r = MatZ::from_str(
            "[[0, 1, 0, 1, 1, 0],\
            [-1, 1, 0, 0, 0, -1],\
            [-1, 0, -1, -1, -1, 0],\
            [-1, 1, 0, 0, 0, 1],\
            [-1, -1, 0, 1, 0, 1],\
            [-1, 0, 0, -1, 0, 1],\
            [0, -1, 0, 0, 0, 0]]",
        )
        .unwrap();

        (params, a, r)
    }

    /// Ensure that the left part of the multiplication to get sa is correctly computed
    #[test]
    fn working_sa_l() {
        let (_, _, r) = get_fixed_trapdoor_for_tag_identity();
        let sa_1 = gen_sa_l(&r);

        let sa_1_cmp = MatZ::from_str(
            "[\
            [1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0],\
            [0, 1, 0, 0, 0, 0, 0, -1, 1, 0, 0, 0, -1],\
            [0, 0, 1, 0, 0, 0, 0, -1, 0, -1, -1, -1, 0],\
            [0, 0, 0, 1, 0, 0, 0, -1, 1, 0, 0, 0, 1],\
            [0, 0, 0, 0, 1, 0, 0, -1, -1, 0, 1, 0, 1],\
            [0, 0, 0, 0, 0, 1, 0, -1, 0, 0, -1, 0, 1],\
            [0, 0, 0, 0, 0, 0, 1, 0, -1, 0, 0, 0, 0],\
            [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],\
            [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],\
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],\
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],\
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],\
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]]",
        )
        .unwrap();
        assert_eq!(sa_1_cmp, sa_1);
    }

    /// Ensure that the right part of the multiplication to get sa is correctly computed
    /// with tag as identity
    #[test]
    fn working_sa_r_identity() {
        let (params, a, _) = get_fixed_trapdoor_for_tag_identity();
        println!("{0}", params.k);
        let tag = MatZq::identity(&params.n, &params.n, &params.q);
        let sa_r = gen_sa_r(&params, &tag, &a);

        println!("{sa_r}");
        let sa_r_cmp = MatZ::from_str(
            "[\
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],\
            [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],\
            [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],\
            [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],\
            [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],\
            [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],\
            [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],\
            [0, 0, 0, 1, 1, 0, 1, 2, 0, 0, 0, 0, 0],\
            [1, 1, 1, 1, 0, 0, 1, -1, 2, 0, 0, 0, 0],\
            [1, 0, 1, 0, 1, 0, 1, 0, -1, 2, 0, 0, 0],\
            [0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 2, 0, 0],\
            [1, 0, 0, 1, 1, 1, 1, 0, 0, 0, -1, 2, 0],\
            [0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, -1, 2]]",
        )
        .unwrap();
        assert_eq!(sa_r_cmp, sa_r);
    }
}

#[cfg(test)]
mod test_compute_s {
    use crate::sample::g_trapdoor::{
        gadget_parameters::GadgetParameters, short_basis_classical::compute_s,
    };
    use qfall_math::{
        integer::{MatZ, Z},
        integer_mod_q::Modulus,
    };
    use std::str::FromStr;

    /// Ensure that the matrix s is computed correctly for a power-of-two modulus
    #[test]
    fn base_2_power_two() {
        let params = GadgetParameters::init_default(2, &Modulus::from(16));

        let s = compute_s(&params);

        let s_cmp = MatZ::from_str(
            "[[2, 0, 0, 0, 0, 0, 0, 0],\
            [-1, 2, 0, 0, 0, 0, 0, 0],\
            [0, -1, 2, 0, 0, 0, 0, 0],\
            [0, 0, -1, 2, 0, 0, 0, 0],\
            [0, 0, 0, 0, 2, 0, 0, 0],\
            [0, 0, 0, 0, -1, 2, 0, 0],\
            [0, 0, 0, 0, 0, -1, 2, 0],\
            [0, 0, 0, 0, 0, 0, -1, 2]]",
        )
        .unwrap();
        assert_eq!(s_cmp, s)
    }

    /// Ensure that the matrix s is computed correctly for a power-of-two modulus
    #[test]
    fn base_2_arbitrary() {
        let modulus = Z::from(0b1100110);
        let params = GadgetParameters::init_default(1, &Modulus::from(&modulus));

        let s = compute_s(&params);

        println!("{s}");
        let s_cmp = MatZ::from_str(
            "[[2, 0, 0, 0, 0, 0, 0],\
            [-1, 2, 0, 0, 0, 0, 1],\
            [0, -1, 2, 0, 0, 0, 1],\
            [0, 0, -1, 2, 0, 0, 0],\
            [0, 0, 0, -1, 2, 0, 0],\
            [0, 0, 0, 0, -1, 2, 1],\
            [0, 0, 0, 0, 0, -1, 1]]",
        )
        .unwrap();
        println!("{}", Z::from(0b1100110));

        assert_eq!(s_cmp, s)
    }

    /// Ensure that the matrix s is computed correctly for a power-of-two modulus
    #[test]
    fn base_5_power_5() {
        let mut params = GadgetParameters::init_default(1, &Modulus::from(625));
        params.k = Z::from(4);
        params.base = Z::from(5);

        let s = compute_s(&params);

        let s_cmp = MatZ::from_str(
            "[[5, 0, 0, 0],\
            [-1, 5, 0, 0],\
            [0, -1, 5, 0],\
            [0, 0, -1, 5]]",
        )
        .unwrap();
        assert_eq!(s_cmp, s)
    }

    /// Ensure that the matrix s is computed correctly for a power-of-two modulus
    #[test]
    fn base_5_arbitrary() {
        let modulus = Z::from_str_b("4123", 5).unwrap();
        let mut params = GadgetParameters::init_default(1, &Modulus::from(&modulus));
        params.k = Z::from(4);
        params.base = Z::from(5);

        let s = compute_s(&params);

        let s_cmp = MatZ::from_str(
            "[[5, 0, 0, 3],\
            [-1, 5, 0, 2],\
            [0, -1, 5, 1],\
            [0, 0, -1, 4]]",
        )
        .unwrap();

        assert_eq!(s_cmp, s)
    }
}

#[cfg(test)]
mod test_compute_w {
    use super::compute_w;
    use crate::sample::g_trapdoor::{
        gadget_classical::gen_gadget_mat, gadget_parameters::GadgetParameters,
    };
    use qfall_math::{
        integer::{MatZ, Z},
        integer_mod_q::{MatZq, Modulus},
        traits::GetNumColumns,
    };
    use std::str::FromStr;

    /// Ensure that `GW = A[I|0] mod q`
    #[test]
    fn working_example_tag_identity() {
        let params = GadgetParameters::init_default(2, &Modulus::from(8));
        let tag = MatZq::identity(2, 2, &params.q);

        let a = MatZq::from_str(
            "[\
            [2, 6, 2, 5, 3, 0, 1, 1, 1, 6, 5, 0, 6],\
            [6, 0, 3, 1, 5, 6, 2, 7, 0, 3, 7, 7, 0]] mod 8",
        )
        .unwrap();

        let w = compute_w(&params, &tag, &a);
        println!("{w}");
        let g = gen_gadget_mat(&params.n, &params.k, &params.base).unwrap();

        let gw = MatZq::from((&(g * w), &params.q));
        let rhs = &a * MatZ::identity(a.get_num_columns(), &params.m_bar);

        assert_eq!(gw, Z::MINUS_ONE * rhs)
    }
}
