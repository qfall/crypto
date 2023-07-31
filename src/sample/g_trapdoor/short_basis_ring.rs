// Copyright © 2023 Marvin Beckmann
//
// This file is part of qFALL-crypto.
//
// qFALL-crypto is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

//! This module contains an implementation to generate a short basis from a ring-based
//! G-Trapdoor and its parity check matrix.

use super::{gadget_parameters::GadgetParametersRing, gadget_ring::find_solution_gadget_ring};
use qfall_math::{
    integer::{MatPolyOverZ, PolyOverZ, Z},
    integer_mod_q::{MatPolynomialRingZq, PolyOverZq, PolynomialRingZq, Zq},
    traits::{
        Concatenate, GetEntry, GetNumColumns, GetNumRows, Pow, SetCoefficient, SetEntry, Tensor,
    },
};

/// Generates a short basis according to [\[1\]](<../index.html#:~:text=[1]>).
/// Also refer to Lemma 5.3 in the eprint version <https://eprint.iacr.org/2011/501.pdf>.
/// Both interpreted in the ring-setting.
///
/// *Note*: At the moment this function does not support tags, this may be added later
/// (and the signature of this function may change).
///
/// The matrix is generated as `[ 1 | 0 | e,  0 | 1 | r, 0 | I ] * [ 0 | I_2, S' | W ]`
/// where `w` is a solution of `g^tw = -A [ I_2 | 0 ] mod q` and `S'` is a
/// reordering of `S` (if `base^k=q` then reversed, otherwise the same as before).
/// This corresponds to an appropriate reordering from
/// [\[1\]](<../index.html#:~:text=[1]>) and Lemma 3.2 from
/// [\[2\]](<../index.html#:~:text=[2]>) in the ring-setting.
///
/// Parameters:
/// - `params`: the gadget parameters with which the trapdoor was generated
/// - `a`: the parity check matrix
/// - `r`: the first part of the trapdoor for `a`
/// - `e`: the second part of the trapdoor for `a`
///
/// Returns a short basis for the lattice `\Lambda^\perp(a)` using the trapdoor `r,e`
///
/// # Examples
/// ```
/// use qfall_crypto::sample::g_trapdoor::gadget_parameters::GadgetParametersRing;
/// use qfall_crypto::sample::g_trapdoor::gadget_ring::gen_trapdoor_ring_lwe;
/// use qfall_crypto::sample::g_trapdoor::short_basis_ring::gen_short_basis_for_trapdoor_ring;
/// use qfall_math::{
///     integer::PolyOverZ,
///     integer_mod_q::{Modulus},
///     rational::Q,
///     traits::{GetNumColumns, GetNumRows},
/// };
///
/// let params = GadgetParametersRing::init_default(8, &Modulus::from(16));
/// let a_bar = PolyOverZ::sample_uniform(&params.n, 0, &params.q).unwrap();
///
/// let (a, r, e) = gen_trapdoor_ring_lwe(&params, &a_bar, &Q::from(5)).unwrap();
///
/// let short_base = gen_short_basis_for_trapdoor_ring(&params, &a, &r, &e);
/// ```
pub fn gen_short_basis_for_trapdoor_ring(
    params: &GadgetParametersRing,
    a: &MatPolynomialRingZq,
    r: &MatPolyOverZ,
    e: &MatPolyOverZ,
) -> MatPolyOverZ {
    let sa_l = gen_sa_l(e, r);
    let sa_r = gen_sa_r(params, a);
    let n = PolyOverZq::from(&params.modulus).get_degree();
    let mut poly_degrees = MatPolyOverZ::new(1, n);
    for i in 0..n {
        let mut x_i = PolyOverZ::default();
        x_i.set_coeff(i, 1).unwrap();
        poly_degrees.set_entry(0, i, x_i).unwrap();
    }
    let basis = poly_degrees.tensor_product(&(sa_l * sa_r));
    // the basis has to be reduced by the modulus to remove high-degrees
    MatPolynomialRingZq::from((&basis, &params.modulus)).get_mat()
}

/// Computes [ 1 | 0 | e,  0 | 1 | r, 0 | I ]
fn gen_sa_l(e: &MatPolyOverZ, r: &MatPolyOverZ) -> MatPolyOverZ {
    let out = e.concat_vertical(r).unwrap();

    let identity_lower_right = MatPolyOverZ::identity(out.get_num_columns(), out.get_num_columns());

    let out = out.concat_vertical(&identity_lower_right).unwrap();
    let identity_left = MatPolyOverZ::identity(out.get_num_rows(), 2);

    identity_left.concat_horizontal(&out).unwrap()
}

/// Computes `[0 | I , S' | W]`
fn gen_sa_r(params: &GadgetParametersRing, a: &MatPolynomialRingZq) -> MatPolyOverZ {
    let ident = MatPolyOverZ::identity(2, 2);
    let right = ident.concat_vertical(&compute_w(params, a)).unwrap();

    let mut s = compute_s(params);
    if params.base.pow(&params.k).unwrap() == Z::from(&params.q) {
        s.reverse_columns();
    }
    let zero = MatPolyOverZ::new(2, s.get_num_columns());
    let left = zero.concat_vertical(&s).unwrap();

    left.concat_horizontal(&right).unwrap()
}

/// Computes `w` with `g^tw = - a[I_2|0] mod qR`
fn compute_w(params: &GadgetParametersRing, a: &MatPolynomialRingZq) -> MatPolyOverZ {
    let minus_one = PolynomialRingZq::from((&PolyOverZ::from(-1), &params.modulus));
    let rhs_1: PolynomialRingZq = a.get_entry(0, 0).unwrap();
    let rhs_2: PolynomialRingZq = a.get_entry(0, 1).unwrap();
    let w_1 =
        find_solution_gadget_ring(&(&minus_one * &rhs_1), &params.k, &params.base).transpose();
    let w_2 =
        find_solution_gadget_ring(&(&minus_one * &rhs_2), &params.k, &params.base).transpose();
    w_1.concat_horizontal(&w_2).unwrap()
}

/// Computes a short basis for the gadget vector.
fn compute_s(params: &GadgetParametersRing) -> MatPolyOverZ {
    let id_k = MatPolyOverZ::identity(&params.k, &params.k);
    let mut sk = &params.base * id_k;
    for i in 0..(sk.get_num_rows() - 1) {
        sk.set_entry(i + 1, i, PolyOverZ::from(-1)).unwrap();
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
            sk.set_entry(i, sk.get_num_columns() - 1, PolyOverZ::from(&q_i))
                .unwrap();
            q = q - q_i;
            q = q.div_exact(&params.base).unwrap();
        }
        sk
    };
    sk
}

#[cfg(test)]
mod test_gen_short_basis_for_trapdoor_ring {
    use super::gen_short_basis_for_trapdoor_ring;
    use crate::sample::g_trapdoor::{
        gadget_parameters::GadgetParametersRing, gadget_ring::gen_trapdoor_ring_lwe,
    };
    use qfall_math::{
        integer::{PolyOverZ, Z},
        integer_mod_q::{MatPolynomialRingZq, Modulus},
        rational::Q,
        traits::{GetEntry, GetNumColumns},
    };

    /// Ensure that every vector within the returned basis is in `\Lambda^\perp(a)`.
    #[test]
    fn is_basis() {
        for n in [5, 10, 12] {
            let params =
                GadgetParametersRing::init_default(n, &Modulus::try_from(&Z::from(16)).unwrap());
            let a_bar = PolyOverZ::sample_uniform(&params.n, 0, &params.q).unwrap();

            let (a, r, e) = gen_trapdoor_ring_lwe(&params, &a_bar, &Q::from(5)).unwrap();

            let short_base = gen_short_basis_for_trapdoor_ring(&params, &a, &r, &e);
            let short_base = MatPolynomialRingZq::from((&short_base, &params.modulus));

            assert_eq!(n * a.get_num_columns(), short_base.get_num_columns());
            let res = a * short_base;
            for i in 0..res.get_num_columns() {
                let entry: PolyOverZ = res.get_entry(0, i).unwrap();
                assert!(entry.is_zero())
            }
        }
    }

    /// Ensure that all entries have a degree of at most n-1.
    #[test]
    fn basis_is_reduced() {
        for n in [5, 10, 12] {
            let params =
                GadgetParametersRing::init_default(n, &Modulus::try_from(&Z::from(16)).unwrap());
            let a_bar = PolyOverZ::sample_uniform(&params.n, 0, &params.q).unwrap();

            let (a, r, e) = gen_trapdoor_ring_lwe(&params, &a_bar, &Q::from(5)).unwrap();

            let short_base = gen_short_basis_for_trapdoor_ring(&params, &a, &r, &e);
            let short_base_reduced =
                MatPolynomialRingZq::from((&short_base, &params.modulus)).get_mat();

            assert_eq!(short_base_reduced, short_base)
        }
    }
}

#[cfg(test)]
mod test_gen_sa {
    use crate::sample::g_trapdoor::{
        gadget_parameters::GadgetParametersRing,
        short_basis_ring::{gen_sa_l, gen_sa_r},
    };
    use qfall_math::{
        integer::MatPolyOverZ,
        integer_mod_q::{MatPolynomialRingZq, Modulus},
    };
    use std::str::FromStr;

    /// Returns a fixed trapdoor and a matrix a for a fixed parameter set.
    fn get_fixed_trapdoor() -> (
        GadgetParametersRing,
        MatPolynomialRingZq,
        MatPolyOverZ,
        MatPolyOverZ,
    ) {
        let params = GadgetParametersRing::init_default(4, &Modulus::from(16));

        let a = MatPolyOverZ::from_str(
            "[[1  1, 4  2 8 8 12, 4  11 10 7 13, 4  9 6 6 12, 4  6 11 1 6, 4  3 10 2 9]]",
        )
        .unwrap();
        let a = MatPolynomialRingZq::from((&a, &params.modulus));

        let r = MatPolyOverZ::from_str("[[4  -1 7 6 -8, 3  0 -2 4, 4  0 3 -4 1, 4  6 4 -1 3]]")
            .unwrap();
        let e =
            MatPolyOverZ::from_str("[[4  -4 8 -3 7, 4  1 -2 2 4, 3  -6 7 -5, 4  -7 10 -12 -15]]")
                .unwrap();

        (params, a, r, e)
    }

    /// Ensure that the left part of the multiplication to get sa is correctly
    /// computed.
    #[test]
    fn working_sa_l() {
        let (_, _, r, e) = get_fixed_trapdoor();
        let sa_l = gen_sa_l(&r, &e);

        let sa_l_cmp = MatPolyOverZ::from_str(
            "[\
                [1  1, 0, 4  -1 7 6 -8, 3  0 -2 4, 4  0 3 -4 1, 4  6 4 -1 3],\
                [0, 1  1, 4  -4 8 -3 7, 4  1 -2 2 4, 3  -6 7 -5, 4  -7 10 -12 -15],\
                [0, 0, 1  1, 0, 0, 0],\
                [0, 0, 0, 1  1, 0, 0],\
                [0, 0, 0, 0, 1  1, 0],\
                [0, 0, 0, 0, 0, 1  1]]",
        )
        .unwrap();

        assert_eq!(sa_l_cmp, sa_l)
    }

    /// Ensure that the right part of the multiplication to get sa is correctly
    /// computed.
    #[test]
    fn working_sa_r() {
        let (params, a, _, _) = get_fixed_trapdoor();
        let sa_r = gen_sa_r(&params, &a);

        let sa_r_cmp = MatPolyOverZ::from_str(
            "[\
                [0, 0, 0, 0, 1  1, 0],\
                [0, 0, 0, 0, 0, 1  1],\
                [0, 0, 0, 1  2, 1  1, 0],\
                [0, 0, 1  2, 1  -1, 1  1, 1  1],\
                [0, 1  2, 1  -1, 0, 1  1, 4  1 0 0 1],\
                [1  2, 1  -1, 0, 0, 1  1, 3  1 1 1]]",
        )
        .unwrap();
        assert_eq!(sa_r_cmp, sa_r);
    }
}

#[cfg(test)]
mod test_compute_s {
    use qfall_math::{
        integer::{MatPolyOverZ, Z},
        integer_mod_q::Modulus,
    };
    use std::str::FromStr;

    use crate::sample::g_trapdoor::{
        gadget_parameters::GadgetParametersRing, short_basis_ring::compute_s,
    };

    /// Ensure that the matrix s is computed correctly for a power-of-two modulus.
    #[test]
    fn base_2_power_two() {
        let params = GadgetParametersRing::init_default(8, &Modulus::from(16));

        let s = compute_s(&params);

        let s_cmp = MatPolyOverZ::from_str(
            "[[1  2, 0, 0, 0],\
                [1  -1, 1  2, 0, 0],\
                [0, 1  -1, 1  2, 0],\
                [0, 0, 1  -1, 1  2]]",
        )
        .unwrap();

        assert_eq!(s_cmp, s)
    }

    /// Ensure that the matrix s is computed correctly for an arbitrary modulus.
    #[test]
    fn base_2_arbitrary() {
        let modulus = Z::from(0b1100110);
        let params = GadgetParametersRing::init_default(1, &Modulus::from(&modulus));

        let s = compute_s(&params);

        let s_cmp = MatPolyOverZ::from_str(
            "[[1  2, 0, 0, 0, 0, 0, 0],\
                [1  -1, 1  2, 0, 0, 0, 0, 1  1],\
                [0, 1  -1, 1  2, 0, 0, 0, 1  1],\
                [0, 0, 1  -1, 1  2, 0, 0, 0],\
                [0, 0, 0, 1  -1, 1  2, 0, 0],\
                [0, 0, 0, 0, 1  -1, 1  2, 1  1],\
                [0, 0, 0, 0, 0, 1  -1, 1  1]]",
        )
        .unwrap();

        assert_eq!(s_cmp, s)
    }

    /// Ensure that the matrix s is computed correctly for a power-of-5 modulus.
    #[test]
    fn base_5_power_5() {
        let mut params = GadgetParametersRing::init_default(1, &Modulus::from(625));
        params.k = Z::from(4);
        params.base = Z::from(5);

        let s = compute_s(&params);

        let s_cmp = MatPolyOverZ::from_str(
            "[[1  5, 0, 0, 0],\
                [1  -1, 1  5, 0, 0],\
                [0, 1  -1, 1  5, 0],\
                [0, 0, 1  -1, 1  5]]",
        )
        .unwrap();

        assert_eq!(s_cmp, s)
    }

    /// Ensure that the matrix s is computed correctly for an arbitrary modulus with
    /// base 5.
    #[test]
    fn base_5_arbitrary() {
        let modulus = Z::from_str_b("4123", 5).unwrap();
        let mut params = GadgetParametersRing::init_default(1, &Modulus::from(&modulus));
        params.k = Z::from(4);
        params.base = Z::from(5);

        let s = compute_s(&params);

        let s_cmp = MatPolyOverZ::from_str(
            "[[1  5, 0, 0, 1  3],\
                [1  -1, 1  5, 0, 1  2],\
                [0, 1  -1, 1  5, 1  1],\
                [0, 0, 1  -1, 1  4]]",
        )
        .unwrap();

        assert_eq!(s_cmp, s)
    }
}

#[cfg(test)]
mod test_compute_w {
    use crate::sample::g_trapdoor::{
        gadget_parameters::GadgetParametersRing,
        gadget_ring::{gen_gadget_ring, gen_trapdoor_ring_lwe},
        short_basis_ring::compute_w,
    };
    use qfall_math::{
        integer::{MatPolyOverZ, PolyOverZ},
        integer_mod_q::{MatPolynomialRingZq, Modulus},
        rational::Q,
        traits::GetNumColumns,
    };

    /// Ensure that `gw = a[I_1|0] mod qR`.
    #[test]
    fn check_w_is_correct_solution() {
        let params = GadgetParametersRing::init_default(8, &Modulus::from(16));
        let a_bar = PolyOverZ::sample_uniform(&params.n, 0, &params.q).unwrap();

        let (a, _, _) = gen_trapdoor_ring_lwe(&params, &a_bar, &Q::from(5)).unwrap();

        let w = compute_w(&params, &a);
        let w = MatPolynomialRingZq::from((&w, &params.modulus));

        let gadget = gen_gadget_ring(&params.k, &params.base).unwrap();
        let gadget = MatPolynomialRingZq::from((&gadget, &params.modulus));

        let gw = gadget.transpose() * w;
        let i0 = -1 * MatPolyOverZ::identity(a.get_num_columns(), 2);
        let i0 = MatPolynomialRingZq::from((&i0, &params.modulus));
        let rhs = a * i0;

        assert_eq!(gw, rhs)
    }
}
