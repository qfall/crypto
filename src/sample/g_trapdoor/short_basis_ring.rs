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
/// Both have to be adopted to the ring-setting here.
///
/// The matrix is generated as:
/// `[ 1 | 0 | e,  0 | 1 | r, 0_{2xn} | I_{kxk} ] * [ 0 | I_2', S' | W' ]`
/// with
/// - `W' := [X^0 | X^1 | ... | X^{n-1}] \otimes W`,
/// - `I_2' := [X^0 | X^1 | ... | X^{n-1}] \otimes I_2` and
/// - `S':= [X^0 | X^1 | ... | X^{n-1}] \otimes S`.
/// Here `W` is a solution of `g^tW = -A [ I_2 | 0 ] mod q`,
/// `S` is a reordered (if `base^k=q` then reversed, otherwise the same as before)
/// short base of `\Lambda^\perp(g^t)`, i.e. `S''` is a reordered short base of `g^t`
/// in the classical case and `S':= [X^0 | X^1 | ... | X^{n-1}] \otimes S''`.
///
/// The appropriate reordering comes from
/// [\[1\]](<../index.html#:~:text=[1]>) and Lemma 3.2 from
/// [\[2\]](<../index.html#:~:text=[2]>).
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
/// };
///
/// let params = GadgetParametersRing::init_default(8, &Modulus::from(16));
/// let a_bar = PolyOverZ::sample_uniform(&params.n, 0, &params.q).unwrap();
///
/// let (a, r, e) = gen_trapdoor_ring_lwe(&params, &a_bar, &Q::from(5)).unwrap();
///
/// let short_base = gen_short_basis_for_trapdoor_ring(&params, &a, &r, &e);
/// ```
///
/// Panics ...
/// - if `e` and `r` are not both of dimensions `1 x k`.
/// - if `a` is not of dimension 1 x 1.
pub fn gen_short_basis_for_trapdoor_ring(
    params: &GadgetParametersRing,
    a: &MatPolynomialRingZq,
    r: &MatPolyOverZ,
    e: &MatPolyOverZ,
) -> MatPolyOverZ {
    let sa_l = gen_sa_l(e, r);
    let sa_r = gen_sa_r(params, a);
    let mut basis = sa_l * sa_r;
    // the basis has to be reduced by the modulus to remove high-degrees
    let ctx_poly = PolyOverZ::from(&PolyOverZq::from(&params.modulus));
    basis.reduce_by_poly(&ctx_poly);
    basis
}

/// Computes [ 1 | 0 | e,  0 | 1 | r, 0_{2xk} | I_{kxk} ]
fn gen_sa_l(e: &MatPolyOverZ, r: &MatPolyOverZ) -> MatPolyOverZ {
    let out = e.concat_vertical(r).unwrap();

    let identity_lower_right = MatPolyOverZ::identity(out.get_num_columns(), out.get_num_columns());

    let out = out.concat_vertical(&identity_lower_right).unwrap();
    let identity_left = MatPolyOverZ::identity(out.get_num_rows(), 2);

    identity_left.concat_horizontal(&out).unwrap()
}

/// Computes `pd \tensor [0_{2xk}, S''] || pd \tensor [I_{2x2}, w]` where
/// `pd := [X^0 | X^1 | ... | X^{n-1}]`.
/// Finally, the sa_r must have `n*m = n*(k+2)` columns.
fn gen_sa_r(params: &GadgetParametersRing, a: &MatPolynomialRingZq) -> MatPolyOverZ {
    let n = params.modulus.get_degree();
    let mut poly_degrees = MatPolyOverZ::new(1, n);
    for i in 0..n {
        let mut x_i = PolyOverZ::default();
        x_i.set_coeff(i, 1).unwrap();
        poly_degrees.set_entry(0, i, x_i).unwrap();
    }

    // compute a short base for `\Lambda^\perp(g^t)` in the classical but interpreted in
    // the ring and by applying the tensor product lift it to a short base for
    // `\Lambda^\perp(g^t)` in the ring.
    let mut s = compute_s(params);
    if params.base.pow(&params.k).unwrap() == Z::from(&params.q) {
        s.reverse_columns();
    }
    let s = poly_degrees.tensor_product(&s);
    let zero = MatPolyOverZ::new(2, &params.k * n);
    let left = zero.concat_vertical(&s).unwrap();

    // compute a solution for `g^tw = - a[I_2|0] mod qR`, but as all `w_i := X^i*w` are
    // also valid solution, we use the tensor product to have a possible solutions for
    // each power.
    let w = compute_w(params, a);
    let ident = MatPolyOverZ::identity(2, 2);
    let right = poly_degrees.tensor_product(&ident.concat_vertical(&w).unwrap());

    left.concat_horizontal(&right).unwrap()
}

/// Computes `w` with `g^tw = - a[I_2|0] mod qR` This is equivalent to finding solutions
/// for `g^tw_0 = -a_0 mod qR` and `g^tw_1 = -a_1 mod qR` and concatenating them after.
fn compute_w(params: &GadgetParametersRing, a: &MatPolynomialRingZq) -> MatPolyOverZ {
    let minus_one = PolynomialRingZq::from((&PolyOverZ::from(-1), &params.modulus));
    let rhs_0: PolynomialRingZq = a.get_entry(0, 0).unwrap();
    let rhs_1: PolynomialRingZq = a.get_entry(0, 1).unwrap();

    let w_0 =
        find_solution_gadget_ring(&(&minus_one * &rhs_0), &params.k, &params.base).transpose();
    let w_1 =
        find_solution_gadget_ring(&(&minus_one * &rhs_1), &params.k, &params.base).transpose();

    w_0.concat_horizontal(&w_1).unwrap()
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
        integer::PolyOverZ,
        integer_mod_q::{MatPolynomialRingZq, Modulus},
        rational::{MatQ, Q},
        traits::{GetEntry, GetNumColumns, GetNumRows, Pow},
    };

    /// Ensure that every vector within the returned basis is in `\Lambda^\perp(a)`.
    #[test]
    fn is_basis() {
        for n in [5, 10, 12] {
            let params = GadgetParametersRing::init_default(n, &Modulus::from(16));
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
            let params = GadgetParametersRing::init_default(n, &Modulus::from(16));
            let a_bar = PolyOverZ::sample_uniform(&params.n, 0, &params.q).unwrap();

            let (a, r, e) = gen_trapdoor_ring_lwe(&params, &a_bar, &Q::from(5)).unwrap();

            let short_base = gen_short_basis_for_trapdoor_ring(&params, &a, &r, &e);
            for i in 0..short_base.get_num_rows() {
                for j in 0..short_base.get_num_columns() {
                    let entry = short_base.get_entry(i, j).unwrap();
                    assert!(entry.get_degree() < n)
                }
            }
        }
    }

    /// Ensure that the orthogonalized short base length is upper bounded by
    /// `(s_1(r) + s_1(e) + 1)*||\tilde S'||`.
    #[test]
    fn ensure_orthogonalized_length_perfect_power() {
        for n in 4..8 {
            let params = GadgetParametersRing::init_default(n, &Modulus::from(32));
            let a_bar = PolyOverZ::sample_uniform(&params.n, 0, &params.q).unwrap();

            let (a, r, e) = gen_trapdoor_ring_lwe(&params, &a_bar, &Q::from(5)).unwrap();

            let short_base = gen_short_basis_for_trapdoor_ring(&params, &a, &r, &e);
            let short_base_embedded = short_base.into_coefficient_embedding_from_matrix(n);

            let orthogonalized_short_basis = MatQ::from(&short_base_embedded).gso();

            // Compute s_1(r) and s_1(e).
            let s1_r = {
                let mut r_max = Q::ZERO;
                let r_embedded = r.into_coefficient_embedding_from_matrix(n);
                for i in 0..r_embedded.get_num_columns() {
                    let r_new = r_embedded
                        .get_column(i)
                        .unwrap()
                        .norm_eucl_sqrd()
                        .unwrap()
                        .sqrt();
                    if r_new > r_max {
                        r_max = r_new
                    }
                }
                r_max
            };
            let s1_e = {
                let mut e_max = Q::ZERO;
                let e_embedded = e.into_coefficient_embedding_from_matrix(n);
                for i in 0..e_embedded.get_num_columns() {
                    let e_new = e_embedded
                        .get_column(i)
                        .unwrap()
                        .norm_eucl_sqrd()
                        .unwrap()
                        .sqrt();
                    if e_new > e_max {
                        e_max = e_new
                    }
                }
                e_max
            };

            // Check that all vectors within the orthogonalized base satisfy the length
            // condition.
            let orth_s_length = 2;
            let upper_bound: Q = (s1_r + s1_e + 1) * orth_s_length;
            for i in 0..orthogonalized_short_basis.get_num_columns() {
                let b_tilde_i = orthogonalized_short_basis.get_column(i).unwrap();

                assert!(b_tilde_i.norm_eucl_sqrd().unwrap() <= upper_bound.pow(2).unwrap())
            }
        }
    }

    /// Ensure that the orthogonalized short base length is upper bounded by
    /// `(s_1(r) + s_1(e) + 1)*||\tilde S'||`.
    #[test]
    fn ensure_orthogonalized_length_not_perfect_power() {
        for n in 4..8 {
            let params = GadgetParametersRing::init_default(n, &Modulus::from(42));
            let a_bar = PolyOverZ::sample_uniform(&params.n, 0, &params.q).unwrap();

            let (a, r, e) = gen_trapdoor_ring_lwe(&params, &a_bar, &Q::from(5)).unwrap();

            let short_base = gen_short_basis_for_trapdoor_ring(&params, &a, &r, &e);
            let short_base_embedded = short_base.into_coefficient_embedding_from_matrix(n);

            let orthogonalized_short_basis = MatQ::from(&short_base_embedded).gso();

            // Compute s_1(r) and s_1(e).
            let s1_r = {
                let mut r_max = Q::ZERO;
                let r_embedded = r.into_coefficient_embedding_from_matrix(n);
                for i in 0..r_embedded.get_num_columns() {
                    let r_new = r_embedded
                        .get_column(i)
                        .unwrap()
                        .norm_eucl_sqrd()
                        .unwrap()
                        .sqrt();
                    if r_new > r_max {
                        r_max = r_new
                    }
                }
                r_max
            };
            let s1_e = {
                let mut e_max = Q::ZERO;
                let e_embedded = e.into_coefficient_embedding_from_matrix(n);
                for i in 0..e_embedded.get_num_columns() {
                    let e_new = e_embedded
                        .get_column(i)
                        .unwrap()
                        .norm_eucl_sqrd()
                        .unwrap()
                        .sqrt();
                    if e_new > e_max {
                        e_max = e_new
                    }
                }
                e_max
            };

            // Check that all vectors within the orthogonalized base satisfy the length
            // condition.
            let orth_s_length = Q::from(5).sqrt();
            let upper_bound: Q = (s1_r + s1_e + 1) * orth_s_length;
            for i in 0..orthogonalized_short_basis.get_num_columns() {
                let b_tilde_i = orthogonalized_short_basis.get_column(i).unwrap();

                assert!(b_tilde_i.norm_eucl_sqrd().unwrap() <= upper_bound.pow(2).unwrap())
            }
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
        integer::{MatPolyOverZ, MatZ, PolyOverZ},
        integer_mod_q::{MatPolynomialRingZq, Modulus, PolyOverZq},
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
        let mut sa_r = gen_sa_r(&params, &a);

        sa_r.reduce_by_poly(&PolyOverZ::from(&PolyOverZq::from(&params.modulus)));

        let sa_r_cmp = MatZ::from_str(
            "[\
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],\
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],\
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],\
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],\
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],\
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],\
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],\
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],\
            [0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],\
            [0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],\
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],\
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 1, 0],\
            [0, 0, 2, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0],\
            [0, 0, 0, 0, 0, 0, 2, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0],\
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0],\
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -1, 0, 0, 0, 0, 0, 0, 1, 1],\
            [0, 2, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, -1, 0, 0, 0, 0],\
            [0, 0, 0, 0, 0, 2, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, -1, 0, 0],\
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, -1],\
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -1, 0, 0, 1, 0, 0, 0, 0, 1, 1],\
            [2, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, -1, 0, -1],\
            [0, 0, 0, 0, 2, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, -1],\
            [0, 0, 0, 0, 0, 0, 0, 0, 2, -1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1, 0, 0],\
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -1, 0, 0, 0, 0, 0, 1, 0, 1, 1, 1]]",
        )
        .unwrap();
        assert_eq!(sa_r_cmp, sa_r.into_coefficient_embedding_from_matrix(4));
    }
}

#[cfg(test)]
mod test_compute_s {
    use crate::sample::g_trapdoor::{
        gadget_parameters::GadgetParametersRing, short_basis_ring::compute_s,
    };
    use qfall_math::{
        integer::{MatPolyOverZ, Z},
        integer_mod_q::Modulus,
    };
    use std::str::FromStr;

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
