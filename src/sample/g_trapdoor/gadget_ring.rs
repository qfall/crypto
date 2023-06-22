// Copyright © 2023 Marvin Beckmann
//
// This file is part of qFALL-crypto.
//
// qFALL-crypto is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

//! This module contains an implementation to generate a gadget trapdoor in a
//! ring-lwe setting.

use super::{gadget_classical::gen_gadget_vec, gadget_parameters::GadgetParametersRing};
use qfall_math::{
    error::MathError,
    integer::{MatPolyOverZ, PolyOverZ, Z},
    integer_mod_q::MatPolynomialRingZq,
    rational::Q,
    traits::{GetEntry, SetCoefficient, SetEntry},
};
use std::str::FromStr;

/// Generates a trapdoor according to Construction 1 in [\[2\]](<../index.html#:~:text=[2]>).
/// Namely:
/// - Generates the gadget matrix: `G`
/// - Samples the trapdoor `R` from the specified distribution in `params`
/// - Outputs
/// `([1 | a_bar | g_1 - (a_bar * r_1 + e_1) | ... | g_k - (a_bar * r_k + e_k) ], r, e)`
/// as a tuple of `(A,R)`, where `R` defines a trapdoor for `A`.
///
/// Parameters:
/// - `params`: all gadget parameters which are required to generate the trapdoor
/// - `a_bar`: the matrix defining the second part of the G-Trapdoor
/// - `tag`: the tag which is hidden within the matrix ´A`
/// - `s`: defining the deviation of the distribution from which `r` and `e` is sampled
///
/// Returns a a parity-check matrix `a` derived from `a_bar` and its gadget-trapdoor
/// `r` under a give tag `h`.
///
/// # Examples
/// ```
/// use qfall_crypto::sample::g_trapdoor::{gadget_parameters::GadgetParametersRing, gadget_ring::gen_trapdoor_ring};
/// use qfall_math::integer::{Z, PolyOverZ};
/// use qfall_math::integer_mod_q::Modulus;
/// use qfall_math::rational::Q;
///
/// let params = GadgetParametersRing::init_default(8, &Modulus::try_from(&Z::from(17)).unwrap());
/// let a_bar = PolyOverZ::sample_uniform(&params.n, &0, &params.q).unwrap();
///
/// let (a, r, e) = gen_trapdoor_ring(&params, &a_bar, &Q::from(10)).unwrap();
/// ```
/// # Errors and Failures
/// - Returns a [`MathError`] of type [`InvalidMatrix`](MathError::InvalidMatrix)
/// or of type [`OutOfBounds`](MathError::OutOfBounds), if `params.k`
/// or `params.n` is either `0`,
/// it is negative or it does not fit into an [`i64`].
pub fn gen_trapdoor_ring(
    params: &GadgetParametersRing,
    a_bar: &PolyOverZ,
    s: &Q,
) -> Result<
    (
        MatPolynomialRingZq,
        MatPolynomialRingZq,
        MatPolynomialRingZq,
    ),
    MathError,
> {
    // Sample `r` and `e` using a provided distribution
    let r = params.distribution.sample(&params.n, &Z::ONE, &params.k, s);
    let e = params.distribution.sample(&params.n, &Z::ONE, &params.k, s);

    // compute the parity check matrix
    // `A = [1 | a | g_1 - (a*r_1 + e_1) | ... | g_k - (a*r_k + e_k)]`
    let mut big_a = MatPolyOverZ::new(1, &(&params.k + Z::from(2)));
    let gadget_vec = gen_gadget_vec(&params.k, &params.base)?;
    let one = PolyOverZ::from_str("1  1")?;
    big_a.set_entry(0, 0, &one)?;
    big_a.set_entry(0, 1, a_bar)?;

    let k = i64::try_from(&params.k)?;

    for i in 0..k {
        let g_i_z = gadget_vec.get_entry(i, 0)?;
        let mut g_i = PolyOverZ::default();
        g_i.set_coeff(0, &g_i_z)?;
        let r_i = r.get_entry(0, i)?;
        let e_i = e.get_entry(0, i)?;

        big_a.set_entry(0, i + 2, g_i - (a_bar * r_i + e_i))?;
    }

    Ok((
        MatPolynomialRingZq::from((&big_a, &params.modulus)),
        MatPolynomialRingZq::from((&r, &params.modulus)),
        MatPolynomialRingZq::from((&e, &params.modulus)),
    ))
}

#[cfg(test)]
mod test_gen_trapdoor_ring {

    use crate::sample::g_trapdoor::{
        gadget_classical::gen_gadget_vec, gadget_parameters::GadgetParametersRing,
        gadget_ring::gen_trapdoor_ring,
    };
    use qfall_math::{
        integer::{PolyOverZ, Z},
        integer_mod_q::{MatPolynomialRingZq, Modulus},
        rational::Q,
        traits::{
            Concatenate, GetCoefficient, GetEntry, GetNumColumns, GetNumRows, Pow, SetCoefficient,
            SetEntry,
        },
    };

    /// assure that the trapdoor `r` returned from [`gen_trapdoor`] is actually a
    /// trapdoor for `a`
    #[test]
    fn is_trapdoor() {
        let modulus = Modulus::try_from(&Z::from(32)).unwrap();
        let params = GadgetParametersRing::init_default(42, &modulus);
        let a_bar = PolyOverZ::sample_uniform(&params.n, &0, &params.q).unwrap();

        // call gen_trapdoor to get matrix a and its 'trapdoor' r
        let (a, r, e) = gen_trapdoor_ring(&params, &a_bar, &Q::from(10)).unwrap();

        let mut unit_vector = MatPolynomialRingZq::new(1, &params.k, &r.get_mod());
        for i in 0..(&params.k).try_into().unwrap() {
            let mut one = PolyOverZ::default();
            one.set_coeff(i, 1).unwrap();
            unit_vector.set_entry(0, i, &one).unwrap();
        }

        // generate the trapdoor for a from r as trapdoor = [[r],[I]]
        let trapdoor = e
            .concat_vertical(&r)
            .unwrap()
            .concat_vertical(&unit_vector)
            .unwrap();

        // ensure G = A*trapdoor (definition of a trapdoor)
        let gadget_vec = gen_gadget_vec(&params.k, &params.base).unwrap();

        let res: MatPolynomialRingZq = &a * &trapdoor;

        assert_eq!(1, res.get_num_columns());
        assert_eq!(1, res.get_num_rows());

        let res_entry: PolyOverZ = res.get_entry(0, 0).unwrap();
        for i in 0..(&params.n).try_into().unwrap() {
            assert_eq!(res_entry.get_coeff(i).unwrap(), params.base.pow(i).unwrap())
        }
    }
}
