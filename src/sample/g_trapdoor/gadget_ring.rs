// Copyright © 2023 Marvin Beckmann
//
// This file is part of qFALL-crypto.
//
// qFALL-crypto is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

//! This module contains an implementation to generate a gadget trapdoor in a
//! ring-lwe setting.

use super::{gadget_classical::find_solution_gadget_mat, gadget_parameters::GadgetParametersRing};
use qfall_math::{
    error::MathError,
    integer::{MatPolyOverZ, PolyOverZ, Z},
    integer_mod_q::{MatPolynomialRingZq, MatZq, PolynomialRingZq},
    rational::Q,
    traits::{Concatenate, GetCoefficient, GetEntry, Pow, SetCoefficient, SetEntry},
};
use std::fmt::Display;

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
/// use qfall_crypto::sample::g_trapdoor::{gadget_parameters::GadgetParametersRing, gadget_ring::gen_trapdoor_ring_lwe};
/// use qfall_math::integer::{Z, PolyOverZ};
/// use qfall_math::integer_mod_q::Modulus;
/// use qfall_math::rational::Q;
///
/// let params = GadgetParametersRing::init_default(8, &Modulus::try_from(&Z::from(17)).unwrap());
/// let a_bar = PolyOverZ::sample_uniform(&params.n, &0, &params.q).unwrap();
///
/// let (a, r, e) = gen_trapdoor_ring_lwe(&params, &a_bar, &Q::from(10)).unwrap();
/// ```
/// # Errors and Failures
/// - Returns a [`MathError`] of type [`InvalidMatrix`](MathError::InvalidMatrix)
/// or of type [`OutOfBounds`](MathError::OutOfBounds), if `params.k`
/// or `params.n` is either `0`,
/// it is negative or it does not fit into an [`i64`].
pub fn gen_trapdoor_ring_lwe(
    params: &GadgetParametersRing,
    a_bar: &PolyOverZ,
    s: &Q,
) -> Result<(MatPolynomialRingZq, MatPolyOverZ, MatPolyOverZ), MathError> {
    // Sample `r` and `e` using a provided distribution
    let r = params.distribution.sample(&params.n, &params.k, s);
    let e = params.distribution.sample(&params.n, &params.k, s);

    // compute the parity check matrix
    // `A = [1 | a | g^t - ar + e]`
    let mut big_a = MatPolyOverZ::new(1, 2);
    big_a.set_entry(0, 0, &PolyOverZ::from(1))?;
    big_a.set_entry(0, 1, a_bar)?;
    let g = gen_gadget_ring(&params.k, &params.base)?;
    big_a = big_a.concat_horizontal(&(g.transpose() - (a_bar * &r + &e)))?;

    Ok((MatPolynomialRingZq::from((&big_a, &params.modulus)), r, e))
}

/// Generates a gadget vector based on its definition in [\[3\]](<../index.html#:~:text=[3]>).
/// This corresponds to a vector `(base ^0, base^1, ..., base^{k-1})` where each entry
/// is a constant polynomial.
///
/// Parameters:
/// - `k`: the size of the gadget vector
/// - `base`: the base with which the entries in the gadget vector are defined
///
/// Returns a gadget vector of length `k` with `base` as its base.
///
/// # Examples
/// ```
/// use qfall_crypto::sample::g_trapdoor::gadget_ring::gen_gadget_ring;
/// use qfall_math::integer::Z;
///
/// let g = gen_gadget_ring(4, &Z::from(2));
/// ```
///
/// # Errors and Failures
/// - Returns a [`MathError`] of type [`InvalidMatrix`](MathError::InvalidMatrix)
/// or of type [`OutOfBounds`](MathError::OutOfBounds), if `k` is either `0`,
/// it is negative or it does not fit into an [`i64`].
pub fn gen_gadget_ring(
    k: impl TryInto<i64> + Display,
    base: &Z,
) -> Result<MatPolyOverZ, MathError> {
    let mut out = MatPolyOverZ::new(k, 1);
    let mut i: i64 = 0;
    while out.set_entry(i, 0, &PolyOverZ::from(base.pow(i)?)).is_ok() {
        i += 1;
    }
    Ok(out)
}

/// Computes an arbitrary solution for `<g^t,x> = value/(Modulus)`.
///
/// Parameters:
/// - `u`: the element for which a solution has to be computed
/// - `k`: the length of a gadget vector
/// - `base`: the base with which the gadget vector is defined
///
/// Returns an arbitrary solution for `<g^t,x> = value/(Modulus)`
///
/// # Examples
/// ```
/// use qfall_math::integer::{Z, MatZ, PolyOverZ};
/// use qfall_math::integer_mod_q::{Zq, Modulus, PolynomialRingZq};
/// use qfall_crypto::sample::g_trapdoor::gadget_ring::gen_gadget_ring;
/// use qfall_crypto::sample::g_trapdoor::gadget_ring::find_solution_gadget_ring;
/// use qfall_crypto::sample::g_trapdoor::gadget_parameters::GadgetParametersRing;
/// use qfall_math::integer_mod_q::MatPolynomialRingZq;
/// use qfall_math::traits::GetEntry;
/// use std::str::FromStr;
///
/// let gp = GadgetParametersRing::init_default(9, &Modulus::from(128));
///
/// let gadget = gen_gadget_ring(&gp.k, &gp.base).unwrap();
/// let gadget = MatPolynomialRingZq::from((&gadget, &gp.modulus));
///
/// let u = PolyOverZ::from_str("10  5 124 12 14 14 1 2 4 1 5").unwrap();
/// let u = PolynomialRingZq::from((&u, &gp.modulus));
///
/// let solution = find_solution_gadget_ring(&u, &gp.k, &gp.base);
/// let solution = MatPolynomialRingZq::from((&solution, &gp.modulus));
/// assert_eq!(u, gadget.dot_product(&solution).unwrap())
/// ```
///
/// # Panics ...
/// - if the modulus of the value is greater than `base^k`.
pub fn find_solution_gadget_ring(u: &PolynomialRingZq, k: &Z, base: &Z) -> MatPolyOverZ {
    let k_i64 = i64::try_from(k).unwrap();
    let ring_poly = u.get_poly();
    let n = ring_poly.get_degree();
    let mut value = MatZq::new(n + 1, 1, u.get_mod().get_q());

    for i in 0..=n {
        let coeff = ring_poly.get_coeff(i).unwrap();
        value.set_entry(i, 0, coeff).unwrap();
    }
    let classical_sol = find_solution_gadget_mat(&value, k, base);

    let mut out = MatPolyOverZ::new(1, k);
    for i in 0..k_i64 {
        let mut poly = PolyOverZ::default();
        for j in 0..=n {
            let entry = classical_sol.get_entry(i + j * k_i64, 0).unwrap();
            poly.set_coeff(j, &entry).unwrap();
        }
        out.set_entry(0, i, &poly).unwrap();
    }
    out
}

#[cfg(test)]
mod test_gen_trapdoor_ring {

    use crate::sample::g_trapdoor::{
        gadget_parameters::GadgetParametersRing, gadget_ring::gen_trapdoor_ring_lwe,
    };
    use qfall_math::{
        integer::{MatPolyOverZ, PolyOverZ, Z},
        integer_mod_q::{MatPolynomialRingZq, Modulus},
        rational::Q,
        traits::{Concatenate, GetCoefficient, GetEntry, GetNumColumns, GetNumRows, Pow},
    };

    /// Computes a trapdoor using the given secrets `(r,e)`
    fn compute_trapdoor(r: &MatPolyOverZ, e: &MatPolyOverZ, k: &Z) -> MatPolyOverZ {
        let i_k = MatPolyOverZ::identity(k, k);

        e.concat_vertical(r).unwrap().concat_vertical(&i_k).unwrap()
    }

    /// Assure that the trapdoor `r` returned from [`gen_trapdoor`] is actually a
    /// trapdoor for `a`
    #[test]
    fn is_trapdoor() {
        let modulus = Modulus::try_from(&Z::from(32)).unwrap();
        let params = GadgetParametersRing::init_default(6, &modulus);
        let a_bar = PolyOverZ::sample_uniform(&params.n, &0, &params.q).unwrap();

        // call gen_trapdoor to get matrix a and its 'trapdoor' r
        let (a, r, e) = gen_trapdoor_ring_lwe(&params, &a_bar, &Q::from(10)).unwrap();

        // generate the trapdoor for a from r as trapdoor = [[e],[r],[I]]
        let trapdoor =
            MatPolynomialRingZq::from((&compute_trapdoor(&r, &e, &params.k), &params.modulus));

        // ensure G = A*trapdoor (definition of a trapdoor)
        let res: MatPolynomialRingZq = &a * &trapdoor;

        assert_eq!(params.k, Z::from(res.get_num_columns()));
        assert_eq!(1, res.get_num_rows());

        for i in 0..(&params.k).try_into().unwrap() {
            let res_entry: PolyOverZ = res.get_entry(0, i).unwrap();
            assert_eq!(res_entry.get_coeff(0).unwrap(), params.base.pow(i).unwrap())
        }
    }
}

#[cfg(test)]
mod test_find_solution_gadget_ring {
    use super::{find_solution_gadget_ring, gen_gadget_ring};
    use crate::sample::g_trapdoor::gadget_parameters::GadgetParametersRing;
    use qfall_math::{
        integer::PolyOverZ,
        integer_mod_q::{MatPolynomialRingZq, Modulus, PolynomialRingZq},
    };
    use std::str::FromStr;

    /// Ensures that the algorithm finds a correct solution such that `<g^t, x> = u`
    #[test]
    fn is_correct_solution() {
        let gp = GadgetParametersRing::init_default(3, &Modulus::from(32));

        let gadget = gen_gadget_ring(&gp.k, &gp.base).unwrap();
        let gadget = MatPolynomialRingZq::from((&gadget, &gp.modulus));

        let u = PolyOverZ::from_str("10  5 124 12 14 14 1 2 4 1 5").unwrap();
        let u = PolynomialRingZq::from((&u, &gp.modulus));

        let solution = find_solution_gadget_ring(&u, &gp.k, &gp.base);
        let solution = MatPolynomialRingZq::from((&solution, &gp.modulus));

        assert_eq!(u, gadget.dot_product(&solution).unwrap())
    }
}
