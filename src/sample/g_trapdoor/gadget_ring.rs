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

pub fn gen_trapdoor(
    params: &GadgetParametersRing,
    s: &Q,
) -> Result<
    (
        MatPolynomialRingZq,
        MatPolynomialRingZq,
        MatPolynomialRingZq,
    ),
    MathError,
> {
    // Sample a polynomial of degree n-1 uniformly at random with coefficients in `[0,q)`
    let a = PolyOverZ::sample_uniform(&params.n, &0, &params.q)?;

    // Sample `r` and `e` using a provided distribution e.g. SampleZ
    let r = params.distribution.sample(&params.n, &Z::ONE, &params.k, s);
    let e = params.distribution.sample(&params.n, &Z::ONE, &params.k, s);

    // compute the parity check matrix
    // `A = [1 | a | g_1 - (a*r_1 + e_1) | ... | g_k - (a*r_k + e_k)]`
    let mut big_a = MatPolyOverZ::new(1, &(&params.k + Z::from(2)))?;
    let gadget_vec = gen_gadget_vec(&params.k, &params.base)?;
    let one = PolyOverZ::from_str("1  1")?;
    big_a.set_entry(0, 0, &one)?;
    big_a.set_entry(0, 1, &a)?;

    let k = i64::try_from(&params.k)?;

    for i in 0..k {
        let g_i_z = gadget_vec.get_entry(i, 0)?;
        let mut g_i = PolyOverZ::default();
        g_i.set_coeff(0, &g_i_z)?;
        let r_i = r.get_entry(0, i)?;
        let e_i = e.get_entry(0, i)?;

        big_a.set_entry(0, i + 2, g_i - (&a * r_i + e_i))?;
    }

    Ok((
        MatPolynomialRingZq::from((&big_a, &params.modulus)),
        MatPolynomialRingZq::from((&r, &params.modulus)),
        MatPolynomialRingZq::from((&e, &params.modulus)),
    ))
}

#[cfg(test)]
mod test {
    use super::gen_trapdoor;
    use crate::sample::g_trapdoor::gadget_parameters::GadgetParametersRing;
    use qfall_math::{
        integer::{PolyOverZ, Z},
        integer_mod_q::Modulus,
        rational::Q,
        traits::GetEntry,
    };

    #[test]
    fn working() {
        let params =
            GadgetParametersRing::init_default(8, &Modulus::try_from(&Z::from(17)).unwrap());

        let (a, r, e) = gen_trapdoor(&params, &Q::from(10)).unwrap();

        for i in 0..7 {
            let a_i: PolyOverZ = a.get_entry(0, i).unwrap();
            println!("a_{i}: {a_i}");
        }
        for i in 0..5 {
            let r_i: PolyOverZ = r.get_entry(0, i).unwrap();
            println!("r_{i}: {r_i}");
        }
        for i in 0..5 {
            let e_i: PolyOverZ = e.get_entry(0, i).unwrap();
            println!("e_{i}: {e_i}");
        }
    }
}