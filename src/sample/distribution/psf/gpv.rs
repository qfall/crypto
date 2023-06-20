// Copyright © 2023 Marvin Beckmann
//
// This file is part of qFALL-crypto.
//
// qFALL-crypto is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

//! Implements a GPV PSF according to: ... using G-Trapdoors (TODO ref etc)

use super::PSF;
use crate::sample::g_trapdoor::{
    gadget_classical::gen_trapdoor, gadget_parameters::GadgetParameters,
    short_basis_classical::gen_short_basis_for_trapdoor,
};
use qfall_math::{
    integer::MatZ,
    integer_mod_q::MatZq,
    rational::{MatQ, Q},
    traits::{GetNumRows, Pow},
};

pub struct PSFGPV {
    pub pp: GadgetParameters,
    pub s: Q,
}

impl PSF<MatZq, MatZ, MatZ, MatZq> for PSFGPV {
    fn trap_gen(&self) -> (MatZq, MatZ) {
        let a_bar = MatZq::sample_uniform(&self.pp.n, &self.pp.m_bar, &self.pp.q);

        let tag = MatZq::identity(&self.pp.n, &self.pp.n, &self.pp.q);

        gen_trapdoor(&self.pp, &a_bar, &tag).unwrap()
    }

    fn samp_d(&self) -> MatZ {
        todo!()
    }

    fn samp_p(&self, a: &MatZq, r: &MatZ, u: &MatZq) -> MatZ {
        let tag = MatZq::identity(&self.pp.n, &self.pp.n, &self.pp.q);

        let short_basis = gen_short_basis_for_trapdoor(&self.pp, &tag, a, r);

        let sol: MatZ = (&a.solve_gaussian_elimination(u).unwrap()).into();

        let center = MatQ::from(&(-1 * &sol));
        let s = &self.s * (Q::from(2) * Q::PI).sqrt();
        let signature: MatZ = MatZ::sample_d(&short_basis, &self.pp.n, &center, &s).unwrap() + sol;

        signature
    }

    fn fa(&self, a: &MatZq, value: &MatZ) -> MatZq {
        a * value
    }

    fn check_dn(&self, sigma: &MatZ) -> bool {
        Q::from(&sigma.norm_eucl_sqrd().unwrap()) <= self.s.pow(2).unwrap() * sigma.get_num_rows()
    }
}
