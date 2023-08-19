// Copyright © 2023 Marvin Beckmann
//
// This file is part of qFALL-crypto.
//
// qFALL-crypto is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

//! Implements a GPV PSF according to [\[1\]](<../index.html#:~:text=[1]>)
//! using G-Trapdoors to generate a short basis and corresponding trapdoor.

use super::PSF;
use crate::sample::g_trapdoor::{
    gadget_classical::gen_trapdoor, gadget_parameters::GadgetParameters,
    short_basis_classical::gen_short_basis_for_trapdoor,
};
use qfall_math::{
    integer::{MatZ, Z},
    integer_mod_q::MatZq,
    rational::{MatQ, Q},
    traits::{GetNumRows, Pow},
};
use serde::{Deserialize, Serialize};

/// A lattice-based implementation of a [`PSF`] according to
/// [\[1\]](<../index.html#:~:text=[1]>) using
/// G-Trapdoors where D_n = \{e \in \Z^m | |e| \leq s \sqrt{m}\}
/// and R_n = Z_q^n.
///
/// Attributes
/// - `gp`: Describes the gadget parameters with which the G-Trapdoor is generated
/// - `s`: The gaussian parameter with which is sampled
///
/// # Examples
/// ```
/// use qfall_crypto::sample::distribution::psf::gpv::PSFGPV;
/// use qfall_crypto::sample::g_trapdoor::gadget_parameters::GadgetParameters;
/// use qfall_math::rational::Q;
/// use qfall_math::integer_mod_q::Modulus;
/// use crate::qfall_crypto::sample::distribution::psf::PSF;
///
/// let modulus = Modulus::from(64);
/// let psf = PSFGPV {
///     gp: GadgetParameters::init_default(8, &modulus),
///     s: Q::from(12),
/// };
///
/// let (a, r) = psf.trap_gen();
/// let domain_sample = psf.samp_d();
/// let range_fa = psf.f_a(&a, &domain_sample);
/// let preimage = psf.samp_p(&a, &r, &range_fa);
///
/// assert!(psf.check_domain(&preimage));
/// ```
#[derive(Serialize, Deserialize)]
pub struct PSFGPV {
    pub gp: GadgetParameters,
    pub s: Q,
}

impl PSF<MatZq, MatZ, MatZ, MatZq> for PSFGPV {
    /// Computes a G-Trapdoor according to the [`GadgetParameters`]
    /// defined in the struct.
    ///
    /// # Examples
    /// ```
    /// use qfall_crypto::sample::distribution::psf::gpv::PSFGPV;
    /// use qfall_crypto::sample::g_trapdoor::gadget_parameters::GadgetParameters;
    /// use qfall_math::rational::Q;
    /// use qfall_math::integer_mod_q::Modulus;
    /// use crate::qfall_crypto::sample::distribution::psf::PSF;
    ///
    /// let modulus = Modulus::from(64);
    /// let psf = PSFGPV {
    ///     gp: GadgetParameters::init_default(8, &modulus),
    ///     s: Q::from(12),
    /// };
    ///
    /// let (a, r) = psf.trap_gen();
    /// ```
    fn trap_gen(&self) -> (MatZq, MatZ) {
        let a_bar = MatZq::sample_uniform(&self.gp.n, &self.gp.m_bar, &self.gp.q);

        let tag = MatZq::identity(&self.gp.n, &self.gp.n, &self.gp.q);

        gen_trapdoor(&self.gp, &a_bar, &tag).unwrap()
    }

    /// Samples in the domain using SampleD with the standard basis and center `0`.
    ///
    /// # Examples
    /// ```
    /// use qfall_crypto::sample::distribution::psf::gpv::PSFGPV;
    /// use qfall_crypto::sample::g_trapdoor::gadget_parameters::GadgetParameters;
    /// use qfall_math::rational::Q;
    /// use qfall_math::integer_mod_q::Modulus;
    /// use crate::qfall_crypto::sample::distribution::psf::PSF;
    ///
    /// let modulus = Modulus::from(64);
    /// let psf = PSFGPV {
    ///     gp: GadgetParameters::init_default(8, &modulus),
    ///     s: Q::from(12),
    /// };
    /// let (a, r) = psf.trap_gen();
    ///
    /// let domain_sample = psf.samp_d();
    /// ```
    fn samp_d(&self) -> MatZ {
        let m = &self.gp.n * &self.gp.k + &self.gp.m_bar;
        MatZ::sample_d_common(&m, &self.gp.n, &self.s).unwrap()
    }

    /// Samples an `e` in the domain using SampleD with a short basis that is generated
    /// from the G-Trapdoor from the conditioned conditioned
    /// discrete gaussian with `f_a(a,e) = u` for a provided syndrome `u`.
    ///
    /// Parameters:
    /// - `a`: The parity-check matrix
    /// - `r`: The G-Trapdoor for `a`
    /// - `u`: The syndrome from the range
    ///
    /// Returns a sample `e` from the domain on the conditioned discrete
    /// gaussian distribution `f_a(a,e) = u`.
    ///
    /// # Examples
    /// ```
    /// use qfall_crypto::sample::distribution::psf::gpv::PSFGPV;
    /// use qfall_crypto::sample::g_trapdoor::gadget_parameters::GadgetParameters;
    /// use qfall_math::rational::Q;
    /// use qfall_math::integer_mod_q::Modulus;
    /// use crate::qfall_crypto::sample::distribution::psf::PSF;
    ///
    /// let modulus = Modulus::from(64);
    /// let psf = PSFGPV {
    ///     gp: GadgetParameters::init_default(8, &modulus),
    ///     s: Q::from(12),
    /// };
    /// let (a, r) = psf.trap_gen();
    /// let domain_sample = psf.samp_d();
    /// let range_fa = psf.f_a(&a, &domain_sample);
    ///
    /// let preimage = psf.samp_p(&a, &r, &range_fa);
    /// assert_eq!(range_fa, psf.f_a(&a, &preimage))
    /// ```
    fn samp_p(&self, a: &MatZq, r: &MatZ, u: &MatZq) -> MatZ {
        let tag = MatZq::identity(&self.gp.n, &self.gp.n, &self.gp.q);
        let short_basis = gen_short_basis_for_trapdoor(&self.gp, &tag, a, r);

        let sol: MatZ = (&a.solve_gaussian_elimination(u).unwrap()).into();

        let center = MatQ::from(&(-1 * &sol));

        sol + MatZ::sample_d(&short_basis, &self.gp.n, &center, &self.s).unwrap()
    }

    /// Implements the efficiently computable function `fa` which here corresponds to
    /// `a*value`
    ///
    /// Parameters:
    /// - `a`: The parity-check matrix of dimensions `n x m`
    /// - `value`: A column vector of length `m`
    ///
    /// Returns `a*value`
    ///
    /// # Examples
    /// ```
    /// use qfall_crypto::sample::distribution::psf::gpv::PSFGPV;
    /// use qfall_crypto::sample::g_trapdoor::gadget_parameters::GadgetParameters;
    /// use qfall_math::rational::Q;
    /// use qfall_math::integer_mod_q::Modulus;
    /// use crate::qfall_crypto::sample::distribution::psf::PSF;
    ///
    /// let modulus = Modulus::from(64);
    /// let psf = PSFGPV {
    ///     gp: GadgetParameters::init_default(8, &modulus),
    ///     s: Q::from(12),
    /// };
    /// let (a, r) = psf.trap_gen();
    /// let domain_sample = psf.samp_d();
    /// let range_fa = psf.f_a(&a, &domain_sample);
    /// ```
    fn f_a(&self, a: &MatZq, value: &MatZ) -> MatZq {
        a * value
    }

    /// Checks whether a value `sigma` is in D_n = \{e \in \Z^m | |e| \leq s \sqrt{m}\}.
    ///
    /// Attributes
    /// - `sigma`: The value for which is checked, if it is in the domain
    ///
    /// Returns true, if `sigma` is in D_n.
    ///
    /// # Examples
    /// ```
    /// use qfall_crypto::sample::distribution::psf::gpv::PSFGPV;
    /// use qfall_crypto::sample::g_trapdoor::gadget_parameters::GadgetParameters;
    /// use qfall_math::rational::Q;
    /// use qfall_math::integer_mod_q::Modulus;
    /// use crate::qfall_crypto::sample::distribution::psf::PSF;
    /// use qfall_math::integer::MatZ;
    /// use qfall_math::traits::GetNumColumns;
    ///
    /// let modulus = Modulus::from(64);
    /// let psf = PSFGPV {
    ///     gp: GadgetParameters::init_default(8, &modulus),
    ///     s: Q::from(12),
    /// };
    /// let (a, r) = psf.trap_gen();
    ///
    /// let vector = psf.samp_d();
    ///
    /// assert!(psf.check_domain(&vector));
    /// ```
    fn check_domain(&self, sigma: &MatZ) -> bool {
        let m = &self.gp.n * &self.gp.k + &self.gp.m_bar;
        Q::from(&sigma.norm_eucl_sqrd().unwrap()) <= self.s.pow(2).unwrap() * &m
            && m == Z::from(sigma.get_num_rows())
    }
}
