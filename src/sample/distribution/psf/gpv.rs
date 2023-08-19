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
/// use crate::qfall_crypto::sample::distribution::psf::PSF;
///
/// let psf = PSFGPV {
///     gp: GadgetParameters::init_default(8, 64),
///     s: Q::from(12),
/// };
///
/// let (a, td) = psf.trap_gen();
/// let domain_sample = psf.samp_d();
/// let range_fa = psf.f_a(&a, &domain_sample);
/// let preimage = psf.samp_p(&a, &td, &range_fa);
///
/// assert!(psf.check_domain(&preimage));
/// ```
#[derive(Serialize, Deserialize)]
pub struct PSFGPV {
    pub gp: GadgetParameters,
    pub s: Q,
}

impl PSF<MatZq, (MatZ, MatQ), MatZ, MatZq> for PSFGPV {
    /// Computes a G-Trapdoor according to the [`GadgetParameters`]
    /// defined in the struct.
    /// It returns a matrix `A` together with a short base and its GSO.
    ///
    /// # Examples
    /// ```
    /// use qfall_crypto::sample::distribution::psf::gpv::PSFGPV;
    /// use qfall_crypto::sample::g_trapdoor::gadget_parameters::GadgetParameters;
    /// use qfall_math::rational::Q;
    /// use crate::qfall_crypto::sample::distribution::psf::PSF;
    ///
    /// let psf = PSFGPV {
    ///     gp: GadgetParameters::init_default(8, 64),
    ///     s: Q::from(12),
    /// };
    ///
    /// let (a, (sh_b, sh_b_gso)) = psf.trap_gen();
    /// ```
    fn trap_gen(&self) -> (MatZq, (MatZ, MatQ)) {
        let a_bar = MatZq::sample_uniform(&self.gp.n, &self.gp.m_bar, &self.gp.q);

        let tag = MatZq::identity(&self.gp.n, &self.gp.n, &self.gp.q);

        let (a, r) = gen_trapdoor(&self.gp, &a_bar, &tag).unwrap();

        let short_base = gen_short_basis_for_trapdoor(&self.gp, &tag, &a, &r);
        let short_base_gso = MatQ::from(&short_base).gso();

        (a, (short_base, short_base_gso))
    }

    /// Samples in the domain using SampleD with the standard basis and center `0`.
    ///
    /// # Examples
    /// ```
    /// use qfall_crypto::sample::distribution::psf::gpv::PSFGPV;
    /// use qfall_crypto::sample::g_trapdoor::gadget_parameters::GadgetParameters;
    /// use qfall_math::rational::Q;
    /// use crate::qfall_crypto::sample::distribution::psf::PSF;
    ///
    /// let psf = PSFGPV {
    ///     gp: GadgetParameters::init_default(8, 64),
    ///     s: Q::from(12),
    /// };
    /// let (a, td) = psf.trap_gen();
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
    /// - `short_base`: The short base for `\Lambda^\perp(A)`
    /// - `short_base_gso`: The precomputed GSO of the short_base
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
    /// use crate::qfall_crypto::sample::distribution::psf::PSF;
    ///
    /// let psf = PSFGPV {
    ///     gp: GadgetParameters::init_default(8, 64),
    ///     s: Q::from(12),
    /// };
    /// let (a, td) = psf.trap_gen();
    /// let domain_sample = psf.samp_d();
    /// let range_fa = psf.f_a(&a, &domain_sample);
    ///
    /// let preimage = psf.samp_p(&a, &td, &range_fa);
    /// assert_eq!(range_fa, psf.f_a(&a, &preimage))
    /// ```
    fn samp_p(&self, a: &MatZq, (short_base, short_base_gso): &(MatZ, MatQ), u: &MatZq) -> MatZ {
        let sol: MatZ = (&a.solve_gaussian_elimination(u).unwrap()).into();

        let center = MatQ::from(&(-1 * &sol));

        sol + MatZ::sample_d_precomputed_gso(
            short_base,
            short_base_gso,
            &self.gp.n,
            &center,
            &self.s,
        )
        .unwrap()
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
    /// use crate::qfall_crypto::sample::distribution::psf::PSF;
    ///
    /// let psf = PSFGPV {
    ///     gp: GadgetParameters::init_default(8, 64),
    ///     s: Q::from(12),
    /// };
    /// let (a, td) = psf.trap_gen();
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
    /// use crate::qfall_crypto::sample::distribution::psf::PSF;
    /// use qfall_math::integer::MatZ;
    /// use qfall_math::traits::GetNumColumns;
    ///
    /// let psf = PSFGPV {
    ///     gp: GadgetParameters::init_default(8, 64),
    ///     s: Q::from(12),
    /// };
    /// let (a, td) = psf.trap_gen();
    ///
    /// let vector = MatZ::new(a.get_num_columns(), 1);
    ///
    /// assert!(psf.check_domain(&vector));
    /// ```
    fn check_domain(&self, sigma: &MatZ) -> bool {
        let m = &self.gp.n * &self.gp.k + &self.gp.m_bar;
        Q::from(&sigma.norm_eucl_sqrd().unwrap()) <= self.s.pow(2).unwrap() * &m
            && m == Z::from(sigma.get_num_rows())
    }
}
