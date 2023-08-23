// Copyright © 2023 Marvin Beckmann
//
// This file is part of qFALL-crypto.
//
// qFALL-crypto is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

//! Implements a GPV PSF over the polynomial ring according to
//! [\[1\]](<../index.html#:~:text=[1]>) and [\[3\]](<../index.html#:~:text=[3]>)
//! using G-Trapdoors to generate a short basis and corresponding trapdoor.

use super::PSF;
use crate::{
    sample::g_trapdoor::{
        gadget_parameters::GadgetParametersRing, gadget_ring::gen_trapdoor_ring_lwe,
        short_basis_ring::gen_short_basis_for_trapdoor_ring,
    },
    utils::rotation_matrix::rot_minus_matrix,
};
use qfall_math::{
    integer::{MatPolyOverZ, MatZ, PolyOverZ, Z},
    integer_mod_q::{MatPolynomialRingZq, MatZq},
    rational::{MatQ, PolyOverQ, Q},
    traits::{FromCoefficientEmbedding, GetNumRows, Pow},
};
use serde::{Deserialize, Serialize};

/// A lattice-based implementation of a [`PSF`] according to
/// [\[1\]](<../index.html#:~:text=[1]>) and [\[3\]](<../index.html#:~:text=[3]>)
/// using G-Trapdoors where D_n = {e ∈ R^m | |iota(e)| <= s sqrt(m*n) }
/// and R_n = R_q.
///
/// Attributes
/// - `gp`: Describes the gadget parameters with which the G-Trapdoor is generated
/// - `s`: The gaussian parameter with which elements from the domain are sampled
/// - `s:td`: The gaussian parameter with which the trapdoor is sampled
///
/// # Examples
/// ```
/// use qfall_crypto::sample::distribution::psf::gpv_ring::PSFGPVRing;
/// use qfall_crypto::sample::g_trapdoor::gadget_parameters::GadgetParametersRing;
/// use qfall_math::rational::Q;
/// use qfall_math::integer_mod_q::Modulus;
/// use crate::qfall_crypto::sample::distribution::psf::PSF;
/// use qfall_math::integer::MatZ;
/// use qfall_math::traits::GetNumColumns;
///
/// let modulus = Modulus::from(512);
/// let psf = PSFGPVRing {
///     gp: GadgetParametersRing::init_default(8, &modulus),
///     s: Q::from(100),
///     s_td: Q::from(1.005_f64),
/// };
///
/// let (a, (r, e)) = psf.trap_gen();
/// let domain_sample = psf.samp_d();
/// let range_fa = psf.f_a(&a, &domain_sample);
/// let preimage = psf.samp_p(&a, &(r,e), &range_fa);
///
/// assert!(psf.check_domain(&preimage));
/// ```
#[derive(Serialize, Deserialize)]
pub struct PSFGPVRing {
    pub gp: GadgetParametersRing,
    pub s: Q,
    pub s_td: Q,
}

impl PSF<MatPolynomialRingZq, (MatPolyOverZ, MatPolyOverZ), MatPolyOverZ, MatPolynomialRingZq>
    for PSFGPVRing
{
    /// Computes a G-Trapdoor according to the [`GadgetParametersRing`]
    /// defined in the struct.
    ///
    /// # Examples
    /// ```
    /// use qfall_crypto::sample::distribution::psf::gpv_ring::PSFGPVRing;
    /// use qfall_crypto::sample::g_trapdoor::gadget_parameters::GadgetParametersRing;
    /// use qfall_math::rational::Q;
    /// use qfall_math::integer_mod_q::Modulus;
    /// use crate::qfall_crypto::sample::distribution::psf::PSF;
    /// use qfall_math::integer::MatZ;
    /// use qfall_math::traits::GetNumColumns;
    ///
    /// let modulus = Modulus::from(512);
    /// let psf = PSFGPVRing {
    ///     gp: GadgetParametersRing::init_default(8, &modulus),
    ///     s: Q::from(100),
    ///     s_td: Q::from(1.005_f64),
    /// };
    /// let (a, (r, e)) = psf.trap_gen();
    /// ```
    fn trap_gen(&self) -> (MatPolynomialRingZq, (MatPolyOverZ, MatPolyOverZ)) {
        let a_bar =
            PolyOverZ::sample_uniform(self.gp.modulus.get_degree() - 1, 0, &self.gp.q).unwrap();
        let (a, r, e) = gen_trapdoor_ring_lwe(&self.gp, &a_bar, &self.s_td).unwrap();

        (a, (r, e))
    }

    /// Samples in the domain using SampleD with the standard basis and center `0`.
    ///
    /// # Examples
    /// ```
    /// use qfall_crypto::sample::distribution::psf::gpv_ring::PSFGPVRing;
    /// use qfall_crypto::sample::g_trapdoor::gadget_parameters::GadgetParametersRing;
    /// use qfall_math::rational::Q;
    /// use qfall_math::integer_mod_q::Modulus;
    /// use crate::qfall_crypto::sample::distribution::psf::PSF;
    /// use qfall_math::integer::MatZ;
    /// use qfall_math::traits::GetNumColumns;
    ///
    /// let modulus = Modulus::from(512);
    /// let psf = PSFGPVRing {
    ///     gp: GadgetParametersRing::init_default(8, &modulus),
    ///     s: Q::from(100),
    ///     s_td: Q::from(1.005_f64),
    /// };
    /// let (a, (r, e)) = psf.trap_gen();
    ///
    /// let domain_sample = psf.samp_d();
    /// ```
    fn samp_d(&self) -> MatPolyOverZ {
        let dimension = self.gp.modulus.get_degree() * (&self.gp.k + 2);
        let sample = MatZ::sample_d_common(dimension, &self.gp.n, &self.s).unwrap();
        MatPolyOverZ::from_coefficient_embedding_to_matrix(&sample, self.gp.modulus.get_degree())
    }

    /// Samples an `e` in the domain using SampleD with a short basis that is generated
    /// from the G-Trapdoor from the conditioned conditioned
    /// discrete gaussian with `f_a(a,e) = u` for a provided syndrome `u`.
    ///
    /// *Note*: the provided parameters `a,r,e,u` must fit together,
    /// otherwise unexpected behavior such as panics may occur.
    ///
    /// Parameters:
    /// - `a`: The parity-check matrix
    /// - `r`: Together with `e` builds a G-Trapdoor for `a`
    /// - `e`: Together with `r` builds a G-Trapdoor for `a`
    /// - `u`: The syndrome from the range
    ///
    /// Returns a sample `e` from the domain on the conditioned discrete
    /// gaussian distribution `f_a(a,e) = u`.
    ///
    /// # Examples
    /// ```
    /// use qfall_crypto::sample::distribution::psf::gpv_ring::PSFGPVRing;
    /// use qfall_crypto::sample::g_trapdoor::gadget_parameters::GadgetParametersRing;
    /// use qfall_math::rational::Q;
    /// use qfall_math::integer_mod_q::Modulus;
    /// use crate::qfall_crypto::sample::distribution::psf::PSF;
    /// use qfall_math::integer::MatZ;
    /// use qfall_math::traits::GetNumColumns;
    ///
    /// let modulus = Modulus::from(512);
    /// let psf = PSFGPVRing {
    ///     gp: GadgetParametersRing::init_default(8, &modulus),
    ///     s: Q::from(100),
    ///     s_td: Q::from(1.005_f64),
    /// };
    /// let (a, (r, e)) = psf.trap_gen();
    ///
    /// let domain_sample = psf.samp_d();
    /// let range_fa = psf.f_a(&a, &domain_sample);
    ///
    /// let preimage = psf.samp_p(&a, &(r,e), &range_fa);
    /// assert_eq!(range_fa, psf.f_a(&a, &preimage))
    /// ```
    fn samp_p(
        &self,
        a: &MatPolynomialRingZq,
        r: &(MatPolyOverZ, MatPolyOverZ),
        u: &MatPolynomialRingZq,
    ) -> MatPolyOverZ {
        // compute solution to `a*x = u`
        // the same as `Rot^-(\iota(a)) \iota(x) = \iota(u)`

        let short_basis = gen_short_basis_for_trapdoor_ring(&self.gp, a, &r.0, &r.1);

        // solve `rot^-(\iota(a)) \iota(x) = \iota(u)` to get solution
        let u_embedded = u
            .get_mat()
            .into_coefficient_embedding_from_matrix(self.gp.modulus.get_degree());
        let a_embedded = a
            .get_mat()
            .into_coefficient_embedding_from_matrix(self.gp.modulus.get_degree());
        let rot_a = rot_minus_matrix(&a_embedded);

        let u_embedded = MatZq::from((&u_embedded, &self.gp.q));
        let rot_a = MatZq::from((&rot_a, &self.gp.q));
        let sol: MatZ = (&rot_a.solve_gaussian_elimination(&u_embedded).unwrap()).into();

        // turn center into a vector of polynomials over Q with maximal degree as the
        // modulus
        let center = MatQ::from(&(-1 * &sol));
        let mut center_embedded = Vec::new();
        for block in 0..(center.get_num_rows() / (self.gp.modulus.get_degree())) {
            let sub_mat = center
                .get_submatrix(
                    block * (self.gp.modulus.get_degree() - 1),
                    (block + 1) * (self.gp.modulus.get_degree() - 1) - 1,
                    0,
                    0,
                )
                .unwrap();
            let embedded_sub_mat = PolyOverQ::from_coefficient_embedding(&sub_mat);
            center_embedded.push(embedded_sub_mat);
        }

        MatPolyOverZ::from_coefficient_embedding_to_matrix(&sol, self.gp.modulus.get_degree())
            + MatPolyOverZ::sample_d(
                &short_basis,
                self.gp.modulus.get_degree(),
                &self.gp.n,
                &center_embedded,
                &self.s,
            )
            .unwrap()
    }

    /// Implements the efficiently computable function `f_a` which here corresponds to
    /// `a*sigma`.
    ///
    /// Parameters:
    /// - `a`: The parity-check matrix of dimensions `n x m`
    /// - `sigma`: A column vector of length `m`
    ///
    /// Returns `a*sigma`
    ///
    /// # Examples
    /// ```
    /// use qfall_crypto::sample::distribution::psf::gpv_ring::PSFGPVRing;
    /// use qfall_crypto::sample::g_trapdoor::gadget_parameters::GadgetParametersRing;
    /// use qfall_math::rational::Q;
    /// use qfall_math::integer_mod_q::Modulus;
    /// use crate::qfall_crypto::sample::distribution::psf::PSF;
    /// use qfall_math::integer::MatZ;
    /// use qfall_math::traits::GetNumColumns;
    ///
    /// let modulus = Modulus::from(512);
    /// let psf = PSFGPVRing {
    ///     gp: GadgetParametersRing::init_default(8, &modulus),
    ///     s: Q::from(100),
    ///     s_td: Q::from(1.005_f64),
    /// };
    /// let (a, (r, e)) = psf.trap_gen();
    ///
    /// let domain_sample = psf.samp_d();
    /// let range_fa = psf.f_a(&a, &domain_sample);
    /// ```
    ///
    /// # Panics ...
    /// - if the number of rows of `sigma` does not match the number of columns of `a`.
    /// - if `sigma` is not a column vector.
    fn f_a(&self, a: &MatPolynomialRingZq, sigma: &MatPolyOverZ) -> MatPolynomialRingZq {
        assert!(
            sigma.is_column_vector(),
            "The input vector is not a column vector."
        );
        let sigma = MatPolynomialRingZq::from((sigma, &a.get_mod()));
        a * sigma
    }

    /// Checks whether a value `sigma` is in
    /// D_n = {e ∈ R^m | |iota(e)| <= s sqrt(m*n) }.
    ///
    /// Attributes
    /// - `sigma`: The value for which is checked, if it is in the domain
    ///
    /// Returns true, if `sigma` is in D_n.
    ///
    /// # Examples
    /// ```
    /// use qfall_crypto::sample::distribution::psf::gpv_ring::PSFGPVRing;
    /// use qfall_crypto::sample::g_trapdoor::gadget_parameters::GadgetParametersRing;
    /// use qfall_math::rational::Q;
    /// use qfall_math::integer_mod_q::Modulus;
    /// use crate::qfall_crypto::sample::distribution::psf::PSF;
    /// use qfall_math::integer::MatZ;
    /// use qfall_math::traits::GetNumColumns;
    ///
    /// let modulus = Modulus::from(512);
    /// let psf = PSFGPVRing {
    ///     gp: GadgetParametersRing::init_default(8, &modulus),
    ///     s: Q::from(100),
    ///     s_td: Q::from(1.005_f64),
    /// };
    /// let (a, (r, e)) = psf.trap_gen();
    ///
    /// let vector = psf.samp_d();
    ///
    /// assert!(psf.check_domain(&vector));
    /// ```
    fn check_domain(&self, sigma: &MatPolyOverZ) -> bool {
        let m = &self.gp.k + 2;
        let nr_coeffs = self.gp.modulus.get_degree();
        let sigma_embedded = sigma.into_coefficient_embedding_from_matrix(nr_coeffs);

        Q::from(&sigma_embedded.norm_eucl_sqrd().unwrap())
            <= self.s.pow(2).unwrap() * sigma_embedded.get_num_rows()
            && m == Z::from(sigma.get_num_rows())
    }
}
