// Copyright © 2023 Marvin Beckmann
//
// This file is part of qFALL-crypto.
//
// qFALL-crypto is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

//! Implements a GPV PSF over the polynomial ring according to
//! [\[1\]](<../index.html#:~:text=[1]>) and [\[2\]](<../index.html#:~:text=[2]>)
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
/// [\[1\]](<index.html#:~:text=[1]>) and [\[2\]](<index.html#:~:text=[2]>)
/// using G-Trapdoors where D_n = {e ∈ R^m | |ι(e)| <= s sqrt(m*n) }
/// and R_n = R_q.
///
/// Attributes
/// - `gp`: Describes the gadget parameters with which the G-Trapdoor is generated
/// - `s`: The Gaussian parameter with which elements from the domain are sampled
/// - `s:td`: The Gaussian parameter with which the trapdoor is sampled
///
/// # Examples
/// ```
/// use qfall_crypto::primitive::psf::PSFGPVRing;
/// use qfall_crypto::sample::g_trapdoor::gadget_parameters::GadgetParametersRing;
/// use qfall_math::rational::Q;
/// use qfall_crypto::primitive::psf::PSF;
///
/// let psf = PSFGPVRing {
///     gp: GadgetParametersRing::init_default(8, 512),
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
    /// use qfall_crypto::primitive::psf::PSFGPVRing;
    /// use qfall_crypto::sample::g_trapdoor::gadget_parameters::GadgetParametersRing;
    /// use qfall_math::rational::Q;
    /// use qfall_crypto::primitive::psf::PSF;
    ///
    /// let psf = PSFGPVRing {
    ///     gp: GadgetParametersRing::init_default(8, 512),
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
    /// use qfall_crypto::primitive::psf::PSFGPVRing;
    /// use qfall_crypto::sample::g_trapdoor::gadget_parameters::GadgetParametersRing;
    /// use qfall_math::rational::Q;
    /// use qfall_crypto::primitive::psf::PSF;
    ///
    /// let psf = PSFGPVRing {
    ///     gp: GadgetParametersRing::init_default(8, 512),
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
    /// discrete Gaussian with `f_a(a,e) = u` for a provided syndrome `u`.
    ///
    /// *Note*: the provided parameters `a, r, e, u` must fit together,
    /// otherwise unexpected behavior such as panics may occur.
    ///
    /// Parameters:
    /// - `a`: The parity-check matrix
    /// - `r`: Together with `e` builds a G-Trapdoor for `a`
    /// - `e`: Together with `r` builds a G-Trapdoor for `a`
    /// - `u`: The syndrome from the range
    ///
    /// Returns a sample `e` from the domain on the conditioned discrete
    /// Gaussian distribution `f_a(a,e) = u`.
    ///
    /// # Examples
    /// ```
    /// use qfall_crypto::primitive::psf::PSFGPVRing;
    /// use qfall_crypto::sample::g_trapdoor::gadget_parameters::GadgetParametersRing;
    /// use qfall_math::rational::Q;
    /// use qfall_crypto::primitive::psf::PSF;
    ///
    /// let psf = PSFGPVRing {
    ///     gp: GadgetParametersRing::init_default(8, 512),
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
        // the same as `Rot^-(ι(a)) ι(x) = ι(u)`

        let short_basis = gen_short_basis_for_trapdoor_ring(&self.gp, a, &r.0, &r.1);

        // solve `rot^-(ι(a)) ι(x) = ι(u)` to get solution
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
                    block * self.gp.modulus.get_degree(),
                    (block + 1) * self.gp.modulus.get_degree() - 1,
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
    /// use qfall_crypto::primitive::psf::PSFGPVRing;
    /// use qfall_crypto::sample::g_trapdoor::gadget_parameters::GadgetParametersRing;
    /// use qfall_math::rational::Q;
    /// use qfall_crypto::primitive::psf::PSF;
    ///
    /// let psf = PSFGPVRing {
    ///     gp: GadgetParametersRing::init_default(8, 512),
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
    /// - if `sigma` is not in the domain.
    fn f_a(&self, a: &MatPolynomialRingZq, sigma: &MatPolyOverZ) -> MatPolynomialRingZq {
        assert!(self.check_domain(sigma));
        let sigma = MatPolynomialRingZq::from((sigma, &a.get_mod()));
        a * sigma
    }

    /// Checks whether a value `sigma` is in D_n = {e ∈ R^m | |ι(e)| <= s sqrt(m*n) }.
    ///
    /// Parameters:
    /// - `sigma`: The value for which is checked, if it is in the domain
    ///
    /// Returns true, if `sigma` is in D_n.
    ///
    /// # Examples
    /// ```
    /// use qfall_crypto::primitive::psf::PSFGPVRing;
    /// use qfall_crypto::sample::g_trapdoor::gadget_parameters::GadgetParametersRing;
    /// use qfall_math::rational::Q;
    /// use qfall_crypto::primitive::psf::PSF;
    ///
    /// let psf = PSFGPVRing {
    ///     gp: GadgetParametersRing::init_default(8, 512),
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

        sigma.is_column_vector()
            && m == Z::from(sigma.get_num_rows())
            && Q::from(&sigma_embedded.norm_eucl_sqrd().unwrap())
                <= self.s.pow(2).unwrap() * sigma_embedded.get_num_rows()
    }
}

#[cfg(test)]
mod test_gpv_psf {
    use super::super::gpv_ring::PSFGPVRing;
    use super::PSF;
    use crate::sample::g_trapdoor::gadget_parameters::GadgetParametersRing;
    use qfall_math::integer::{MatPolyOverZ, PolyOverZ};
    use qfall_math::integer_mod_q::MatPolynomialRingZq;
    use qfall_math::rational::Q;
    use qfall_math::traits::{GetNumColumns, GetNumRows, SetEntry};

    fn compute_s(n: i64) -> Q {
        ((2 * 2 * Q::from(1.005_f64) * Q::from(n).sqrt() + 1) * 2) * 4
    }

    /// Ensures that `samp_d` actually computes values that are in D_n.
    #[test]
    fn samp_d_samples_from_dn() {
        let (n, modulus) = (5, 123456789);
        let psf = PSFGPVRing {
            gp: GadgetParametersRing::init_default(n, modulus),
            s: Q::from(1000),
            s_td: Q::from(1.005_f64),
        };

        for _ in 0..5 {
            assert!(psf.check_domain(&psf.samp_d()));
        }
    }

    /// Ensures that `samp_p` actually computes preimages that are also in the correct
    /// domain.
    #[test]
    fn samp_p_preimage_and_domain() {
        for (n, modulus) in [(5, i32::MAX - 57), (6, i32::MAX)] {
            let psf = PSFGPVRing {
                gp: GadgetParametersRing::init_default(n, modulus),
                s: compute_s(n),
                s_td: Q::from(1.005_f64),
            };
            let (a, r) = psf.trap_gen();
            let domain_sample = psf.samp_d();
            let range_fa = psf.f_a(&a, &domain_sample);

            let preimage = psf.samp_p(&a, &r, &range_fa);

            assert_eq!(range_fa, psf.f_a(&a, &preimage));
            assert!(psf.check_domain(&preimage));
        }
    }

    /// Ensures that `f_a` returns `a*sigma`.
    #[test]
    fn f_a_works_as_expected() {
        for (n, modulus) in [(5, 256), (6, 128)] {
            let psf = PSFGPVRing {
                gp: GadgetParametersRing::init_default(n, modulus),
                s: compute_s(n),
                s_td: Q::from(1.005_f64),
            };
            let (a, _) = psf.trap_gen();
            let domain_sample = psf.samp_d();

            let domain_sample_2 = MatPolynomialRingZq::from((&domain_sample, &a.get_mod()));
            assert_eq!(&a * &domain_sample_2, psf.f_a(&a, &domain_sample));
        }
    }

    /// Ensures that `f_a` panics if a value is provided, that is not within the domain.
    /// Sigma is not a vector.
    #[test]
    #[should_panic]
    fn f_a_sigma_not_in_domain_matrix() {
        let psf = PSFGPVRing {
            gp: GadgetParametersRing::init_default(8, 1024),
            s: compute_s(8),
            s_td: Q::from(1.005_f64),
        };
        let (a, _) = psf.trap_gen();
        let not_in_domain = MatPolyOverZ::new(a.get_num_columns(), 2);

        let _ = psf.f_a(&a, &not_in_domain);
    }

    /// Ensures that `f_a` panics if a value is provided, that is not within the domain.
    /// Sigma is not of the correct length.
    #[test]
    #[should_panic]
    fn f_a_sigma_not_in_domain_incorrect_length() {
        let psf = PSFGPVRing {
            gp: GadgetParametersRing::init_default(8, 1024),
            s: compute_s(8),
            s_td: Q::from(1.005_f64),
        };
        let (a, _) = psf.trap_gen();
        let not_in_domain = MatPolyOverZ::new(a.get_num_columns() - 1, 1);

        let _ = psf.f_a(&a, &not_in_domain);
    }

    /// Ensures that `f_a` panics if a value is provided, that is not within the domain.
    /// Sigma is too long.
    #[test]
    #[should_panic]
    fn f_a_sigma_not_in_domain_too_long() {
        let psf = PSFGPVRing {
            gp: GadgetParametersRing::init_default(8, 1024),
            s: compute_s(8),
            s_td: Q::from(1.005_f64),
        };
        let (a, _) = psf.trap_gen();
        let not_in_domain = psf.s.round()
            * a.get_num_columns()
            * 8
            * MatPolyOverZ::identity(a.get_num_columns(), 1);

        let _ = psf.f_a(&a, &not_in_domain);
    }

    /// Ensures that `check_domain` works for vectors with the correct length.
    #[test]
    fn check_domain_as_expected() {
        let psf = PSFGPVRing {
            gp: GadgetParametersRing::init_default(9, 1024),
            s: compute_s(9),
            s_td: Q::from(1.005_f64),
        };
        let (a, _) = psf.trap_gen();
        let value = PolyOverZ::from(psf.s.round() * 3);
        let mut in_domain = MatPolyOverZ::new(a.get_num_columns(), 1);
        for i in 0..in_domain.get_num_rows() {
            in_domain.set_entry(i, 0, &value).unwrap();
        }

        assert!(psf.check_domain(&MatPolyOverZ::new(a.get_num_columns(), 1)));
        assert!(psf.check_domain(&in_domain));
    }

    /// Ensures that `check_domain` returns false for values that are not in the domain.
    #[test]
    fn check_domain_not_in_dn() {
        let psf = PSFGPVRing {
            gp: GadgetParametersRing::init_default(8, 1024),
            s: compute_s(8),
            s_td: Q::from(1.005_f64),
        };
        let (a, _) = psf.trap_gen();

        let matrix = MatPolyOverZ::new(a.get_num_columns(), 2);
        let too_short = MatPolyOverZ::new(a.get_num_columns() - 1, 1);
        let too_long = MatPolyOverZ::new(a.get_num_columns() + 1, 1);
        let entry_too_large = psf.s.round()
            * a.get_num_columns()
            * 8
            * MatPolyOverZ::identity(a.get_num_columns(), 1);

        assert!(!psf.check_domain(&matrix));
        assert!(!psf.check_domain(&too_long));
        assert!(!psf.check_domain(&too_short));
        assert!(!psf.check_domain(&entry_too_large));
    }
}
