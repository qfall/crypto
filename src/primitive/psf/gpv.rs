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
/// G-Trapdoors where D_n = {e ∈ Z^m | |e| <= s sqrt(m)}
/// and R_n = Z_q^n.
///
/// Attributes
/// - `gp`: Describes the gadget parameters with which the G-Trapdoor is generated
/// - `s`: The gaussian parameter with which is sampled
///
/// # Examples
/// ```
/// use qfall_crypto::primitive::psf::gpv::PSFGPV;
/// use qfall_crypto::sample::g_trapdoor::gadget_parameters::GadgetParameters;
/// use qfall_math::rational::Q;
/// use qfall_crypto::primitive::psf::PSF;
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
    /// use qfall_crypto::primitive::psf::gpv::PSFGPV;
    /// use qfall_crypto::sample::g_trapdoor::gadget_parameters::GadgetParameters;
    /// use qfall_math::rational::Q;
    /// use qfall_crypto::primitive::psf::PSF;
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
    /// use qfall_crypto::primitive::psf::gpv::PSFGPV;
    /// use qfall_crypto::sample::g_trapdoor::gadget_parameters::GadgetParameters;
    /// use qfall_math::rational::Q;
    /// use qfall_crypto::primitive::psf::PSF;
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
    /// *Note*: the provided parameters `a,r,u` must fit together,
    /// otherwise unexpected behavior such as panics may occur.
    ///
    /// Parameters:
    /// - `a`: The parity-check matrix
    /// - `short_base`: The short base for `Λ^⟂(A)`
    /// - `short_base_gso`: The precomputed GSO of the short_base
    /// - `u`: The syndrome from the range
    ///
    /// Returns a sample `e` from the domain on the conditioned discrete
    /// gaussian distribution `f_a(a,e) = u`.
    ///
    /// # Examples
    /// ```
    /// use qfall_crypto::primitive::psf::gpv::PSFGPV;
    /// use qfall_crypto::sample::g_trapdoor::gadget_parameters::GadgetParameters;
    /// use qfall_math::rational::Q;
    /// use qfall_crypto::primitive::psf::PSF;
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

    /// Implements the efficiently computable function `f_a` which here corresponds to
    /// `a*sigma`. The sigma must be from the domain, i.e. D_n.
    ///
    /// Parameters:
    /// - `a`: The parity-check matrix of dimensions `n x m`
    /// - `sigma`: A column vector of length `m`
    ///
    /// Returns `a*sigma`
    ///
    /// # Examples
    /// ```
    /// use qfall_crypto::primitive::psf::gpv::PSFGPV;
    /// use qfall_crypto::sample::g_trapdoor::gadget_parameters::GadgetParameters;
    /// use qfall_math::rational::Q;
    /// use qfall_crypto::primitive::psf::PSF;
    ///
    /// let psf = PSFGPV {
    ///     gp: GadgetParameters::init_default(8, 64),
    ///     s: Q::from(12),
    /// };
    /// let (a, td) = psf.trap_gen();
    /// let domain_sample = psf.samp_d();
    /// let range_fa = psf.f_a(&a, &domain_sample);
    /// ```
    ///
    /// # Panics ...
    /// - if `sigma` is not in the domain.
    fn f_a(&self, a: &MatZq, sigma: &MatZ) -> MatZq {
        assert!(self.check_domain(sigma));
        a * sigma
    }

    /// Checks whether a value `sigma` is in D_n = {e ∈ Z^m | |e| <= s sqrt(m)}.
    ///
    /// Parameters:
    /// - `sigma`: The value for which is checked, if it is in the domain
    ///
    /// Returns true, if `sigma` is in D_n.
    ///
    /// # Examples
    /// ```
    /// use qfall_crypto::primitive::psf::PSF;
    /// use qfall_crypto::primitive::psf::gpv::PSFGPV;
    /// use qfall_crypto::sample::g_trapdoor::gadget_parameters::GadgetParameters;
    /// use qfall_math::rational::Q;
    ///
    /// let psf = PSFGPV {
    ///     gp: GadgetParameters::init_default(8, 64),
    ///     s: Q::from(12),
    /// };
    /// let (a, td) = psf.trap_gen();
    ///
    /// let vector = psf.samp_d();
    ///
    /// assert!(psf.check_domain(&vector));
    /// ```
    fn check_domain(&self, sigma: &MatZ) -> bool {
        let m = &self.gp.n * &self.gp.k + &self.gp.m_bar;
        sigma.is_column_vector()
            && m == Z::from(sigma.get_num_rows())
            && Q::from(&sigma.norm_eucl_sqrd().unwrap()) <= self.s.pow(2).unwrap() * &m
    }
}

#[cfg(test)]
mod test_gpv_psf {
    use super::super::gpv::PSFGPV;
    use super::PSF;
    use crate::sample::g_trapdoor::gadget_parameters::GadgetParameters;
    use qfall_math::integer::MatZ;
    use qfall_math::rational::Q;
    use qfall_math::traits::{GetNumColumns, GetNumRows, SetEntry};

    /// Ensures that `samp_d` actually computes values that are in D_n.
    #[test]
    fn samp_d_samples_from_dn() {
        for (n, modulus) in [(5, 256), (10, 128), (15, 157)] {
            let psf = PSFGPV {
                gp: GadgetParameters::init_default(n, modulus),
                s: Q::from(10),
            };

            for _ in 0..5 {
                assert!(psf.check_domain(&psf.samp_d()));
            }
        }
    }

    /// Ensures that `samp_p` actually computes preimages that are also in the correct
    /// domain.
    #[test]
    fn samp_p_preimage_and_domain() {
        for (n, modulus) in [(5, 256), (6, 128)] {
            let psf = PSFGPV {
                gp: GadgetParameters::init_default(n, modulus),
                s: Q::from(10),
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
            let psf = PSFGPV {
                gp: GadgetParameters::init_default(n, modulus),
                s: Q::from(10),
            };
            let (a, _) = psf.trap_gen();
            let domain_sample = psf.samp_d();

            assert_eq!(&a * &domain_sample, psf.f_a(&a, &domain_sample));
        }
    }

    /// Ensures that `f_a` panics if a value is provided, that is not within the domain.
    /// Sigma is not a vector.
    #[test]
    #[should_panic]
    fn f_a_sigma_not_in_domain_matrix() {
        let psf = PSFGPV {
            gp: GadgetParameters::init_default(8, 128),
            s: Q::from(10),
        };
        let (a, _) = psf.trap_gen();
        let not_in_domain = MatZ::new(a.get_num_columns(), 2);

        let _ = psf.f_a(&a, &not_in_domain);
    }

    /// Ensures that `f_a` panics if a value is provided, that is not within the domain.
    /// Sigma is not of the correct length.
    #[test]
    #[should_panic]
    fn f_a_sigma_not_in_domain_incorrect_length() {
        let psf = PSFGPV {
            gp: GadgetParameters::init_default(8, 128),
            s: Q::from(10),
        };
        let (a, _) = psf.trap_gen();
        let not_in_domain = MatZ::new(a.get_num_columns() - 1, 1);

        let _ = psf.f_a(&a, &not_in_domain);
    }

    /// Ensures that `f_a` panics if a value is provided, that is not within the domain.
    /// Sigma is too long.
    #[test]
    #[should_panic]
    fn f_a_sigma_not_in_domain_too_long() {
        let psf = PSFGPV {
            gp: GadgetParameters::init_default(8, 128),
            s: Q::from(10),
        };
        let (a, _) = psf.trap_gen();
        let not_in_domain =
            psf.s.round() * a.get_num_columns() * MatZ::identity(a.get_num_columns(), 1);

        let _ = psf.f_a(&a, &not_in_domain);
    }

    /// Ensures that `check_domain` works for vectors with the correct length.
    #[test]
    fn check_domain_as_expected() {
        let psf = PSFGPV {
            gp: GadgetParameters::init_default(8, 128),
            s: Q::from(10),
        };
        let (a, _) = psf.trap_gen();
        let value = psf.s.round();
        let mut in_domain = MatZ::new(a.get_num_columns(), 1);
        for i in 0..in_domain.get_num_rows() {
            in_domain.set_entry(i, 0, &value).unwrap();
        }

        assert!(psf.check_domain(&MatZ::new(a.get_num_columns(), 1)));
        assert!(psf.check_domain(&in_domain));
    }

    /// Ensures that `check_domain` returns false for values that are not in the domain.
    #[test]
    fn check_domain_not_in_dn() {
        let psf = PSFGPV {
            gp: GadgetParameters::init_default(8, 128),
            s: Q::from(10),
        };
        let (a, _) = psf.trap_gen();

        let matrix = MatZ::new(a.get_num_columns(), 2);
        let too_short = MatZ::new(a.get_num_columns() - 1, 1);
        let too_long = MatZ::new(a.get_num_columns() + 1, 1);
        let entry_too_large =
            psf.s.round() * a.get_num_columns() * MatZ::identity(a.get_num_columns(), 1);

        assert!(!psf.check_domain(&matrix));
        assert!(!psf.check_domain(&too_long));
        assert!(!psf.check_domain(&too_short));
        assert!(!psf.check_domain(&entry_too_large));
    }
}
