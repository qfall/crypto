// Copyright © 2023 Marvin Beckmann
//
// This file is part of qFALL-crypto.
//
// qFALL-crypto is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

//! This module contains the [`GadgetParameters`] object,
//! which contains all parameters which are needed to generate a classical G-Trapdoor.

use super::trapdoor_distribution::{PlusMinusOneZero, TrapdoorDistribution};
use qfall_math::{integer::Z, integer_mod_q::Modulus, traits::Pow};
use serde::{Deserialize, Serialize};

/// Collects all parameters which are necessary to compute a g-trapdoor.
/// You can either use [`GadgetParameters::init_default`] or set all values
/// and distributions yourself.
///
/// Attributes:
/// - `n`: the security parameter
/// - `k`: the size of the gadget vector: mostly taken as `log_base(q)`
/// - `m_bar`: has to be chose appropriately for regularity and for
/// the distribution to be subgaussian
/// - `base`: the base with which the gadget-vector and matrix are generated
/// - `q`: the modulus
/// - `distribution`: the distribution from which the matrix `\bar A` is sampled
///
/// # Examples
/// ```
/// use qfall_crypto::sample::g_trapdoor::gadget_parameters::GadgetParameters;
/// use qfall_math::integer::Z;
/// use qfall_math::integer_mod_q::Modulus;
///
/// let params = GadgetParameters::init_default(42, &Modulus::from(42));
/// ```
#[derive(Serialize, Deserialize)]
pub struct GadgetParameters {
    pub n: Z,
    pub k: Z,
    pub m_bar: Z,
    pub base: Z,
    pub q: Modulus,
    pub distribution: Box<dyn TrapdoorDistribution>,
}

impl GadgetParameters {
    /// Initializes default values for [`GadgetParameters`] to create a classical
    /// G-trapdoor.
    ///
    /// - `base = 2` is taken from [\[1\]](<../index.html#:~:text=[1]>): Theorem 1.
    /// - `k = log_2_ceil(q)` is taken from [\[1\]](<../index.html#:~:text=[1]>):
    /// Theorem 1.
    /// - `w = n * k`: As it is required to match the dimension of the gadget matrix,
    /// hence it has to  equal to `n * size of gadget_vec`
    /// - `m_bar = n * k + \log(n)^2`: is taken from [\[1\]](<../index.html#:~:text=[1]>)
    /// as a function satisfying `m_bar = n \log q + \omega(\log n)`
    /// - the distribution is taken as [`PlusMinusOneZero`],
    /// see the example from [\[1\]](<../index.html#:~:text=[1]>): after statistical instantiation in section 3.2
    ///
    /// Parameters:
    /// - `n`: the security parameter for the generation
    /// - `modulus`: the modulus over which the TrapGen operates
    ///
    /// Returns an instantiation of default GadgetParameters based on the references mentioned above
    ///
    /// # Examples
    /// ```
    /// use qfall_crypto::sample::g_trapdoor::gadget_parameters::GadgetParameters;
    /// use qfall_math::integer::Z;
    /// use qfall_math::integer_mod_q::Modulus;
    ///
    /// let params = GadgetParameters::init_default(42, &Modulus::from(42));
    /// ```
    ///
    /// # Panics...
    /// - ... if the security parameter `n` is not in `\[1, i64::MAX\]`.
    pub fn init_default(n: impl Into<Z>, modulus: &Modulus) -> Self {
        // panic if n < 1 (security parameter must be positive) and not larger than
        // [`i64`] because downstream matrices can be at most that size
        let n = n.into();
        assert!(n >= Z::ONE && n <= Z::from(i64::MAX));

        let base = Z::from(2);
        let log_q = Z::from(modulus).log_ceil(&base).unwrap();
        let n_log_q = &n * &log_q;
        let log_n = n.log_ceil(&base).unwrap();
        let m_bar = &n_log_q + &log_n.pow(2).unwrap();
        Self {
            n,
            k: log_q,
            m_bar,
            base,
            q: modulus.clone(),
            distribution: Box::new(PlusMinusOneZero),
        }
    }
}
