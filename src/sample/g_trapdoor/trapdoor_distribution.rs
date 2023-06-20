// Copyright © 2023 Marvin Beckmann
//
// This file is part of qFALL-crypto.
//
// qFALL-crypto is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

//! This module contains all implementation of TrapdoorDistributions from
//! which the matrix `\bar A` is sampled in the trapdoor generation algorithm.

use qfall_math::integer::MatZ;
use serde::{Deserialize, Serialize};

/// This trait should be implemented by all distributions which should be
/// used to generate a trapdoor.
/// Within the [`GadgetParameters`](super::GadgetParameters) there is distribution which
/// can be defined, from which the matrix `\bar A` is sampled.
///
/// Parameters:
/// - `m_bar`: number of rows of the matrix that is sampled
/// - `w`: number of columns of the matrix that is sampled
///
/// Returns a matrix which is sampled according to the defined distribution
#[typetag::serde]
pub trait TrapdoorDistribution {
    fn sample(&self, m_bar: i64, w: i64) -> MatZ;
}

/// A distribution which samples a matrix of type [`MatZ`] with entries in `\{-1,0,1\}`
/// with probability `1/4` for `-1` and `1` an probability `1/2` for `0`
#[derive(Serialize, Deserialize)]
pub struct PlusMinusOneZero;

#[typetag::serde]
impl TrapdoorDistribution for PlusMinusOneZero {
    /// Sample a matrix from distribution with probability `1/2` for `0`
    /// and `1/4` each for `+/-1`.
    ///
    /// Parameters:
    /// - `m_bar`: number of columns of the matrix
    /// - `w`: number of rows of the matrix
    ///
    /// Returns a matrix where each entry is sampled independently with probability
    /// `1/2` for `0` and `1/4` each for `+/-1`.
    ///
    /// # Examples
    /// ```
    /// use qfall_crypto::sample::g_trapdoor::trapdoor_distribution::{PlusMinusOneZero, TrapdoorDistribution};
    /// use qfall_math::integer::Z;
    /// use qfall_math::integer_mod_q::Modulus;
    ///
    /// let mat = PlusMinusOneZero.sample(42, 24);
    /// ```
    fn sample(&self, m_bar: i64, w: i64) -> MatZ {
        let mat_1 = MatZ::sample_uniform(m_bar, w, &0, &2).unwrap();
        let mat_2 = MatZ::sample_uniform(m_bar, w, &0, &2).unwrap();
        mat_1 - mat_2
    }
}

#[cfg(test)]
mod test_pm_one_zero {
    use super::PlusMinusOneZero;
    use super::TrapdoorDistribution;
    use qfall_math::integer::Z;
    use qfall_math::traits::GetEntry;

    /// ensure that the distribution samples in its correct range
    #[test]
    fn correct_range() {
        let sample = PlusMinusOneZero.sample(10, 5);

        for i in 0..10 {
            for j in 0..5 {
                assert!(
                    Z::MINUS_ONE <= sample.get_entry(i, j).unwrap()
                        && Z::ONE >= sample.get_entry(i, j).unwrap()
                );
            }
        }
    }
}
