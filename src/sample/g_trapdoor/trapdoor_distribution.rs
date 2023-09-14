// Copyright © 2023 Marvin Beckmann
//
// This file is part of qFALL-crypto.
//
// qFALL-crypto is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

//! This module contains all implementation of TrapdoorDistributions from
//! which the matrix `A_bar` is sampled in the trapdoor generation algorithm.

use qfall_math::{
    integer::{MatPolyOverZ, MatZ, PolyOverZ, Z},
    rational::Q,
    traits::SetEntry,
};
use serde::{Deserialize, Serialize};

/// This trait should be implemented by all distributions which should be
/// used to generate a trapdoor.
#[typetag::serde]
pub trait TrapdoorDistribution {
    /// Sample from a matrix according to a predefined distribution.
    ///
    /// Parameters:
    /// - `m_bar`: number of rows of the matrix that is sampled
    /// - `w`: number of columns of the matrix that is sampled
    ///
    /// Returns a matrix which is sampled according to the defined distribution.
    fn sample(&self, m_bar: &Z, w: &Z) -> MatZ;
}

/// This trait should be implemented by all distributions which should be
/// used to generate a trapdoor over a ring.
#[typetag::serde]
pub trait TrapdoorDistributionRing {
    /// Sample a matrix of polynomials of length `n` with entries sampled
    /// using a predefined distribution.
    ///
    /// Parameters:
    /// - `n`: length of the polynomial
    /// - `nr_cols`: number of columns of the matrix
    /// - `s`: the Gaussian parameter used for SampleZ
    ///
    /// Returns a matrix where each entry is a polynomials of length `n`, sampled
    /// using the defined distribution.
    fn sample(&self, n: &Z, nr_cols: &Z, s: &Q) -> MatPolyOverZ;
}

/// A distribution which samples a matrix of type [`MatZ`] with entries in `\{-1,0,1\}`
/// with probability `1/4` for `-1` and `1` an probability `1/2` for `0`
#[derive(Serialize, Deserialize)]
pub struct PlusMinusOneZero;

/// A distribution which samples a row vector of type [`MatPolyOverZ`] where each
/// coefficient is a polynomial of degree `n-1` and each coefficient of the polynomial
/// is sampled using [`Z::sample_discrete_gauss`]
#[derive(Serialize, Deserialize)]
pub struct SampleZ;

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
    ///
    /// let mat = PlusMinusOneZero.sample(&42.into(), &24.into());
    /// ```
    ///
    /// # Panics...
    /// - if `m_bar` or `w` does not fit into in `i64` or is smaller than `1`.
    fn sample(&self, m_bar: &Z, w: &Z) -> MatZ {
        let mat_1 = MatZ::sample_uniform(m_bar, w, 0, 2).unwrap();
        let mat_2 = MatZ::sample_uniform(m_bar, w, 0, 2).unwrap();
        mat_1 - mat_2
    }
}

#[typetag::serde]
impl TrapdoorDistributionRing for SampleZ {
    /// Sample a matrix of polynomials of length `n` with entries sampled
    /// using [`Z::sample_discrete_gauss`].
    ///
    /// Parameters:
    /// - `n`: length of the polynomial
    /// - `nr_cols`: number of columns of the matrix
    /// - `s`: the Gaussian parameter used for SampleZ
    ///
    /// Returns a matrix where each entry is a polynomials of length `n`, sampled
    /// using [`Z::sample_discrete_gauss`].
    ///
    /// # Examples
    /// ```
    /// use qfall_crypto::sample::g_trapdoor::trapdoor_distribution::{SampleZ, TrapdoorDistributionRing};
    ///
    /// let mat = SampleZ.sample(&42.into(), &24.into(), &3.into());
    /// ```
    ///
    /// # Panics...
    /// - if `n`, `nr_rows` or `nr_cols` does not fit into in `i64`
    /// or is smaller than `1`.
    fn sample(&self, n: &Z, nr_cols: &Z, s: &Q) -> MatPolyOverZ {
        let n = i64::try_from(n).unwrap();
        let nr_cols = i64::try_from(nr_cols).unwrap();
        let mut out_mat = MatPolyOverZ::new(1, nr_cols);
        for j in 0..nr_cols {
            let sample = PolyOverZ::sample_discrete_gauss(n - 1, n, 0, s).unwrap();
            out_mat.set_entry(0, j, &sample).unwrap();
        }

        out_mat
    }
}

#[cfg(test)]
mod test_pm_one_zero {
    use super::PlusMinusOneZero;
    use super::TrapdoorDistribution;
    use qfall_math::integer::Z;
    use qfall_math::traits::GetEntry;

    /// Ensure that the distribution samples in its correct range.
    #[test]
    fn correct_range() {
        let sample = PlusMinusOneZero.sample(&10.into(), &5.into());

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

#[cfg(test)]
mod test_sample_z {
    use super::{SampleZ, TrapdoorDistributionRing};
    use qfall_math::{
        rational::Q,
        traits::{GetNumRows, Pow},
    };

    /// Ensure that the distribution samples are in the correct range,
    #[test]
    fn correct_range_high_prob() {
        for _ in 0..20 {
            let s = 5.into();
            let sample = SampleZ.sample(&10.into(), &15.into(), &s);

            // it should be the same as sampling a vector with 10*15 entries
            let coeff_embedding = sample
                .transpose()
                .into_coefficient_embedding_from_matrix(10);

            // test for concentration bound
            assert!(
                Q::from(coeff_embedding.norm_eucl_sqrd().unwrap())
                    <= s.pow(2).unwrap() * coeff_embedding.get_num_rows()
            );
        }
    }
}
