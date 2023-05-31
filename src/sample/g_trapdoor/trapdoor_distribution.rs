// Copyright © 2023 Marvin Beckmann
//
// This file is part of qFALL-crypto.
//
// qFALL-crypto is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

//! This module contains all implementation of TrapdoorDistributions from
//! which the matrix `\bar A` is sampled in the trapdoor generation algorithm.

use qfall_math::{
    integer::{MatPolyOverZ, MatZ, PolyOverZ, Z},
    rational::Q,
    traits::{SetCoefficient, SetEntry},
};

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
pub trait TrapdoorDistribution {
    fn sample(&self, m_bar: &Z, w: &Z) -> MatZ;
}

pub trait TrapdoorDistributionRing {
    fn sample(&self, n: &Z, nr_rows: &Z, nr_cols: &Z, s: &Q) -> MatPolyOverZ;
}

/// A distribution which samples a matrix of type [`MatZ`] with entries in `\{-1,0,1\}`
/// with probability `1/4` for `-1` and `1` an probability `1/2` for `0`
pub struct PlusMinusOneZero;

pub struct SampleZ;

impl TrapdoorDistribution for PlusMinusOneZero {
    /// Sample a matrix from distribution with probability `1/2` for `0`
    /// and `1/4` each for `+/-1`
    ///
    /// Parameters:
    /// - `m_bar`: number of columns of the matrix
    /// - `w`: number of rows of the matrix
    ///
    /// # Examples
    /// ```
    /// use qfall_crypto::sample::g_trapdoor::trapdoor_distribution::{PlusMinusOneZero, TrapdoorDistribution};
    /// use qfall_math::integer::Z;
    /// use qfall_math::integer_mod_q::Modulus;
    ///
    /// let mat = PlusMinusOneZero.sample(&42.into(), &24.into());
    /// ```
    /// # Panics...
    /// - ... `m_bar` or `w` does not fit into in `i64` or is smaller than `1`.
    fn sample(&self, m_bar: &Z, w: &Z) -> MatZ {
        let mat_1 = MatZ::sample_uniform(m_bar, w, &0, &2).unwrap();
        let mat_2 = MatZ::sample_uniform(m_bar, w, &0, &2).unwrap();
        mat_1 - mat_2
    }
}

impl TrapdoorDistributionRing for SampleZ {
    fn sample(&self, n: &Z, nr_rows: &Z, nr_cols: &Z, s: &Q) -> MatPolyOverZ {
        let n = i64::try_from(n).unwrap();
        let nr_rows = i64::try_from(nr_rows).unwrap();
        let nr_cols = i64::try_from(nr_cols).unwrap();
        let mut out_mat = MatPolyOverZ::new(nr_rows, nr_cols).unwrap();
        for i in 0..nr_rows {
            for j in 0..nr_cols {
                let mut sample = PolyOverZ::default();
                for k in 0..n {
                    let sample_z = Z::sample_discrete_gauss(&n, &Z::ZERO, s).unwrap();
                    sample.set_coeff(k, &sample_z).unwrap();
                }
                out_mat.set_entry(i, j, &sample).unwrap();
            }
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

    /// ensure that the distribution samples in its correct range
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