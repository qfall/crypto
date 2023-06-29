// Copyright © 2023 Niklas Siemer
//
// This file is part of qFALL-crypto.
//
// qFALL-crypto is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

//! This module contains an implementation of the collision-resistant
//! SIS-based hash function.
//!
//! The main references are listed in the following:
//! - \[1\] Peikert, Chris (2016).
//! A decade of lattice cryptography.
//! In: Theoretical Computer Science 10.4.
//! <https://web.eecs.umich.edu/~cpeikert/pubs/lattice-survey.pdf>

use qfall_math::{error::MathError, integer::Z, integer_mod_q::MatZq};
use serde::{Deserialize, Serialize};

#[derive(Debug, Serialize, Deserialize)]
pub struct SISHash {
    key: MatZq, // implicitely contains n = nrows, m = ncols, q = modulus
}

impl SISHash {
    /// Generates a new secret key for an [`SISHash`] instance.
    ///
    /// Parameters:
    /// - `n`: specifies the security parameter and number of rows
    ///   of the uniform at random instantiated matrix `A`
    /// - `m`: specifies the number of columns of matrix `A`
    /// - `q`: specifies the modulus
    ///
    /// Returns a new instance of a [`SISHash`] function with freshly
    /// chosen secret key `A` of type [`MatZq`], dimensions `n x m`,
    /// and modulus `q`. Otherwise, a [`MathError`] is returned, if `n <= 0`
    /// or `m < n log q`, or `q <= ⌈sqrt(n log q)⌉` as then collision-resistance is not ensured.
    ///
    /// # Examples
    /// ```
    /// use qfall_crypto::construction::hash::SISHash;
    ///
    /// let hash = SISHash::gen(5, 18, 11).unwrap();
    /// ```
    ///
    /// # Errors and Failures
    /// - Returns a [`MathError`] of type [`InvalidIntegerInput`](MathError::InvalidIntegerInput)
    /// if `n <= 0` or `m < n log q`, or `q <= ⌈sqrt(n log q)⌉`
    /// as then collision-resistance is not ensured.
    pub fn gen(n: impl Into<Z>, m: impl Into<Z>, q: impl Into<Z>) -> Result<Self, MathError> {
        let n: Z = n.into();
        let m: Z = m.into();
        let q: Z = q.into();

        if n < Z::ONE {
            return Err(MathError::InvalidIntegerInput(String::from(
                "n must be chosen bigger than 0.",
            )));
        }

        // computed according to bullet point 3 of section 4.1.1 in Decade
        let m_bar = (&n * q.log(&2).unwrap()).ceil();

        // m >= m_bar according to bullet point 3 of section 4.1.1 in Decade
        if m < m_bar {
            return Err(MathError::InvalidIntegerInput(String::from(
                "m was chosen smaller than n log q, but it must be larger to satisfy the pigeonhole principle.",
            )));
        }
        // q > ⌈sqrt(m_bar)⌉ according to bullet point 3 + 1 of section 4.1.1 in Decade
        if q <= m_bar.sqrt().ceil() {
            return Err(MathError::InvalidIntegerInput(String::from(
                "q was chosen smaller than ⌈sqrt(n log q)⌉, but it must be larger to satisfy the pigeonhole principle.",
            )));
        }

        let mat_a = MatZq::sample_uniform(&n, &m, q);

        Ok(Self { key: mat_a })
    }

    /// Applies f_A to `value`, i.e. computes `A * value`.
    ///
    /// Parameters:
    /// - `value`: specifies an element from the domain,
    ///   i.e. a column vector of length `m` with modulus `q`
    ///
    /// Returns the hash digest of dimension `n`.
    ///
    /// # Examples
    /// ```
    /// use qfall_crypto::construction::hash::SISHash;
    /// use qfall_math::integer_mod_q::MatZq;
    /// use std::str::FromStr;
    /// let hash = SISHash::gen(1, 3, 7).unwrap();
    /// let value = MatZq::from_str("[[1],[2],[3]] mod 7").unwrap();
    ///
    /// hash.hash(&value);
    /// ```
    ///
    /// # Panics ...
    /// - if `value` isn't a column vector.
    /// - if `value` has a different modulus than the [`SISHash`] instance.
    /// - if `value` has mismatching dimensions, i.e. the vector length isn't `m`.
    pub fn hash(&self, value: &MatZq) -> MatZq {
        if !value.is_column_vector() {
            panic!("The hashed value has to be a column vector!");
        }

        &self.key * value
    }
}

#[cfg(test)]
mod test_gen {
    use super::{SISHash, Z};
    use qfall_math::traits::{GetNumColumns, GetNumRows};

    /// Checks whether too small chosen `n` results in an error.
    #[test]
    fn invalid_n() {
        let res_0 = SISHash::gen(0, 2, 2);
        let res_1 = SISHash::gen(-1, 2, 2);
        let res_2 = SISHash::gen(i64::MIN, 2, 2);

        assert!(res_0.is_err());
        assert!(res_1.is_err());
        assert!(res_2.is_err());
    }

    /// Checks whether too small chosen `m` results in an error.
    #[test]
    fn invalid_m() {
        let res_0 = SISHash::gen(1, 0, 2);
        let res_1 = SISHash::gen(2, 2, 2);
        let res_2 = SISHash::gen(4, 5, i64::MAX);

        assert!(res_0.is_err());
        assert!(res_1.is_err());
        assert!(res_2.is_err());
    }

    /// Checks whether too small chosen `q` results in an error.
    #[test]
    fn invalid_q() {
        let res_0 = SISHash::gen(10, 50, 6);
        let res_1 = SISHash::gen(5, 50, 4);

        assert!(res_0.is_err());
        assert!(res_1.is_err());
    }

    /// Ensures that a working example returns a proper instance.
    #[test]
    fn working_example() {
        let hash = SISHash::gen(5, 18, 11).unwrap();

        assert_eq!(5, hash.key.get_num_rows());
        assert_eq!(18, hash.key.get_num_columns());
        assert_eq!(Z::from(11), Z::from(hash.key.get_mod()));
    }

    /// Ensures that the expected availability is provided.
    #[test]
    fn availability() {
        let _ = SISHash::gen(4i8, 4i8, 4i8);
        let _ = SISHash::gen(4i8, 4i16, 4i32);
        let _ = SISHash::gen(4u8, 4i64, 4u16);
        let _ = SISHash::gen(4u64, 4u32, 4);
        let _ = SISHash::gen(Z::ONE, 4i64, 4u16);
        let _ = SISHash::gen(Z::ONE, Z::from(2), Z::from(2));
    }
}

#[cfg(test)]
mod test_hash {
    use super::{MatZq, SISHash};

    /// Ensures that non-column-vectors result in a panic.
    #[should_panic]
    #[test]
    fn not_column_vec() {
        let hash = SISHash::gen(1, 3, 7).unwrap();
        let value = MatZq::new(1, 3, 7);

        hash.hash(&value);
    }

    /// Ensures that mismatching dimensions result in a panic.
    #[should_panic]
    #[test]
    fn mismatching_dimensions() {
        let hash = SISHash::gen(1, 3, 7).unwrap();
        let value = MatZq::new(4, 1, 7);

        hash.hash(&value);
    }

    /// Ensures that mismatching moduli result in a panic.
    #[should_panic]
    #[test]
    fn mismatching_moduli() {
        let hash = SISHash::gen(1, 3, 7).unwrap();
        let value = MatZq::new(3, 1, 8);

        hash.hash(&value);
    }

    /// Ensures that a working example returns a proper instance.
    #[test]
    fn working_example() {
        let hash = SISHash::gen(1, 3, 7).unwrap();
        let value = MatZq::new(3, 1, 7);

        hash.hash(&value);
    }
}
