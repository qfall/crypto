// Copyright © 2023 Phil Milewski
//
// This file is part of qFALL-crypto.
//
// qFALL-crypto is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

//! This module contains hashes into different domains.

use qfall_math::{
    integer::Z,
    integer_mod_q::{MatZq, Modulus, Zq},
    traits::SetEntry,
};
use sha2::{Digest, Sha256};
use std::fmt::Display;

/// Computes the sha256 hash value of a given String literal.
///
/// Parameters:
/// - `string`: specifies the value that is hashed.
///
/// # Examples
/// ```
/// use qfall_crypto::primitive::hash::sha256;
///
/// let string = "Hello World!";
/// let hash = sha256(string);
/// assert_eq!("7f83b1657ff1fc53b92dc18148a1d65dfc2d4b1fa3d677284addd200126d9069", hash);
/// ```
pub fn sha256(string: &str) -> String {
    let mut hasher = Sha256::new();
    hasher.update(string);
    let result = hasher.finalize();
    format!("{:x}", result)
}

/// Hashes a given String literal into a `Zq`.
///
/// Parameters:
/// - `string`: specifies the value that is hashed.
/// - `modulus`: specifies the modulus of the returned `Zq` value
///
/// # Examples
/// ```
/// use qfall_crypto::primitive::hash::hash_to_zq_sha256;
/// use qfall_math::{integer::Z, integer_mod_q::{Modulus, Zq}};
/// use std::str::FromStr;
///
/// let string = "Hello World!";
/// let modulus = &Modulus::from_str("7").unwrap();
///
/// let hash: Zq = hash_to_zq_sha256("Hello World!", &modulus);
/// assert_eq!(Zq::try_from((2,7)).unwrap(), hash)
/// ```
pub fn hash_to_zq_sha256(string: &str, modulus: &Modulus) -> Zq {
    let modulus_new: Z = modulus.into();
    let bitsize = (&modulus_new).bits();
    let mut hex = "".to_string();
    let string2 = modulus_new.to_string() + " " + string;

    for i in 0..(bitsize / 128 + 1)
    // hashing into e.g. Zq with 256 bit length of q from 256 bit will result in
    // lower values to be up to two times as likely as higher values.
    // Doubling the bit size of the hashed number will
    // reduce this difference to 1/2^n which is negligible.
    {
        hex = hex + &sha256(&(i.to_string() + " " + &string2));
    }

    Zq::try_from((Z::from_str_b(&hex, 16).unwrap(), modulus)).unwrap()
}

/// Hashes a given String literal into a [`MatZq`] .
///
/// Parameters:
/// - `string`: specifies the value that is hashed
/// - `num_rows`: specifies the number of rows of the result
/// - `num_cols`: specifies the number of columns of the result
/// - `modulus`: specifies the modulus of the returned [`MatZq`] value
///
/// # Examples
/// ```
/// use qfall_crypto::primitive::hash::hash_to_mat_zq_sha256;
/// use qfall_math::{integer::Z, integer_mod_q::{Modulus, MatZq}};
/// use std::str::FromStr;
///
/// let string = "Hello World!";
/// let modulus = &Modulus::from_str("7").unwrap();
///
/// let hash: MatZq = hash_to_mat_zq_sha256(string, 2, 2, &modulus);
/// assert_eq!(MatZq::from_str("[[6, 3],[5, 2]] mod 7").unwrap(), hash);
/// ```
///
/// # Panics ...
/// - ... if the number of rows or columns is less or equal to `0`.
pub fn hash_to_mat_zq_sha256(
    string: &str,
    num_rows: impl TryInto<i64> + Display + Into<i64>,
    num_cols: impl TryInto<i64> + Display + Into<i64>,
    modulus: &Modulus,
) -> MatZq {
    let num_rows_new: i64 = num_rows.into();
    let num_cols_new: i64 = num_cols.into();
    if num_cols_new <= 0 || num_rows_new <= 0 {
        panic!("The number of rows and number of columns must be at least one.");
    }
    let mut matrix = MatZq::new(num_rows_new, num_cols_new, modulus).unwrap();
    let new_string = format!(" {num_rows_new} {num_cols_new} {string}");
    for i in 0..num_rows_new {
        for j in 0..num_cols_new {
            matrix
                .set_entry(
                    i,
                    j,
                    hash_to_zq_sha256(&format!("{i} {j}{new_string}"), modulus),
                )
                .unwrap();
        }
    }
    matrix
}

#[cfg(test)]
mod tests_sha {
    use crate::primitive::hash::{hash_to_mat_zq_sha256, hash_to_zq_sha256, sha256, Z};
    use qfall_math::{
        integer_mod_q::{MatZq, Modulus, Zq},
        traits::Pow,
    };
    use std::str::FromStr;

    // ensure sha256 works
    #[test]
    fn test_sha256() {
        let str1 = "Hello World!";
        let str2 = "qfall";

        let hash1 = sha256(str1);
        let hash2 = sha256(str2);

        assert_eq!(
            "7f83b1657ff1fc53b92dc18148a1d65dfc2d4b1fa3d677284addd200126d9069",
            hash1
        );
        assert_eq!(
            "eb6ed1369a670050bd04b24036e8c29144b0f6b10166dc9c8b4987a6026c715f",
            hash2
        );
    }

    // ensure hashing into [`Zq`] works as intended
    #[test]
    fn test_hash_to_zq_sha256() {
        let str1 = "Hello World!";
        let str2 = "qfall";

        let hash1 = hash_to_zq_sha256(str1, &Modulus::from_str("256").unwrap());
        let hash2 = hash_to_zq_sha256(str2, &Modulus::from_str("16").unwrap());

        assert_eq!(Zq::try_from((150, 256)).unwrap(), hash1);
        assert_eq!(Zq::try_from((12, 16)).unwrap(), hash2);
    }

    // ensure hashing into [`Zq`] hits the whole domain not just the first 256 bit
    #[test]
    fn test_hash_to_zq_sha256_large() {
        let str1 = "Hello World!";

        let mut large = false;
        for i in 0..5 {
            if Z::from(hash_to_zq_sha256(
                &(i.to_string() + str1),
                &Modulus::try_from_z(&Z::from(271).pow(100).unwrap()).unwrap(),
            )) > Z::from(u64::MAX)
            {
                large = true;
            }
        }

        assert!(large);
    }

    // ensure hashing into [`MatZq`] works as intended
    #[test]
    fn test_hash_to_mat_zq_sha256() {
        let str1 = "Hello World!";
        let str2 = "qfall";

        let hash1 = hash_to_mat_zq_sha256(str1, 2, 2, &Modulus::from_str("256").unwrap());
        let hash2 = hash_to_mat_zq_sha256(str2, 2, 2, &Modulus::from_str("16").unwrap());

        assert_eq!(
            MatZq::from_str("[[159, 26],[249, 141]] mod 256").unwrap(),
            hash1
        );
        assert_eq!(MatZq::from_str("[[3, 12],[9, 12]] mod 16").unwrap(), hash2);
    }

    // ensure hashing into [`MatZq`] works as intended
    #[test]
    #[should_panic]
    fn test_hash_to_mat_zq_sha256_negative_dimensions() {
        let str1 = "Hello World!";

        let _ = hash_to_mat_zq_sha256(str1, 0, 0, &Modulus::from_str("256").unwrap());
    }
}
