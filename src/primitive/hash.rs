// Copyright © 2023 Phil Milewski
//
// This file is part of qFALL-crypto.
//
// qFALL-crypto is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

//! This module contains hashes into different domains.

use qfall_math::utils::index::evaluate_indices;
use qfall_math::{
    integer::{MatPolyOverZ, Z},
    integer_mod_q::{MatPolynomialRingZq, MatZq, Modulus, ModulusPolynomialRingZq, Zq},
    traits::SetEntry,
};
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::fmt::Display;

/// Computes the sha256 hash value of a given String literal.
///
/// Parameters:
/// - `string`: specifies the value that is hashed.
///
/// Returns the sha256 value of the given string as a hex string.
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
/// Returns a [`Zq`] as a hash value for the given string.
///
///  # Examples
/// ```
/// use qfall_crypto::primitive::hash::hash_to_zq_sha256;
/// use qfall_math::{integer::Z, integer_mod_q::{Modulus, Zq}};
///
/// let string = "Hello World!";
/// let modulus = Modulus::from(7);
///
/// let hash: Zq = hash_to_zq_sha256("Hello World!", &modulus);
/// assert_eq!(Zq::from((2, 7)), hash)
/// ```
pub fn hash_to_zq_sha256(string: &str, modulus: impl Into<Modulus>) -> Zq {
    let modulus = modulus.into();
    let modulus_new: Z = (&modulus).into();
    let bitsize = modulus_new.bits();
    let mut hex = "".to_string();
    let string2 = format!("{modulus_new} {string}");

    for i in 0..=bitsize / 128
    // hashing into e.g. Zq with 256 bit length of q from 256 bit will result in
    // lower values to be up to two times as likely as higher values.
    // Doubling the bit size of the hashed number will
    // reduce this difference to 1/2^n which is negligible.
    // https://crypto.stackexchange.com/questions/37305/how-can-i-instantiate-a-generalized-hash-function
    {
        hex = hex + &sha256(&format!("{i} {string2}"));
    }

    Zq::from((Z::from_str_b(&hex, 16).unwrap(), modulus))
}

/// Hashes a given String literal into a [`MatZq`] .
///
/// Parameters:
/// - `string`: specifies the value that is hashed
/// - `num_rows`: specifies the number of rows of the result
/// - `num_cols`: specifies the number of columns of the result
/// - `modulus`: specifies the modulus of the returned [`MatZq`] value
///
/// Returns a [`MatZq`] as a hash for the given string.
///
/// # Examples
/// ```
/// use qfall_crypto::primitive::hash::hash_to_mat_zq_sha256;
/// use qfall_math::{integer::Z, integer_mod_q::{Modulus, MatZq}};
/// use std::str::FromStr;
///
/// let string = "Hello World!";
/// let modulus = Modulus::from(7);
///
/// let hash: MatZq = hash_to_mat_zq_sha256(string, 2, 2, &modulus);
/// assert_eq!(MatZq::from_str("[[6, 3],[5, 2]] mod 7").unwrap(), hash);
/// ```
///
/// # Panics ...
/// - if the number of rows or columns is less or equal to `0` or does not fit into an [`i64`].
pub fn hash_to_mat_zq_sha256(
    string: &str,
    num_rows: impl TryInto<i64> + Display,
    num_cols: impl TryInto<i64> + Display,
    modulus: impl Into<Modulus>,
) -> MatZq {
    let modulus = modulus.into();
    let (num_rows_new, num_cols_new) = evaluate_indices(num_rows, num_cols).unwrap();
    let mut matrix = MatZq::new(num_rows_new, num_cols_new, modulus.clone());

    let new_string = format!("{num_rows_new} {num_cols_new} {string}");
    for i in 0..num_rows_new {
        for j in 0..num_cols_new {
            matrix
                .set_entry(
                    i,
                    j,
                    hash_to_zq_sha256(&format!("{i} {j} {new_string}"), &modulus),
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
        traits::{Distance, Pow},
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

        let hash1 = hash_to_zq_sha256(str1, &Modulus::from(256));
        let hash2 = hash_to_zq_sha256(str2, &Modulus::from(16));

        assert_eq!(Zq::from((150, 256)), hash1);
        assert_eq!(Zq::from((12, 16)), hash2);
    }

    // ensure hashing into [`Zq`] hits the whole domain not just the first 256 bit
    #[test]
    fn test_hash_to_zq_sha256_large() {
        let str1 = "Hello World!";

        let mut large = false;
        for i in 0..5 {
            if hash_to_zq_sha256(
                &(i.to_string() + str1),
                &Modulus::from(&Z::from(271).pow(100).unwrap()),
            )
            .distance(Z::ZERO)
                > Z::from(u64::MAX)
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

        let hash1 = hash_to_mat_zq_sha256(str1, 2, 2, &Modulus::from(256));
        let hash2 = hash_to_mat_zq_sha256(str2, 2, 2, &Modulus::from(16));

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

        let _ = hash_to_mat_zq_sha256(str1, 0, 0, &Modulus::from(16));
    }
}

/// Hash object to hash a String into a [`MatZq`].
/// The object fixes the modulus and the corresponding dimensions.
///
/// Parameters:
/// - `modulus`: Defines the range in which each entry is hashed
/// - `rows`: Defines the number of rows of the hash value
/// - `cols`: Defines the number of columns of the hash value
///
/// Returns a [`MatZq`] as a hash for the given string.
///
/// # Examples
/// ```
/// use qfall_crypto::primitive::hash::{HashMatZq, HashInto};
/// use qfall_crypto::primitive::hash::hash_to_mat_zq_sha256;
/// use qfall_math::{integer::Z, integer_mod_q::{Modulus, MatZq}};
/// use std::str::FromStr;
///
/// let modulus = Modulus::from(7);
///
/// let hasher = HashMatZq {
///     modulus: modulus,
///     rows: 17,
///     cols: 3,
/// };
/// let hash_val = hasher.hash("Hello");
/// ```
#[derive(Serialize, Deserialize)]
pub struct HashMatZq {
    pub modulus: Modulus,
    pub rows: i64,
    pub cols: i64,
}
pub trait HashInto<Domain> {
    fn hash(&self, m: &str) -> Domain;
}
impl HashInto<MatZq> for HashMatZq {
    /// Hashes a given String literal into a [`MatZq`].
    /// The dimensions and the modulus is fixed by the hash object.
    ///
    /// Parameters:
    /// - `string`: specifies the value that is hashed
    ///
    /// Returns a [`MatZq`] as a hash for the given string.
    ///
    /// # Examples
    /// ```
    /// use qfall_crypto::primitive::hash::{HashMatZq, HashInto};
    /// use qfall_crypto::primitive::hash::hash_to_mat_zq_sha256;
    /// use qfall_math::{integer::Z, integer_mod_q::{Modulus, MatZq}};
    /// use std::str::FromStr;
    ///
    /// let modulus = Modulus::from(7);
    ///
    /// let hasher = HashMatZq {
    ///     modulus: modulus,
    ///     rows: 17,
    ///     cols: 3,
    /// };
    /// let hash_val = hasher.hash("Hello");
    /// ```
    fn hash(&self, m: &str) -> MatZq {
        hash_to_mat_zq_sha256(m, self.rows, self.cols, &self.modulus)
    }
}

/// Hash object to hash a String into a [`MatPolynomialRingZq`].
/// The object fixes the modulus and the corresponding dimensions.
///
/// Parameters:
/// - `modulus`: Defines the range in which each entry is hashed
/// - `rows`: Defines the number of rows of the hash value
/// - `cols`: Defines the number of columns of the hash value
///
/// Returns a [`MatZq`] as a hash for the given string.
///
/// # Examples
/// ```
/// use qfall_crypto::primitive::hash::{HashMatPolynomialRingZq, HashInto};
/// use qfall_crypto::sample::g_trapdoor::gadget_parameters::GadgetParametersRing;
/// use qfall_math::{
///     integer_mod_q::Modulus,
///     traits::{GetNumColumns, GetNumRows},
/// };
///
/// let gp = GadgetParametersRing::init_default(10, &Modulus::from(99));
///
/// let hasher = HashMatPolynomialRingZq {
///     modulus: gp.modulus,
///     rows: 17,
///     cols: 3,
/// };
/// let hash_val = hasher.hash("Hello");
/// ```
#[derive(Serialize, Deserialize)]
pub struct HashMatPolynomialRingZq {
    pub modulus: ModulusPolynomialRingZq,
    pub rows: i64,
    pub cols: i64,
}
impl HashInto<MatPolynomialRingZq> for HashMatPolynomialRingZq {
    /// Hashes a given String literal into a [`MatPolynomialRingZq`].
    /// The dimensions and the modulus is fixed by the hash object.
    ///
    /// Parameters:
    /// - `string`: specifies the value that is hashed
    ///
    /// Returns a [`MatPolynomialRingZq`] as a hash for the given string.
    ///
    /// # Examples
    /// ```
    /// use qfall_crypto::primitive::hash::{HashMatPolynomialRingZq, HashInto};
    /// use qfall_crypto::sample::g_trapdoor::gadget_parameters::GadgetParametersRing;
    /// use qfall_math::{
    ///     integer_mod_q::Modulus,
    ///     traits::{GetNumColumns, GetNumRows},
    /// };
    ///
    /// let gp = GadgetParametersRing::init_default(10, &Modulus::from(99));
    ///
    /// let hasher = HashMatPolynomialRingZq {
    ///     modulus: gp.modulus,
    ///     rows: 17,
    ///     cols: 3,
    /// };
    /// let hash_val = hasher.hash("Hello");
    /// ```
    fn hash(&self, m: &str) -> MatPolynomialRingZq {
        let highest_deg = self.modulus.get_degree();
        let embedding =
            (&hash_to_mat_zq_sha256(m, self.rows * highest_deg, self.cols, &self.modulus.get_q()))
                .into();
        let poly_mat = MatPolyOverZ::from_coefficient_embedding_to_matrix(&embedding, highest_deg);
        MatPolynomialRingZq::from((&poly_mat, &self.modulus))
    }
}

#[cfg(test)]
mod hash_into_mat_polynomial_ring_zq {
    use super::{HashInto, HashMatPolynomialRingZq};
    use crate::sample::g_trapdoor::gadget_parameters::GadgetParametersRing;
    use qfall_math::{
        integer::PolyOverZ,
        integer_mod_q::Modulus,
        traits::{GetEntry, GetNumColumns, GetNumRows},
    };

    /// Ensure that the hash function maps into the correct dimension and it is also
    /// static, i.e. the same value is returned, when the same value  is hashed.
    #[test]
    fn correct_dimensions() {
        let gp = GadgetParametersRing::init_default(10, &Modulus::from(99));

        let hasher = HashMatPolynomialRingZq {
            modulus: gp.modulus,
            rows: 17,
            cols: 3,
        };
        let hash_val = hasher.hash("Hello");
        let hash_val_2 = hasher.hash("Hello");
        let entry: PolyOverZ = hash_val.get_entry(0, 0).unwrap();

        assert_eq!(hasher.rows, hash_val.get_num_rows());
        assert_eq!(hasher.cols, hash_val.get_num_columns());
        assert_eq!(hash_val, hash_val_2);
        assert_eq!(9, entry.get_degree())
    }
}
