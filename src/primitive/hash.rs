// Copyright © 2023 Phil Milewski
//
// This file is part of qFALL-crypto.
//
// qFALL-crypto is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

//! This module contains hashes into different domains.

use std::str::FromStr;

use qfall_math::{
    integer::Z,
    integer_mod_q::{MatZq, Zq},
};
use sha2::{Digest, Sha256};

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
/// ```
#[allow(dead_code)]
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
/// use qfall_math::{integer::Z, integer_mod_q::Zq};
///
/// let string = "Hello World!";
/// let modulus = Z::from(11);
/// let hash: Zq = hash_to_zq_sha256(string, &modulus);
/// ```
#[allow(dead_code)]
pub fn hash_to_zq_sha256(string: &str, modulus: &Z) -> Zq {
    let bitsize = modulus.bits();
    let mut hex = "".to_string();
    let string2 = modulus.to_string() + string;
    for i in 0..(bitsize / 128 + 1)
    // hashing into e.g. Zq with 256 bit length of q from 256 bit will result in
    // lower values to be up to two times as likely as higher values.
    // Doubling the bit size of the hashed number will
    // reduce this difference to 1/2^n which is negligible.
    {
        hex = hex + &sha256(&(i.to_string() + " " + &string2));
    }
    Zq::try_from_z_z(&Z::from_str_b(&hex, 16).unwrap(), &modulus).unwrap()
}

/// Hashes a given String literal into a `MatZq` column vector.
///
/// Parameters:
/// - `string`: specifies the value that is hashed.
/// - `modulus`: specifies the modulus of the returned `MatZq` value
/// - `dimension`: specifies the number of entries of the vector.
///
/// # Examples
/// ```
/// use qfall_crypto::primitive::hash::hash_to_vec_zq_sha256;
/// use qfall_math::{integer::Z, integer_mod_q::MatZq};
///
/// let string = "Hello World!";
/// let modulus = Z::from(11);
/// let hash: MatZq = hash_to_vec_zq_sha256(string, &modulus, 10);
/// ```
#[allow(dead_code)]
pub fn hash_to_vec_zq_sha256(string: &str, modulus: &Z, dimension: u64) -> MatZq {
    let mut matrix = "[[".to_string();
    for i in 0..dimension - 1 {
        matrix = matrix
            + &Z::from(hash_to_zq_sha256(
                &(i.to_string() + " " + &dimension.to_string() + " " + string),
                &modulus,
            ))
            .to_string()
            + "],[";
    }
    matrix = matrix
        + &Z::from(hash_to_zq_sha256(
            &(dimension.to_string() + " " + string),
            &modulus,
        ))
        .to_string()
        + "]]";

    MatZq::from_str(&format!("{} mod {}", matrix, modulus)).unwrap()
}

#[cfg(test)]
mod tests_sha {
    use crate::primitive::hash::sha256;

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
}
