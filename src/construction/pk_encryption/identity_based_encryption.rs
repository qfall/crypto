// Copyright © 2023 Phil Milewski
//
// This file is part of qFALL-crypto.
//
// qFALL-crypto is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

//! This module contains an implementation of the IND-CPA secure
//! identity based public key encryption scheme.
//!
//! The main references are listed in the following:
//! - \[1\] Gentry, Craig and Peikert, Chris and Vaikuntanathan, Vinod (2008).
//! Trapdoors for hard lattices and new cryptographic constructions.
//! In: Proceedings of the fortieth annual ACM symposium on Theory of computing.
//! <https://dl.acm.org/doi/pdf/10.1145/1374376.1374407>

use super::PKEncryption;
use serde::{Deserialize, Serialize};

/// This struct manages and stores the public parameters of a [`IBE`]
/// public key encryption instance.
///
/// Attributes:
/// - `n`: specifies the security parameter, which is not equal to the bit-security level
/// - `m`: defines the dimension of the underlying lattice
/// - `q`: specifies the modulus over which the encryption is computed
/// - `r`: specifies the gaussian parameter used for SampleD,
///   i.e. used for encryption
/// - `alpha`: specifies the gaussian parameter used for independent
///   sampling from χ, i.e. for multiple discrete Gaussian samples used
///   for key generation
///todo
/// # Examples
/// ```
/// use qfall_crypto::construction::pk_encryption::{DualRegev, PKEncryption};
/// use qfall_math::integer::Z;
/// // setup public parameters and key pair
/// let dual_regev = DualRegev::default();
/// let (pk, sk) = dual_regev.gen();
///
/// // encrypt a bit
/// let msg = Z::ZERO; // must be a bit, i.e. msg = 0 or 1
/// let cipher = dual_regev.enc(&pk, &msg);
///
/// // decrypt
/// let m = dual_regev.dec(&sk, &cipher);
///
/// assert_eq!(msg, m);
/// ```
#[derive(Debug, Serialize, Deserialize)]
pub struct IBE {
    n: Z,       // security parameter
    m: Z,       // number of rows of matrix A
    q: Modulus, // modulus
    r: Q,       // gaussian parameter for sampleD
    alpha: Q,   // gaussian parameter for sampleZ

    storage: HashMap<String, Domain>,
}
