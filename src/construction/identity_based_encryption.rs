// Copyright © 2023 Phil Milewski
//
// This file is part of qFALL-crypto.
//
// qFALL-crypto is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

//! This module provides the trait a struct should implement if it is an
//! instance of a identity based public key encryption scheme. Furthermore,
//! it contains cryptographic schemes implementing the `IBE` trait.
//! The main references are listed in the following
//! and will be further referenced in submodules by these numbers:
//! - \[1\] Gentry, Craig and Peikert, Chris and Vaikuntanathan, Vinod (2008).
//! Trapdoors for hard lattices and new cryptographic constructions.
//! In: Proceedings of the fortieth annual ACM symposium on Theory of computing.
//! <https://dl.acm.org/doi/pdf/10.1145/1374376.1374407>
//! - \[2\] Regev, Oded (2009).
//! On lattices, learning with errors, random linear codes, and cryptography.
//! In: Journal of the ACM 6.
//! <https://dl.acm.org/doi/pdf/10.1145/1568318.1568324>

mod dual_regev_ibe;
pub use dual_regev_ibe::DualRegevIBE;

use qfall_math::integer::Z;

pub trait IBE {
    type MasterPublicKey;
    type MasterSecretKey;
    type SecretKey;
    type Cipher;
    type Identity;

    // algorithm to initialize an [`IBE`] scheme
    fn setup(&self) -> (Self::MasterPublicKey, Self::MasterSecretKey);

    // algorithm used do extract an identity specific secret key
    fn extract(
        &mut self,
        pk: &Self::MasterPublicKey,
        sk: &Self::MasterSecretKey,
        identity: &Self::Identity,
    ) -> Self::SecretKey;

    /// Encrypts the provided `message` using the master public key `mpk` and `identity` of the recipient.
    ///
    /// Parameters:
    /// - `mpk`: specifies the master public key used for this IBE
    /// - `identity`: specifies the recipient that should be able to decrypt the encrypted message
    /// - `message`: specifies the message to be encrypted
    ///
    /// Returns the encryption of `message` as a [`Self::Cipher`] instance.
    fn enc(
        &self,
        mpk: &Self::MasterPublicKey,
        identity: &Self::Identity,
        message: impl Into<Z>,
    ) -> Self::Cipher;

    /// Decrypts the provided `cipher` using the extracted secret key `sk`.
    ///
    /// Parameters:
    /// - `sk`: specifies the extracted secret key used for decryption
    /// - `cipher`: specifies the ciphertext to be decrypted
    ///
    /// Returns the decryption of `cipher` as a [`Z`] instance.
    fn dec(&self, sk: &Self::SecretKey, cipher: &Self::Cipher) -> Z;
}
