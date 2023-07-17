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

mod dual_regev_ibe;
pub use dual_regev_ibe::DualRegevIBE;

use qfall_math::integer::Z;

pub trait IBE {
    type MasterPublicKey;
    type MasterSecretKey;
    type SecretKey;
    type Cipher;
    type Identity;

    fn gen(&self) -> (Self::MasterPublicKey, Self::MasterSecretKey);
    fn extract(
        &mut self,
        pk: &Self::MasterPublicKey,
        sk: &Self::MasterSecretKey,
        identity: &Self::Identity,
    ) -> Self::SecretKey;
    fn enc(
        &self,
        pk: &Self::MasterPublicKey,
        identity: &Self::Identity,
        message: impl Into<Z>,
    ) -> Self::Cipher;
    fn dec(&self, sk: &Self::SecretKey, cipher: &Self::Cipher) -> Z;
}
