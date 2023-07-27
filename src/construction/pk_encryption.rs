// Copyright © 2023 Niklas Siemer
//
// This file is part of qFALL-crypto.
//
// qFALL-crypto is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

//! This module provides the trait a struct should implement if it is an
//! instance of a public key encryption scheme. Furthermore, it contains
//! cryptographic schemes implementing the `PKEncryption` trait.

mod dual_regev_discrete_gauss;
mod regev;
pub use dual_regev_discrete_gauss::DualRegevWithDiscreteGaussianRegularity;
pub use regev::Regev;

use qfall_math::integer::Z;

pub trait PKEncryption {
    type PublicKey;
    type SecretKey;
    type Cipher;

    fn gen(&self) -> (Self::PublicKey, Self::SecretKey);
    fn enc(&self, pk: &Self::PublicKey, message: impl Into<Z>) -> Self::Cipher;
    fn dec(&self, sk: &Self::SecretKey, cipher: &Self::Cipher) -> Z;
}
