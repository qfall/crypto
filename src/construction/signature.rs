// Copyright © 2023 Marvin Beckmann
//
// This file is part of qFALL-crypto.
//
// qFALL-crypto is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

//! This module provides the trait a struct should implement if it is an
//! instance of a signature scheme. Furthermore, it contains cryptographic signatures
//! implementing the `Signature` trait.

pub mod fdh;
mod pfdh;

/// This trait captures the essential functionalities each signature scheme has to support.
/// These include
/// - `gen`: to generate a public key and secret key pair
/// - `sign`: which allows to create a signature for a message
/// - `vfy`: which allows to check if a signature is valid for a certain message
///
/// Note: The gen does not take in the parameter `1^n`, as this is a public parameter,
/// which shall be defined by the struct implementing this trait.
pub trait SignatureScheme {
    type SecretKey;
    type PublicKey;
    type Signature;

    fn gen(&mut self) -> (Self::PublicKey, Self::SecretKey);
    fn sign(&mut self, m: String, sk: &Self::SecretKey, pk: &Self::PublicKey) -> Self::Signature;
    fn vfy(&self, m: String, sigma: &Self::Signature, pk: &Self::PublicKey) -> bool;
}
