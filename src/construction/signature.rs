// Copyright © 2023 Marcel Luca Schmidt, Marvin Beckmann
//
// This file is part of qFALL-crypto.
//
// qFALL-crypto is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

//! This module provides the trait a struct should implement if it is an
//! instance of a signature scheme. Furthermore, it contains cryptographic signatures
//! implementing the `Signature` trait.
//! - \[1\] Gentry, Craig, Chris Peikert, and Vinod Vaikuntanathan.
//! "Trapdoors for hard lattices and new cryptographic constructions."
//! Proceedings of the fortieth annual ACM symposium on Theory of computing. 2008.
//! <https://doi.org/10.1145/1374376.1374407>

mod fdh;
mod pfdh;

pub use fdh::FDH;
pub use pfdh::PFDH;

/// This trait should be implemented by every signature scheme.
/// It captures the essential functionalities each signature scheme has to support.
///
/// Note: The gen does not take in the parameter `1^n`, as this is a public parameter,
/// which shall be defined by the struct implementing this trait.
pub trait SignatureScheme {
    /// The type of the secret key.
    type SecretKey;
    /// The type of the public key.
    type PublicKey;
    /// The type of the signature.
    type Signature;

    /// Generates a public key and a secret key from the attributes the
    /// struct has, which implements this trait.
    ///
    /// Returns the public key and the secret key.
    fn gen(&mut self) -> (Self::PublicKey, Self::SecretKey);

    /// Signs a message using the secret key (and potentially the public key).
    ///
    /// Returns the resulting signature.
    fn sign(&mut self, m: String, sk: &Self::SecretKey, pk: &Self::PublicKey) -> Self::Signature;

    /// Verifies that a signature is valid for a message by using the public key.
    ///
    /// Returns the result of the verification as a boolean.
    fn vfy(&self, m: String, sigma: &Self::Signature, pk: &Self::PublicKey) -> bool;
}
