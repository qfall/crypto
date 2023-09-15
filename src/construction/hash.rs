// Copyright © 2023 Niklas Siemer
//
// This file is part of qFALL-crypto.
//
// qFALL-crypto is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

//! This module contains implementations of hash functions.

pub mod sha256;
mod sis;

pub use sis::SISHash;

/// This trait should be implemented by hashes with domain [`str`].
pub trait HashInto<DigestSpace> {
    /// Hashes a given String literal.
    ///
    /// Paramters:
    /// - `m`: specifies the string message to be hashed
    ///
    /// Returns a hash of type Domain.
    fn hash(&self, m: &str) -> DigestSpace;
}
