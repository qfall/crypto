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
//!
//! - \[1\] Peikert, Chris (2016).
//! A decade of lattice cryptography.
//! In: Theoretical Computer Science 10.4.
//! <https://web.eecs.umich.edu/~cpeikert/pubs/lattice-survey.pdf>
//! - \[3\] Regev, Oded (2009).
//! On lattices, learning with errors, random linear codes, and cryptography.
//! In: Journal of the ACM 6.
//! <https://dl.acm.org/doi/pdf/10.1145/1568318.1568324>

mod dual_regev;
mod regev;
pub use dual_regev::DualRegev;
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
