// Copyright © 2023 Niklas Siemer, Marvin Beckmann
//
// This file is part of qFALL-crypto.
//
// qFALL-crypto is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

//! This module contains fundamental cryptographic constructions, on which other
//! constructions can be build on.
//! Among others these include encryption schemes and signature schemes.
//! A construction is build the same way:
//!
//! 1. A trait that combines the common feature, e.g.
//! [`public key encryption`](pk_encryption::PKEncryption).
//! 2. Explicit implementations of the trait, e.g.
//! [`RingLPR`](pk_encryption::RingLPR).

pub mod hash;
pub mod identity_based_encryption;
pub mod pk_encryption;
pub mod signature;
