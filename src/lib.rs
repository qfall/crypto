// Copyright © 2023 Niklas Siemer, Marvin Beckmann
//
// This file is part of qFALL-crypto.
//
// qFALL-crypto is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

//! # What is qFALL-crypto?
//! qFall-crypto provides cryptographic basics such as mathematical primitives,
//! fundamental lattice-based cryptographic constructions, and samplable distributions/
//! possibilities to sample instances of lattice problems to prototype
//! lattice-based cryptographic constructions and more.
//!
//! Currently qFALL-crypto supports 3 main construction types:
//! - [Identity-Based Encryptions](construction::identity_based_encryption::IBEScheme)
//! - [Public-Key Encryptions](construction::pk_encryption::PKEncryptionScheme)
//! - [Signatures](construction::signature::SignatureScheme)
//!
//! These are identified by traits and then implemented for specific constructions, e.g.
//! [`RingLPR`](construction::pk_encryption::RingLPR).
//! Our library has further primitives useful for prototyping such as
//! [`PSFs`](primitive::psf::PSF) that can be used to implement constructions.
//!
//! qFALL-crypto is free software: you can redistribute it and/or modify it under
//! the terms of the Mozilla Public License Version 2.0 as published by the
//! Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.
//!
//! ## Tutorial + Website
//! You can find a dedicated [tutorial](https://qfall.github.io/book/index.html) to qFALL-crypto on our [website](https://qfall.github.io/).
//! The tutorial explains the basic steps starting from installation and
//! continues with basic usage.
//! qFALL-crypto is co-developed together with qFALL-math which provides the basic
//! foundation that is used to implement the cryptographic constructions.

pub mod construction;
pub mod primitive;
pub mod sample;
pub mod utils;
