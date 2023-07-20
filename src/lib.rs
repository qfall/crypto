// Copyright © 2023 Niklas Siemer
//
// This file is part of qFALL-crypto.
//
// qFALL-crypto is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

//! This crate provides basics such as mathematical primitives, fundamental
//! lattice-based cryptographic constructions, and samplable distributions/
//! possibilities to sample instances of lattice problems to prototype
//! lattice-based cryptographic constructions and more.

pub mod construction;
pub mod primitive;
pub mod sample;
pub mod utils;
