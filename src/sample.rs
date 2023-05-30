// Copyright © 2023 Niklas Siemer
//
// This file is part of qFALL-crypto.
//
// qFALL-crypto is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

//! This module contains anything that should be easily samplable for lattice-based
//! cryptography. This includes distributions like `DiscreteGauss`, instances of
//! lattice problems like `LWE`, and trapdoors.

pub mod distribution;
pub mod g_trapdoor;
pub mod lattice_problem;
