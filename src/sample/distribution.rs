// Copyright © 2023 Niklas Siemer
//
// This file is part of qFALL-crypto.
//
// qFALL-crypto is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

//! This module collects different distributions and implementations to draw
//! samples according to a specified distribution.
//! - \[1\] Gentry, Craig, Chris Peikert, and Vinod Vaikuntanathan.
//! "Trapdoors for hard lattices and new cryptographic constructions."
//! Proceedings of the fortieth annual ACM symposium on Theory of computing. 2008.
//! <https://doi.org/10.1145/1374376.1374407>

pub mod binomial;
pub mod discrete_gauss;
pub mod psf;
pub mod uniform;
