// Copyright © 2023 Sven Moog
//
// This file is part of qFALL-crypto.
//
// qFALL-crypto is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.
//! This file collects the benchmarks from other files.

use criterion::criterion_main;

pub mod pfdh;
pub mod regev;

criterion_main! {regev::benches, pfdh::benches}
