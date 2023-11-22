// Copyright © 2023 Sven Moog
//
// This file is part of qFALL-crypto.
//
// qFALL-crypto is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.
//! This file collects the benchmarks from other files.

use criterion::criterion_main;

/// This `enum` is used to decide whether a small or large
/// security parameter `n` is used for the benchmark.
#[derive(PartialEq)]
pub(crate) enum SizeN {
    Small,
    Large,
}

mod pfdh;
mod pk_encryption;

criterion_main! {pk_encryption::regev::benches, pk_encryption::ring_lpr::benches, pfdh::benches}
