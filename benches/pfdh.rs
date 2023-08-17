// Copyright © 2023 Marvin Beckmann
//
// This file is part of qFALL-crypto.
//
// qFALL-crypto is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

use criterion::{criterion_group, Criterion};
use qfall_crypto::construction::signature::{pfdh::Pfdh, SignatureScheme};

/// Performs a full instantiation with an additional signing and verifying of a signature.
fn pfdh_cycle(n: i64) {
    let mut pfdh = Pfdh::init_gpv(n, 113, 17, 128);

    let m = "Hello World!";

    let (pk, sk) = pfdh.gen();
    let sigma = pfdh.sign(m.to_owned(), &sk, &pk);

    pfdh.vfy(m.to_owned(), &sigma, &pk);
}

/// Benchmark [pfdh_cycle] with `n = 8`.
///
/// This benchmark can be run with for example:
/// - `cargo criterion PFDH\ n=8`
/// - `cargo bench --bench benchmarks PFDH\ n=8`
/// - `cargo flamegraph --bench benchmarks -- --bench PFDH\ n=8`
///
/// Shorter variants or regex expressions can also be used to specify the
/// benchmark name. The `\ ` is used to escape the space, alternatively,
/// quotation marks can be used.
fn bench_pfdh_cycle(c: &mut Criterion) {
    c.bench_function("PFDH n=8", |b| b.iter(|| pfdh_cycle(8)));
}

criterion_group!(benches, bench_pfdh_cycle);
