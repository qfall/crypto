// Copyright © 2023 Sven Moog
//
// This file is part of qFALL-crypto.
//
// qFALL-crypto is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

use criterion::*;
use qfall_crypto::construction::pk_encryption::PKEncryption;
use qfall_crypto::construction::pk_encryption::Regev;
use qfall_math::integer::Z;

/// Performs a full-cycle of gen, enc, dec with regev.
fn regev_cycle(n: i64) {
    let msg = Z::ONE;
    let regev = Regev::new_from_n(n).unwrap();

    let (pk, sk) = regev.gen();
    let cipher = regev.enc(&pk, &msg);
    let _ = regev.dec(&sk, &cipher);
}

/// Benchmark [regev_cycle] with `n = 50`.
///
/// This benchmark can be run with for example:
/// - `cargo criterion Regev\ n=50`
/// - `cargo bench --bench benchmarks Regev\ n=50`
/// - `cargo flamegraph --bench benchmarks -- --bench Regev\ n=50`
///
/// Shorter variants or regex expressions can also be used to specify the
/// benchmark name. The `\ ` is used to escape the space, alternatively,
/// quotation marks can be used.
fn bench_regev_cycle(c: &mut Criterion) {
    c.bench_function("Regev n=50", |b| b.iter(|| regev_cycle(50)));
}

/// Benchmark [regev_cycle] with `n = 10, 20, 30, 40, 50, 60`
///
/// This benchmark can be run with for example:
/// - `cargo criterion "Regev\ n\ sweep"`
/// - `cargo criterion Regev\ n\ sweep/n=20` (only run the n=20 benchmark).
/// - `cargo criterion 'Regev.*n=20'` (only run the n=20 benchmark).
/// - `cargo bench --bench benchmarks Regev\ n\ sweep`
///
/// Shorter variants or regex expressions can also be used to specify the
/// benchmark name. The `\ ` is used to escape the space, alternatively,
/// quotation marks can be used.
fn bench_regev_cycle_n_sweep(c: &mut Criterion) {
    let mut group = c.benchmark_group("Regev n sweep");

    for n in [10, 20, 30, 40, 50, 60].iter() {
        group.bench_function(format!("n={n}"), |b| b.iter(|| regev_cycle(*n)));
    }

    group.finish();
}

criterion_group!(benches, bench_regev_cycle, bench_regev_cycle_n_sweep);
