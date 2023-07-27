// Copyright © 2023 Marvin Beckmann
//
// This file is part of qFALL-crypto.
//
// qFALL-crypto is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

//! The main references are listed in the following
//! and will be further referenced in submodules by these numbers:
//! - \[1\] Micciancio, D., Peikert, C. (2012).
//! Trapdoors for Lattices: Simpler, Tighter, Faster, Smaller.
//! In: Pointcheval, D., Johansson, T. (eds) Advances in Cryptology – EUROCRYPT 2012.
//! EUROCRYPT 2012. Lecture Notes in Computer Science, vol 7237.
//! Springer, Berlin, Heidelberg. <https://doi.org/10.1007/978-3-642-29011-4_41>
//! - \[2\] Cash, D., Hofheinz, D., Kiltz, E., & Peikert, C. (2012).
//! Bonsai trees, or how to delegate a lattice basis. Journal of cryptology, 25, 601-639.
//! <https://doi.org/10.1007/s00145-011-9105-2>

use self::{gadget_classical::gen_trapdoor, gadget_parameters::GadgetParameters};
use qfall_math::{
    integer::{MatZ, Z},
    integer_mod_q::{MatZq, Modulus},
};

pub mod gadget_classical;
pub mod gadget_parameters;
pub mod short_basis_classical;
pub mod trapdoor_distribution;

/// Computes a trapdoor with default values.
///
/// - `params` is computed using [`GadgetParameters::init_default`].
/// - `tag = I_n` is taken from [\[1\]](<index.html#:~:text=[1]>): Algorithm 1
///
/// Parameters:
/// - `n`: the security parameter
/// - `modulus`: the modulus for the trapdoor
///
/// Returns a matrix `a` and its gadget-trapdoor `r` as in [\[1\]](<index.html#:~:text=[1]>): Algorithm 1 for some fixed set of parameters [`GadgetParameters::init_default`].
///
/// # Examples
/// ```
/// use qfall_crypto::sample::g_trapdoor::gen_trapdoor_default;
/// use qfall_math::integer::Z;
/// use qfall_math::integer_mod_q::Modulus;
///
/// let (a,r) = gen_trapdoor_default(42, &Modulus::from(101));
/// ```
///
/// # Panics ...
/// - if the security parameter `n` is not in `\[1, i64::MAX\]`.
pub fn gen_trapdoor_default(n: impl Into<Z>, modulus: &Modulus) -> (MatZq, MatZ) {
    // panic if n < 1 (security parameter must be positive)
    let n = n.into();
    assert!(n >= Z::ONE);

    let params = GadgetParameters::init_default(n, modulus);

    // a_bar <-$ Z_q^{n * m_bar}
    let a_bar = MatZq::sample_uniform(&params.n, &params.m_bar, &params.q);

    // tag = I_n
    let tag = MatZq::identity(&params.n, &params.n, &params.q);

    // we can unwrap, as we compute the parameters on our own and
    // they should always work
    gen_trapdoor(&params, &a_bar, &tag).unwrap()
}

#[cfg(test)]
mod test_gen_trapdoor_default {
    use super::gen_trapdoor_default;
    use qfall_math::integer_mod_q::Modulus;

    #[test]
    fn working() {
        let (_, _) = gen_trapdoor_default(100, &Modulus::from(32));
    }
}
