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
//! - \[2\] El Bansarkhani, R., Buchmann, J. (2014). Improvement and Efficient
//! Implementation of a Lattice-Based Signature Scheme. In: Lange, T., Lauter, K.,
//! Lisoněk, P. (eds) Selected Areas in Cryptography -- SAC 2013. SAC 2013. Lecture Notes
//! in Computer Science(), vol 8282. Springer, Berlin, Heidelberg.
//! <https://doi.org/10.1007/978-3-662-43414-7_3>
//! - \[3\] Gür, K.D., Polyakov, Y., Rohloff, K., Ryan, G.W. and Savas, E., 2018,
//! January. Implementation and evaluation of improved gaussian sampling for lattice
//!  trapdoors. In Proceedings of the 6th Workshop on Encrypted Computing & Applied
//! Homomorphic Cryptography (pp. 61-71). <https://dl.acm.org/doi/pdf/10.1145/3267973.3267975>

use crate::sample::g_trapdoor::{
    gadget_parameters::GadgetParametersRing, gadget_ring::gen_trapdoor_ring_lwe,
};

use self::{gadget_classical::gen_trapdoor, gadget_parameters::GadgetParameters};
use qfall_math::{
    integer::{MatZ, PolyOverZ, Z},
    integer_mod_q::{MatPolynomialRingZq, MatZq, Modulus},
    rational::Q,
};

pub mod gadget_classical;
pub mod gadget_parameters;
pub mod gadget_ring;
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

/// Computes a trapdoor with default values.
///
/// - `params` is computed using [`GadgetParametersRing::init_default`].
///
/// Parameters:
/// - `n`: the security parameter
/// - `modulus`: the modulus for the trapdoor
///
/// Returns a matrix `a` and its gadget-trapdoor `(r,e)` as in [\[2\]](<index.html#:~:text=[2]>):
/// Construction 1 for some fixed set of parameters [`GadgetParameters::init_default`].
///
/// # Examples
/// ```
/// use qfall_crypto::sample::g_trapdoor::gen_trapdoor_ring_default;
/// use qfall_math::integer::Z;
/// use qfall_math::integer_mod_q::Modulus;
///
/// let (a,r, e) = gen_trapdoor_ring_default(100, &Modulus::try_from(&Z::from(29)).unwrap(), 10);;
/// ```
///
/// # Panics...
/// - ... if the security parameter `n` is not in `\[1, i64::MAX\]`.
pub fn gen_trapdoor_ring_default(
    n: impl Into<Z>,
    modulus: &Modulus,
    s: impl Into<Q>,
) -> (
    MatPolynomialRingZq,
    MatPolynomialRingZq,
    MatPolynomialRingZq,
) {
    // panic if n < 1 (security parameter must be positive)
    let n = n.into();
    assert!(n >= Z::ONE);
    let s = s.into();

    let params = GadgetParametersRing::init_default(n, modulus);

    // a_bar <-$ Zq[X]^n
    let a_bar = PolyOverZ::sample_uniform(&params.n, &0, &params.q).unwrap();

    // we can unwrap, as we compute the parameters on our own and
    // they should always work
    gen_trapdoor_ring_lwe(&params, &a_bar, &s).unwrap()
}

#[cfg(test)]
mod test_gen_trapdoor_default {
    use super::gen_trapdoor_default;
    use qfall_math::integer_mod_q::Modulus;

    /// Ensure that computing a a trapdoor is working when called with the default
    /// parameters
    #[test]
    fn working() {
        let (_, _) = gen_trapdoor_default(100, &Modulus::from(32));
    }
}

#[cfg(test)]
mod test_gen_trapdoor_ring_default {
    use super::gen_trapdoor_ring_default;
    use qfall_math::{integer::Z, integer_mod_q::Modulus};

    /// Ensure that computing a a trapdoor is working when called with the default
    /// parameters
    #[test]
    fn working() {
        let (_, _, _) =
            gen_trapdoor_ring_default(100, &Modulus::try_from(&Z::from(29)).unwrap(), 10);
    }
}
