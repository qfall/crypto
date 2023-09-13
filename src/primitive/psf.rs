// Copyright © 2023 Marvin Beckmann
//
// This file is part of qFALL-crypto.
//
// qFALL-crypto is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

//! Contains the [`PSF`] trait, which is then subsequently implemented by
//! explicit constructions such as [`PSFGPV`].
//!
//! The main references are listed in the following
//! and will be further referenced in submodules by these numbers:
//! - \[1\] Micciancio, D., Peikert, C. (2012).
//! Trapdoors for Lattices: Simpler, Tighter, Faster, Smaller.
//! In: Pointcheval, D., Johansson, T. (eds) Advances in Cryptology – EUROCRYPT 2012.
//! EUROCRYPT 2012. Lecture Notes in Computer Science, vol 7237.
//! Springer, Berlin, Heidelberg. <https://doi.org/10.1007/978-3-642-29011-4_41>
//! - \[2\] Gür, K.D., Polyakov, Y., Rohloff, K., Ryan, G.W. and Savas, E., 2018,
//! January. Implementation and evaluation of improved Gaussian sampling for lattice
//!  trapdoors. In Proceedings of the 6th Workshop on Encrypted Computing & Applied
//! Homomorphic Cryptography (pp. 61-71). <https://dl.acm.org/doi/pdf/10.1145/3267973.3267975>

mod gpv;
mod gpv_ring;

pub use gpv::PSFGPV;
pub use gpv_ring::PSFGPVRing;

/// This trait should be implemented by all constructions that are
/// actual implementations of a preimage sampleable function.
/// A formal definition for these PSFs can be found in
/// [\[1\]](<index.html#:~:text=[1]>)
pub trait PSF<A, Trapdoor, Domain, Range> {
    /// Samples a parity-check matrix and a trapdoor for that matrix
    fn trap_gen(&self) -> (A, Trapdoor);
    /// Samples an element in the domain according to a specified distribution
    fn samp_d(&self) -> Domain;
    /// Samples an element `e` in the domain according to a specified distribution
    /// conditioned on `f_a(a, e) = u`
    fn samp_p(&self, a: &A, r: &Trapdoor, u: &Range) -> Domain;
    /// Implements the efficiently computable function `fa`
    /// which is uniquely classified by `a`
    fn f_a(&self, a: &A, sigma: &Domain) -> Range;
    /// Checks whether an element is in the correct domain (and not just the correct type)
    fn check_domain(&self, sigma: &Domain) -> bool;
}
