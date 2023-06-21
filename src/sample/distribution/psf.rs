// Copyright © 2023 Marvin Beckmann
//
// This file is part of qFALL-crypto.
//
// qFALL-crypto is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

//! Contains the [`PSF`] trait, which is then subsequently implemented by
//! explicit constructions such as [`PSFGPV`](gpv::PSFGPV).

pub mod gpv;

/// This trait should be implemented by all constructions that are
/// actual implementations of a [`PSF`].
/// A formal definition for these PSFs can be found in
/// [\[1\]](<../index.html#:~:text=[1]>)
pub trait PSF<A, Trapdoor, Domain, Range> {
    /// Samples a parity-check matrix and a trapdoor for that matrix
    fn trap_gen(&self) -> (A, Trapdoor);
    /// Samples an element in the domain according to a specified distribution
    fn samp_d(&self) -> Domain;
    /// Samples an element `e` in the domain according to a specified distribution
    /// conditioned on `fa(a, e) = u`
    fn samp_p(&self, a: &A, r: &Trapdoor, u: &Range) -> Domain;
    /// Implements the efficiently computable function `fa`
    /// which is uniquely classified by `a`
    fn fa(&self, a: &A, sigma: &Domain) -> Range;
    /// Checks whether an element is in the correct domain (and not just the correct type)
    fn check_dn(&self, sigma: &Domain) -> bool;
}
