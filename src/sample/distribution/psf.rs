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

pub trait PSF<A, Trapdoor, Domain, Range> {
    fn trap_gen(&self) -> (A, Trapdoor);
    fn samp_d(&self) -> Domain;
    fn samp_p(&self, a: &A, r: &Trapdoor, u: &Range) -> Domain;
    fn fa(&self, a: &A, sigma: &Domain) -> Range;
    fn check_dn(&self, sigma: &Domain) -> bool;
}
