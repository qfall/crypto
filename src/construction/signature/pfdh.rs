// Copyright © 2023 Phil Milewski
//
// This file is part of qFALL-crypto.
//
// qFALL-crypto is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

//! This Module contains a general implementation of the [`Pfdh`] scheme.
//!
//! Implementation of a [`Pfdh`]-signature scheme are thereby fairly easy,
//! see [`Pfdh::init_gpv`](crate::construction::signature::pfdh::gpv) that
//! works with every PSF and a corresponding hash function.

use super::SignatureScheme;
use crate::{primitive::hash::HashInto, sample::distribution::psf::PSF};
use qfall_math::{integer::Z, traits::Pow};
use serde::{Deserialize, Serialize};
use std::marker::PhantomData;

pub mod gpv;
pub mod serialize;

#[derive(Serialize)]
pub struct Pfdh<
    A,
    Trapdoor,
    Domain: Serialize + for<'a> Deserialize<'a>,
    Range,
    T: PSF<A, Trapdoor, Domain, Range> + Serialize + for<'a> Deserialize<'a>,
    Hash: HashInto<Range> + Serialize + for<'a> Deserialize<'a>,
    Randomness: Into<Z>,
> {
    pub psf: Box<T>,
    pub hash: Box<Hash>,
    pub randomness_length: Z,

    // The parameters below can be ignored, they are just there for generic usage
    #[serde(skip_serializing)]
    pub _a_type: PhantomData<A>,
    #[serde(skip_serializing)]
    pub _trapdoor_type: PhantomData<Trapdoor>,
    #[serde(skip_serializing)]
    pub _domain_type: PhantomData<Domain>,
    #[serde(skip_serializing)]
    pub _range_type: PhantomData<Range>,
    #[serde(skip_serializing)]
    pub _randomness_length: PhantomData<Randomness>,
}

impl<A, Trapdoor, Domain, Range, T, Hash, Randomness> SignatureScheme
    for Pfdh<A, Trapdoor, Domain, Range, T, Hash, Randomness>
where
    Domain: Clone + Serialize + for<'a> Deserialize<'a>,
    Range: PartialEq<Range>,
    T: PSF<A, Trapdoor, Domain, Range> + Serialize + for<'a> Deserialize<'a>,
    Hash: HashInto<Range> + Serialize + for<'a> Deserialize<'a>,
    Randomness: Into<Z> + Clone,
{
    type SecretKey = Trapdoor;
    type PublicKey = A;
    type Signature = (Domain, Z);

    /// Generates a trapdoor by calling the `trap_gen` of the psf
    fn gen(&mut self) -> (Self::PublicKey, Self::SecretKey) {
        self.psf.trap_gen()
    }

    /// Firstly generate randomness
    /// It hashes the message and randomness into the domain and then computes a signature using
    /// `samp_p` from the psf with the trapdoor.
    fn sign(&mut self, m: String, sk: &Self::SecretKey, pk: &Self::PublicKey) -> Self::Signature {
        let randomness =
            Z::sample_uniform(0, Z::from(2).pow(&self.randomness_length).unwrap()).unwrap();
        let u = (self.hash).hash(&format!("{m} {randomness} {}", &self.randomness_length));
        let signature_part1 = self.psf.samp_p(pk, sk, &u);

        (signature_part1, randomness)
    }

    /// Checks if a signature is firstly within D_n, and then checks if
    /// the signature is actually a valid preimage under `fa` of `hash(m||r)`.
    fn vfy(&self, m: String, sigma: &Self::Signature, pk: &Self::PublicKey) -> bool {
        if !self.psf.check_domain(&sigma.0) {
            return false;
        }

        let u = (self.hash).hash(&format!("{m} {} {}", sigma.1, &self.randomness_length));

        self.psf.f_a(pk, &sigma.0) == u
    }
}
