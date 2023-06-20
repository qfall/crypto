// Copyright © 2023 Marvin Beckmann
//
// This file is part of qFALL-crypto.
//
// qFALL-crypto is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

//! This Module contains a general implementation of the [`Fdh`] scheme,
//! which only has to be instantiated with a corresponding PSF, a storage and
//! a corresponding hash function.
//!
//! Explicit implementation of a [`Fdh`]-signature scheme are thereby fairly easy,
//! see [`Fdh::init_gpv`].

use super::SignatureScheme;
use crate::{primitive::hash::HashInto, sample::distribution::psf::PSF};
use serde::{Deserialize, Serialize};
use std::{collections::HashMap, marker::PhantomData};

pub mod gpv;
pub mod serialize;

#[derive(Serialize)]
pub struct Fdh<
    A,
    Trapdoor,
    Domain: Serialize + for<'a> Deserialize<'a>,
    Range,
    T: PSF<A, Trapdoor, Domain, Range> + Serialize + for<'a> Deserialize<'a>,
    Hash: HashInto<Range> + Serialize + for<'a> Deserialize<'a>,
> {
    pub psf: Box<T>,
    pub storage: HashMap<String, Domain>,
    pub hash: Box<Hash>,

    // The parameters below can be ignored, they are just there for generic usage
    #[serde(skip_serializing)]
    pub _a_type: PhantomData<A>,
    #[serde(skip_serializing)]
    pub _trapdoor_type: PhantomData<Trapdoor>,
    #[serde(skip_serializing)]
    pub _range_type: PhantomData<Range>,
}

impl<A, Trapdoor, Domain, Range, T, Hash> SignatureScheme
    for Fdh<A, Trapdoor, Domain, Range, T, Hash>
where
    Domain: Clone + Serialize + for<'a> Deserialize<'a>,
    Range: PartialEq<Range>,
    T: PSF<A, Trapdoor, Domain, Range> + Serialize + for<'a> Deserialize<'a>,
    Hash: HashInto<Range> + Serialize + for<'a> Deserialize<'a>,
{
    type SecretKey = Trapdoor;
    type PublicKey = A;
    type Signature = Domain;

    fn gen(&mut self) -> (Self::PublicKey, Self::SecretKey) {
        self.psf.trap_gen()
    }

    fn sign(&mut self, m: String, sk: &Self::SecretKey, pk: &Self::PublicKey) -> Self::Signature {
        // check if it is in the HashMap
        if let Some(sigma) = self.storage.get(&m) {
            return sigma.clone();
        }

        let u = (self.hash).hash(&m);
        let signature = self.psf.samp_p(pk, sk, &u);

        // insert signature in HashMap
        self.storage.insert(m, signature.clone());
        signature
    }

    fn vfy(&self, m: String, sigma: &Self::Signature, pk: &Self::PublicKey) -> bool {
        if !self.psf.check_dn(sigma) {
            return false;
        }

        let u = (self.hash).hash(&m);

        self.psf.fa(pk, sigma) == u
    }
}
