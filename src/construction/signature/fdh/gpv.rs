// Copyright © 2023 Marvin Beckmann
//
// This file is part of qFALL-crypto.
//
// qFALL-crypto is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

//! A classical implementation of the [`Fdh`] scheme using the [`PSFGPV`]
//! according to [\[1\]](<../index.html#:~:text=[1]>).

use super::Fdh;
use crate::{
    primitive::hash::HashMatZq,
    sample::{distribution::psf::gpv::PSFGPV, g_trapdoor::gadget_parameters::GadgetParameters},
};
use qfall_math::{
    integer::{MatZ, Z},
    integer_mod_q::{MatZq, Modulus},
    rational::{MatQ, Q},
};
use std::{collections::HashMap, marker::PhantomData};

impl Fdh<MatZq, (MatZ, MatQ), MatZ, MatZq, PSFGPV, HashMatZq> {
    /// Initializes an FDH signature scheme from a [`PSFGPV`].
    ///
    /// This function corresponds to an implementation of an FDH-signature
    /// scheme with the explicit PSF [`PSFGPV`] which is generated using
    /// the default of [`GadgetParameters`].
    ///
    /// Parameters:
    /// - `n`: The security parameter
    /// - `modulus`: The modulus used for the G-Trapdoors
    /// - `s`: The gaussian parameter with which is sampled
    ///
    /// Returns an explicit implementation of a FDH-signature scheme.
    ///
    /// # Example
    /// ```
    /// use qfall_crypto::construction::signature::{fdh::Fdh, SignatureScheme};
    ///
    /// let m = "Hello World!";
    ///
    /// let mut fdh = Fdh::init_gpv(4, 113, 17);
    /// let (pk, sk) = fdh.gen();
    ///
    /// let sigma = fdh.sign(m.to_string(), &sk, &pk);
    ///
    /// assert!(fdh.vfy(m.to_string(), &sigma, &pk))
    /// ```
    ///
    /// # Panics ...
    /// - if `modulus <= 1`.
    pub fn init_gpv(n: impl Into<Z>, modulus: impl Into<Modulus>, s: impl Into<Q>) -> Self {
        let n = n.into();
        let n_i64 = i64::try_from(&n).unwrap();
        let modulus = modulus.into();
        let psf = PSFGPV {
            gp: GadgetParameters::init_default(&n, &modulus),
            s: s.into(),
        };
        Self {
            psf: Box::new(psf),
            storage: HashMap::new(),
            hash: Box::new(HashMatZq {
                modulus,
                rows: n_i64,
                cols: 1,
            }),
            _a_type: PhantomData,
            _trapdoor_type: PhantomData,
            _range_type: PhantomData,
        }
    }
}

#[cfg(test)]
mod text_fdh {
    use super::Fdh;
    use crate::{
        construction::signature::SignatureScheme, primitive::hash::HashMatZq,
        sample::distribution::psf::gpv::PSFGPV,
    };
    use qfall_math::{
        integer::{MatZ, Z},
        integer_mod_q::MatZq,
        rational::{MatQ, Q},
        traits::Pow,
    };

    /// Ensure that the generated signature is valid.
    #[test]
    fn ensure_valid_signature_is_generated() {
        let n = Z::from(4);
        let k = Z::from(6);
        // `s >= ||\tilde short_base|| * omega(\sqrt{\log m})`,
        // here `\log(2*n*k) = omega(\sqrt{\log m}))` (Theorem 4.1 - GPV08)
        let s: Q = ((&n * &k).sqrt() + 1) * Q::from(2) * (Z::from(2) * &n * &k).log(2).unwrap();
        let modulus = Z::from(2).pow(&k).unwrap();

        let mut fdh = Fdh::init_gpv(n, &modulus, &s);
        let (pk, sk) = fdh.gen();

        for i in 0..10 {
            let m = format!("Hello World! {}", i);

            let sigma = fdh.sign(m.to_owned(), &sk, &pk);

            assert_eq!(&sigma, &fdh.sign(m.to_owned(), &sk, &pk));
            assert!(fdh.vfy(m.to_owned(), &sigma, &pk))
        }
    }

    /// Ensure that an entry is actually added to the local storage.
    #[test]
    fn storage_filled() {
        let mut fdh = Fdh::init_gpv(5, 1024, 10);

        let m = "Hello World!";
        let (pk, sk) = fdh.gen();
        let _ = fdh.sign(m.to_owned(), &sk, &pk);

        assert!(fdh.storage.contains_key(m))
    }

    /// Ensure that after deserialization the HashMap still contains all entries.
    #[test]
    fn reload_hashmap() {
        let mut fdh = Fdh::init_gpv(5, 1024, 10);

        // fill one entry in the HashMap
        let m = "Hello World!";
        let (pk, sk) = fdh.gen();
        let _ = fdh.sign(m.to_owned(), &sk, &pk);

        let fdh_string = serde_json::to_string(&fdh).expect("Unable to create a json object");
        let fdh_2: Fdh<MatZq, (MatZ, MatQ), MatZ, MatZq, PSFGPV, HashMatZq> =
            serde_json::from_str(&fdh_string).unwrap();

        assert_eq!(fdh.storage, fdh_2.storage);
    }
}
