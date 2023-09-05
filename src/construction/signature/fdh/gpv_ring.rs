// Copyright © 2023 Marvin Beckmann
//
// This file is part of qFALL-crypto.
//
// qFALL-crypto is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

//! A ring implementation of the [`Fdh`] scheme using the [`PSFGPVRing`]
//! according to [\[1\]](<../index.html#:~:text=[1]>).

use super::Fdh;
use crate::{
    construction::hash::sha256::HashMatPolynomialRingZq, primitive::psf::gpv_ring::PSFGPVRing,
    sample::g_trapdoor::gadget_parameters::GadgetParametersRing,
};
use qfall_math::{
    integer::{MatPolyOverZ, Z},
    integer_mod_q::{MatPolynomialRingZq, Modulus},
    rational::Q,
};
use std::{collections::HashMap, marker::PhantomData};

impl
    Fdh<
        MatPolynomialRingZq,
        (MatPolyOverZ, MatPolyOverZ),
        MatPolyOverZ,
        MatPolynomialRingZq,
        PSFGPVRing,
        HashMatPolynomialRingZq,
    >
{
    /// Initializes an FDH signature scheme from a [`PSFGPVRing`].
    /// The trapdoor is sampled with a gaussian parameter of 1.005
    /// as done in [\[3\]](<index.html#:~:text=[3]>) who derived it from
    /// [\[5\]](<index.html#:~:text=[5]>).
    ///
    /// This function corresponds to an implementation of an FDH-signature
    /// scheme with the explicit PSF [`PSFGPVRing`] which is generated using
    /// the default of [`GadgetParametersRing`].
    ///
    /// Parameters:
    /// - `n`: The security parameter
    /// - `modulus`: The modulus used for the G-Trapdoors
    /// - `s`: The gaussian parameter with which is sampled
    ///
    /// Returns an explicit implementation of an FDH-signature scheme.
    ///
    /// # Example
    /// ```
    /// use qfall_crypto::construction::signature::{fdh::Fdh, SignatureScheme};
    ///
    /// let mut fdh = Fdh::init_gpv_ring(8, 512, 100);
    /// let (pk, sk) = fdh.gen();
    ///
    /// let m = &format!("Hello World!");
    ///
    /// let sigma = fdh.sign(m.to_owned(), &sk, &pk);
    /// assert!(fdh.vfy(m.to_owned(), &sigma, &pk));
    /// ```
    ///
    /// # Panics ...
    /// - if `modulus <= 1`.
    pub fn init_gpv_ring(n: impl Into<Z>, modulus: impl Into<Modulus>, s: impl Into<Q>) -> Self {
        let n = n.into();
        let modulus = modulus.into();
        let s = s.into();
        let psf = PSFGPVRing {
            gp: GadgetParametersRing::init_default(&n, &modulus),
            s,
            s_td: Q::from(1.005_f64),
        };
        let modulus = psf.gp.modulus.clone();
        Self {
            psf: Box::new(psf),
            storage: HashMap::new(),
            hash: Box::new(HashMatPolynomialRingZq {
                modulus,
                rows: 1,
                cols: 1,
            }),
            _a_type: PhantomData,
            _trapdoor_type: PhantomData,
            _range_type: PhantomData,
        }
    }
}

#[cfg(test)]
mod test_fdh {
    use super::{Fdh, PSFGPVRing};
    use crate::{
        construction::hash::sha256::HashMatPolynomialRingZq,
        construction::signature::SignatureScheme,
    };
    use qfall_math::{integer::MatPolyOverZ, integer_mod_q::MatPolynomialRingZq, rational::Q};

    const MODULUS: i64 = 512;
    const N: i64 = 8;
    fn compute_s() -> Q {
        ((2 * 2 * Q::from(1.005_f64) * Q::from(N).sqrt() + 1) * 2) * 4
    }

    /// Ensure that the generated signature is valid.
    #[test]
    fn ensure_valid_signature_is_generated() {
        let mut fdh = Fdh::init_gpv_ring(N, MODULUS, &compute_s());
        let (pk, sk) = fdh.gen();

        for i in 0..10 {
            let m = &format!("Hello World! {i}");

            let sigma = fdh.sign(m.to_owned(), &sk, &pk);

            assert!(
                fdh.vfy(m.to_owned(), &sigma, &pk),
                "This is a probabilistic test and may fail with negligible probability. \
                As n is rather small here, try to rerun the test and check whether the \
                test fails again."
            )
        }
    }

    /// Ensure that an entry is actually added to the local storage.
    #[test]
    fn storage_filled() {
        let mut fdh = Fdh::init_gpv_ring(N, MODULUS, &compute_s());

        let m = "Hello World!";
        let (pk, sk) = fdh.gen();
        let sign_1 = fdh.sign(m.to_owned(), &sk, &pk);
        let sign_2 = fdh.sign(m.to_owned(), &sk, &pk);

        assert!(fdh.storage.contains_key(m));
        assert_eq!(sign_1, sign_2);
    }

    /// Ensure that after deserialization the HashMap still contains all entries.
    #[test]
    fn reload_hashmap() {
        let mut fdh = Fdh::init_gpv_ring(N, MODULUS, &compute_s());

        // fill one entry in the HashMap
        let m = "Hello World!";
        let (pk, sk) = fdh.gen();
        let _ = fdh.sign(m.to_owned(), &sk, &pk);

        let fdh_string = serde_json::to_string(&fdh).expect("Unable to create a json object");

        let fdh_2: Fdh<
            MatPolynomialRingZq,
            (MatPolyOverZ, MatPolyOverZ),
            MatPolyOverZ,
            MatPolynomialRingZq,
            PSFGPVRing,
            HashMatPolynomialRingZq,
        > = serde_json::from_str(&fdh_string).unwrap();

        assert_eq!(fdh.storage, fdh_2.storage);
    }
}
