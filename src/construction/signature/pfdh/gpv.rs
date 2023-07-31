// Copyright © 2023 Phil Milewski
//
// This file is part of qFALL-crypto.
//
// qFALL-crypto is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

//! A classical implementation of the [`Pfdh`] scheme using the [`PSFGPV`]
//! according to [\[1\]](<../index.html#:~:text=[1]>).

use super::Pfdh;
use crate::{
    primitive::hash::HashMatZq,
    sample::{distribution::psf::gpv::PSFGPV, g_trapdoor::gadget_parameters::GadgetParameters},
};
use qfall_math::{
    integer::{MatZ, Z},
    integer_mod_q::{MatZq, Modulus},
    rational::Q,
};
use std::marker::PhantomData;

impl Pfdh<MatZq, MatZ, MatZ, MatZq, PSFGPV, HashMatZq> {
    /// Initializes an PFDH signature scheme from a [`PSFGPV`].
    ///
    /// This function corresponds to an implementation of an PFDH-signature
    /// scheme with the explicit PSF [`PSFGPV`] which is generated using
    /// the default of [`GadgetParameters`].
    ///
    /// Parameters:
    /// - `n`: The security parameter
    /// - `modulus`: The modulus used for the G-Trapdoors
    /// - `s`: The gaussian parameter with which is sampled
    /// - `randomness_length`: the number of bits used for the randomness
    ///
    /// Returns an explicit implementation of a PFDH-signature scheme.
    ///
    /// # Example
    /// ```
    /// use qfall_crypto::construction::signature::pfdh::Pfdh;
    /// use qfall_math::integer::Z;
    /// use qfall_math::integer_mod_q::Modulus;
    /// use qfall_math::rational::Q;
    /// use crate::qfall_crypto::construction::signature::SignatureScheme;
    ///
    /// let s = Q::from(17);
    /// let n = Z::from(4);
    /// let modulus = Modulus::try_from(&Z::from(113)).unwrap();
    ///
    /// let mut pfdh = Pfdh::init_gpv(n, &modulus, &s, 128);
    ///
    /// let m = "Hello World!";
    ///
    /// let (pk, sk) = pfdh.gen();
    /// let sigma = pfdh.sign(m.to_owned(), &sk, &pk);
    ///
    /// assert!(pfdh.vfy(m.to_owned(), &sigma, &pk))
    /// ```
    pub fn init_gpv<Integer: Into<Z>>(
        n: impl Into<Z>,
        modulus: &Modulus,
        s: &Q,
        randomness_length: Integer,
    ) -> Self {
        let n = n.into();
        let r: Z = randomness_length.into();
        let psf = PSFGPV {
            gp: GadgetParameters::init_default(&n, modulus),
            s: s.clone(),
        };
        let n = i64::try_from(&n).unwrap();
        Self {
            psf: Box::new(psf),
            hash: Box::new(HashMatZq {
                modulus: modulus.clone(),
                rows: n,
                cols: 1,
            }),
            randomness_length: r,
            _a_type: PhantomData,
            _trapdoor_type: PhantomData,
            _domain_type: PhantomData,
            _range_type: PhantomData,
        }
    }
}

#[cfg(test)]
mod text_fdh {
    use super::Pfdh;
    use crate::construction::signature::SignatureScheme;
    use qfall_math::{integer::Z, integer_mod_q::Modulus, rational::Q, traits::Pow};

    /// Ensure that the generated signature is valid
    #[test]
    fn ensure_valid_signature_is_generated() {
        let n = Z::from(4);
        let k = Z::from(6);
        // `s >= ||\tilde short_base|| * omega(\sqrt{\log m})`,
        // here `\log(2*n*k) = omega(\sqrt{\log m}))` (Theorem 4.1 - GPV08)
        let s: Q = ((&n * &k).sqrt() + 1) * Q::from(2) * (Z::from(2) * &n * &k).log(2).unwrap();
        let modulus = Modulus::try_from(&Z::from(2).pow(&k).unwrap()).unwrap();

        let mut pfdh = Pfdh::init_gpv(n, &modulus, &s, 128);
        let (pk, sk) = pfdh.gen();

        for i in 0..10 {
            let m = format!("Hello World! {}", i);

            let sigma = pfdh.sign(m.to_owned(), &sk, &pk);

            assert!(pfdh.vfy(m.to_owned(), &sigma, &pk))
        }
    }
}
