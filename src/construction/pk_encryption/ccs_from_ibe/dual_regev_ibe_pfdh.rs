// Copyright © 2023 Niklas Siemer
//
// This file is part of qFALL-crypto.
//
// qFALL-crypto is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

//! A classical implementation of the [`CCSfromIBE`] scheme using
//! the [`DualRegevIBE`] and [`Pfdh`].

use super::CCSfromIBE;
use crate::{
    construction::{
        hash::sha256::HashMatZq, identity_based_encryption::DualRegevIBE, signature::pfdh::Pfdh,
    },
    primitive::psf::gpv::PSFGPV,
};
use qfall_math::{
    integer::{MatZ, Z},
    integer_mod_q::{MatZq, Modulus},
    rational::{MatQ, Q},
};

impl CCSfromIBE<DualRegevIBE, Pfdh<MatZq, (MatZ, MatQ), MatZ, MatZq, PSFGPV, HashMatZq>> {
    /// Initializes a [`CCSfromIBE`] PK encryption scheme from a [`DualRegevIBE`] and a [`Pfdh`] signature.
    ///
    /// Parameters:
    /// - `n`: specifies the security parameter
    /// - `modulus`: specifies the modulus
    /// - `r`: specifies the Gaussian parameter used for the [`PSFGPV`] for the [`Pfdh`]
    /// - `randomness_length`: specifies the number of bits added to the message before signing
    /// - `alpha`: specifies the Gaussian parameter used for encryption in
    /// [`DualRegev`](crate::construction::pk_encryption::DualRegev) in the [`DualRegevIBE`]
    ///
    /// Returns an explicit implementation of an IND-CCA-secure public key encryption scheme.
    ///
    /// # Example
    /// ```
    /// use qfall_crypto::construction::pk_encryption::CCSfromIBE;
    ///
    /// let mut scheme = CCSfromIBE::init_dr_pfdh(4, 13933, 4, 10.77, 0.0021);
    /// ```
    ///
    /// # Panics ...
    /// - if `modulus <= 1`.
    /// - if `n < 1`.
    pub fn init_dr_pfdh(
        n: impl Into<Z>, // security parameter
        modulus: impl Into<Modulus>,
        randomness_length: impl Into<Z>, // added for to the message before signing
        r: impl Into<Q>,                 // gaussian parameter for PSF
        alpha: impl Into<Q>,             // gaussian parameter for Dual Regev Encryption
    ) -> Self {
        let n = n.into();
        let modulus = modulus.into();
        let r = r.into();

        let dr_ibe = DualRegevIBE::new(&n, &modulus, &r, alpha);
        let pfdh = Pfdh::init_gpv(n, modulus, r, randomness_length);

        Self {
            ibe: dr_ibe,
            signature: pfdh,
        }
    }

    /// Initializes a [`CCSfromIBE`] PK encryption scheme from a [`DualRegevIBE`] and a [`Pfdh`] signature
    /// from a given `n > 0`.
    ///
    /// Parameters:
    /// - `n`: specifies the security parameter
    ///
    /// Returns an explicit implementation of an IND-CCA-secure public key encryption scheme
    /// chosen with appropriate parameters for given `n`..
    ///
    /// # Example
    /// ```
    /// use qfall_crypto::construction::pk_encryption::CCSfromIBE;
    ///
    /// let mut scheme = CCSfromIBE::init_dr_pfdh_from_n(4);
    /// ```
    ///
    /// # Panics ...
    /// - if `n < 4`.
    pub fn init_dr_pfdh_from_n(n: impl Into<Z>) -> Self {
        let n = n.into();
        assert!(
            n > Z::from(3),
            "n needs to be chosen larger than 3 for this function to work properly."
        );

        let ibe = DualRegevIBE::new_from_n(&n);
        let pfdh = Pfdh::init_gpv(&n, &ibe.dual_regev.q, &ibe.psf.s, &n);

        Self {
            ibe,
            signature: pfdh,
        }
    }
}

#[cfg(test)]
mod test_ccs_from_ibe {
    use super::CCSfromIBE;
    use crate::construction::pk_encryption::PKEncryptionMut;
    use qfall_math::integer::Z;

    /// Checks whether the full-cycle of gen, enc, dec works properly
    /// for message 0 and small n.
    #[test]
    fn cycle_zero() {
        let msg = Z::ZERO;
        let mut scheme = CCSfromIBE::init_dr_pfdh_from_n(4);

        let (pk, sk) = scheme.gen();
        let cipher = scheme.enc(&pk, &msg);
        let m = scheme.dec(&sk, &cipher);
        assert_eq!(msg, m);
    }

    /// Checks whether the full-cycle of gen, enc, dec works properly
    /// for message 1 and small n.
    #[test]
    fn cycle_one() {
        let msg = Z::ONE;
        let mut scheme = CCSfromIBE::init_dr_pfdh_from_n(4);

        let (pk, sk) = scheme.gen();
        let cipher = scheme.enc(&pk, &msg);
        let m = scheme.dec(&sk, &cipher);
        assert_eq!(msg, m);
    }
}
