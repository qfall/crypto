// Copyright © 2023 Niklas Siemer
//
// This file is part of qFALL-crypto.
//
// qFALL-crypto is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

//! This module contains a general implementation of an IND-CCA secure
//! public key encryption scheme constructed
//! via an [`IdentityBasedEncryption`] scheme and a [`SignatureScheme`].

use super::PKEncryptionMut;
use crate::construction::{
    identity_based_encryption::IdentityBasedEncryption, signature::SignatureScheme,
};
use qfall_math::integer::Z;
use serde::{Deserialize, Serialize};

pub mod dual_regev_ibe_pfdh;

/// This struct manages and stores the public parameters of an [`CCSfromIBE`]
/// public key encryption construction based on [\[5\]](<index.html#:~:text=[5]>).
///
/// Attributes:
/// - `ibe`: specifies the IBE scheme used in this construction
/// - `signature`: specifies the signature scheme used in this construction
///
/// # Examples
/// ```
/// todo!();
/// ```
#[derive(Serialize, Deserialize, Clone)]
pub struct CCSfromIBE<IBE: IdentityBasedEncryption, Signature: SignatureScheme>
where
    IBE::Cipher: ToString,
{
    pub ibe: IBE,
    pub signature: Signature,
}

impl<IBE, Signature> PKEncryptionMut for CCSfromIBE<IBE, Signature>
where
    IBE: IdentityBasedEncryption,
    Signature: SignatureScheme,
    IBE::Cipher: ToString,
    IBE::MasterPublicKey: Clone,
    Signature::PublicKey: Into<IBE::Identity> + Clone,
{
    type Cipher = (Signature::PublicKey, IBE::Cipher, Signature::Signature);
    type PublicKey = IBE::MasterPublicKey;
    type SecretKey = (IBE::MasterPublicKey, IBE::MasterSecretKey);

    /// Generates a (pk, sk) pair for the CCS construction
    /// by following these steps:
    /// - (mpk, msk) = ibe.setup()
    ///
    /// Then, `pk = mpk` and `sk = (mpk, msk)` are returned.
    ///
    /// # Examples
    /// ```
    /// todo!();
    /// ```
    fn gen(&mut self) -> (Self::PublicKey, Self::SecretKey) {
        let (pk, sk) = self.ibe.setup();
        (pk.clone(), (pk, sk))
    }

    /// Generates an encryption of `message mod 2` for the provided public key by following these steps:
    /// - (vrfy_key, sign_key) = signature.gen()
    /// - c = ibe.enc(mpk, vrfy_key, message), i.e. encrypt `message mod 2` with respect to identity `vrfy_key`
    /// - sigma = signature.sign(c, sign_key, vrfy_key), i.e. sign message `c`
    ///
    /// Then, the ciphertext `(vrfy_key, c, sigma)` is returned.
    ///
    /// Parameters:
    /// - `pk`: specifies the public key `pk = A`
    /// - `message`: specifies the message that should be encrypted
    ///
    /// Returns a cipher consisting of a tuple `cipher = (vrfy_key, c, sigma)`.
    ///
    /// # Examples
    /// ```
    /// todo!();
    /// ```
    fn enc(&mut self, pk: &Self::PublicKey, message: impl Into<Z>) -> Self::Cipher {
        let message = message.into().modulo(2);

        let (vrfy_key, sign_key) = self.signature.gen();

        let c = self.ibe.enc(pk, &vrfy_key.clone().into(), message);
        let sigma = self.signature.sign(c.to_string(), &sign_key, &vrfy_key);
        (vrfy_key, c, sigma)
    }

    /// Decrypts the provided `cipher` using the secret key `sk` by following these steps:
    /// - if signature.vrfy(c, sigma, vrfy_key) is not successful, output -1, otherwise proceed
    /// - secret_key = ibe.extract(mpk, msk, vrfy_key), i.e. extract the secret key for identity `vrfy_key`
    /// - ibe.dec(secret_key, c)
    ///
    /// Then, the resulting decryption is returned.
    ///
    /// Parameters:
    /// - `sk`: specifies the secret key `sk = (mpk, msk)`
    /// - `cipher`: specifies the cipher containing `cipher = (vrfy_key, c, sigma)`
    ///
    /// Returns the decryption of `cipher` as a [`Z`] instance.
    ///
    /// # Examples
    /// ```
    /// todo!();
    /// ```
    fn dec(&mut self, sk: &Self::SecretKey, cipher: &Self::Cipher) -> Z {
        if !self
            .signature
            .vfy(cipher.1.to_string(), &cipher.2, &cipher.0)
        {
            return Z::MINUS_ONE;
        }

        // as this secret key is never computed or seen by anybody else than the owner of the master secret key,
        // it does not hurt the security premises of the scheme to ignore the storage of the IBE's key extraction,
        // which is required in order to not make `self` mutable accross all public key encryption schemes
        let secret = self.ibe.extract(&sk.0, &sk.1, &cipher.0.clone().into());
        self.ibe.dec(&secret, &cipher.1)
    }
}
