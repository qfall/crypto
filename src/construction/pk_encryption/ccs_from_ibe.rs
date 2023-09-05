// Copyright © 2023 Niklas Siemer
//
// This file is part of qFALL-crypto.
//
// qFALL-crypto is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

//! This module contains an implementation of the IND-CCA secure
//! public key encryption scheme constructed via an IBE scheme and a signature.

use super::PKEncryption;
use crate::construction::{
    identity_based_encryption::IdentityBasedEncryption, signature::SignatureScheme,
};
use qfall_math::integer::Z;
use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize)]
pub struct CCSfromIBE<IBE: IdentityBasedEncryption, Signature: SignatureScheme> {
    ibe: IBE,
    signature: Signature,
}

impl<IBE, Signature> PKEncryption for CCSfromIBE<IBE, Signature>
where
    IBE: IdentityBasedEncryption,
    Signature: SignatureScheme,
{
    type Cipher = (
        SignatureScheme::PublicKey,
        IdentityBasedEncryption::Cipher,
        SignatureScheme::Signature,
    );
    type PublicKey = IdentityBasedEncryption::MasterPublicKey;
    type SecretKey = (
        IdentityBasedEncryption::MasterPublicKey,
        IdentityBasedEncryption::MasterSecretKey,
    );

    fn gen(&self) -> (Self::PublicKey, Self::SecretKey) {
        let (pk, sk) = self.ibe.setup();
        (pk, (pk, sk))
    }

    fn enc(&self, pk: &Self::PublicKey, message: impl Into<Z>) -> Self::Cipher {
        let message = message.into().modulo(2).unwrap();

        let (sign_key, vrfy_key) = self.signature.gen();

        let c = self.ibe.enc(pk, vrfy_key, message);
        let sigma = self.signature.sign(c.to_string(), sk, pk);
        (vrfy_key, c, sigma)
    }

    fn dec(&self, sk: &Self::SecretKey, cipher: &Self::Cipher) -> Z {
        if !self.signature.vfy(cipher.1, cipher.2, cipher.0) {
            return Z::MINUS_ONE;
        }

        let secret = self.ibe.extract(sk.0, sk.1, cipher.0);
        self.ibe.dec(secret, cipher.1)
    }
}
