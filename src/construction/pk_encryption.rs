// Copyright © 2023 Niklas Siemer
//
// This file is part of qFALL-crypto.
//
// qFALL-crypto is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

//! This module provides the trait a struct should implement if it is an
//! instance of a public key encryption scheme. Furthermore, it contains
//! cryptographic schemes implementing the `PKEncryption` trait.
//!
//! The main references are listed in the following:
//! - \[1\] Peikert, Chris (2016).
//! A decade of lattice cryptography.
//! In: Theoretical Computer Science 10.4.
//! <https://web.eecs.umich.edu/~cpeikert/pubs/lattice-survey.pdf>
//! - \[2\] Gentry, Craig and Peikert, Chris and Vaikuntanathan, Vinod (2008).
//! Trapdoors for hard lattices and new cryptographic constructions.
//! In: Proceedings of the fortieth annual ACM symposium on Theory of computing.
//! <https://dl.acm.org/doi/pdf/10.1145/1374376.1374407>
//! - \[3\] Regev, Oded (2009).
//! On lattices, learning with errors, random linear codes, and cryptography.
//! In: Journal of the ACM 6.
//! <https://dl.acm.org/doi/pdf/10.1145/1568318.1568324>
//! \[4\] Lindner, R., and C. Peikert (2011).
//! Better key sizes (and attacks) for LWE-based encryption.
//! In: Topics in Cryptology -  RSA Conference 2011, Springer.
//! <https://eprint.iacr.org/2010/613.pdf>

mod dual_regev;
mod dual_regev_discrete_gauss;
mod lpr;
mod regev;
mod regev_discrete_gauss;
mod ring_lpr;
pub use dual_regev::DualRegev;
pub use dual_regev_discrete_gauss::DualRegevWithDiscreteGaussianRegularity;
pub use lpr::LPR;
pub use regev::Regev;
pub use regev_discrete_gauss::RegevWithDiscreteGaussianRegularity;
pub use ring_lpr::RingLPR;

use qfall_math::integer::Z;

/// This trait should be implemented by every public key encryption scheme.
/// It offers a simple interface to use and implement PKEs.
pub trait PKEncryption {
    type PublicKey;
    type SecretKey;
    type Cipher;

    /// Generates a public key pair `(pk, sk)` suitable for the specific scheme.
    ///
    /// Returns a tuple `(pk, sk)` consisting of [`Self::PublicKey`] and [`Self::SecretKey`].
    fn gen(&self) -> (Self::PublicKey, Self::SecretKey);

    /// Encrypts the provided `message` using the public key `pk`.
    ///
    /// Parameters:
    /// - `pk`: specifies the public key used for encryption
    /// - `message`: specifies the message to be encrypted
    ///
    /// Returns the encryption of `message` as a [`Self::Cipher`] instance.
    fn enc(&self, pk: &Self::PublicKey, message: impl Into<Z>) -> Self::Cipher;

    /// Decrypts the provided `cipher` using the secret key `sk`.
    ///
    /// Parameters:
    /// - `sk`: specifies the secret key used for decryption
    /// - `cipher`: specifies the ciphertext to be decrypted
    ///
    /// Returns the decryption of `cipher` as a [`Z`] instance.
    fn dec(&self, sk: &Self::SecretKey, cipher: &Self::Cipher) -> Z;
}

/// This trait generically implements multi-bit encryption
/// for any scheme implementing the [`PKEncryption`] trait.
///
/// It splits the given ciphertext up into its bits and
/// stores the individual encrypted bits as a vector of ciphertexts.
pub trait GenericMultiBitEncryption: PKEncryption {
    /// Encrypts multiple bits by appending several encryptions of single bits.
    /// The order of single ciphers is `[c0, c1, ..., cn]`, where `c0` is the least significant bit.
    /// Negative values are not allowed. Hence, the absolute value is being encrypted.
    ///
    /// Parameters:
    /// - `pk`: specifies the public key
    /// - `message`: specifies the message that should be encryted
    ///
    /// Returns a cipher of type [`Vec`] containing [`PKEncryption::Cipher`].
    fn enc_multiple_bits(&self, pk: &Self::PublicKey, message: impl Into<Z>) -> Vec<Self::Cipher> {
        let message: Z = message.into().abs();

        let bits = message.to_bits();
        let mut out = vec![];
        for bit in bits {
            if bit {
                out.push(self.enc(pk, Z::ONE));
            } else {
                out.push(self.enc(pk, Z::ZERO));
            }
        }

        out
    }

    /// Decrypts a multiple bit ciphertext.
    ///
    /// Parameters:
    /// - `sk`: specifies the secret key used for decryption
    /// - `cipher`: specifies a slice of ciphers containing several [`PKEncryption::Cipher`] instances
    /// to be decrypted
    ///
    /// Returns the decryption of `cipher` as a [`Z`] instance.
    fn dec_multiple_bits(&self, sk: &Self::SecretKey, cipher: &[Self::Cipher]) -> Z {
        let mut bits = vec![];

        for item in cipher {
            if self.dec(sk, item) == Z::ZERO {
                bits.push(false);
            } else {
                bits.push(true);
            }
        }

        Z::from_bits(&bits)
    }
}
