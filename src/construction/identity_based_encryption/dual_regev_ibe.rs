// Copyright © 2023 Phil Milewski
//
// This file is part of qFALL-crypto.
//
// qFALL-crypto is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

//! This module contains an implementation of the IND-CPA secure
//! identity based public key encryption scheme. The encryption scheme is based
//! on [`DualRegev`].
//!
//! The main references are listed in the following:
//! - \[1\] Gentry, Craig and Peikert, Chris and Vaikuntanathan, Vinod (2008).
//! Trapdoors for hard lattices and new cryptographic constructions.
//! In: Proceedings of the fortieth annual ACM symposium on Theory of computing.
//! <https://dl.acm.org/doi/pdf/10.1145/1374376.1374407>

use super::IBE;
use crate::{
    construction::pk_encryption::{DualRegev, PKEncryption},
    primitive::hash::hash_to_mat_zq_sha256,
    sample::{
        distribution::psf::{gpv::PSFGPV, PSF},
        g_trapdoor::gadget_parameters::GadgetParameters,
    },
};
use qfall_math::{
    integer::{MatZ, Z},
    integer_mod_q::{MatZq, Modulus, Zq},
    rational::Q,
    traits::{Concatenate, Distance, GetEntry, GetNumColumns, GetNumRows, SetEntry},
};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// This struct manages and stores the public parameters of a [`IBE`]
/// public key encryption instance.
///
/// Attributes:
/// - `n`: specifies the security parameter, which is not equal to the bit-security level
/// - `m`: defines the dimension of the underlying lattice
/// - `q`: specifies the modulus over which the encryption is computed
/// - `r`: specifies the gaussian parameter used for SampleD,
///   i.e. used for encryption
/// - `alpha`: specifies the gaussian parameter used for independent
///   sampling from χ, i.e. for multiple discrete Gaussian samples used
///   for key generation
///todo
/// # Examples
/// ```
/// ```
#[derive(Serialize, Deserialize)]
pub struct DualRegevIBE {
    n: Z,       // security parameter
    m: Z,       // number of rows of matrix A
    q: Modulus, // modulus
    r: Q,       // gaussian parameter for sampleD
    alpha: Q,   // gaussian parameter for sampleZ
    psf: PSFGPV,
    storage: HashMap<String, MatZ>,
}

impl IBE for DualRegevIBE {
    type Cipher = MatZq;
    type MasterPublicKey = MatZq;
    type MasterSecretKey = MatZ;
    type SecretKey = MatZ;
    type Identity = String;

    /// Generates a (pk, sk) pair for the Dual Regev public key encryption scheme
    /// by following these steps:
    /// - s <- Z_q^n
    /// - A <- Z_q^{n x m}
    /// - x <- χ^m
    /// - p = A^t * s + x
    ///
    /// Then, `pk = (A, p)` and `sk = s` is output.
    ///
    /// # Examples
    /// ```
    /// use qfall_crypto::construction::pk_encryption::{PKEncryption, DualRegev};
    /// let dual_regev = DualRegev::default();
    ///
    /// let (pk, sk) = dual_regev.gen();
    /// ```
    fn gen(&self) -> (Self::MasterPublicKey, Self::MasterSecretKey) {
        self.psf.trap_gen()
    }

    /// Given an identity it extracts a corresponding secret key by using samp_p
    /// of the given [`PSF`].
    ///
    /// Parameters:
    /// - `master_pk`: The master public key for the encryption scheme
    /// - `master_sk`: Zhe master secret key of the encryption scheme, namely
    /// the trapdoor for the [`PSF`]
    /// - `identity`: The identity, for which the corresponding secret key
    ///  should be returned
    ///
    /// Returns the corresponding secret key of `identity` under public key
    /// `master_pk`.
    ///
    /// # Examples
    /// ```
    /// use qfall_crypto::construction::identity_based_encryption::{IBE, DualRegevIBE};
    /// let mut ibe = DualRegevIBE::default();
    /// let (master_pk, master_sk) = ibe.gen();
    /// let id = String::from("Hello World!");
    ///
    /// let sk = ibe.extract(&master_pk, &master_sk, &id);
    /// ```
    fn extract(
        &mut self,
        master_pk: &Self::MasterPublicKey,
        master_sk: &Self::MasterSecretKey,
        identity: &Self::Identity,
    ) -> Self::SecretKey {
        println!("Begin Extract");
        // check if it is in the HashMap
        if let Some(value) = self.storage.get(identity) {
            return value.clone();
        }

        let u = hash_to_mat_zq_sha256(&identity, &self.n, 1, &self.q);
        let secret_key = self.psf.samp_p(&master_pk, &master_sk, &u);

        // insert secret key in HashMap
        self.storage.insert(identity.clone(), secret_key.clone()); //todo insert for different pk and sk
        println!("End extract");
        secret_key
    }

    /// Generates an encryption of `message mod 2` for the provided public key by following these steps:
    /// e <- SampleD over lattice Z^m, center 0 with gaussian parameter r
    /// - u = A * e
    /// - c = p^t * e + message *  ⌊q/2⌋
    ///
    /// Then, `cipher = (u, c)` is output.
    ///
    /// Parameters:
    /// - `pk`: specifies the public key, which contains two matrices `pk = (A, p)`
    /// - `message`: specifies the message that should be encryted
    ///
    /// Returns a cipher of the form `cipher = (u, c)` for [`MatZq`] `u` and [`Zq`] `c`.
    ///
    /// # Examples
    /// ```
    /// use qfall_crypto::construction::pk_encryption::{PKEncryption, DualRegev};
    /// let dual_regev = DualRegev::default();
    /// let (pk, sk) = dual_regev.gen();
    ///
    /// let cipher = dual_regev.enc(&pk, 1);
    /// ```
    fn enc(
        &self,
        master_pk: &Self::MasterPublicKey,
        identity: &Self::Identity,
        message: impl Into<Z>,
    ) -> Self::Cipher {
        println!("Begin Enc");

        let identity_based_pk = hash_to_mat_zq_sha256(&identity, &self.n, 1, &self.q);
        println!("{}\n identibased_pk: {}", master_pk, identity_based_pk);

        let message = Zq::from((message, 2));
        let message = message.get_value();

        // s <- Z_q^n
        let vec_s_t = MatZq::sample_uniform(1, &self.n, &self.q);
        // e <- χ^(m+1)
        let vec_e_t = MatZq::sample_discrete_gauss(
            1,
            (master_pk.get_num_columns() + 1),
            &self.q,
            &self.n,
            0,
            &self.alpha * Z::from(&self.q),
        )
        .unwrap();
        println!("\n{vec_e_t}\n\n");
        let mut msg_q_half = MatZq::identity(1, 1, &self.q);
        msg_q_half
            .set_entry(0, 0, message * Z::from(&self.q).div_floor(&Z::from(2)))
            .unwrap();
        let message_entry = (&vec_s_t * identity_based_pk + msg_q_half);

        // c^t = s^t * A + e^t + [0^{1xn} | msg *  ⌊q/2⌋]
        let mut c = ((vec_s_t * master_pk)
            .concat_horizontal(&message_entry)
            .unwrap()
            + vec_e_t)
            .transpose();
        println!("\nmatrix done\n\n");
        // hide message in last entry
        // compute msg * ⌊q/2⌋

        // set last entry of c = last_entry + msg * ⌊q/2⌋

        c
    }

    /// Decrypts the provided `cipher` using the secret key `sk` by following these steps:
    /// - x = c - s^t * u
    /// - if x mod q is closer to ⌊q/2⌋ than to 0, output 1. Otherwise, output 0.
    ///
    /// Parameters:
    /// - `sk`: specifies the secret key `sk = s`
    /// - `cipher`: specifies the cipher containing `cipher = (u, c)`
    ///
    /// Returns the decryption of `cipher` as a [`Z`] instance.
    ///
    /// # Examples
    /// ```
    /// use qfall_crypto::construction::pk_encryption::{PKEncryption, DualRegev};
    /// use qfall_math::integer::Z;
    /// let dual_regev = DualRegev::default();
    /// let (pk, sk) = dual_regev.gen();
    /// let cipher = dual_regev.enc(&pk, 1);
    ///
    /// let m = dual_regev.dec(&sk, &cipher);
    ///
    /// assert_eq!(Z::ONE, m);
    /// ```
    fn dec(&self, sk: &Self::SecretKey, cipher: &Self::Cipher) -> Z {
        println!("Begin Dec");

        let tmp = (Z::MINUS_ONE * sk)
            .concat_vertical(&MatZ::identity(1, 1))
            .unwrap();
        let result: Zq = (cipher.transpose() * tmp).get_entry(0, 0).unwrap();

        let q_half = Z::from(&self.q).div_floor(&Z::from(2));

        if result.distance(Z::ZERO) > result.distance(q_half) {
            Z::ONE
        } else {
            Z::ZERO
        }
    }
}

impl Default for DualRegevIBE {
    /// Initializes a [`DualRegevIBE`] struct with parameters generated by `DualRegev::new_from_n(2)`.
    /// This parameter choice is not secure as the dimension of the lattice is too small,
    /// but it provides an efficient working example.
    ///
    /// # Examples
    /// ```
    /// use qfall_crypto::construction::pk_encryption::DualRegev;
    ///
    /// let dual_regev = DualRegev::default();
    /// ```
    fn default() -> Self {
        let n = Z::from(2);
        let m = Z::from(16);
        let q = Modulus::from(443);
        let r = Q::from(4);
        let alpha = Q::from((1, 64));

        Self::new(n, m, q, r, alpha) //todo
    }
}

impl DualRegevIBE {
    pub fn new(
        n: Z,       // security parameter
        m: Z,       // number of rows of matrix A
        q: Modulus, // modulus
        r: Q,       // gaussian parameter for sampleD
        alpha: Q,   // gaussian parameter for sampleZ
    ) -> Self {
        let gadget = GadgetParameters::init_default(&n, &q);
        let psf = PSFGPV {
            gp: gadget,
            s: Q::from(20),
        }; //todo standard deviation
        let gadget = GadgetParameters::init_default(&n, &q);
        Self {
            n: n,
            m: gadget.m_bar,
            q: q,
            r: r,
            alpha: alpha,
            psf: psf,
            storage: HashMap::new(),
        }
    }
}

#[cfg(test)]
mod test_dual_regev_ibe {
    use super::DualRegevIBE;
    use crate::construction::identity_based_encryption::IBE;
    use qfall_math::integer::Z;

    /// Checks whether the full-cycle of gen, enc, dec works properly
    /// for message 0 and small n.
    #[test]
    fn cycle_zero_small_n() {
        let msg = Z::ZERO;
        let id = String::from("Hello World!");
        let mut cryptosystem = DualRegevIBE::default();

        let (pk, sk) = cryptosystem.gen();
        let id_sk = cryptosystem.extract(&pk, &sk, &id);
        let cipher = cryptosystem.enc(&pk, &id, &msg);
        let m = cryptosystem.dec(&id_sk, &cipher);
        assert_eq!(msg, m);
    }

    /// for message 0 and small n.
    #[test]
    fn cycle_zero_small_n12() {
        let msg = Z::ZERO;
        let id = String::from("Hel213lo World!");
        let mut cryptosystem = DualRegevIBE::default();

        let (pk, sk) = cryptosystem.gen();
        let id_sk = cryptosystem.extract(&pk, &sk, &id);
        let cipher = cryptosystem.enc(&pk, &id, &msg);
        let m = cryptosystem.dec(&id_sk, &cipher);
        assert_eq!(msg, m);
    }
}
