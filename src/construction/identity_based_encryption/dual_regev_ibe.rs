// Copyright © 2023 Phil Milewski
//
// This file is part of qFALL-crypto.
//
// qFALL-crypto is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

//! This module contains an implementation of the IND-CPA secure
//! identity based public key encryption scheme. The encryption scheme is based
//! on [`DualRegevIBE`].

use super::IBE;
use crate::{
    primitive::hash::hash_to_mat_zq_sha256,
    sample::{
        distribution::psf::{gpv::PSFGPV, PSF},
        g_trapdoor::gadget_parameters::GadgetParameters,
    },
};
use qfall_math::{
    error::MathError,
    integer::{MatZ, Z},
    integer_mod_q::{MatZq, Modulus, Zq},
    rational::Q,
    traits::{Concatenate, Distance, GetEntry, Pow, SetEntry},
};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// This struct manages and stores the public parameters of a [`IBE`]
/// public key encryption instance based on [\[1\]](<index.html#:~:text=[1]>).
///
/// Attributes:
/// - `n`: specifies the security parameter, which is not equal to the bit-security level
/// - `m`: defines the dimension of the underlying lattice
/// - `q`: specifies the modulus over which the encryption is computed
/// - `r`: specifies the gaussian parameter used by the [`PSF`],
///   i.e. used for encryption
/// - `alpha`: specifies the gaussian parameter used for independent
///   sampling from χ, i.e. for multiple discrete Gaussian samples used
///   for key generation
///
/// # Examples
/// ```
/// use qfall_crypto::construction::identity_based_encryption::{DualRegevIBE, IBE};
/// use qfall_math::integer::Z;
/// // setup public parameters and key pair
/// let mut ibe = DualRegevIBE::default();
/// let (pk, sk) = ibe.gen();
///
/// // extract a identity based secret key
/// let identity = String::from("identity");
/// let id_sk = ibe.extract(&pk, &sk, &identity);
///
/// // encrypt a bit
/// let msg = Z::ZERO; // must be a bit, i.e. msg = 0 or 1
/// let cipher = ibe.enc(&pk, &identity, &msg);
///
/// // decrypt
/// let m = ibe.dec(&id_sk, &cipher);
///
/// assert_eq!(msg, m)
/// ```
#[derive(Serialize, Deserialize)]
pub struct DualRegevIBE {
    n: Z,       // security parameter
    m: Z,       // number of rows of matrix A
    q: Modulus, // modulus
    r: Q,       // gaussian parameter the [`PSF`]
    alpha: Q,   // gaussian parameter the encryption
    psf: PSFGPV,
    storage: HashMap<String, MatZ>,
}

impl DualRegevIBE {
    /// Initializes a [`DualRegevIBE`] struct with parameters generated by
    /// `DualRegev::new(n, q, r, alpha)`
    ///
    /// Returns an [`DualRegevIBE`] instance.
    ///
    /// # Examples
    /// ```
    /// use qfall_crypto::construction::identity_based_encryption::DualRegevIBE;
    ///
    /// let ibe = DualRegevIBE::new(4, 54983, 14, 0.0025);
    /// ```
    pub fn new(
        n: impl Into<Z>,       // security parameter
        q: impl Into<Modulus>, // modulus
        r: impl Into<Q>,       // gaussian parameter for sampleD
        alpha: impl Into<Q>,   // gaussian parameter for sampleZ
    ) -> Self {
        let n = n.into();
        let q = q.into();
        let r = r.into();
        let alpha = alpha.into();

        let gadget = GadgetParameters::init_default(&n, &q);

        let log_q = Z::from(&q).log_ceil(2).unwrap();
        let n_log_q = &n * &log_q;
        let m = &gadget.m_bar + n_log_q;

        let psf = PSFGPV {
            gp: gadget,
            s: r.clone(),
        };
        Self {
            n: n,
            m: m,
            q: q,
            r: r,
            alpha: alpha,
            psf: psf,
            storage: HashMap::new(),
        }
    }

    /// Initializes a [`DualRegevIBE`] struct with parameters generated by `DualRegev::new_from_n(n)`.
    ///
    /// **WARNING:** Due to the [`PSF`] this schemes extract algorithm is slow for n > 5.
    ///
    /// Returns an [`DualRegevIBE`] instance.
    ///
    /// # Examples
    /// ```
    /// use qfall_crypto::construction::identity_based_encryption::DualRegevIBE;
    ///
    /// let dual_regev = DualRegevIBE::new_from_n(4);
    /// ```
    pub fn new_from_n(n: impl Into<Z>) -> Self {
        let n: Z = n.into();
        if n < Z::from(2) {
            panic!("Security parameter n has to be larger than 1");
        }

        let n_i64 = i64::try_from(&n).unwrap();
        // these powers are chosen according to experience s.t. at least every
        // fifth generation of public parameters outputs a valid pair
        // the exponent is only tested for n < 8
        let power = match n_i64 {
            2..=3 => 10,
            4 => 7,
            5..=7 => 6,
            _ => 5,
        };

        // generate prime q in [n^power / 2, n^power]
        let upper_bound: Z = n.pow(power).unwrap();
        let lower_bound = upper_bound.div_ceil(&Z::from(2));
        // prime used due to guide from GPV08 after Proposition 8.1
        // on how to choose appropriate parameters, but prime is not
        // necessarily needed for this scheme to be correct or secure
        let q = Modulus::from(Z::sample_prime_uniform(&lower_bound, &upper_bound).unwrap());

        let gadget = GadgetParameters::init_default(&n, &q);
        let log_q = Z::from(&q).log_ceil(2).unwrap();
        let n_log_q = &n * &log_q;
        let m = &gadget.m_bar + n_log_q;
        let r: Q = m.sqrt();
        let alpha = 1 / (&r * 2 * (&m + Z::ONE).sqrt() * (&n).log(2).unwrap());

        let psf = PSFGPV {
            gp: gadget,
            s: r.clone(),
        };
        Self {
            n: n.into(),
            m: m,
            q: q,
            r: r,
            alpha: alpha,
            psf: psf,
            storage: HashMap::new(),
        }
    }

    /// Checks the public parameters for security according to Theorem 1.1
    /// and Lemma 5.4 of [\[2\]](<index.html#:~:text=[2]>), as well as
    /// the requirements of [\[1\]](<index.html#:~:text=[1]>)`s eprint version
    /// at Section 7.1 of [GPV08 - eprint](https://eprint.iacr.org/2007/432.pdf).
    ///
    /// The required properties are:
    /// - q >= 5 * r * (m + 1)
    /// - r >= sqrt(m)
    /// - m > (n + 1) * log(q)
    ///
    /// Returns an empty result if the public parameters guarantees security w.r.t. `n`
    /// or a [`MathError`] if the instance would not be secure.
    ///
    /// # Examples
    /// ```
    /// use qfall_crypto::construction::identity_based_encryption::DualRegevIBE;
    /// let ibe = DualRegevIBE::default();
    ///
    /// assert!(ibe.check_security().is_ok());
    /// ```
    ///
    /// # Errors and Failures
    /// - Returns a [`MathError`] of type [`InvalidIntegerInput`](MathError::InvalidIntegerInput)
    /// if at least one parameter was not chosen appropriately for a
    /// secure Dual Regev public key encryption instance.
    pub fn check_security(&self) -> Result<(), MathError> {
        let q = Q::from(&self.q);

        // Security requirements
        // q >= 5 * r * (m + 1)
        if &q < &((5 * &self.r) * (&self.m + Q::ONE)) {
            return Err(MathError::InvalidIntegerInput(String::from(
                "Security is not guaranteed as q < 5 * r * (m + 1), but q >= 5 * r * (m + 1) is required.",
            )));
        }

        // r >= sqrt(m)
        if &self.r < &self.m.sqrt() {
            return Err(MathError::InvalidIntegerInput(String::from(
                "Security is not guaranteed as r < sqrt(m), but r >= sqrt(m) is required.",
            )));
        }

        // m >= (n + 1) * log(q)
        if &Q::from(&self.m) <= &((&self.n + 1) * &q.log(2).unwrap()) {
            return Err(MathError::InvalidIntegerInput(String::from(
        "Security is not guaranteed as m <= (n + 1) * log(q), but m > (n + 1) * log(q) is required.",
    )));
        }

        Ok(())
    }

    /// Checks the public parameters for
    /// correctness according to Lemma 5.1 of [\[2\]](<index.html#:~:text=[2]>).
    ///
    /// The required properties are:
    /// - α <= 1/(r *sqrt(m)* log(n)
    ///
    /// **WARNING:** Some requirements are missing to ensure overwhelming correctness of the scheme.
    ///
    /// Returns an empty result if the public parameters guarantee correctness
    /// with overwhelming probability or a [`MathError`] if the instance would
    /// not be correct.
    ///
    /// # Examples
    /// ```
    /// use qfall_crypto::construction::identity_based_encryption::DualRegevIBE;
    /// let ibe = DualRegevIBE::default();
    ///
    /// assert!(ibe.check_correctness().is_ok());
    /// ```
    ///
    /// # Errors and Failures
    /// - Returns a [`MathError`] of type [`InvalidIntegerInput`](MathError::InvalidIntegerInput)
    /// if at least one parameter was not chosen appropriately for a
    /// correct Dual Regev IBE public key encryption instance.
    pub fn check_correctness(&self) -> Result<(), MathError> {
        if self.n <= Z::ONE {
            return Err(MathError::InvalidIntegerInput(String::from(
                "n must be chosen bigger than 1.",
            )));
        }

        // α <= 1/(r *sqrt(m)* log(n))
        if &self.alpha > &(1 / (2 * &self.r * (&self.m + Z::ONE).sqrt()) * self.n.log(2).unwrap()) {
            return Err(MathError::InvalidIntegerInput(String::from(
                "Correctness is not guaranteed as α > 1/(r *sqrt(m)* log(n)), but α <= 1/(2 * r *sqrt(m)* log(n)) is required.",
            )));
        }

        Ok(())
    }
}

impl Default for DualRegevIBE {
    /// Initializes a [`DualRegevIBE`] struct with parameters generated by `DualRegevIBE::new_from_n(4)`.
    /// This parameter choice is not secure as the dimension of the lattice is too small,
    /// but it provides an efficient working example.
    ///
    /// Returns an [`DualRegevIBE`] instance.
    ///
    /// # Examples
    /// ```
    /// use qfall_crypto::construction::identity_based_encryption::DualRegevIBE;
    ///
    /// let ibe = DualRegevIBE::default();
    /// ```
    fn default() -> Self {
        DualRegevIBE::new_from_n(4)
    }
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
    /// use qfall_crypto::construction::identity_based_encryption::{DualRegevIBE, IBE};
    /// let ibe = DualRegevIBE::default();
    ///
    /// let (pk, sk) = ibe.gen();
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
    ///
    /// let id = String::from("identity");
    /// let sk = ibe.extract(&master_pk, &master_sk, &id);
    /// ```
    fn extract(
        &mut self,
        master_pk: &Self::MasterPublicKey,
        master_sk: &Self::MasterSecretKey,
        identity: &Self::Identity,
    ) -> Self::SecretKey {
        // check if it is in the HashMap
        if let Some(value) = self.storage.get(identity) {
            return value.clone();
        }

        let u = hash_to_mat_zq_sha256(&identity, &self.n, 1, &self.q);
        let secret_key = self.psf.samp_p(&master_pk, &master_sk, &u);

        // insert secret key in HashMap
        self.storage.insert(identity.clone(), secret_key.clone()); //todo insert for different pk and sk

        secret_key
    }

    /// Generates an encryption of `message mod 2` for the provided public key
    /// and identity by following these steps:
    /// e = <- χ^(m+1)
    /// p = H(id)
    /// - u = A * e
    /// - c = p^t * e + message *  ⌊q/2⌋
    ///
    /// Then, `cipher = [u | c]` is output.
    ///
    /// Parameters:
    /// - `master_pk`: specifies the public key, which contains two matrices `pk = (A, p)`
    /// - `identity`: specifies the identity used for encryption
    /// - `message`: specifies the message that should be encrypted
    ///
    /// Returns a cipher of the form `cipher = [u | c]` for [`MatZq`] `u`.
    ///
    /// # Examples
    /// ```
    /// use qfall_crypto::construction::identity_based_encryption::{DualRegevIBE, IBE};
    /// let ibe = DualRegevIBE::default();
    /// let (pk, sk) = ibe.gen();
    ///
    /// let id = String::from("identity");
    /// let cipher = ibe.enc(&pk, &id, 1);
    /// ```
    fn enc(
        &self,
        master_pk: &Self::MasterPublicKey,
        identity: &Self::Identity,
        message: impl Into<Z>,
    ) -> Self::Cipher {
        let identity_based_pk = hash_to_mat_zq_sha256(&identity, &self.n, 1, &self.q);

        let message = Zq::from((message, 2));
        let message = message.get_value();

        // s <- Z_q^n
        let vec_s_t = MatZq::sample_uniform(1, &self.n, &self.q);
        // e <- χ^(m+1)
        let vec_e_t = MatZq::sample_discrete_gauss(
            1,
            &self.m + 1,
            &self.q,
            &self.n,
            0,
            &self.alpha * Z::from(&self.q),
        )
        .unwrap();

        // p^t * e + message *  ⌊q/2⌋
        let mut msg_q_half = MatZq::identity(1, 1, &self.q);
        msg_q_half
            .set_entry(0, 0, message * Z::from(&self.q).div_floor(&Z::from(2)))
            .unwrap();
        let message_entry = &vec_s_t * identity_based_pk + msg_q_half;

        // [A * e | c]`
        let c = ((vec_s_t * master_pk)
            .concat_horizontal(&message_entry)
            .unwrap()
            + vec_e_t)
            .transpose();

        c
    }

    /// Decrypts the provided `cipher` using the secret key `sk` by following these steps:
    /// - u = H(id)
    /// - x = c - s^t * u
    /// - if x mod q is closer to ⌊q/2⌋ than to 0, output 1. Otherwise, output 0.
    ///
    /// Parameters:
    /// - `sk_id`: specifies the secret key `sk = s` obtained by extract
    /// - `cipher`: specifies the cipher containing `cipher = [u | c]`
    ///
    /// Returns the decryption of `cipher` as a [`Z`] instance.
    ///
    /// # Examples
    /// ```
    /// use qfall_crypto::construction::identity_based_encryption::{DualRegevIBE, IBE};
    /// use qfall_math::integer::Z;
    /// // setup public parameters and key pair
    /// let mut ibe = DualRegevIBE::default();
    /// let (pk, sk) = ibe.gen();
    ///
    /// // extract a identity based secret key
    /// let identity = String::from("identity");
    /// let id_sk = ibe.extract(&pk, &sk, &identity);
    ///
    /// // encrypt a bit
    /// let msg = Z::ZERO; // must be a bit, i.e. msg = 0 or 1
    /// let cipher = ibe.enc(&pk, &identity, &msg);
    ///
    /// // decrypt
    /// let m = ibe.dec(&id_sk, &cipher);
    ///
    /// assert_eq!(msg, m)
    /// ```
    fn dec(&self, sk_id: &Self::SecretKey, cipher: &Self::Cipher) -> Z {
        let tmp = (Z::MINUS_ONE * sk_id)
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

#[cfg(test)]
mod test_dual_regev_ibe {
    use super::DualRegevIBE;
    use crate::construction::identity_based_encryption::IBE;
    use qfall_math::integer::Z;

    /// Checks whether `new` is available for types implementing [`Into<Z>`].
    #[test]
    fn new_availability() {
        let _ = DualRegevIBE::new(2u8, 2u16, 2u32, 2u64);
        let _ = DualRegevIBE::new(2u16, 2u64, 2i32, 2i64);
        let _ = DualRegevIBE::new(2i16, 2i64, 2u32, 2u8);
        let _ = DualRegevIBE::new(Z::from(2), &Z::from(2), 2u8, 2i8);
    }

    /// Ensures that `new_from_n` is available for types implementing [`Into<Z>`].
    #[test]
    fn availability() {
        let _ = DualRegevIBE::new_from_n(4u8);
        let _ = DualRegevIBE::new_from_n(4u16);
        let _ = DualRegevIBE::new_from_n(4u32);
        let _ = DualRegevIBE::new_from_n(4u64);
        let _ = DualRegevIBE::new_from_n(4i8);
        let _ = DualRegevIBE::new_from_n(4i16);
        let _ = DualRegevIBE::new_from_n(4i32);
        let _ = DualRegevIBE::new_from_n(4i64);
        let _ = DualRegevIBE::new_from_n(Z::from(4));
        let _ = DualRegevIBE::new_from_n(&Z::from(4));
    }

    /// Checks whether `new_from_n` returns an error for invalid input n.
    #[test]
    #[should_panic]
    fn invalid_n() {
        DualRegevIBE::new_from_n(1);
    }

    /// Checks whether the full-cycle of gen, extract, enc, dec works properly
    /// for message 0 and the default.
    #[test]
    fn cycle_zero_default() {
        let msg = Z::ZERO;
        let id = String::from("Hello World!");
        let mut cryptosystem = DualRegevIBE::default();

        let (pk, sk) = cryptosystem.gen();
        let id_sk = cryptosystem.extract(&pk, &sk, &id);
        let cipher = cryptosystem.enc(&pk, &id, &msg);
        let m = cryptosystem.dec(&id_sk, &cipher);

        assert_eq!(msg, m)
    }

    /// Checks whether the full-cycle of gen, extract, enc, dec works properly
    /// for message 0 and small n.
    #[test]
    fn cycle_zero_small_n() {
        let msg = Z::ZERO;
        let id = String::from("Hel213lo World!");
        let mut cryptosystem = DualRegevIBE::new_from_n(5);

        let (pk, sk) = cryptosystem.gen();
        let id_sk = cryptosystem.extract(&pk, &sk, &id);
        let cipher = cryptosystem.enc(&pk, &id, &msg);
        let m = cryptosystem.dec(&id_sk, &cipher);
        assert_eq!(msg, m);
    }

    /// for message 0 and small n.
    #[test]
    fn new_from_n() {
        for i in 1..=5 {
            let msg = Z::ONE;
            let id = String::from(format!("Hello World!{i}"));
            let mut cryptosystem = DualRegevIBE::default();

            cryptosystem.check_security().unwrap();
            cryptosystem.check_correctness().unwrap();

            let (pk, sk) = cryptosystem.gen();

            let id_sk = cryptosystem.extract(&pk, &sk, &id);
            for _j in 1..=100 {
                let cipher = cryptosystem.enc(&pk, &id, &msg);
                let m = cryptosystem.dec(&id_sk, &cipher);

                assert_eq!(msg, m);
            }
        }
    }
}
