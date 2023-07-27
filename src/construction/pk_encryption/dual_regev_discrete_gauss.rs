// Copyright © 2023 Niklas Siemer
//
// This file is part of qFALL-crypto.
//
// qFALL-crypto is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

//! This module contains an implementation of the IND-CPA secure
//! public key Dual Regev encryption scheme with an instantiation of the regularity lemma
//! via a discrete Gaussian distribution.
//!
//! The main references are listed in the following:
//! - \[1\] Gentry, Craig and Peikert, Chris and Vaikuntanathan, Vinod (2008).
//! Trapdoors for hard lattices and new cryptographic constructions.
//! In: Proceedings of the fortieth annual ACM symposium on Theory of computing.
//! <https://dl.acm.org/doi/pdf/10.1145/1374376.1374407>

use super::PKEncryption;
use qfall_math::{
    error::MathError,
    integer::Z,
    integer_mod_q::{MatZq, Modulus, Zq},
    rational::Q,
    traits::{Distance, Pow},
};
use serde::{Deserialize, Serialize};

/// This struct manages and stores the public parameters of a [`DualRegevWithDiscreteGaussianRegularity`]
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
///
/// # Examples
/// ```
/// use qfall_crypto::construction::pk_encryption::{DualRegevWithDiscreteGaussianRegularity, PKEncryption};
/// use qfall_math::integer::Z;
/// // setup public parameters and key pair
/// let dual_regev = DualRegevWithDiscreteGaussianRegularity::default();
/// let (pk, sk) = dual_regev.gen();
///
/// // encrypt a bit
/// let msg = Z::ZERO; // must be a bit, i.e. msg = 0 or 1
/// let cipher = dual_regev.enc(&pk, &msg);
///
/// // decrypt
/// let m = dual_regev.dec(&sk, &cipher);
///
/// assert_eq!(msg, m);
/// ```
#[derive(Debug, Serialize, Deserialize)]
pub struct DualRegevWithDiscreteGaussianRegularity {
    n: Z,       // security parameter
    m: Z,       // number of rows of matrix A
    q: Modulus, // modulus
    r: Q,       // gaussian parameter for sampleD
    alpha: Q,   // gaussian parameter for sampleZ
}

impl DualRegevWithDiscreteGaussianRegularity {
    /// Instantiates a [`DualRegevWithDiscreteGaussianRegularity`] PK encryption instance with the
    /// specified parameters if they ensure a secure and correct instance.
    ///
    /// **WARNING:** The given parameters are not checked for security nor
    /// correctness of the scheme.
    /// If you want to check your parameters for provable security and correctness,
    /// use [`DualRegevWithDiscreteGaussianRegularity::check_correctness`] and [`DualRegevWithDiscreteGaussianRegularity::check_security`].
    /// Or use [`DualRegevWithDiscreteGaussianRegularity::new_from_n`] for generating secure and correct
    /// public parameters for [`DualRegevWithDiscreteGaussianRegularity`] according to your choice of `n`.
    ///
    /// Parameters:
    /// - `n`: specifies the security parameter and number of rows
    ///   of the uniform at random instantiated matrix `A`
    /// - `m`: specifies the number of columns of matrix `A`
    /// - `q`: specifies the modulus
    /// - `r`: specifies the gaussian parameter used for SampleD,
    ///   i.e. used for encryption
    /// - `alpha`: specifies the gaussian parameter used for independent
    ///   sampling from χ, i.e. for multiple discrete Gaussian samples used
    ///   for key generation
    ///
    /// Returns a correct and secure [`DualRegevWithDiscreteGaussianRegularity`] PK encryption instance or
    /// a [`MathError`] if the instance would not be correct or secure.
    ///
    /// # Examples
    /// ```
    /// use qfall_crypto::construction::pk_encryption::DualRegevWithDiscreteGaussianRegularity;
    ///
    /// let dual_regev = DualRegevWithDiscreteGaussianRegularity::new(2, 16, 443, 4, 0.15625);
    /// ```
    ///
    /// # Panics ...
    /// - if the given modulus `q <= 1`.
    pub fn new(
        n: impl Into<Z>,
        m: impl Into<Z>,
        q: impl Into<Z>,
        r: impl Into<Q>,
        alpha: impl Into<Q>,
    ) -> Self {
        let n: Z = n.into();
        let m: Z = m.into();
        let q: Z = q.into();
        let r: Q = r.into();
        let alpha: Q = alpha.into();

        let q = Modulus::from(&q);

        Self { n, m, q, r, alpha }
    }

    /// Generates a new [`DualRegevWithDiscreteGaussianRegularity`] instance, i.e. a new set of suitable
    /// (provably secure and correct) public parameters,
    /// given the security parameter `n`.
    ///
    /// Parameters:
    /// - `n`: specifies the security parameter and number of rows
    ///   of the uniform at random instantiated matrix `A`
    ///
    /// Returns a correct and secure [`DualRegevWithDiscreteGaussianRegularity`] PK encryption instance or
    /// a [`MathError`] if the given `n <= 1`.
    ///
    /// # Examples
    /// ```
    /// use qfall_crypto::construction::pk_encryption::DualRegevWithDiscreteGaussianRegularity;
    ///
    /// let dual_regev = DualRegevWithDiscreteGaussianRegularity::new_from_n(2).unwrap();
    /// ```
    ///
    /// # Errors and Failures
    /// - Returns a [`MathError`] of type [`InvalidIntegerInput`](MathError::InvalidIntegerInput)
    /// if `n <= 1`.
    pub fn new_from_n(n: impl Into<Z>) -> Result<Self, MathError> {
        let n = n.into();
        if n <= Z::ONE {
            return Err(MathError::InvalidIntegerInput(String::from(
                "n must be chosen bigger than 1.",
            )));
        }

        let mut m: Z;
        let mut q: Modulus;
        let mut r: Q;
        let mut alpha: Q;
        (m, q, r, alpha) = Self::gen_new_public_parameters(&n);
        let mut out = Self {
            n: n.clone(),
            m,
            q,
            r,
            alpha,
        };
        while out.check_correctness().is_err() || out.check_security().is_err() {
            (m, q, r, alpha) = Self::gen_new_public_parameters(&n);
            out = Self {
                n: n.clone(),
                m,
                q,
                r,
                alpha,
            };
        }

        Ok(out)
    }

    /// Generates new public parameters, which must not be secure or correct
    /// depending on the random choice of `q`. At least every fifth execution
    /// of this function should output a valid set of public parameters,
    /// ensuring a secure and correct PK encryption scheme.
    ///
    /// Parameters:
    /// - `n`: specifies the security parameter and number of rows
    ///   of the uniform at random instantiated matrix `A`
    ///
    /// Returns a set of public parameters `(m, q, r, alpha)` chosen according to
    /// the provided `n`.
    ///
    /// # Examples
    /// ```compile_fail
    /// use qfall_crypto::construction::pk_encryption::DualRegevWithDiscreteGaussianRegularity;
    /// use qfall_math::integer::Z;
    /// let n = Z::from(2);
    ///
    /// let (m, q, r, alpha) = DualRegevWithDiscreteGaussianRegularity::gen_new_public_parameters(&n);
    /// ```
    fn gen_new_public_parameters(n: &Z) -> (Z, Modulus, Q, Q) {
        let n_i64 = i64::try_from(n).unwrap();
        // these powers are chosen according to experience s.t. at least every
        // fifth generation of public parameters outputs a valid pair
        let power = match n_i64 {
            2 => 9,
            3 => 8,
            4..=5 => 7,
            6..=8 => 6,
            9..=12 => 5,
            13..=30 => 4,
            _ => 3,
        };

        // generate prime q in [n^power / 2, n^power]
        let upper_bound: Z = n.pow(power).unwrap();
        let lower_bound = upper_bound.div_ceil(&Z::from(2));
        let q = Z::sample_prime_uniform(&lower_bound, &upper_bound).unwrap();

        // choose m = 2 (n+1) lg q
        let m = (Z::from(2) * (n + Z::ONE) * q.log(10).unwrap()).ceil();

        // choose r = log m
        let r = m.log(2).unwrap();

        // alpha = 1/(sqrt(m) * log^2 m)
        let alpha = 1 / (m.sqrt() * m.log(2).unwrap().pow(2).unwrap());

        let q = Modulus::from(&q);

        (m, q, r, alpha)
    }

    /// Checks a provided set of public parameters according to their validity
    /// regarding correctness and completeness according to
    /// Theorem 7.1 and Lemma 8.2 of [\[1\]](<index.html#:~:text=[1]>).
    ///
    /// Returns an empty result or a [`MathError`] if the instance would
    /// not be correct.
    ///
    /// # Examples
    /// ```
    /// use qfall_crypto::construction::pk_encryption::DualRegevWithDiscreteGaussianRegularity;
    /// let dr = DualRegevWithDiscreteGaussianRegularity::default();
    ///
    /// let is_valid = dr.check_correctness().is_ok();
    /// ```
    ///
    /// # Errors and Failures
    /// - Returns a [`MathError`] of type [`InvalidIntegerInput`](MathError::InvalidIntegerInput)
    /// if at least one parameter was not chosen appropriately for a
    /// correct DualRegevWithDiscreteGaussianRegularity public key encryption instance.
    pub fn check_correctness(&self) -> Result<(), MathError> {
        let q: Z = Z::from(&self.q);

        if self.n <= Z::ONE {
            return Err(MathError::InvalidIntegerInput(String::from(
                "n must be chosen bigger than 1.",
            )));
        }
        if !q.is_prime() {
            return Err(MathError::InvalidIntegerInput(String::from(
                "q must be prime.",
            )));
        }

        // Completeness requirements
        // q >= 5 * r * (m+1)
        if Q::from(q) < 5 * &self.r * (&self.m + Z::ONE) {
            return Err(MathError::InvalidIntegerInput(String::from(
                "Completeness is not guaranteed as q < 5rm, but q >= 5rm is required.",
            )));
        }
        // α <= 1/(r * sqrt(m+1) * ω(sqrt(log n))
        if self.alpha > 1 / (&self.r * (&self.m + Z::ONE).sqrt() * self.n.log(2).unwrap().sqrt()) {
            return Err(MathError::InvalidIntegerInput(String::from(
                "Completeness is not guaranteed as α > 1/(r*sqrt(m)*ω(sqrt(log n)), but α <= 1/(r*sqrt(m)*ω(sqrt(log n)) is required.",
            )));
        }

        Ok(())
    }

    /// Checks a provided set of public parameters according to their validity
    /// regarding security according to
    /// Theorem 7.1 and Lemma 8.4 of [\[1\]](<index.html#:~:text=[1]>).
    ///
    /// Returns an empty result or a [`MathError`] if the instance would
    /// not be secure.
    ///
    /// # Examples
    /// ```
    /// use qfall_crypto::construction::pk_encryption::DualRegevWithDiscreteGaussianRegularity;
    /// let dr = DualRegevWithDiscreteGaussianRegularity::default();
    ///
    /// let is_valid = dr.check_security().is_ok();
    /// ```
    ///
    /// # Errors and Failures
    /// - Returns a [`MathError`] of type [`InvalidIntegerInput`](MathError::InvalidIntegerInput)
    /// if at least one parameter was not chosen appropriately for a
    /// secure DualRegevWithDiscreteGaussianRegularity public key encryption instance.
    pub fn check_security(&self) -> Result<(), MathError> {
        let q: Z = Z::from(&self.q);

        // Security requirements
        // q * α >= n
        if &q * &self.alpha < Q::from(&self.n) {
            return Err(MathError::InvalidIntegerInput(String::from(
                "Security is not guaranteed as q * α < n, but q * α >= n is required.",
            )));
        }
        // m >= 2n lg (q)
        if Q::from(&self.m) < 2 * &self.n * q.log(10).unwrap() {
            return Err(MathError::InvalidIntegerInput(String::from(
                "Security is not guaranteed as m < 2(n + 1) lg (q), but m >= 2(n + 1) lg (q) is required.",
            )));
        }
        // r >= ω( sqrt( log m ) )
        if self.r < self.m.log(2).unwrap().sqrt() {
            return Err(MathError::InvalidIntegerInput(String::from(
                "Security is not guaranteed as r < sqrt( log m ) and r >= ω(sqrt(log m)) is required."
            )));
        }

        Ok(())
    }

    /// This function instantiates a 128-bit secure [`DualRegevWithDiscreteGaussianRegularity`] scheme.
    ///
    /// The public parameters used for this scheme were generated
    /// via `DualRegevWithDiscreteGaussianRegularity::new_from_n(350)`
    /// and its bit-security determined via the [lattice estimator](https://github.com/malb/lattice-estimator).
    pub fn secure128() -> Self {
        Self::new(350, 5248, 29892991, 12.357, 0.00009)
    }
}

impl Default for DualRegevWithDiscreteGaussianRegularity {
    /// Initializes a [`DualRegevWithDiscreteGaussianRegularity`] struct with parameters
    /// generated by `DualRegevWithDiscreteGaussianRegularity::new_from_n(2)`.
    /// This parameter choice is not secure as the dimension of the lattice is too small,
    /// but it provides an efficient working example.
    ///
    /// # Examples
    /// ```
    /// use qfall_crypto::construction::pk_encryption::DualRegevWithDiscreteGaussianRegularity;
    ///
    /// let dual_regev = DualRegevWithDiscreteGaussianRegularity::default();
    /// ```
    fn default() -> Self {
        let n = Z::from(2);
        let m = Z::from(16);
        let q = Modulus::from(443);
        let r = Q::from(4);
        let alpha = Q::from((1, 64));

        Self { n, m, q, r, alpha }
    }
}

impl PKEncryption for DualRegevWithDiscreteGaussianRegularity {
    type Cipher = (MatZq, Zq);
    type PublicKey = (MatZq, MatZq);
    type SecretKey = MatZq;

    /// Generates a (pk, sk) pair for the Dual Regev public key encryption scheme
    /// by following these steps:
    /// - e <- SampleD over lattice Z^m, center 0 with gaussian parameter r
    /// - A <- Z_q^{n x m}
    /// - p = A * e
    ///
    /// Then, `pk = (A, u)` and `sk = e` is output.
    ///
    /// # Examples
    /// ```
    /// use qfall_crypto::construction::pk_encryption::{PKEncryption, DualRegevWithDiscreteGaussianRegularity};
    /// let dual_regev = DualRegevWithDiscreteGaussianRegularity::default();
    ///
    /// let (pk, sk) = dual_regev.gen();
    /// ```
    fn gen(&self) -> (Self::PublicKey, Self::SecretKey) {
        // e <- SampleD over lattice Z^m, center 0 with gaussian parameter r
        let vec_e = MatZq::sample_d_common(&self.m, &self.q, &self.n, &self.r).unwrap();
        // A <- Z_q^{n x m}
        let mat_a = MatZq::sample_uniform(&self.n, &self.m, &self.q);

        // u = A * e
        let vec_u = &mat_a * &vec_e;

        // pk = (A, u), sk = e
        ((mat_a, vec_u), vec_e)
    }

    /// Generates an encryption of `message mod 2` for the provided public key by following these steps:
    /// - vec_x <- χ^m, x <- χ
    /// - p = A^t * s + vec_x
    /// - c = u^t * s + x + message *  ⌊q/2⌋
    ///
    /// Then, `cipher = (p, c)` is output.
    ///
    /// Parameters:
    /// - `pk`: specifies the public key, which contains two matrices `pk = (A, u)`
    /// - `message`: specifies the message that should be encryted
    ///
    /// Returns a cipher of the form `cipher = (p, c)` for [`MatZq`] `u` and [`Zq`] `c`.
    ///
    /// # Examples
    /// ```
    /// use qfall_crypto::construction::pk_encryption::{PKEncryption, DualRegevWithDiscreteGaussianRegularity};
    /// let dual_regev = DualRegevWithDiscreteGaussianRegularity::default();
    /// let (pk, sk) = dual_regev.gen();
    ///
    /// let cipher = dual_regev.enc(&pk, 1);
    /// ```
    fn enc(&self, pk: &Self::PublicKey, message: impl Into<Z>) -> Self::Cipher {
        // generate message = message mod 2
        let message: Z = message.into();
        let message = Zq::from((&message, 2));
        let message = message.get_value();

        // s <- Z_q^n
        let vec_s = MatZq::sample_uniform(&self.n, 1, &self.q);
        // vec_x <- χ^m
        let vec_x = MatZq::sample_discrete_gauss(
            &self.m,
            1,
            &self.q,
            &self.n,
            0,
            &(&self.alpha * Z::from(&self.q)),
        )
        .unwrap();

        // x <- χ
        let x = Z::sample_discrete_gauss(&self.n, 0, &(&self.alpha * Z::from(&self.q))).unwrap();

        // p = u^t * s + vec_x
        let vec_p = &pk.0.transpose() * &vec_s + vec_x;
        // c = u^t * s + x + msg *  ⌊q/2⌋
        let q_half = Z::from(&self.q).div_floor(&Z::from(2));
        let c = pk.1.dot_product(&vec_s).unwrap() + x + message * q_half;

        (vec_p, c)
    }

    /// Decrypts the provided `cipher` using the secret key `sk` by following these steps:
    /// - x = c - e^t * p
    /// - if x mod q is closer to ⌊q/2⌋ than to 0, output 1. Otherwise, output 0.
    ///
    /// Parameters:
    /// - `sk`: specifies the secret key `sk = e`
    /// - `cipher`: specifies the cipher containing `cipher = (p, c)`
    ///
    /// Returns the decryption of `cipher` as a [`Z`] instance.
    ///
    /// # Examples
    /// ```
    /// use qfall_crypto::construction::pk_encryption::{PKEncryption, DualRegevWithDiscreteGaussianRegularity};
    /// use qfall_math::integer::Z;
    /// let dual_regev = DualRegevWithDiscreteGaussianRegularity::default();
    /// let (pk, sk) = dual_regev.gen();
    /// let cipher = dual_regev.enc(&pk, 1);
    ///
    /// let m = dual_regev.dec(&sk, &cipher);
    ///
    /// assert_eq!(Z::ONE, m);
    /// ```
    fn dec(&self, sk: &Self::SecretKey, cipher: &Self::Cipher) -> Z {
        let result = &cipher.1 - sk.dot_product(&cipher.0).unwrap();

        let q_half = Z::from(&self.q).div_floor(&Z::from(2));

        if result.distance(Z::ZERO) > result.distance(q_half) {
            Z::ONE
        } else {
            Z::ZERO
        }
    }
}

#[cfg(test)]
mod test_pp_generation {
    use super::DualRegevWithDiscreteGaussianRegularity;
    use super::Z;

    /// Checks whether `new` is available for types implementing [`Into<Z>`].
    #[test]
    fn new_availability() {
        let _ = DualRegevWithDiscreteGaussianRegularity::new(2u8, 2u16, 2u32, 2u64, 2i8);
        let _ = DualRegevWithDiscreteGaussianRegularity::new(2u16, 2u64, 2i32, 2i64, 2i16);
        let _ = DualRegevWithDiscreteGaussianRegularity::new(2i16, 2i64, 2u32, 2u8, 2u16);
        let _ =
            DualRegevWithDiscreteGaussianRegularity::new(Z::from(2), &Z::from(2), 2u8, 2i8, 2u32);
    }

    /// Checks whether `new_from_n` works properly for different choices of n.
    #[test]
    fn suitable_security_params() {
        let n_choices = [
            2, 3, 4, 5, 6, 7, 8, 9, 10, 25, 50, 75, 100, 250, 500, 750, 1000, 2500, 5000, 5001,
            10000,
        ];

        for n in n_choices {
            assert!(DualRegevWithDiscreteGaussianRegularity::new_from_n(n).is_ok());
        }
    }

    /// Checks whether the [`Default`] parameter choice is suitable.
    #[test]
    fn default_suitable() {
        let dr = DualRegevWithDiscreteGaussianRegularity::default();

        assert!(dr.check_correctness().is_ok());
        assert!(dr.check_security().is_ok());
    }

    /// Checks whether the generated public parameters from `new_from_n` are
    /// valid choices according to security and correctness of the scheme.
    #[test]
    fn choice_valid() {
        let n_choices = [
            2, 3, 4, 5, 6, 7, 8, 9, 10, 25, 50, 75, 100, 250, 500, 750, 1000, 2500, 5000, 5001,
            10000,
        ];

        for n in n_choices {
            let dr = DualRegevWithDiscreteGaussianRegularity::new_from_n(n).unwrap();
            assert!(dr.check_correctness().is_ok());
            assert!(dr.check_security().is_ok());
        }
    }

    /// Ensures that `new_from_n` is available for types implementing [`Into<Z>`].
    #[test]
    fn new_from_n_availability() {
        let _ = DualRegevWithDiscreteGaussianRegularity::new_from_n(2u8);
        let _ = DualRegevWithDiscreteGaussianRegularity::new_from_n(2u16);
        let _ = DualRegevWithDiscreteGaussianRegularity::new_from_n(2u32);
        let _ = DualRegevWithDiscreteGaussianRegularity::new_from_n(2u64);
        let _ = DualRegevWithDiscreteGaussianRegularity::new_from_n(2i8);
        let _ = DualRegevWithDiscreteGaussianRegularity::new_from_n(2i16);
        let _ = DualRegevWithDiscreteGaussianRegularity::new_from_n(2i32);
        let _ = DualRegevWithDiscreteGaussianRegularity::new_from_n(2i64);
        let _ = DualRegevWithDiscreteGaussianRegularity::new_from_n(Z::from(2));
        let _ = DualRegevWithDiscreteGaussianRegularity::new_from_n(&Z::from(2));
    }

    /// Checks whether `new_from_n` returns an error for invalid input n.
    #[test]
    fn invalid_n() {
        assert!(DualRegevWithDiscreteGaussianRegularity::new_from_n(1).is_err());
        assert!(DualRegevWithDiscreteGaussianRegularity::new_from_n(0).is_err());
        assert!(DualRegevWithDiscreteGaussianRegularity::new_from_n(-1).is_err());
    }

    /// Checks whether `secure128` outputs a new instance with correct and secure parameters.
    #[test]
    fn secure128_validity() {
        let dr = DualRegevWithDiscreteGaussianRegularity::secure128();

        assert!(dr.check_correctness().is_ok());
        assert!(dr.check_security().is_ok());
    }
}

#[cfg(test)]
mod test_dual_regev {
    use super::DualRegevWithDiscreteGaussianRegularity;
    use crate::construction::pk_encryption::PKEncryption;
    use qfall_math::integer::Z;

    /// Checks whether the full-cycle of gen, enc, dec works properly
    /// for message 0 and small n.
    #[test]
    fn cycle_zero_small_n() {
        let msg = Z::ZERO;
        let dr = DualRegevWithDiscreteGaussianRegularity::default();

        let (pk, sk) = dr.gen();
        let cipher = dr.enc(&pk, &msg);
        let m = dr.dec(&sk, &cipher);
        assert_eq!(msg, m);
    }

    /// Checks whether the full-cycle of gen, enc, dec works properly
    /// for message 1 and small n.
    #[test]
    fn cycle_one_small_n() {
        let msg = Z::ONE;
        let dr = DualRegevWithDiscreteGaussianRegularity::default();

        let (pk, sk) = dr.gen();
        let cipher = dr.enc(&pk, &msg);
        let m = dr.dec(&sk, &cipher);
        assert_eq!(msg, m);
    }

    /// Checks whether the full-cycle of gen, enc, dec works properly
    /// for message 0 and larger n.
    #[test]
    fn cycle_zero_large_n() {
        let msg = Z::ZERO;
        let dr = DualRegevWithDiscreteGaussianRegularity::new_from_n(30).unwrap();

        let (pk, sk) = dr.gen();
        let cipher = dr.enc(&pk, &msg);
        let m = dr.dec(&sk, &cipher);
        assert_eq!(msg, m);
    }

    /// Checks whether the full-cycle of gen, enc, dec works properly
    /// for message 1 and larger n.
    #[test]
    fn cycle_one_large_n() {
        let msg = Z::ONE;
        let dr = DualRegevWithDiscreteGaussianRegularity::new_from_n(30).unwrap();

        let (pk, sk) = dr.gen();
        let cipher = dr.enc(&pk, &msg);
        let m = dr.dec(&sk, &cipher);
        assert_eq!(msg, m);
    }
}
