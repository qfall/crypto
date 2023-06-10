// Copyright © 2023 Niklas Siemer
//
// This file is part of qFALL-crypto.
//
// qFALL-crypto is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

//! This module contains an implementation of the IND-CCA secure
//! public key Dual Regev encryption scheme.
//!
//! The main references are listed in the following:
//! - \[1\] Gentry, Craig and Peikert, Chris and Vaikuntanathan, Vinod (2008).
//! Trapdoors for hard lattices and new cryptographic constructions.
//! In: Proceedings of the fortieth annual ACM symposium on Theory of computing.
//! <https://dl.acm.org/doi/pdf/10.1145/1374376.1374407>

use super::PKEncryption;
use qfall_math::{
    error::MathError,
    integer::{MatZ, Z},
    integer_mod_q::{MatZq, Modulus, Zq},
    rational::{MatQ, Q},
    traits::{Distance, Pow},
};
use serde::{Deserialize, Serialize};

#[derive(Debug, Serialize, Deserialize)]
pub struct DualRegev {
    n: Z,       // security parameter
    m: Z,       // number of rows of matrix A
    q: Modulus, // modulus
    r: Q,       // gaussian parameter for sampleD
    alpha: Q,   // gaussian parameter for sampleZ
}

impl DualRegev {
    /// Instantiates a [`DualRegev`] PK encryption instance with the
    /// specified parameters if they ensure a secure and correct instance.
    ///
    /// Parameters:
    /// - `n`: specifies the security parameter and number of rows
    ///   of the uniform at random instantiated matrix `A`
    /// - `m`: specifies the number of columns of matrix `A`
    /// - `q`: specifies the modulus
    /// - `r`: specifies the gaussian parameter used for SampleD,
    ///   i.e. used for encryption
    /// - `alpha`: specifies the gaussian parameter used for independent
    ///   sampling from chi, i.e. for multiple discrete Gaussian samples used
    ///   for key generation
    ///
    /// Returns a correct and secure [`DualRegev`] PK encryption instance or
    /// a [`MathError`] if the instance would not be correct or secure.
    ///
    /// # Example
    /// ```
    /// use qfall_crypto::construction::pk_encryption::DualRegev;
    ///
    /// let dual_regev = DualRegev::new(2, 16, 401, 4, 8).unwrap();
    /// ```
    ///
    /// # Errors and Failures
    /// - Returns a [`MathError`] of type [`InvalidIntegerInput`](MathError::InvalidIntegerInput)
    /// if at least one parameter was not chosen appropriately for a
    /// correct and secure DualRegev public key encryption instance.
    pub fn new(
        n: impl Into<Z>,
        m: impl Into<Z>,
        q: impl Into<Z>,
        r: impl Into<Q>,
        alpha: impl Into<Q>,
    ) -> Result<Self, MathError> {
        let n: Z = n.into();
        let m: Z = m.into();
        let q: Z = q.into();
        let r: Q = r.into();
        let alpha: Q = alpha.into();

        Self::check_params(&n, &m, &q, &r, &alpha)?;

        let q = Modulus::try_from(&q).unwrap();

        Ok(Self { n, m, q, r, alpha })
    }

    /// Generates a new [`DualRegev`] instance, i.e. a new set of suitable public parameters,
    /// given the security parameter `n`.
    ///
    /// Parameters:
    /// - `n`: specifies the security parameter and number of rows
    ///   of the uniform at random instantiated matrix `A`
    ///
    /// Returns a correct and secure [`DualRegev`] PK encryption instance or
    /// a [`MathError`] if the given `n <= 1`.
    ///
    /// # Example
    /// ```
    /// use qfall_crypto::construction::pk_encryption::DualRegev;
    ///
    /// let dual_regev = DualRegev::new_from_n(2).unwrap();
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
        let mut q: Z;
        let mut r: Q;
        let mut alpha: Q;
        (m, q, r, alpha) = Self::gen_new_public_parameters(&n);
        while Self::check_params(&n, &m, &q, &r, &alpha).is_err() {
            (m, q, r, alpha) = Self::gen_new_public_parameters(&n);
        }

        let q = Modulus::try_from_z(&q).unwrap();

        Ok(Self { n, m, q, r, alpha })
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
    /// # Example
    /// ```compile_fail
    /// use qfall_crypto::construction::pk_encryption::DualRegev;
    /// use qfall_math::integer::Z;
    /// let n = Z::from(2);
    ///
    /// let (m, q, r, alpha) = DualRegev::gen_new_public_parameters(&n);
    /// ```
    fn gen_new_public_parameters(n: &Z) -> (Z, Z, Q, Q) {
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
            31..=5000 => 3,
            _ => 2,
        };

        // generate prime q in [n^power / 2, n^power]
        // TODO: Replace by q = Z::sample_prime_uniform once implemented
        let upper_bound: Z = n.pow(power).unwrap();
        let lower_bound = upper_bound.div_ceil(&Z::from(2));
        let mut q = Z::sample_uniform(&lower_bound, &upper_bound).unwrap();
        while !q.is_prime() {
            q = Z::sample_uniform(&lower_bound, &upper_bound).unwrap();
        }

        // choose m = 2 (n+1) lg q
        let m = (Z::from(2) * (n + Z::ONE) * q.log(&10).unwrap()).round();

        // choose r = log m
        let r = m.log(&2).unwrap();

        // alpha = 1/(sqrt(m) * log^2 m)
        // TODO: remove ceil, when log is applicable to Qs
        let alpha = 1 / m.sqrt() * (m.log(&2).unwrap()).ceil().log(&2).unwrap();

        (m, q, r, alpha)
    }

    /// Checks a provided set of public parameters according to their validity
    /// regarding security and completeness according to
    /// Lemma 8.2 and 8.4 of [\[1\]](<index.html#:~:text=[1]>).
    ///
    /// Parameters:
    /// - `n`: specifies the security parameter and number of rows
    ///   of the uniform at random instantiated matrix `A`
    /// - `m`: specifies the number of columns of matrix `A`
    /// - `q`: specifies the modulus
    /// - `r`: specifies the gaussian parameter used for SampleD,
    ///   i.e. used for encryption
    /// - `alpha`: specifies the gaussian parameter used for independent
    ///   sampling from chi, i.e. for [`MatZ::sample_discrete_gauss`] used
    ///   for key generation
    ///
    /// Returns an empty result or a [`MathError`] if the instance would
    /// not be correct or secure.
    ///
    /// # Example
    /// ```compile_fail
    /// use qfall_crypto::construction::pk_encryption::DualRegev;
    /// use qfall_math::{integer::Z, rational::Q};
    /// let n = Z::from(2);
    /// let m = Z::from(16);
    /// let q = Z::from(401);
    /// let r = Q::from(4);
    /// let alpha = Q::from(8);
    ///
    /// let is_valid = DualRegev::check_params(&n, &m, &q, &r, &alpha).is_err();
    /// ```
    ///
    /// # Errors and Failures
    /// - Returns a [`MathError`] of type [`InvalidIntegerInput`](MathError::InvalidIntegerInput)
    /// if at least one parameter was not chosen appropriately for a
    /// correct and secure DualRegev public key encryption instance.
    fn check_params(n: &Z, m: &Z, q: &Z, r: &Q, alpha: &Q) -> Result<(), MathError> {
        if n <= &Z::ONE {
            return Err(MathError::InvalidIntegerInput(String::from(
                "n must be chosen bigger than 1.",
            )));
        }
        if !q.is_prime() {
            return Err(MathError::InvalidIntegerInput(String::from(
                "q must be prime.",
            )));
        }

        // Security requirements
        // q * α >= n
        if q * alpha < Q::from(n) {
            return Err(MathError::InvalidIntegerInput(String::from(
                "Security is not guaranteed as q * α < n, but q * α >= n is required.",
            )));
        }
        // m >= 2(n + 1) lg (q)
        if Q::from(m) < 2 * (n + 1) * q.log(&10).unwrap() {
            return Err(MathError::InvalidIntegerInput(String::from(
                "Security is not guaranteed as m < 2(n + 1) lg (q), but m >= 2(n + 1) lg (q) is required."
            )));
        }
        // r >= ω( sqrt( log m ) )
        if r < &m.log(&2).unwrap().sqrt() {
            return Err(MathError::InvalidIntegerInput(String::from(
                "Security is not guaranteed as r < sqrt( log m ) and r >= ω(sqrt(log m)) is required."
            )));
        }

        // Completeness requirements
        // q >= 5 * r * m
        if Q::from(q) < 5 * r * m {
            return Err(MathError::InvalidIntegerInput(String::from(
                "Completeness is not guaranteed as q < 5rm, but q >= 5rm is required.",
            )));
        }
        // α <= 1/(r * sqrt(m) * ω(sqrt(log n))
        if alpha > &(1 / (r * m.sqrt() * n.log(&2).unwrap().sqrt())) {
            return Err(MathError::InvalidIntegerInput(String::from(
                "Completeness is not guaranteed as α > 1/(r*sqrt(m)*ω(sqrt(log n)), but α <= 1/(r*sqrt(m)*ω(sqrt(log n)) is required."
            )));
        }

        Ok(())
    }
}

impl Default for DualRegev {
    /// Initializes a [`DualRegev`] struct with parameters generated by `DualRegev::new_from_n(2)`.
    /// This parameter choice is not secure as the dimension of the lattice is too small,
    /// but it provides an efficient working example.
    ///
    /// # Example
    /// ```
    /// use qfall_crypto::construction::pk_encryption::DualRegev;
    ///
    /// let dual_regev = DualRegev::default();
    /// ```
    fn default() -> Self {
        let n = Z::from(2);
        let m = Z::from(16);
        let q = Modulus::try_from(&Z::from(401)).unwrap();
        let r = Q::from(4);
        let alpha = Q::from(8);

        Self { n, m, q, r, alpha }
    }
}

impl PKEncryption for DualRegev {
    type Cipher = (MatZq, Zq);
    type PublicKey = (MatZq, MatZq);
    type SecretKey = MatZq;

    /// Generates a (pk, sk) pair for the Dual Regev public key encryption scheme
    /// by following these steps:
    /// - s <- Z_q^n
    /// - A <- Z_q^{n x m}
    /// - x <- χ^m
    /// - p = A^t * s + x
    ///
    /// Then, `pk = (A, p)` and `sk = s` is output.
    ///
    /// # Example
    /// ```
    /// use qfall_crypto::construction::pk_encryption::{PKEncryption, DualRegev};
    /// let dual_regev = DualRegev::default();
    ///
    /// let (pk, sk) = dual_regev.gen();
    /// ```
    fn gen(&self) -> (Self::PublicKey, Self::SecretKey) {
        // s <- Z_q^n
        let vec_s = MatZq::sample_uniform(&self.n, 1, &self.q).unwrap();

        // A <- Z_q^{n x m}
        let mat_a = MatZq::sample_uniform(&self.n, &self.m, &self.q).unwrap();
        // x <- χ^m
        let vec_x =
            MatZq::sample_discrete_gauss(&self.m, 1, &self.q, &self.n, 0, &self.alpha).unwrap();
        // p = A^t * s + x
        let vec_p = mat_a.transpose() * &vec_s + vec_x;

        // pk = (A, p), sk = s
        ((mat_a, vec_p), vec_s)
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
    /// # Example
    /// ```
    /// use qfall_crypto::construction::pk_encryption::{PKEncryption, DualRegev};
    /// let dual_regev = DualRegev::default();
    /// let (pk, sk) = dual_regev.gen();
    ///
    /// let cipher = dual_regev.enc(&pk, 1);
    /// ```
    fn enc(&self, pk: &Self::PublicKey, message: impl Into<Z>) -> Self::Cipher {
        // generate message = message mod 2
        let message: Z = message.into();
        let message = Zq::try_from_z_z(&message, &Z::from(2)).unwrap();
        let message = message.get_value();

        // e <- SampleD over lattice Z^m, center 0 with gaussian parameter r
        let basis = MatZ::identity(&self.m, &self.m).unwrap();
        let center = MatQ::new(&self.m, 1).unwrap();
        // TODO: Replace by MatZq::sample_d once available
        let vec_e = MatZ::sample_d(&basis, &self.n, &center, &self.r).unwrap();
        let vec_e = MatZq::from((&vec_e, &self.q));

        // u = A * e
        let vec_u = &pk.0 * &vec_e;
        // c = p^t * e + msg *  ⌊q/2⌋
        let q_half = Z::from(&self.q).div_floor(&Z::from(2));
        let c = pk.1.dot_product(&vec_e).unwrap() + message * q_half;

        (vec_u, c)
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
    /// # Example
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
    use super::DualRegev;
    use super::Z;

    /// Checks whether `new` works properly for correct and secure parameter choices
    #[test]
    fn new_suitable() {
        assert!(DualRegev::new(2, 16, 401, 4, 8).is_ok());
        assert!(DualRegev::new(6, 62, 24781, 6, 21).is_ok());
        assert!(DualRegev::new(20, 210, 99823, 8, 43).is_ok());
    }

    /// Checks whether `new` returns an error for insecure or not complete public parameters
    #[test]
    fn new_unsuitable() {
        assert!(DualRegev::new(2, 16, 401, 4, 20).is_err());
        assert!(DualRegev::new(2, 16, 401, 6, 8).is_err());
        assert!(DualRegev::new(2, 16, 500, 4, 8).is_err());
        assert!(DualRegev::new(2, 16, 399, 4, 8).is_err());
        assert!(DualRegev::new(2, 30, 401, 4, 8).is_err());
        assert!(DualRegev::new(1, 16, 401, 4, 20).is_err());
    }

    /// Checks whether `new_from_n` works properly for different choices of n
    #[test]
    fn suitable_security_params() {
        let n_choices = [
            2, 3, 4, 5, 6, 8, 9, 12, 13, 30, 21, 50, 100, 250, 500, 1000, 2500, 5000, 5001, 10000,
        ];

        for n in n_choices {
            assert!(DualRegev::new_from_n(n).is_ok());
        }
    }

    /// Checks whether the [`Default`] parameter choice is suitable
    #[test]
    fn default_suitable() {
        let dr = DualRegev::default();

        assert!(DualRegev::check_params(&dr.n, &dr.m, &Z::from(&dr.q), &dr.r, &dr.alpha).is_ok());
    }

    /// Checks whether the generated public parameters from `new_from_n` are
    /// valid choices according to security and correctness of the scheme
    #[test]
    fn choice_valid() {
        let n_choices = [2, 3, 5, 8, 10, 14, 25, 50, 125, 300, 600, 1200, 4000, 6000];

        for n in n_choices {
            let dr = DualRegev::new_from_n(n).unwrap();
            assert!(
                DualRegev::check_params(&dr.n, &dr.m, &Z::from(&dr.q), &dr.r, &dr.alpha).is_ok()
            );
        }
    }

    /// Ensures that `new_from_n` is available for types implementing [`Into<Z>`]
    #[test]
    fn availability() {
        let _ = DualRegev::new_from_n(2u8);
        let _ = DualRegev::new_from_n(2u16);
        let _ = DualRegev::new_from_n(2u32);
        let _ = DualRegev::new_from_n(2u64);
        let _ = DualRegev::new_from_n(2i8);
        let _ = DualRegev::new_from_n(2i16);
        let _ = DualRegev::new_from_n(2i32);
        let _ = DualRegev::new_from_n(2i64);
        let _ = DualRegev::new_from_n(Z::from(2));
        let _ = DualRegev::new_from_n(&Z::from(2));
    }

    /// Checks whether `new_from_n` returns an error for invalid input n
    #[test]
    fn invalid_n() {
        assert!(DualRegev::new_from_n(1).is_err());
        assert!(DualRegev::new_from_n(0).is_err());
        assert!(DualRegev::new_from_n(-1).is_err());
    }
}

#[cfg(test)]
mod test_dual_regev {
    use super::DualRegev;
    use crate::construction::pk_encryption::PKEncryption;
    use qfall_math::integer::Z;

    /// Checks whether the full-cycle of gen, enc, dec works properly
    /// for message 0 and small n
    #[test]
    fn cycle_zero_small_n() {
        let msg = Z::ZERO;
        let dr = DualRegev::new_from_n(5).unwrap();

        let (pk, sk) = dr.gen();
        let cipher = dr.enc(&pk, &msg);
        let m = dr.dec(&sk, &cipher);
        assert_eq!(msg, m);
    }

    /// Checks whether the full-cycle of gen, enc, dec works properly
    /// for message 1 and small n
    #[test]
    fn cycle_one_small_n() {
        let msg = Z::ONE;
        let dr = DualRegev::new_from_n(5).unwrap();

        let (pk, sk) = dr.gen();
        let cipher = dr.enc(&pk, &msg);
        let m = dr.dec(&sk, &cipher);
        assert_eq!(msg, m);
    }

    /// Checks whether the full-cycle of gen, enc, dec works properly
    /// for message 0 and larger n
    #[test]
    fn cycle_zero_large_n() {
        let msg = Z::ZERO;
        let dr = DualRegev::new_from_n(50).unwrap();

        let (pk, sk) = dr.gen();
        let cipher = dr.enc(&pk, &msg);
        let m = dr.dec(&sk, &cipher);
        assert_eq!(msg, m);
    }

    /// Checks whether the full-cycle of gen, enc, dec works properly
    /// for message 1 and larger n
    #[test]
    fn cycle_one_large_n() {
        let msg = Z::ONE;
        let dr = DualRegev::new_from_n(50).unwrap();

        let (pk, sk) = dr.gen();
        let cipher = dr.enc(&pk, &msg);
        let m = dr.dec(&sk, &cipher);
        assert_eq!(msg, m);
    }
}
