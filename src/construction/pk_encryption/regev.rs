// Copyright © 2023 Niklas Siemer
//
// This file is part of qFALL-crypto.
//
// qFALL-crypto is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

//! This module contains an implementation of the IND-CPA secure
//! public key Regev encryption scheme.
//!
//! The main references are listed in the following:
//! - \[1\] Regev, Oded (2009).
//! On lattices, learning with errors, random linear codes, and cryptography.
//! In: Journal of the ACM 6.
//! <https://dl.acm.org/doi/pdf/10.1145/1568318.1568324>
//! - \[2\] Peikert, Chris (2016).
//! A decade of lattice cryptography.
//! In: Theoretical Computer Science 10.4.
//! <https://web.eecs.umich.edu/~cpeikert/pubs/lattice-survey.pdf>
use super::PKEncryption;
use qfall_math::{
    error::MathError,
    integer::{MatZ, Z},
    integer_mod_q::{MatZq, Modulus, Zq},
    rational::Q,
    traits::{Concatenate, Distance, GetEntry, GetNumRows, Pow, SetEntry},
};
use serde::{Deserialize, Serialize};

/// This struct manages and stores the public parameters of a [`Regev`]
/// public key encryption instance.
///
/// Attributes:
/// - `n`: specifies the security parameter, which is not equal to the bit-security level
/// - `m`: defines the dimension of the underlying lattice
/// - `q`: specifies the modulus over which the encryption is computed
/// - `alpha`: specifies the gaussian parameter used for independent
///   sampling from χ, i.e. for multiple discrete Gaussian samples used
///   for key generation
///
/// # Examples
/// ```
/// use qfall_crypto::construction::pk_encryption::{Regev, PKEncryption};
/// use qfall_math::integer::Z;
/// // setup public parameters and key pair
/// let regev = Regev::default();
/// let (pk, sk) = regev.gen();
///
/// // encrypt a bit
/// let msg = Z::ZERO; // must be a bit, i.e. msg = 0 or 1
/// let cipher = regev.enc(&pk, &msg);
///
/// // decrypt
/// let m = regev.dec(&sk, &cipher);
///
/// assert_eq!(msg, m);
/// ```
#[derive(Debug, Serialize, Deserialize)]
pub struct Regev {
    n: Z,       // security parameter
    m: Z,       // number of rows of matrix A
    q: Modulus, // modulus
    alpha: Q,   // gaussian parameter for sampleZ
}

impl Regev {
    /// Instantiates a [`Regev`] PK encryption instance with the
    /// specified parameters if they ensure a secure and correct instance.
    ///
    /// **WARNING:** The given parameters are not checked for security nor
    /// correctness of the scheme.
    /// If you want to check your parameters for provable security and correctness,
    /// use [`Regev::check_correctness`] and [`Regev::check_security`].
    /// Or use [`Regev::new_from_n`] for generating secure and correct
    /// public parameters for [`Regev`] according to your choice of `n`.
    ///
    /// Parameters:
    /// - `n`: specifies the security parameter and number of rows
    ///   of the uniform at random instantiated matrix `A`
    /// - `m`: specifies the number of columns of matrix `A`
    /// - `q`: specifies the modulus
    /// - `alpha`: specifies the gaussian parameter used for independent
    ///   sampling from chi, i.e. for multiple discrete Gaussian samples used
    ///   for key generation
    ///
    /// Returns a correct and secure [`Regev`] PK encryption instance or
    /// a [`MathError`] if the instance would not be correct or secure.
    ///
    /// # Examples
    /// ```
    /// use qfall_crypto::construction::pk_encryption::Regev;
    ///
    /// let regev = Regev::new(3, 16, 13, 2).unwrap();
    /// ```
    ///
    /// # Panics ...
    /// - if the given modulus `q <= 1`.
    pub fn new(
        n: impl Into<Z>,
        m: impl Into<Z>,
        q: impl Into<Z>,
        alpha: impl Into<Q>,
    ) -> Result<Self, MathError> {
        let n: Z = n.into();
        let m: Z = m.into();
        let q: Z = q.into();
        let alpha: Q = alpha.into();

        let q = Modulus::from(&q);

        Ok(Self { n, m, q, alpha })
    }

    /// Generates a new [`Regev`] instance, i.e. a new set of suitable
    /// (provably secure and correct) public parameters,
    /// given the security parameter `n` for `n >= 10`.
    ///
    /// Parameters:
    /// - `n`: specifies the security parameter and number of rows
    ///   of the uniform at random instantiated matrix `A`
    ///
    /// Returns a correct and secure [`Regev`] PK encryption instance or
    /// a [`MathError`] if the given `n < 10`.
    ///
    /// # Examples
    /// ```
    /// use qfall_crypto::construction::pk_encryption::Regev;
    ///
    /// let regev = Regev::new_from_n(15).unwrap();
    /// ```
    ///
    /// # Errors and Failures
    /// - Returns a [`MathError`] of type [`InvalidIntegerInput`](MathError::InvalidIntegerInput)
    /// if `n <= 1`.
    ///
    /// Panics...
    /// - if `n` does not fit into an [`i64`].
    pub fn new_from_n(n: impl Into<Z>) -> Result<Self, MathError> {
        let n = n.into();
        if n <= Z::from(9) {
            return Err(MathError::InvalidIntegerInput(String::from(
                "Choose n >= 10 as this function does not return parameters ensuring proper correctness of the scheme otherwise.",
            )));
        }

        let mut m: Z;
        let mut q: Modulus;
        let mut alpha: Q;
        (m, q, alpha) = Self::gen_new_public_parameters(&n);
        let mut out = Self {
            n: n.clone(),
            m,
            q,
            alpha,
        };
        while out.check_correctness().is_err() || out.check_security().is_err() {
            (m, q, alpha) = Self::gen_new_public_parameters(&n);
            out = Self {
                n: n.clone(),
                m,
                q,
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
    /// Returns a set of public parameters `(m, q, alpha)` chosen according to
    /// the provided `n`.
    ///
    /// # Examples
    /// ```compile_fail
    /// use qfall_crypto::construction::pk_encryption::Regev;
    /// use qfall_math::integer::Z;
    /// let n = Z::from(2);
    ///
    /// let (m, q, alpha) = Regev::gen_new_public_parameters(&n);
    /// ```
    ///
    /// Panics...
    /// - if `n` does not fit into an [`i64`].
    fn gen_new_public_parameters(n: &Z) -> (Z, Modulus, Q) {
        let n_i64 = i64::try_from(n).unwrap();
        // these powers are chosen according to experience s.t. at least every
        // fifth generation of public parameters outputs a valid pair
        let power = match n_i64 {
            2..=4 => 5,
            5 => 4,
            _ => 3,
        };

        // generate prime q in [n^power / 2, n^power]
        let upper_bound: Z = n.pow(power).unwrap();
        let lower_bound = upper_bound.div_ceil(&Z::from(2));
        let q = Z::sample_prime_uniform(&lower_bound, &upper_bound).unwrap();

        // choose m = (n+1) log q
        let m = (n + Z::ONE) * q.log(2).unwrap().ceil();

        // alpha = 1/(sqrt(n) * log^2 n)
        let alpha = 1 / (n.sqrt() * n.log(2).unwrap().pow(2).unwrap());

        let q = Modulus::from(q);

        (m, q, alpha)
    }

    /// Checks a provided set of public parameters according to their
    /// completeness according to Lemma 5.1 of [\[1\]](<index.html#:~:text=[1]>).
    ///
    /// **WARNING:** Some requirements are missing to ensure overwhelming correctness of the scheme.
    ///
    /// Returns an empty result or a [`MathError`] if the instance would
    /// not be correct.
    ///
    /// # Examples
    /// ```
    /// use qfall_crypto::construction::pk_encryption::Regev;
    /// let regev = Regev::default();
    ///
    /// let is_valid = regev.check_correctness().is_ok();
    /// ```
    ///
    /// # Errors and Failures
    /// - Returns a [`MathError`] of type [`InvalidIntegerInput`](MathError::InvalidIntegerInput)
    /// if at least one parameter was not chosen appropriately for a
    /// correct Regev public key encryption instance.
    pub fn check_correctness(&self) -> Result<(), MathError> {
        let q = Z::from(&self.q);

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

        // Correctness requirements
        // α = o (1 / ( sqrt(n) * log n ) )
        if self.alpha > 1 / (self.n.sqrt() * self.n.log(2).unwrap()) {
            return Err(MathError::InvalidIntegerInput(String::from(
                "Completeness is not guaranteed as α >= 1 / (sqrt(n) * log n), but α < 1 / (sqrt(n) * log n) is required."
            )));
        }

        if 20 * self.m.sqrt() * &self.alpha > Q::from(q) {
            return Err(MathError::InvalidIntegerInput(String::from(
                "Completeness is not guaranteed as 5 * sqrt(m) * α > q/4, but 5 * sqrt(m) * α <= q/4 is required."
            )));
        }

        Ok(())
    }

    /// Checks a provided set of public parameters according to their validity
    /// regarding security according to Theorem 1.1
    /// and Lemma 5.4 of [\[1\]](<index.html#:~:text=[1]>).
    ///
    /// Returns an empty result or a [`MathError`] if the instance would
    /// not be secure.
    ///
    /// # Examples
    /// ```
    /// use qfall_crypto::construction::pk_encryption::Regev;
    /// let regev = Regev::default();
    ///
    /// let is_valid = regev.check_security().is_ok();
    /// ```
    ///
    /// # Errors and Failures
    /// - Returns a [`MathError`] of type [`InvalidIntegerInput`](MathError::InvalidIntegerInput)
    /// if at least one parameter was not chosen appropriately for a
    /// secure Regev public key encryption instance.
    pub fn check_security(&self) -> Result<(), MathError> {
        let q = Z::from(&self.q);

        // Security requirements
        // q * α >= 2 sqrt(n)
        if &q * &self.alpha < 2 * self.n.sqrt() {
            return Err(MathError::InvalidIntegerInput(String::from(
                "Security is not guaranteed as q * α < 2 sqrt(n), but q * α >= 2 sqrt(n) is required.",
            )));
        }
        // m > (n + 1) log q
        if self.m <= ((&self.n + Z::ONE) * q.log(2).unwrap()).ceil() {
            return Err(MathError::InvalidIntegerInput(String::from(
                "Security is not guaranteed as m <= (n + 1) log q,
                but m > (n + 1) log q is required.",
            )));
        }

        Ok(())
    }

    /// This function instantiates a 128-bit secure [`Regev`] scheme.
    ///
    /// The public parameters used for this scheme were generated via `Regev::new_from_n(350)`
    /// and its bit-security determined via the [lattice estimator](https://github.com/malb/lattice-estimator).
    pub fn secure128() -> Self {
        Self::new(230, 5313, 7764299, 0.0011).unwrap()
    }
}

impl Default for Regev {
    /// Initializes a [`Regev`] struct with parameters generated by `Regev::new_from_n(3)`.
    /// This parameter choice is not secure as the dimension of the lattice is too small,
    /// but it provides an efficient working example.
    ///
    /// # Examples
    /// ```
    /// use qfall_crypto::construction::pk_encryption::Regev;
    ///
    /// let regev = Regev::default();
    /// ```
    fn default() -> Self {
        let n = Z::from(13);
        let m = Z::from(154);
        let q = Modulus::from(1427);
        let alpha = Q::from(0.02);

        Self { n, m, q, alpha }
    }
}

impl PKEncryption for Regev {
    type Cipher = MatZq;
    type PublicKey = MatZq;
    type SecretKey = MatZq;

    /// Generates a (pk, sk) pair for the Regev public key encryption scheme
    /// by following these steps:
    /// - A <- Z_q^{n x m}
    /// - s <- Z_q^n
    /// - e^t <- χ^m
    /// - b^t = s^t * A + e^t
    /// - A = [A^t | b]^t
    ///
    /// Then, `pk = A` and `sk = s` is output.
    ///
    /// # Examples
    /// ```
    /// use qfall_crypto::construction::pk_encryption::{PKEncryption, Regev};
    /// let regev = Regev::default();
    ///
    /// let (pk, sk) = regev.gen();
    /// ```
    fn gen(&self) -> (Self::PublicKey, Self::SecretKey) {
        // A <- Z_q^{n x m}
        let mat_a = MatZq::sample_uniform(&self.n, &self.m, &self.q);
        // s <- Z_q^n
        let vec_s = MatZq::sample_uniform(&self.n, 1, &self.q);
        // e^t <- χ^m
        let vec_e_t = MatZq::sample_discrete_gauss(
            1,
            &self.m,
            &self.q,
            &self.n,
            0,
            &self.alpha * Z::from(&self.q),
        )
        .unwrap();

        // b^t = s^t * A + e^t
        let vec_b_t = vec_s.transpose() * &mat_a + vec_e_t;

        // A = [A^t | b]^t
        let mat_a = mat_a.concat_vertical(&vec_b_t).unwrap();

        // pk = A, sk = s
        (mat_a, vec_s)
    }

    /// Generates an encryption of `message mod 2` for the provided public key by following these steps:
    /// - x <- Z_2^m
    /// - c = A * x + [0^{1 x n} | msg *  ⌊q/2⌋]^t
    ///
    /// Then, cipher `c` is output.
    ///
    /// Parameters:
    /// - `pk`: specifies the public key `pk = A`
    /// - `message`: specifies the message that should be encryted
    ///
    /// Returns a cipher `c` of type [`MatZq`].
    ///
    /// # Examples
    /// ```
    /// use qfall_crypto::construction::pk_encryption::{PKEncryption, Regev};
    /// let regev = Regev::default();
    /// let (pk, sk) = regev.gen();
    ///
    /// let cipher = regev.enc(&pk, 1);
    /// ```
    fn enc(&self, pk: &Self::PublicKey, message: impl Into<Z>) -> Self::Cipher {
        // generate message = message mod 2
        let message: Z = message.into();
        let message = Zq::from((message, 2));
        let message = message.get_value();

        // x <- Z_2^m
        let vec_x = MatZ::sample_uniform(&self.m, 1, 0, 2).unwrap();

        // c = A * x + [0^{1xn} | msg *  ⌊q/2⌋]^t
        let mut c = pk * vec_x;

        // hide message in last entry
        // compute msg * ⌊q/2⌋
        let msg_q_half = message * Z::from(&self.q).div_floor(&Z::from(2));
        // set last entry of c = last_entry + msg * ⌊q/2⌋
        let last_entry: Zq = c.get_entry(c.get_num_rows() - 1, 0).unwrap();
        c.set_entry(c.get_num_rows() - 1, 0, last_entry + msg_q_half)
            .unwrap();

        c
    }

    /// Decrypts the provided `cipher` using the secret key `sk` by following these steps:
    /// - x = [-sk^t | 1] * c
    /// - if x mod q is closer to ⌊q/2⌋ than to 0, output 1. Otherwise, output 0.
    ///
    /// Parameters:
    /// - `sk`: specifies the secret key `sk = s`
    /// - `cipher`: specifies the cipher containing `cipher = c`
    ///
    /// Returns the decryption of `cipher` as a [`Z`] instance.
    ///
    /// # Examples
    /// ```
    /// use qfall_crypto::construction::pk_encryption::{PKEncryption, Regev};
    /// use qfall_math::integer::Z;
    /// let regev = Regev::default();
    /// let (pk, sk) = regev.gen();
    /// let cipher = regev.enc(&pk, 1);
    ///
    /// let m = regev.dec(&sk, &cipher);
    ///
    /// assert_eq!(Z::ONE, m);
    /// ```
    fn dec(&self, sk: &Self::SecretKey, cipher: &Self::Cipher) -> Z {
        let result = (-1i8 * sk)
            .concat_vertical(&MatZq::identity(1, 1, &self.q))
            .unwrap()
            .dot_product(cipher)
            .unwrap();

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
    use super::Regev;
    use super::Z;

    /// Checks whether `new` is available for types implementing [`Into<Z>`].
    #[test]
    fn new_availability() {
        let _ = Regev::new(2u8, 2u16, 2u32, 2u64);
        let _ = Regev::new(2u16, 2u64, 2i32, 2i64);
        let _ = Regev::new(2i16, 2i64, 2u32, 2u8);
        let _ = Regev::new(Z::from(2), &Z::from(2), 2u8, 2i8);
    }

    /// Checks whether `new_from_n` works properly for different choices of n.
    #[test]
    fn suitable_security_params() {
        let n_choices = [
            10, 11, 12, 13, 14, 25, 50, 100, 250, 500, 1000, 2500, 5000, 5001, 10000,
        ];

        for n in n_choices {
            assert!(Regev::new_from_n(n).is_ok());
        }
    }

    /// Checks whether the [`Default`] parameter choice is suitable.
    #[test]
    fn default_suitable() {
        let regev = Regev::default();

        assert!(regev.check_correctness().is_ok());
        assert!(regev.check_security().is_ok());
    }

    /// Checks whether the generated public parameters from `new_from_n` are
    /// valid choices according to security and correctness of the scheme.
    #[test]
    fn choice_valid() {
        let n_choices = [10, 14, 25, 50, 125, 300, 600, 1200, 4000, 6000];

        for n in n_choices {
            let regev = Regev::new_from_n(n).unwrap();
            assert!(regev.check_correctness().is_ok());
            assert!(regev.check_security().is_ok());
        }
    }

    /// Ensures that `new_from_n` is available for types implementing [`Into<Z>`].
    #[test]
    fn availability() {
        let _ = Regev::new_from_n(10u8);
        let _ = Regev::new_from_n(10u16);
        let _ = Regev::new_from_n(10u32);
        let _ = Regev::new_from_n(10u64);
        let _ = Regev::new_from_n(10i8);
        let _ = Regev::new_from_n(10i16);
        let _ = Regev::new_from_n(10i32);
        let _ = Regev::new_from_n(10i64);
        let _ = Regev::new_from_n(Z::from(10));
        let _ = Regev::new_from_n(&Z::from(10));
    }

    /// Checks whether `new_from_n` returns an error for invalid input n.
    #[test]
    fn invalid_n() {
        assert!(Regev::new_from_n(1).is_err());
        assert!(Regev::new_from_n(0).is_err());
        assert!(Regev::new_from_n(-1).is_err());
    }

    /// Checks whether `secure128` outputs a new instance with correct and secure parameters.
    #[test]
    fn secure128_validity() {
        let regev = Regev::secure128();

        assert!(regev.check_correctness().is_ok());
        assert!(regev.check_security().is_ok());
    }
}

#[cfg(test)]
mod test_regev {
    use super::Regev;
    use crate::construction::pk_encryption::PKEncryption;
    use qfall_math::integer::Z;

    /// Checks whether the full-cycle of gen, enc, dec works properly
    /// for message 0 and small n.
    #[test]
    fn cycle_zero_small_n() {
        let msg = Z::ZERO;
        let regev = Regev::default();

        let (pk, sk) = regev.gen();
        let cipher = regev.enc(&pk, &msg);
        let m = regev.dec(&sk, &cipher);
        assert_eq!(msg, m);
    }

    /// Checks whether the full-cycle of gen, enc, dec works properly
    /// for message 1 and small n.
    #[test]
    fn cycle_one_small_n() {
        let msg = Z::ONE;
        let regev = Regev::default();

        let (pk, sk) = regev.gen();
        let cipher = regev.enc(&pk, &msg);
        let m = regev.dec(&sk, &cipher);
        assert_eq!(msg, m);
    }

    /// Checks whether the full-cycle of gen, enc, dec works properly
    /// for message 0 and larger n.
    #[test]
    fn cycle_zero_large_n() {
        let msg = Z::ZERO;
        let regev = Regev::new_from_n(50).unwrap();

        let (pk, sk) = regev.gen();
        let cipher = regev.enc(&pk, &msg);
        let m = regev.dec(&sk, &cipher);
        assert_eq!(msg, m);
    }

    /// Checks whether the full-cycle of gen, enc, dec works properly
    /// for message 1 and larger n.
    #[test]
    fn cycle_one_large_n() {
        let msg = Z::ONE;
        let regev = Regev::new_from_n(50).unwrap();

        let (pk, sk) = regev.gen();
        let cipher = regev.enc(&pk, &msg);
        let m = regev.dec(&sk, &cipher);
        assert_eq!(msg, m);
    }
}
