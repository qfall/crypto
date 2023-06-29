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
    /// # Example
    /// ```
    /// use qfall_crypto::construction::pk_encryption::Regev;
    ///
    /// let regev = Regev::new(3, 16, 13, 2).unwrap();
    /// ```
    ///
    /// # Errors and Failures
    /// - Returns a [`MathError`] of type [`InvalidIntegerInput`](MathError::InvalidIntegerInput)
    /// if at least one parameter was not chosen appropriately for a
    /// correct and secure Regev public key encryption instance.
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

        Self::check_params(&n, &m, &q, &alpha)?;

        let q = Modulus::from(&q);

        Ok(Self { n, m, q, alpha })
    }

    /// Generates a new [`Regev`] instance, i.e. a new set of suitable public parameters,
    /// given the security parameter `n`.
    ///
    /// Parameters:
    /// - `n`: specifies the security parameter and number of rows
    ///   of the uniform at random instantiated matrix `A`
    ///
    /// Returns a correct and secure [`Regev`] PK encryption instance or
    /// a [`MathError`] if the given `n <= 1`.
    ///
    /// # Example
    /// ```
    /// use qfall_crypto::construction::pk_encryption::Regev;
    ///
    /// let regev = Regev::new_from_n(4).unwrap();
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
        let mut alpha: Q;
        (m, q, alpha) = Self::gen_new_public_parameters(&n);
        while Self::check_params(&n, &m, &q, &alpha).is_err() {
            (m, q, alpha) = Self::gen_new_public_parameters(&n);
        }

        let q = Modulus::from(&q);

        Ok(Self { n, m, q, alpha })
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
    /// use qfall_crypto::construction::pk_encryption::Regev;
    /// use qfall_math::integer::Z;
    /// let n = Z::from(2);
    ///
    /// let (m, q, alpha) = Regev::gen_new_public_parameters(&n);
    /// ```
    fn gen_new_public_parameters(n: &Z) -> (Z, Z, Q) {
        // generate prime q in [n^2, 2 n^2]
        let lower_bound: Z = n.pow(2).unwrap();
        let upper_bound = 2 * &lower_bound;
        let q = Z::sample_prime_uniform(&lower_bound, &upper_bound).unwrap();

        // choose m = 1.05 * (n+1) log q
        let m = ((n + Z::ONE) * q.log(&2).unwrap() * 1.05f32).round();

        // alpha = 1/(sqrt(n) * log^2 n)
        // TODO: remove ceil, when log is applicable to Qs
        let alpha = 1 / n.sqrt() * (n.log(&2).unwrap()).ceil().log(&2).unwrap();

        (m, q, alpha)
    }

    /// Checks a provided set of public parameters according to their validity
    /// regarding security and completeness according to Theorem 1.1,
    /// Lemma 5.1 and 5.4 of [\[1\]](<index.html#:~:text=[1]>).
    ///
    /// Parameters:
    /// - `n`: specifies the security parameter and number of rows
    ///   of the uniform at random instantiated matrix `A`
    /// - `m`: specifies the number of columns of matrix `A`
    /// - `q`: specifies the modulus
    /// - `alpha`: specifies the gaussian parameter used for independent
    ///   sampling from chi, i.e. for [`MatZ::sample_discrete_gauss`] used
    ///   for key generation
    ///
    /// Returns an empty result or a [`MathError`] if the instance would
    /// not be correct or secure.
    ///
    /// # Example
    /// ```compile_fail
    /// use qfall_crypto::construction::pk_encryption::Regev;
    /// use qfall_math::{integer::Z, rational::Q};
    /// let n = Z::from(2);
    /// let m = Z::from(16);
    /// let q = Z::from(401);
    /// let alpha = Q::from(8);
    ///
    /// let is_valid = Regev::check_params(&n, &m, &q, &alpha).is_err();
    /// ```
    ///
    /// # Errors and Failures
    /// - Returns a [`MathError`] of type [`InvalidIntegerInput`](MathError::InvalidIntegerInput)
    /// if at least one parameter was not chosen appropriately for a
    /// correct and secure Regev public key encryption instance.
    fn check_params(n: &Z, m: &Z, q: &Z, alpha: &Q) -> Result<(), MathError> {
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
        // q * α >= 2 sqrt(n)
        if q * alpha < 2 * n.sqrt() {
            return Err(MathError::InvalidIntegerInput(String::from(
                "Security is not guaranteed as q * α < n, but q * α >= n is required.",
            )));
        }
        // m >= (n + 1) log (q)
        if Q::from(m) <= (n + 1) * q.log(&2).unwrap() {
            return Err(MathError::InvalidIntegerInput(String::from(
                "Security is not guaranteed as m < 2(n + 1) lg (q),
                but m >= 2(n + 1) lg (q) is required.",
            )));
        }

        // Completeness requirements
        // q >= n^2
        if q < &n.pow(2).unwrap() {
            return Err(MathError::InvalidIntegerInput(String::from(
                "Completeness is not guaranteed as q < n^2, but q >= n^2 is required.",
            )));
        }
        // q <= 2 n^2
        if q > &(2 * n.pow(2).unwrap()) {
            return Err(MathError::InvalidIntegerInput(String::from(
                "Completeness is not guaranteed as q > 2 * n^2, but q <= 2 * n^2 is required.",
            )));
        }
        // α = o (1 / ( sqrt(n) * log n ) )
        if alpha >= &(1 / (n.sqrt() * n.log(&2).unwrap())) {
            return Err(MathError::InvalidIntegerInput(String::from(
                "Completeness is not guaranteed as α >= 1 / (sqrt(n) * log n), but α < 1 / (sqrt(n) * log n) is required."
            )));
        }

        Ok(())
    }

    /// This function instantiates a 128-bit secure [`Regev`] scheme.
    ///
    /// **WARNING:** Arora-GB seems to be very effective due to the uniform error-distribution.
    /// Here, the bit-security could be significantly lower - down to 2^36 operations.
    /// Although we believe that this value is so low due to wrongly chosen parameters on our side.
    /// TODO: Further investigate whether parameter choice for lattice estimator suited the interface.
    /// Current input to lattice estimator:
    /// ```python
    /// x = LWEParameters(n=617, q=2927, Xe=ND.Uniform(0, 1), Xs=ND.DiscreteGaussian(7.29))
    /// LWE.estimate(x)
    /// ```
    ///
    /// The public parameters used for this scheme were generated via `Regev::new_from_n(50)`
    /// and its bit-security determined via the [lattice estimator](https://github.com/malb/lattice-estimator).
    pub fn secure128() -> Self {
        Self::new(50, 617, 2927, 18.1).unwrap()
    }
}

impl Default for Regev {
    /// Initializes a [`Regev`] struct with parameters generated by `Regev::new_from_n(2)`.
    /// This parameter choice is not secure as the dimension of the lattice is too small,
    /// but it provides an efficient working example.
    ///
    /// # Example
    /// ```
    /// use qfall_crypto::construction::pk_encryption::Regev;
    ///
    /// let regev = Regev::default();
    /// ```
    fn default() -> Self {
        let n = Z::from(3);
        let m = Z::from(16);
        let q = Modulus::from(13);
        let alpha = Q::from(2);

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
    /// # Example
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
        let vec_e_t =
            MatZq::sample_discrete_gauss(1, &self.m, &self.q, &self.n, 0, &self.alpha).unwrap();

        // b^t = s^t * A + e^t
        let vec_b_t = vec_s.transpose() * &mat_a + vec_e_t;

        // A = [A^t | b]^t
        let mat_a = mat_a.concat_vertical(&vec_b_t).unwrap();

        // pk = A, sk = s
        (mat_a, vec_s)
    }

    /// Generates an encryption of `message mod 2` for the provided public key by following these steps:
    /// e <- SampleD over lattice Z^m, center 0 with gaussian parameter r
    /// - x <- Z_2^m
    /// - c = A * x + [0^{1xn} | msg *  ⌊q/2⌋]^t
    ///
    /// Then, cipher `c` is output.
    ///
    /// Parameters:
    /// - `pk`: specifies the public key, which contains two matrices `pk = (A, p)`
    /// - `message`: specifies the message that should be encryted
    ///
    /// Returns a cipher `c` of type [`MatZq`].
    ///
    /// # Example
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
        let vec_x = MatZ::sample_uniform(&self.m, 1, &0, &2).unwrap();

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
    /// # Example
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

    /// Checks whether `new` works properly for correct and secure parameter choices.
    #[test]
    fn new_suitable() {
        assert!(Regev::new(3, 17, 17, 2).is_ok());
        assert!(Regev::new(6, 42, 53, 4).is_ok());
        assert!(Regev::new(40, 499, 3109, 16).is_ok());
    }

    /// Checks whether `new` returns an error for insecure or not complete public parameters.
    #[test]
    fn new_unsuitable() {
        assert!(Regev::new(3, 16, 13, 20).is_err());
        assert!(Regev::new(3, 16, 14, 2).is_err());
        assert!(Regev::new(3, 4, 13, 2).is_err());
        assert!(Regev::new(1, 16, 13, 20).is_err());
    }

    /// Checks whether `new_from_n` works properly for different choices of n.
    #[test]
    fn suitable_security_params() {
        let n_choices = [
            3, 4, 5, 6, 8, 9, 12, 13, 30, 21, 50, 100, 250, 500, 1000, 2500, 5000, 5001, 10000,
        ];

        for n in n_choices {
            assert!(Regev::new_from_n(n).is_ok());
        }
    }

    /// Checks whether the [`Default`] parameter choice is suitable.
    #[test]
    fn default_suitable() {
        let regev = Regev::default();

        assert!(Regev::check_params(&regev.n, &regev.m, &Z::from(&regev.q), &regev.alpha).is_ok());
    }

    /// Checks whether the generated public parameters from `new_from_n` are
    /// valid choices according to security and correctness of the scheme.
    #[test]
    fn choice_valid() {
        let n_choices = [3, 5, 8, 10, 14, 25, 50, 125, 300, 600, 1200, 4000, 6000];

        for n in n_choices {
            let regev = Regev::new_from_n(n).unwrap();
            assert!(
                Regev::check_params(&regev.n, &regev.m, &Z::from(&regev.q), &regev.alpha).is_ok()
            );
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

        let res = Regev::check_params(&regev.n, &regev.m, &Z::from(&regev.q), &regev.alpha);
        assert!(res.is_ok());
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
        let regev = Regev::new_from_n(5).unwrap();

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
        let regev = Regev::new_from_n(5).unwrap();

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
