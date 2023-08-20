// Copyright © 2023 Niklas Siemer
//
// This file is part of qFALL-crypto.
//
// qFALL-crypto is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

//! This module contains an implementation of the IND-CPA secure
//! public key Dual Regev encryption scheme.

use super::{GenericMultiBitEncryption, PKEncryption};
use qfall_math::{
    error::MathError,
    integer::{MatZ, Z},
    integer_mod_q::{MatZq, Modulus, Zq},
    rational::Q,
    traits::{Concatenate, Distance, GetEntry, Pow, SetEntry},
};
use serde::{Deserialize, Serialize};

/// This struct manages and stores the public parameters of a [`DualRegev`]
/// public key encryption instance.
///
/// Attributes:
/// - `n`: specifies the security parameter, which is not equal to the bit-security level
/// - `m`: defines the dimension of the underlying lattice
/// - `q`: specifies the modulus over which the encryption is computed
/// - `alpha`:  specifies the gaussian parameter used for independent
/// sampling from the discrete Gaussian distribution
///
/// # Examples
/// ```
/// use qfall_crypto::construction::pk_encryption::{DualRegev, PKEncryption};
/// use qfall_math::integer::Z;
/// // setup public parameters and key pair
/// let dual_regev = DualRegev::default();
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
pub struct DualRegev {
    n: Z,       // security parameter
    m: Z,       // number of rows of matrix A
    q: Modulus, // modulus
    alpha: Q,   // gaussian parameter for sampleZ
}

impl DualRegev {
    /// Instantiates a [`DualRegev`] PK encryption instance with the
    /// specified parameters.
    ///
    /// **WARNING:** The given parameters are not checked for security nor
    /// correctness of the scheme.
    /// If you want to check your parameters for provable security and correctness,
    /// use [`DualRegev::check_correctness`] and [`DualRegev::check_security`].
    /// Or use [`DualRegev::new_from_n`] for generating secure and correct
    /// public parameters for [`DualRegev`] according to your choice of `n`.
    ///
    /// Parameters:
    /// - `n`: specifies the security parameter and number of rows
    ///   of the uniform at random instantiated matrix `A`
    /// - `m`: specifies the number of columns of matrix `A`
    /// - `q`: specifies the modulus
    /// - `alpha`:  specifies the gaussian parameter used for independent
    /// sampling from the discrete Gaussian distribution
    ///
    /// Returns [`DualRegev`] PK encryption instance.
    ///
    /// # Examples
    /// ```
    /// use qfall_crypto::construction::pk_encryption::DualRegev;
    ///
    /// let dual_regev = DualRegev::new(3, 16, 13, 2);
    /// ```
    ///
    /// # Panics ...
    /// - if the given modulus `q <= 1`.
    pub fn new(
        n: impl Into<Z>,
        m: impl Into<Z>,
        q: impl Into<Modulus>,
        alpha: impl Into<Q>,
    ) -> Self {
        let n: Z = n.into();
        let m: Z = m.into();
        let q: Modulus = q.into();
        let alpha: Q = alpha.into();

        Self { n, m, q, alpha }
    }

    /// Generates a new [`DualRegev`] instance, i.e. a new set of suitable
    /// (provably secure and correct) public parameters,
    /// given the security parameter `n` for `n >= 10`.
    ///
    /// Parameters:
    /// - `n`: specifies the security parameter and number of rows
    ///   of the uniform at random instantiated matrix `A`
    ///
    /// Returns a correct and secure [`DualRegev`] PK encryption instance or
    /// a [`MathError`] if the given `n < 10`.
    ///
    /// # Examples
    /// ```
    /// use qfall_crypto::construction::pk_encryption::DualRegev;
    ///
    /// let dual_regev = DualRegev::new_from_n(15);
    /// ```
    ///
    /// Panics...
    /// - if `n < 10`.
    /// - if `n` does not fit into an [`i64`].
    pub fn new_from_n(n: impl Into<Z>) -> Self {
        let n = n.into();
        if n < Z::from(10) {
            panic!("Choose n >= 10 as this function does not return parameters ensuring proper correctness of the scheme otherwise.");
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

        out
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
    /// use qfall_crypto::construction::pk_encryption::DualRegev;
    /// use qfall_math::integer::Z;
    /// let n = Z::from(2);
    ///
    /// let (m, q, alpha) = DualRegev::gen_new_public_parameters(&n);
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
        // prime used due to guide from GPV08 after Proposition 8.1
        // on how to choose appropriate parameters, but prime is not
        // necessarily needed for this scheme to be correct or secure
        let q = Z::sample_prime_uniform(&lower_bound, &upper_bound).unwrap();

        // choose m = (n+1) log q
        let m = (n + Z::ONE) * q.log(2).unwrap().ceil();

        // α = 1/(2 * sqrt(n) * log^2 n)
        let alpha = 1 / (2 * n.sqrt() * n.log(2).unwrap().pow(2).unwrap());

        let q = Modulus::from(q);

        (m, q, alpha)
    }

    /// Checks the public parameters for
    /// correctness according to Lemma 5.1 of [\[3\]](<index.html#:~:text=[3]>).
    ///
    /// The required properties are:
    /// - α = o (1 / ( sqrt(n) * log n ) )
    /// - concentration bound with r=5: r * sqrt(m) * α > q/4
    ///
    /// **WARNING:** Some requirements are missing to ensure overwhelming correctness of the scheme.
    ///
    /// Returns an empty result if the public parameters guarantee correctness
    /// with overwhelming probability or a [`MathError`] if the instance would
    /// not be correct.
    ///
    /// # Examples
    /// ```
    /// use qfall_crypto::construction::pk_encryption::DualRegev;
    /// let dual_regev = DualRegev::default();
    ///
    /// let is_valid = dual_regev.check_correctness().is_ok();
    /// ```
    ///
    /// # Errors and Failures
    /// - Returns a [`MathError`] of type [`InvalidIntegerInput`](MathError::InvalidIntegerInput)
    /// if at least one parameter was not chosen appropriately for a
    /// correct Dual Regev public key encryption instance.
    pub fn check_correctness(&self) -> Result<(), MathError> {
        let q = Z::from(&self.q);

        if self.n <= Z::ONE {
            return Err(MathError::InvalidIntegerInput(String::from(
                "n must be chosen bigger than 1.",
            )));
        }

        // Correctness requirements
        // α = o (1 / ( sqrt(n) * log n ) )
        if self.alpha > 1 / (self.n.sqrt() * self.n.log(2).unwrap()) {
            return Err(MathError::InvalidIntegerInput(String::from(
                "Correctness is not guaranteed as α >= 1 / (sqrt(n) * log n), but α < 1 / (sqrt(n) * log n) is required."
            )));
        }
        // concentration bound with r=5 -> r * sqrt(m) * α > q/4
        if 20 * self.m.sqrt() * &self.alpha > Q::from(q) {
            return Err(MathError::InvalidIntegerInput(String::from(
                "Correctness is not guaranteed as 5 * sqrt(m) * α > q/4, but 5 * sqrt(m) * α <= q/4 is required."
            )));
        }

        Ok(())
    }

    /// Checks the public parameters for security according to Theorem 1.1
    /// and Lemma 5.4 of [\[3\]](<index.html#:~:text=[3]>).
    ///
    /// The required properties are:
    /// - q * α >= 2 sqrt(n)
    /// - m > (n + 1) log q
    ///
    /// Returns an empty result if the public parameters guarantees security w.r.t. `n`
    /// or a [`MathError`] if the instance would not be secure.
    ///
    /// # Examples
    /// ```
    /// use qfall_crypto::construction::pk_encryption::DualRegev;
    /// let dual_regev = DualRegev::default();
    ///
    /// let is_valid = dual_regev.check_security().is_ok();
    /// ```
    ///
    /// # Errors and Failures
    /// - Returns a [`MathError`] of type [`InvalidIntegerInput`](MathError::InvalidIntegerInput)
    /// if at least one parameter was not chosen appropriately for a
    /// secure Dual Regev public key encryption instance.
    pub fn check_security(&self) -> Result<(), MathError> {
        let q = Z::from(&self.q);

        // Security requirements
        // q * α >= 2 sqrt(n)
        if &q * &self.alpha < 2 * self.n.sqrt() {
            return Err(MathError::InvalidIntegerInput(String::from(
                "Security is not guaranteed as q * α < 2 * sqrt(n), but q * α >= 2 * sqrt(n) is required.",
            )));
        }
        // m > (n + 1) log q
        if self.m <= ((&self.n + Z::ONE) * q.log(2).unwrap()).ceil() {
            return Err(MathError::InvalidIntegerInput(String::from(
                "Security is not guaranteed as m <= (n + 1) log q, but m > (n + 1) log q is required.",
            )));
        }

        Ok(())
    }

    /// This function instantiates a 128-bit secure [`DualRegev`] scheme.
    ///
    /// The public parameters used for this scheme were generated via `DualRegev::new_from_n(350)`
    /// and its bit-security determined via the [lattice estimator](https://github.com/malb/lattice-estimator).
    pub fn secure128() -> Self {
        Self::new(230, 5313, 7764299, 0.0011)
    }
}

impl Default for DualRegev {
    /// Initializes a [`DualRegev`] struct with parameters generated by `DualRegev::new_from_n(13)`.
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
        let n = Z::from(13);
        let m = Z::from(154);
        let q = Modulus::from(1427);
        let alpha = Q::from(0.01);

        Self { n, m, q, alpha }
    }
}

impl PKEncryption for DualRegev {
    type Cipher = MatZq;
    type PublicKey = MatZq;
    type SecretKey = MatZ;

    /// Generates a (pk, sk) pair for the Dual Regev public key encryption scheme
    /// by following these steps:
    /// - A <- Z_q^{n x m}
    /// - x <- {0,1}^m
    /// - u = A * x
    /// - A = [A | u]
    ///
    /// Then, `pk = A` and `sk = x` of type [`MatZq`] are returned.
    ///
    /// # Examples
    /// ```
    /// use qfall_crypto::construction::pk_encryption::{PKEncryption, DualRegev};
    /// let dual_regev = DualRegev::default();
    ///
    /// let (pk, sk) = dual_regev.gen();
    /// ```
    fn gen(&self) -> (Self::PublicKey, Self::SecretKey) {
        // A <- Z_q^{n x m}
        let mat_a = MatZq::sample_uniform(&self.n, &self.m, &self.q);
        // x <- Z_2^m
        let vec_x = MatZ::sample_uniform(&self.m, 1, 0, 2).unwrap();

        // u = A * x
        let vec_u = &mat_a * &vec_x;

        // A = [A | u]
        let mat_a = mat_a.concat_horizontal(&vec_u).unwrap();

        // pk = A, sk = x
        (mat_a, vec_x)
    }

    /// Generates an encryption of `message mod 2` for the provided public key by following these steps:
    /// - s <- Z_q^n
    /// - e <- χ^(m+1)
    /// - c^t = s^t * A + e^t + [0^{1xn} | msg *  ⌊q/2⌋]
    /// where χ is discrete Gaussian distributed with center 0 and Gaussian parameter q * α.
    ///
    /// Then, the ciphertext `c` is returned as a vector of type [`MatZq`].
    ///
    /// Parameters:
    /// - `pk`: specifies the public key `pk = A`
    /// - `message`: specifies the message that should be encryted
    ///
    /// Returns a cipher `c` of type [`MatZq`].
    ///
    /// # Examples
    /// ```
    /// use qfall_crypto::construction::pk_encryption::{PKEncryption, DualRegev};
    /// let dual_regev = DualRegev::default();
    /// let (pk, sk) = dual_regev.gen();
    ///
    /// let cipher = dual_regev.enc(&pk, 1);
    /// ```
    fn enc(&self, pk: &Self::PublicKey, message: impl Into<Z>) -> Self::Cipher {
        // generate message = message mod 2
        let message: Z = message.into().modulo(2);

        // s <- Z_q^n
        let vec_s_t = MatZq::sample_uniform(1, &self.n, &self.q);
        // e <- χ^(m+1)
        let vec_e_t = MatZq::sample_discrete_gauss(
            1,
            &(&self.m + 1),
            &self.q,
            &self.n,
            0,
            &self.alpha * Z::from(&self.q),
        )
        .unwrap();

        // c^t = s^t * A + e^t + [0^{1xn} | msg *  ⌊q/2⌋]
        let mut c = (vec_s_t * pk + vec_e_t).transpose();

        // hide message in last entry
        // compute msg * ⌊q/2⌋
        let msg_q_half = message * Z::from(&self.q).div_floor(&Z::from(2));
        // set last entry of c = last_entry + msg * ⌊q/2⌋
        let last_entry: Zq = c.get_entry(-1, 0).unwrap();
        c.set_entry(-1, 0, last_entry + msg_q_half).unwrap();

        c
    }

    /// Decrypts the provided `cipher` using the secret key `sk` by following these steps:
    /// - x = c^t * [-sk^t | 1]^t
    /// - if x mod q is closer to ⌊q/2⌋ than to 0, output 1. Otherwise, output 0.
    ///
    /// Parameters:
    /// - `sk`: specifies the secret key `sk = x`
    /// - `cipher`: specifies the cipher containing `cipher = c`
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

// adds generic multi-bit encryption to this scheme
impl GenericMultiBitEncryption for DualRegev {}

#[cfg(test)]
mod test_pp_generation {
    use super::DualRegev;
    use super::Z;

    /// Checks whether `new` is available for types implementing [`Into<Z>`].
    #[test]
    fn new_availability() {
        let _ = DualRegev::new(2u8, 2u16, 2u32, 2u64);
        let _ = DualRegev::new(2u16, 2u64, 2i32, 2i64);
        let _ = DualRegev::new(2i16, 2i64, 2u32, 2u8);
        let _ = DualRegev::new(Z::from(2), &Z::from(2), 2u8, 2i8);
    }

    /// Checks whether `new_from_n` works properly for different choices of n.
    #[test]
    fn suitable_security_params() {
        let n_choices = [
            10, 11, 12, 13, 14, 25, 50, 100, 250, 500, 1000, 2500, 5000, 5001, 10000,
        ];

        for n in n_choices {
            let _ = DualRegev::new_from_n(n);
        }
    }

    /// Checks whether the [`Default`] parameter choice is suitable.
    #[test]
    fn default_suitable() {
        let dr = DualRegev::default();

        assert!(dr.check_correctness().is_ok());
        assert!(dr.check_security().is_ok());
    }

    /// Checks whether the generated public parameters from `new_from_n` are
    /// valid choices according to security and correctness of the scheme.
    #[test]
    fn choice_valid() {
        let n_choices = [10, 14, 25, 50, 125, 300, 600, 1200, 4000, 6000];

        for n in n_choices {
            let dr = DualRegev::new_from_n(n);
            assert!(dr.check_correctness().is_ok());
            assert!(dr.check_security().is_ok());
        }
    }

    /// Ensures that `new_from_n` is available for types implementing [`Into<Z>`].
    #[test]
    fn availability() {
        let _ = DualRegev::new_from_n(10u8);
        let _ = DualRegev::new_from_n(10u16);
        let _ = DualRegev::new_from_n(10u32);
        let _ = DualRegev::new_from_n(10u64);
        let _ = DualRegev::new_from_n(10i8);
        let _ = DualRegev::new_from_n(10i16);
        let _ = DualRegev::new_from_n(10i32);
        let _ = DualRegev::new_from_n(10i64);
        let _ = DualRegev::new_from_n(Z::from(10));
        let _ = DualRegev::new_from_n(&Z::from(10));
    }

    /// Checks whether `new_from_n` returns an error for invalid input n.
    #[test]
    #[should_panic]
    fn invalid_n() {
        DualRegev::new_from_n(9);
    }

    /// Checks whether `secure128` outputs a new instance with correct and secure parameters.
    #[test]
    fn secure128_validity() {
        let dr = DualRegev::secure128();

        assert!(dr.check_correctness().is_ok());
        assert!(dr.check_security().is_ok());
    }
}

#[cfg(test)]
mod test_dual_regev {
    use super::DualRegev;
    use crate::construction::pk_encryption::PKEncryption;
    use qfall_math::integer::Z;

    /// Checks whether the full-cycle of gen, enc, dec works properly
    /// for message 0 and small n.
    #[test]
    fn cycle_zero_small_n() {
        let msg = Z::ZERO;
        let dr = DualRegev::default();

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
        let dr = DualRegev::default();

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
        let dr = DualRegev::new_from_n(50);

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
        let dr = DualRegev::new_from_n(50);

        let (pk, sk) = dr.gen();
        let cipher = dr.enc(&pk, &msg);
        let m = dr.dec(&sk, &cipher);
        assert_eq!(msg, m);
    }

    /// Checks that modulus 2 is applied correctly.
    #[test]
    fn modulus_application() {
        let messages = [2, 3, i64::MAX, i64::MIN];
        let dr = DualRegev::default();
        let (pk, sk) = dr.gen();

        for msg in messages {
            let msg_mod = Z::from(msg.rem_euclid(2));

            let cipher = dr.enc(&pk, &msg);
            let m = dr.dec(&sk, &cipher);

            assert_eq!(msg_mod, m);
        }
    }
}

#[cfg(test)]
mod test_multi_bits {
    use super::{DualRegev, GenericMultiBitEncryption, PKEncryption};
    use qfall_math::integer::Z;

    /// Checks whether the multi-bit encryption cycle works properly
    /// for small and large positive values.
    #[test]
    fn positive() {
        let values = [3, 13, 23, 230, 501, 1024, i64::MAX];

        for value in values {
            let msg = Z::from(value);
            let scheme = DualRegev::default();

            let (pk, sk) = scheme.gen();
            let cipher = scheme.enc_multiple_bits(&pk, &msg);
            let m = scheme.dec_multiple_bits(&sk, &cipher);

            assert_eq!(msg, m);
        }
    }

    /// Checks whether the multi-bit encryption cycle works properly
    /// for zero.
    #[test]
    fn zero() {
        let msg = Z::ZERO;
        let scheme = DualRegev::default();

        let (pk, sk) = scheme.gen();
        let cipher = scheme.enc_multiple_bits(&pk, &msg);
        let m = scheme.dec_multiple_bits(&sk, &cipher);

        assert_eq!(msg, m);
    }

    /// Checks whether the multi-bit encryption cycle works properly
    /// for small and large negative values, which are not encrypted itself,
    /// but their absolute value.
    #[test]
    fn negative() {
        let values = [-3, -13, -23, -230, -501, -1024, i64::MIN];

        for value in values {
            let msg = Z::from(value);
            let scheme = DualRegev::default();

            let (pk, sk) = scheme.gen();
            let cipher = scheme.enc_multiple_bits(&pk, &msg);
            let m = scheme.dec_multiple_bits(&sk, &cipher);

            assert_eq!(msg.abs(), m);
        }
    }
}
