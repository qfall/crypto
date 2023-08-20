// Copyright © 2023 Niklas Siemer
//
// This file is part of qFALL-crypto.
//
// qFALL-crypto is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

//! This module contains an implementation of the IND-CPA secure
//! public key Ring-LPR encryption scheme.

use super::PKEncryption;
use qfall_math::{
    error::MathError,
    integer::{PolyOverZ, Z},
    integer_mod_q::{Modulus, ModulusPolynomialRingZq, PolyOverZq, PolynomialRingZq, Zq},
    rational::Q,
    traits::{Distance, GetCoefficient, Pow, SetCoefficient},
};
use serde::{Deserialize, Serialize};

/// This struct manages and stores the public parameters of a [`RingLPR`]
/// public key encryption instance.
///
/// Attributes:
/// - `n`: specifies the security parameter, which is not equal to the bit-security level
/// - `q`: specifies the modulus over which the encryption is computed
/// - `alpha`: specifies the gaussian parameter used for independent
/// sampling from the discrete Gaussian distribution
///
/// # Examples
/// ```
/// use qfall_crypto::construction::pk_encryption::{RingLPR, PKEncryption};
/// use qfall_math::integer::Z;
/// // setup public parameters and key pair
/// let lpr = RingLPR::default();
/// let (pk, sk) = lpr.gen();
///
/// // encrypt a bit
/// let msg = Z::ZERO; // must be a bit, i.e. msg = 0 or 1
/// let cipher = lpr.enc(&pk, &msg);
///
/// // decrypt
/// let m = lpr.dec(&sk, &cipher);
///
/// assert_eq!(msg, m);
/// ```
#[derive(Debug, Serialize, Deserialize)]
pub struct RingLPR {
    n: Z,                       // security parameter
    q: ModulusPolynomialRingZq, // modulus
    alpha: Q,                   // gaussian parameter for sampleZ
}

impl RingLPR {
    /// Instantiates a [`RingLPR`] PK encryption instance with the
    /// specified parameters.
    ///
    /// **WARNING:** The given parameters are not checked for security nor
    /// correctness of the scheme.
    /// If you want to check your parameters for provable security and correctness,
    /// use [`RingLPR::check_correctness`] and [`RingLPR::check_security`].
    /// Or use [`RingLPR::new_from_n`] for generating secure and correct
    /// public parameters for [`RingLPR`] according to your choice of `n`.
    ///
    /// Parameters:
    /// - `n`: specifies the security parameter and number of rows
    ///   of the uniform at random instantiated matrix `A`
    /// - `q`: specifies the modulus
    /// - `alpha`: specifies the gaussian parameter used for independent
    /// sampling from the discrete Gaussian distribution
    ///
    /// Returns a correct and secure [`RingLPR`] PK encryption instance or
    /// a [`MathError`] if the instance would not be correct or secure.
    ///
    /// # Examples
    /// ```
    /// use qfall_crypto::construction::pk_encryption::RingLPR;
    ///
    /// let lpr = RingLPR::new(3, 13, 2);
    /// ```
    ///
    /// # Panics ...
    /// - if the given modulus `q <= 1`.
    pub fn new(n: impl Into<Z>, q: impl Into<Modulus>, alpha: impl Into<Q>) -> Self {
        let n: Z = n.into();

        // mod = (X^n + 1) mod q
        let mut q = PolyOverZq::from((1, q));
        q.set_coeff(&n, 1).unwrap();
        let q = ModulusPolynomialRingZq::from(&q);

        let alpha: Q = alpha.into();

        Self { n, q, alpha }
    }

    /// Generates a new [`RingLPR`] instance, i.e. a new set of suitable
    /// (provably secure and correct) public parameters,
    /// given the security parameter `n` for `n >= 10`.
    ///
    /// Parameters:
    /// - `n`: specifies the security parameter and number of rows
    ///   of the uniform at random instantiated matrix `A`
    ///
    /// Returns a correct and secure [`RingLPR`] PK encryption instance or
    /// a [`MathError`] if the given `n < 10`.
    ///
    /// # Examples
    /// ```
    /// use qfall_crypto::construction::pk_encryption::RingLPR;
    ///
    /// let lpr = RingLPR::new_from_n(16);
    /// ```
    ///
    /// Panics...
    /// - if `n < 10`
    /// - if `n` does not fit into an [`i64`].
    pub fn new_from_n(n: impl Into<Z>) -> Self {
        let n = n.into();
        assert!(
            n >= Z::from(10),
            "Choose n >= 10 as this function does not return parameters ensuring proper correctness of the scheme otherwise."
        );

        let mut q: ModulusPolynomialRingZq;
        let mut alpha: Q;
        (q, alpha) = Self::gen_new_public_parameters(&n);
        let mut out = Self {
            n: n.clone(),
            q,
            alpha,
        };
        while out.check_correctness().is_err() || out.check_security().is_err() {
            (q, alpha) = Self::gen_new_public_parameters(&n);
            out = Self {
                n: n.clone(),
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
    /// Returns a set of public parameters `(q, alpha)` chosen according to
    /// the provided `n`.
    ///
    /// # Examples
    /// ```compile_fail
    /// use qfall_crypto::construction::pk_encryption::RingLPR;
    /// use qfall_math::integer::Z;
    /// let n = Z::from(2);
    ///
    /// let (q, alpha) = RingLPR::gen_new_public_parameters(&n);
    /// ```
    ///
    /// Panics...
    /// - if `n` does not fit into an [`i64`].
    fn gen_new_public_parameters(n: &Z) -> (ModulusPolynomialRingZq, Q) {
        let n_i64 = i64::try_from(n).unwrap();

        // generate prime q in [n^3 / 2, n^3]
        let upper_bound: Z = n.pow(3).unwrap();
        let lower_bound = upper_bound.div_ceil(&Z::from(2));
        // prime used due to guide from GPV08 after Proposition 8.1
        // on how to choose appropriate parameters, but prime is not
        // necessarily needed for this scheme to be correct or secure
        let q = Z::sample_prime_uniform(&lower_bound, &upper_bound).unwrap();

        // Found out by experience as the bound is not tight enough to ensure correctness for large n.
        // Hence, a small factor roughly of max(log n - 4, 1) has to be applied.
        // Checked for 100 parameter sets with 10 cyles each and no mismatching decryptions occurred.
        let factor = match n_i64 {
            1..=20 => 1,
            21..=40 => 2,
            41..=80 => 3,
            81..=160 => 4,
            _ => 5,
        };
        // α = 1/(sqrt(n) * log^2 n)
        let alpha = 1 / (factor * n.sqrt() * n.log(2).unwrap().pow(3).unwrap());

        // mod = (X^n + 1) mod q
        let mut q = PolyOverZq::from((1, q));
        q.set_coeff(n, 1).unwrap();
        let q = ModulusPolynomialRingZq::from(&q);

        (q, alpha)
    }

    /// Checks the public parameters for
    /// correctness according to Lemma 3.1 of [\[4\]](<index.html#:~:text=[4]>).
    ///
    /// The required properties are:
    /// - α = o (1 / (sqrt(n) * log^3 n))
    ///
    /// **WARNING**: This bound is not tight. Hence, we added a small factor
    /// loosely corresponding to max(log n - 4, 1) below to ensure correctness
    /// with overwhelming proability.
    ///
    /// Returns an empty result if the public parameters guarantee correctness
    /// with overwhelming probability or a [`MathError`] if the instance would
    /// not be correct.
    ///
    /// # Examples
    /// ```
    /// use qfall_crypto::construction::pk_encryption::RingLPR;
    /// let lpr = RingLPR::default();
    ///
    /// let is_valid = lpr.check_correctness().is_ok();
    /// ```
    ///
    /// # Errors and Failures
    /// - Returns a [`MathError`] of type [`InvalidIntegerInput`](MathError::InvalidIntegerInput)
    /// if at least one parameter was not chosen appropriately for a
    /// correct RingLPR public key encryption instance.
    /// - Returns a [`MathError`] of type [`ConversionError`](MathError::ConversionError)
    /// if the value does not fit into an [`i64`]
    pub fn check_correctness(&self) -> Result<(), MathError> {
        let n_i64 = i64::try_from(&self.n)?;

        if self.n <= Z::ONE {
            return Err(MathError::InvalidIntegerInput(String::from(
                "n must be chosen bigger than 1.",
            )));
        }

        // TODO: Add check whether n is power of two

        // Found out by experience as the bound is not tight enough to ensure correctness for large n.
        // Hence, a small factor roughly of max(log n - 4, 1) has to be applied.
        // Checked for 100 parameter sets with 10 cyles each and no mismatching decryptions occurred.
        let factor = match n_i64 {
            1..=20 => 1,
            21..=40 => 2,
            41..=80 => 3,
            81..=160 => 4,
            _ => 5,
        };
        // α = o (1 / sqrt(n) * log^3 n ))
        if self.alpha > 1 / (factor * self.n.sqrt() * self.n.log(2).unwrap().pow(3).unwrap()) {
            return Err(MathError::InvalidIntegerInput(String::from(
                "Correctness is not guaranteed as α >= 1 / (sqrt(n) * log^3 n),\
                but α < 1 / (sqrt(n) * log^3 n) is required. Please check the documentation!",
            )));
        }

        Ok(())
    }

    /// Checks the public parameters for security according to Section 2.2
    /// and Lemma 3.2 of [\[4\]](<index.html#:~:text=[4]>).
    ///
    /// The required properties are:
    /// - q * α >= 2 sqrt(n)
    ///
    /// Returns an empty result if the public parameters guarantee security
    /// w.r.t. `n` or a [`MathError`] if the instance would
    /// not be secure.
    ///
    /// # Examples
    /// ```
    /// use qfall_crypto::construction::pk_encryption::RingLPR;
    /// let lpr = RingLPR::default();
    ///
    /// let is_valid = lpr.check_security().is_ok();
    /// ```
    ///
    /// # Errors and Failures
    /// - Returns a [`MathError`] of type [`InvalidIntegerInput`](MathError::InvalidIntegerInput)
    /// if at least one parameter was not chosen appropriately for a
    /// secure RingLPR public key encryption instance.
    pub fn check_security(&self) -> Result<(), MathError> {
        let q = Z::from(&self.q.get_q());

        // Security requirements
        // q * α >= 2 sqrt(n)
        if &q * &self.alpha < 2 * self.n.sqrt() {
            return Err(MathError::InvalidIntegerInput(String::from(
                "Security is not guaranteed as q * α < 2 * sqrt(n), but q * α >= 2 * sqrt(n) is required.",
            )));
        }

        Ok(())
    }

    /// This function instantiates a 128-bit secure [`RingLPR`] scheme.
    ///
    /// The public parameters used for this scheme were generated via `RingLPR::new_from_n(350)`
    /// and its bit-security determined via the [lattice estimator](https://github.com/malb/lattice-estimator).
    pub fn secure128() -> Self {
        Self::new(512, 92897729, 0.000005)
    }
}

impl Default for RingLPR {
    /// Initializes a [`RingLPR`] struct with parameters generated by `RingLPR::new_from_n(3)`.
    /// This parameter choice is not secure as the dimension of the lattice is too small,
    /// but it provides an efficient working example.
    ///
    /// # Examples
    /// ```
    /// use qfall_crypto::construction::pk_encryption::RingLPR;
    ///
    /// let lpr = RingLPR::default();
    /// ```
    fn default() -> Self {
        Self::new(16, 2399, 0.0039)
    }
}

impl PKEncryption for RingLPR {
    type Cipher = (PolynomialRingZq, PolynomialRingZq);
    type PublicKey = (PolynomialRingZq, PolynomialRingZq);
    type SecretKey = PolynomialRingZq;

    /// Generates a (pk, sk) pair for the RingLPR public key encryption scheme
    /// by following these steps:
    /// - a <- R_q
    /// - s <- χ
    /// - e <- χ
    /// - b = s * a + e
    /// where χ is discrete Gaussian distributed with center 0 and Gaussian parameter q * α.
    ///
    /// Then, `pk = (a, b)` and `sk = s` are returned.
    ///
    /// # Examples
    /// ```
    /// use qfall_crypto::construction::pk_encryption::{PKEncryption, RingLPR};
    /// let lpr = RingLPR::default();
    ///
    /// let (pk, sk) = lpr.gen();
    /// ```
    fn gen(&self) -> (Self::PublicKey, Self::SecretKey) {
        // a <- R_q
        let a = PolynomialRingZq::sample_uniform(&self.q);
        // s <- χ
        let s = PolynomialRingZq::sample_discrete_gauss(
            &self.q,
            &self.n,
            0,
            &self.alpha * &self.q.get_q(),
        )
        .unwrap();
        // e <- χ
        let e = PolynomialRingZq::sample_discrete_gauss(
            &self.q,
            &self.n,
            0,
            &self.alpha * &self.q.get_q(),
        )
        .unwrap();

        // b = s * a + e
        let b = &a * &s + e;

        // pk = (a, b), sk = s
        ((a, b), s)
    }

    /// Generates an encryption of `message mod 2^(n-2)` for the provided public key by following these steps:
    /// - r <- χ
    /// - e1 <- χ
    /// - e2 <- χ
    /// - u = a * r + e1
    /// - v = b * r + e2 + mu * q/2
    /// - c = (u, v)
    /// where χ is discrete Gaussian distributed with center 0 and Gaussian parameter q * α.
    ///
    /// Then, cipher `c = (u, v)` as a polynomial of type [`PolynomialRingZq`] is returned.
    ///
    /// Parameters:
    /// - `pk`: specifies the public key `pk = (a, b)`
    /// - `message`: specifies the message that should be encryted
    ///
    /// Returns a cipher `c = (u, v)` of types [`PolynomialRingZq`].
    ///
    /// # Examples
    /// ```
    /// use qfall_crypto::construction::pk_encryption::{PKEncryption, RingLPR};
    /// let lpr = RingLPR::default();
    /// let (pk, sk) = lpr.gen();
    ///
    /// let cipher = lpr.enc(&pk, 15);
    /// ```
    fn enc(&self, pk: &Self::PublicKey, message: impl Into<Z>) -> Self::Cipher {
        // ensure mu has at most n bits
        let message: Z = message.into().abs();
        let mu = message.modulo(Z::from(2).pow(&self.n - 2).unwrap());
        // set mu_q_half to polynomial with n {0,1} coefficients
        let bits = mu.to_bits();
        let mut mu_q_half = PolynomialRingZq::from((&PolyOverZ::default(), &self.q));
        let q_half = self.q.get_q().div_floor(&Z::from(2));
        for (i, bit) in bits.iter().enumerate() {
            if *bit {
                mu_q_half.set_coeff(i, &q_half).unwrap();
            }
        }

        // r <- χ
        let r = PolynomialRingZq::sample_discrete_gauss(
            &self.q,
            &self.n,
            0,
            &self.alpha * &self.q.get_q(),
        )
        .unwrap();
        // e1 <- χ
        let e1 = PolynomialRingZq::sample_discrete_gauss(
            &self.q,
            &self.n,
            0,
            &self.alpha * &self.q.get_q(),
        )
        .unwrap();
        // e2 <- χ
        let e2 = PolynomialRingZq::sample_discrete_gauss(
            &self.q,
            &self.n,
            0,
            &self.alpha * &self.q.get_q(),
        )
        .unwrap();

        // u = a * r + e1
        let u = &pk.0 * &r + e1;
        // v = b * r + e2 + mu * q/2
        let v = &pk.1 * &r + e2 + mu_q_half;

        // c = (u, v)
        (u, v)
    }

    /// Decrypts the provided `cipher` using the secret key `sk` by following these steps:
    /// - v - s * u
    /// - result = 0
    /// - for each coefficient of v - s * u:
    ///   - check whether the coefficient mod q is closer to ⌊q/2⌋ than to 0.
    ///     If so, add 2^coefficient to result.
    /// - return result
    ///
    /// Parameters:
    /// - `sk`: specifies the secret key `sk = s`
    /// - `cipher`: specifies the cipher containing `cipher = (u, v)`
    ///
    /// Returns the decryption of `cipher` as a [`Z`] instance.
    ///
    /// # Examples
    /// ```
    /// use qfall_crypto::construction::pk_encryption::{PKEncryption, RingLPR};
    /// use qfall_math::integer::Z;
    /// let lpr = RingLPR::default();
    /// let (pk, sk) = lpr.gen();
    /// let cipher = lpr.enc(&pk, 212);
    ///
    /// let m = lpr.dec(&sk, &cipher);
    ///
    /// assert_eq!(Z::from(212), m);
    /// ```
    fn dec(&self, sk: &Self::SecretKey, cipher: &Self::Cipher) -> Z {
        // res = v - s * u
        let result = &cipher.1 - sk * &cipher.0;

        let q_half = self.q.get_q().div_floor(&Z::from(2));

        // check for each coefficient whether it's closer to 0 or q/2
        // if closer to q/2 -> add 2^i to result
        let mut vec = vec![];
        for i in 0..result.get_degree() {
            let coeff = Zq::from((result.get_coeff(i).unwrap(), self.q.get_q()));
            if coeff.distance(&q_half) < coeff.distance(Z::ZERO) {
                vec.push(true);
            } else {
                vec.push(false);
            }
        }

        Z::from_bits(&vec)
    }
}

#[cfg(test)]
mod test_pp_generation {
    use super::RingLPR;
    use super::Z;

    /// Checks whether `new` is available for types implementing [`Into<Z>`].
    #[test]
    fn new_availability() {
        let _ = RingLPR::new(2u8, 2u32, 2u64);
        let _ = RingLPR::new(2u16, 2i32, 2i64);
        let _ = RingLPR::new(2i16, 2u32, 2u8);
        let _ = RingLPR::new(Z::from(2), 2u8, 2i8);
    }

    /// Checks whether `new_from_n` works properly for different choices of n.
    #[test]
    fn suitable_security_params() {
        let n_choices = [16, 32, 64, 128, 256, 512, 1024];

        for n in n_choices {
            let _ = RingLPR::new_from_n(n);
        }
    }

    /// Checks whether the [`Default`] parameter choice is suitable.
    #[test]
    fn default_suitable() {
        let scheme = RingLPR::default();

        assert!(scheme.check_correctness().is_ok());
        assert!(scheme.check_security().is_ok());
    }

    /// Checks whether the generated public parameters from `new_from_n` are
    /// valid choices according to security and correctness of the scheme.
    #[test]
    fn choice_valid() {
        let n_choices = [16, 32, 64, 128, 256, 512, 1024];

        for n in n_choices {
            let scheme = RingLPR::new_from_n(n);

            assert!(scheme.check_correctness().is_ok());
            assert!(scheme.check_security().is_ok());
        }
    }

    /// Ensures that `new_from_n` is available for types implementing [`Into<Z>`].
    #[test]
    fn availability() {
        let _ = RingLPR::new_from_n(16u8);
        let _ = RingLPR::new_from_n(16u16);
        let _ = RingLPR::new_from_n(16u32);
        let _ = RingLPR::new_from_n(16u64);
        let _ = RingLPR::new_from_n(16i8);
        let _ = RingLPR::new_from_n(16i16);
        let _ = RingLPR::new_from_n(16i32);
        let _ = RingLPR::new_from_n(16i64);
        let _ = RingLPR::new_from_n(Z::from(16));
        let _ = RingLPR::new_from_n(&Z::from(16));
    }

    /// Checks whether `new_from_n` returns an error for invalid input n.
    #[test]
    #[should_panic]
    fn invalid_n() {
        RingLPR::new_from_n(9);
    }

    /// Checks whether `secure128` outputs a new instance with correct and secure parameters.
    #[test]
    fn secure128_validity() {
        let scheme = RingLPR::secure128();

        assert!(scheme.check_correctness().is_ok());
        assert!(scheme.check_security().is_ok());
    }
}

#[cfg(test)]
mod test_ring_lpr {
    use super::RingLPR;
    use crate::construction::pk_encryption::PKEncryption;
    use qfall_math::integer::Z;

    /// Checks whether the full-cycle of gen, enc, dec works properly
    /// for several messages and small n.
    #[test]
    fn cycle_small_n() {
        let scheme = RingLPR::default();
        let (pk, sk) = scheme.gen();
        let messages = [0, 1, 2, 15, 70, 256, 580, 1000, 4000, 8000];

        for message in messages {
            let cipher = scheme.enc(&pk, message);
            let m = scheme.dec(&sk, &cipher);

            assert_eq!(Z::from(message), m);
        }
    }

    /// Checks whether the full-cycle of gen, enc, dec works properly
    /// for several messages and larger n.
    #[test]
    fn cycle_large_n() {
        let scheme = RingLPR::new_from_n(64);
        let (pk, sk) = scheme.gen();
        let messages = [
            0,
            1,
            2,
            15,
            70,
            256,
            580,
            1_000,
            4_000,
            8_000,
            20_000,
            80_000,
            240_000,
            4_000_000,
            100_000_000,
        ];

        for message in messages {
            let cipher = scheme.enc(&pk, message);
            let m = scheme.dec(&sk, &cipher);

            assert_eq!(Z::from(message), m);
        }
    }

    /// Checks that modulus 2^(n-1) is applied correctly.
    #[test]
    fn modulus_application() {
        let messages = [32768];
        let scheme = RingLPR::default();
        let (pk, sk) = scheme.gen();

        for msg in messages {
            let cipher = scheme.enc(&pk, &msg);
            let m = scheme.dec(&sk, &cipher);

            assert_eq!(Z::ZERO, m);
        }
    }
}
