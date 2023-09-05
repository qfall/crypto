// Copyright © 2023 Niklas Siemer
//
// This file is part of qFALL-crypto.
//
// qFALL-crypto is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

//! This module contains functions to quickly instantiate
//! common moduli for ring-based lattice cryptography.

use qfall_math::{
    error::MathError,
    integer_mod_q::{Modulus, ModulusPolynomialRingZq, PolyOverZq},
    traits::SetCoefficient,
};
use std::fmt::Display;

/// Outputs a [`ModulusPolynomialRingZq`] of the form `X^n + 1 mod modulus`.
///
/// Parameters:
/// - `n`: specifies the degree of the modulus polynomial
/// - `modulus`: specifies the modulus of the modulus polynomial
///
/// Returns a [`ModulusPolynomialRingZq`] of the form `X^n + 1 mod modulus` or
/// a [`MathError`] if `modulus <= 1`, `n < 0`, or `n` does not into an [`i64`].
///
/// # Examples
/// ```
/// use qfall_crypto::utils::common_moduli::new_anticyclic;
///
/// let poly_mod = new_anticyclic(8, 17);
/// ```
///
/// # Errors and Failures
/// - Returns a [`MathError`] of type [`OutOfBounds`](MathError::OutOfBounds) if
/// the `n` is negative or it does not fit into an [`i64`].
///
/// # Panics ...
/// - if the `modulus` is not larger than `1`.
pub fn new_anticyclic(
    n: impl TryInto<i64> + Display,
    modulus: impl Into<Modulus>,
) -> Result<ModulusPolynomialRingZq, MathError> {
    let mut poly = PolyOverZq::from((1, modulus));
    poly.set_coeff(n, 1)?;
    Ok(ModulusPolynomialRingZq::from(&poly))
}

/// Outputs a [`ModulusPolynomialRingZq`] of the form `X^n - 1 mod modulus`.
///
/// Parameters:
/// - `n`: specifies the degree of the modulus polynomial
/// - `modulus`: specifies the modulus of the modulus polynomial
///
/// Returns a [`ModulusPolynomialRingZq`] of the form `X^n - 1 mod modulus` or
/// a [`MathError`] if `modulus <= 1`, `n < 0`, or `n` does not into an [`i64`].
///
/// # Examples
/// ```
/// use qfall_crypto::utils::common_moduli::new_cyclic;
///
/// let poly_mod = new_cyclic(8, 17);
/// ```
///
/// # Errors and Failures
/// - Returns a [`MathError`] of type [`OutOfBounds`](MathError::OutOfBounds) if
/// the `n` is negative or it does not fit into an [`i64`].
///
/// # Panics ...
/// - if the `modulus` is not larger than `1`.
pub fn new_cyclic(
    n: impl TryInto<i64> + Display,
    modulus: impl Into<Modulus>,
) -> Result<ModulusPolynomialRingZq, MathError> {
    let mut poly = PolyOverZq::from((-1, modulus));
    poly.set_coeff(n, 1)?;
    Ok(ModulusPolynomialRingZq::from(&poly))
}

#[cfg(test)]
mod test_new_anticyclic {
    use super::new_anticyclic;
    use qfall_math::{integer::Z, integer_mod_q::PolyOverZq, traits::GetCoefficient};

    /// Ensure that the modulus polynomial has the specified degree.
    #[test]
    fn degree() {
        let degrees = [1, 4, 7, 16, 32, 128];
        for degree in degrees {
            let poly_mod = new_anticyclic(degree, 7).unwrap();

            assert_eq!(degree, poly_mod.get_degree());
        }
    }

    /// Check whether the method outputs the correct polynomial.
    #[test]
    fn correct_polynomial() {
        let degrees = [1, 4, 7, 16, 32, 128];
        for degree in degrees {
            let poly_mod = new_anticyclic(degree, 7).unwrap();
            let poly_zq = PolyOverZq::from(&poly_mod);

            assert_eq!(Z::ONE, poly_zq.get_coeff(degree).unwrap());
            assert_eq!(Z::ONE, poly_zq.get_coeff(0).unwrap());
            for i in 1..degree {
                assert_eq!(Z::ZERO, poly_zq.get_coeff(i).unwrap());
            }
        }
    }

    /// Ensures that the correct modulus is set as
    /// the integer modulus of the output modulus polynomial.
    #[test]
    fn correct_modulus() {
        let moduli = [7, 10, i64::MAX];
        for modulus in moduli {
            let poly_mod = new_anticyclic(2, modulus).unwrap();

            assert_eq!(Z::from(modulus), poly_mod.get_q());
        }
    }

    /// Ensures that an invalid degree for the modulus polynomial results in an error.
    #[test]
    fn invalid_n() {
        let res = new_anticyclic(-1, 7);

        assert!(res.is_err());
    }

    /// Ensures that an invalid modulus for the modulus polynomial results in a panic.
    #[test]
    #[should_panic]
    fn invalid_modulus() {
        let _ = new_anticyclic(2, 0);
    }
}

#[cfg(test)]
mod test_new_cyclic {
    use super::new_cyclic;
    use qfall_math::{integer::Z, integer_mod_q::PolyOverZq, traits::GetCoefficient};

    /// Ensure that the modulus polynomial has the specified degree.
    #[test]
    fn degree() {
        let degrees = [1, 4, 7, 16, 32, 128];
        for degree in degrees {
            let poly_mod = new_cyclic(degree, 7).unwrap();

            assert_eq!(degree, poly_mod.get_degree());
        }
    }

    /// Check whether the method outputs the correct polynomial.
    #[test]
    fn correct_polynomial() {
        let degrees = [1, 4, 7, 16, 32, 128];
        for degree in degrees {
            let poly_mod = new_cyclic(degree, 7).unwrap();
            let poly_zq = PolyOverZq::from(&poly_mod);

            assert_eq!(Z::ONE, poly_zq.get_coeff(degree).unwrap());
            assert_eq!(Z::from(6), poly_zq.get_coeff(0).unwrap());
            for i in 1..degree {
                assert_eq!(Z::ZERO, poly_zq.get_coeff(i).unwrap());
            }
        }
    }

    /// Ensures that the correct modulus is set as
    /// the integer modulus of the output modulus polynomial.
    #[test]
    fn correct_modulus() {
        let moduli = [7, 10, i64::MAX];
        for modulus in moduli {
            let poly_mod = new_cyclic(2, modulus).unwrap();

            assert_eq!(Z::from(modulus), poly_mod.get_q());
        }
    }

    /// Ensures that an invalid degree for the modulus polynomial results in an error.
    #[test]
    fn invalid_n() {
        let res = new_cyclic(-1, 7);

        assert!(res.is_err());
    }

    /// Ensures that an invalid modulus for the modulus polynomial results in a panic.
    #[test]
    #[should_panic]
    fn invalid_modulus() {
        let _ = new_cyclic(2, 0);
    }
}
