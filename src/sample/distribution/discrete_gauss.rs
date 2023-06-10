// Copyright © 2023 Marvin Beckmann
//
// This file is part of qFALL-crypto.
//
// qFALL-crypto is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

//! This module contains an implementation to sample from a lattice using a discrete gaussian distribution.

use qfall_math::{
    integer::{MatZ, Z},
    integer_mod_q::{MatZq, Modulus, Zq},
    rational::Q,
    traits::{Concatenate, GetEntry, GetNumColumns, GetNumRows, Pow, SetEntry},
};

/// Efficiently samples a preimage for a syndrome `u` with a trapdoor gadget-trapdoor `R`
/// according to Algorithm 3 in [\[1\]](<../index.html#:~:text=[1]>).
///
/// *Note*: This Function does not provide security as the perturbation vector is not yet
/// implemented correctly.
/// At the moment this function only works for tags, that correspond to an identity
/// matrix, this will be changed in the future. If its not the identity, the function
/// will panic.
///
/// TODO: Main issue: the length of the final vector should depend on `s`.
/// Formally `s` is only important in the presence of the perturbation vector.
/// Otherwise the `s` does not have any influence on the length of the final vector.
/// THe vector currently returned therefore does not have a length associated in any way
/// with `s`.
///
/// Namely:
/// - Samples perturbation: `p`
/// - Computes the inverse of the tag
/// - Computes the shifted syndrome `v = H^{-1}(u - (a*p))`
/// - Sample a preimage of the syndrome with the trapdoor.
/// - Concatenate all vectors to a new vector `z`
/// - Return `p + [R^t | I]^t z`
///
/// Parameters:
/// - `a`: the parity check matrix
/// - `trapdoor`: the g-trapdoor for the matrix `a`
/// - `u`: the syndrome for which sample_d has to compute a preimage
/// - `k`: the length of the gadget vector
/// - `s`: the standard deviation with which the element shall be sampled
/// - `base`: the base in which the gadget vector is represented
///
/// Returns a preimage `x` satisfying `a*x = u` using the trapdoor.
///
/// # Examples
/// ```
/// use qfall_crypto::sample::distribution::discrete_gauss::sample_d_efficient;
/// use qfall_crypto::sample::g_trapdoor::gadget_parameters::GadgetParameters;
/// use qfall_crypto::sample::g_trapdoor::gen_trapdoor_default;
/// use qfall_math::{integer::Z, integer_mod_q::{Modulus, MatZq}, rational::Q};
///
/// let base = Z::from(2);
/// let s = Q::from(1);
/// let q = Modulus::try_from_z(&Z::from(1024)).unwrap();
/// let gadget_parameters = GadgetParameters::init_default(4, &q);
///
/// let (a, r) = gen_trapdoor_default(4, &q);
///
/// let tag = MatZq::identity(&gadget_parameters.n, &gadget_parameters.n, &q).unwrap();
///
/// let u = MatZq::sample_uniform(&gadget_parameters.n, 1, &q).unwrap();
///
/// let sample = sample_d_efficient(
///     4,
///     &a,
///     &r,
///     &u,
///     &tag,
///     (&gadget_parameters.k).try_into().unwrap(),
///     &s,
///     &base,
/// );
///
/// assert_eq!(u, a * MatZq::from((&sample, &q)));
/// ```
#[allow(clippy::too_many_arguments)]
pub fn sample_d_efficient(
    n: impl Into<Z>,
    a: &MatZq,
    trapdoor: &MatZ,
    u: &MatZq,
    tag: &MatZq,
    k: i64,
    s: &Q,
    base: &Z,
) -> MatZ {
    let n = n.into();
    // efficient sampleD samples with standard deviation `s \omega(\sqrt{\log n})`,
    // hence we have to rescale the provided standard deviation accordingly
    let _s = s / (n.log_ceil(&Z::from(2)).unwrap());

    // TODO: sample perturbation
    let p = compute_perturbation(trapdoor.get_num_rows() + trapdoor.get_num_columns());

    // TODO: invert tag, change this once is_identity is implemented
    if &MatZq::identity(tag.get_num_rows(), tag.get_num_columns(), tag.get_mod()).unwrap() != tag {
        panic!("Tags that are not the identity matrix are not supported at the moment.")
    }
    let tag_inv = tag;

    // compute new syndrome
    let v = tag_inv * (u - (a * MatZq::from((&p, &a.get_mod()))));

    let z: MatZ = if base.pow(k).unwrap() == a.get_mod().into() {
        // use bucket sampling for syndrome if `q = base^k

        // std_deviation = `base * \omega(\sqrt{\log n})`
        let std_deviation = base * n.log(&Z::from(2)).unwrap();
        let mut z = bucket_sampling(
            &v.get_entry(0, 0).unwrap(),
            &a.get_num_rows().into(),
            &std_deviation,
            k,
            base,
        );
        for i in 1..v.get_num_rows() {
            let v_i = v.get_entry(i, 0).unwrap();
            z = z
                .concat_vertical(&bucket_sampling(
                    &v_i,
                    &a.get_num_rows().into(),
                    &std_deviation,
                    k,
                    base,
                ))
                .unwrap();
        }
        z
    } else {
        // use (classical) sampleD for syndrome if `q != base ^k`
        todo!();
    };

    // shift to correct coset and add perturbation
    p + trapdoor
        .concat_vertical(
            &MatZ::identity(trapdoor.get_num_columns(), trapdoor.get_num_columns()).unwrap(),
        )
        .unwrap()
        * z
}

/// This function computes the non-spherical Gaussian perturbation which is added to the syndrome.
/// WIP(TODO): As this function is not yet implemented properly, this function is not
/// explained further for now.
///
/// Parameters:
/// - `m`: the size of the perturbation vector
///
/// Returns a perturbation vector of size `m`
fn compute_perturbation(m: i64) -> MatZ {
    MatZ::new(m, 1).unwrap()
}

/// Implements bucket-sampling according to [\[1\]](<../index.html#:~:text=[1]>).
/// Note: only works for `q`'s that are a power of `base`.
///
/// Parameters:
/// - `u`: the syndrome for which we want a vector `x`, such that `g^t * x = u`
/// - `n`: the security parameter
/// - `s`: the standard deviation with which the element is sampled
/// - `k`: the length of the gadget vector and the vector `x` we want to compute
/// - `base`: the amount of buckets and the base of the gadget vector
///
/// Returns a vector `x` which if multiplied by `g^t` maps to `u`.
///
/// # Examples
/// ```
/// use qfall_crypto::sample::distribution::discrete_gauss::bucket_sampling;
/// use qfall_crypto::sample::g_trapdoor::gadget_classical::gen_gadget_vec;
/// use qfall_math::{integer::Z, rational::Q, traits::GetEntry};
///
/// let base = Z::from(4);
/// let n = Z::from(5);
/// let k = 20;
/// let s = Q::from(2);
///
/// let u = Z::from(1234);
///
/// let sample = bucket_sampling(&u, &n, &s, k, &base);
/// let gadget_vec = gen_gadget_vec(k, &base).unwrap();
/// assert_eq!(
///     u,
///     (gadget_vec.transpose() * sample).get_entry(0, 0).unwrap()
/// );
/// ```
pub fn bucket_sampling(u: &Z, n: &Z, s: &Q, k: i64, base: &Z) -> MatZ {
    let mut buckets = vec![Vec::new(); i64::try_from(base).unwrap() as usize];
    let base_modulus = Modulus::try_from_z(base).unwrap();

    compute_out_vec(n, s, k, &base_modulus, &mut buckets, u)
}

/// Helper function for `bucket_sampling`.
/// This function samples a vector `out` such that
/// `g^t * out = u`.
///
/// Parameters:
/// - `n`: the security parameter
/// - `s`: the standard deviation with which the element is sampled
/// - `k`: the size of the vector
/// - `base_modulus`: the amount of buckets
/// - `buckets`: the buckets in which the samples are saved.
/// - `u`: the value which restraints the sample
///
/// Returns a vector such that `g^t * out = u`
fn compute_out_vec(
    n: &Z,
    s: &Q,
    k: i64,
    base_modulus: &Modulus,
    buckets: &mut [Vec<Z>],
    u: &Z,
) -> MatZ {
    let mut out = MatZ::new(k, 1).unwrap();
    let mut u = u.clone();
    for i in 0..k {
        let remainder = Zq::from_z_modulus(&u, base_modulus);
        let remainder_usize = i64::try_from(&remainder.get_value()).unwrap() as usize;
        fill_buckets(n, s, remainder_usize, base_modulus, buckets);
        let sample = buckets[remainder_usize].pop().unwrap();
        out.set_entry(i, 0, &sample).unwrap();

        u = (u - &sample).div_exact(&base_modulus.into()).unwrap();
    }
    out
}

/// Helper function for `compute_out_vec`.
/// This function checks for a certain remainder `s`, whether
/// there are elements in the related bucket, if not, it
/// samples until there is a sample in the correct bucket.
///
/// Parameters:
/// - `n`: the security parameter
/// - `s`: the standard deviation with which the element is sampled
/// - `matched_remainder`: the bucket for which we look for a sample
/// - `base_modulus`: the amount of buckets
/// - `buckets`: the buckets in which the samples are saved.
///
/// Ensures that the bucket indexed with `matched_remainder` contains a sample
fn fill_buckets(
    n: &Z,
    s: &Q,
    matched_remainder: usize,
    base_modulus: &Modulus,
    buckets: &mut [Vec<Z>],
) {
    // this function takes in the standard deviation, but calls `Z::sample_discrete_gauss`
    // hence we have to rescale the `s`
    let gaussian_param = s * (Q::from(2) * Q::PI).sqrt();
    while buckets[matched_remainder].is_empty() {
        let sample = Z::sample_discrete_gauss(n, &Z::ZERO, &gaussian_param).unwrap();
        let remainder = Zq::from_z_modulus(&sample, base_modulus);
        let remainder_usize = i64::try_from(&remainder.get_value()).unwrap() as usize;
        buckets[remainder_usize].push(sample);
    }
}

#[cfg(test)]
mod test_sample_d_efficient {
    use super::sample_d_efficient;
    use crate::sample::g_trapdoor::{
        gadget_classical::gen_trapdoor, gadget_parameters::GadgetParameters, gen_trapdoor_default,
    };
    use qfall_math::{
        integer::Z,
        integer_mod_q::{MatZq, Modulus},
        rational::Q,
        traits::SetEntry,
    };

    /// returns a triangular matrix with `1` on the diagonal which is thereby invertible.
    fn generate_invertible_tag(n: u64, q: &Modulus) -> MatZq {
        let mut out = MatZq::identity(n, n, q).unwrap();

        for i in 0..n {
            for j in i + 1..n {
                out.set_entry(i, j, Z::sample_uniform(&0, &q).unwrap())
                    .unwrap()
            }
        }
        out
    }

    /// ensure that the sampled vector is actually a preimage of `f_a`, i.e.
    /// `Ax = u` for q being a power of two
    #[test]
    fn working_with_tag_identity_q_power_two() {
        let base = Z::from(2);
        let s = Q::from(10);
        let q = Modulus::try_from_z(&Z::from(16777216)).unwrap();
        let gadget_parameters = GadgetParameters::init_default(3, &q);

        let (a, r) = gen_trapdoor_default(3, &q);

        let tag = MatZq::identity(&gadget_parameters.n, &gadget_parameters.n, &q).unwrap();

        let u = MatZq::sample_uniform(&gadget_parameters.n, 1, &q).unwrap();

        let sample = sample_d_efficient(
            3,
            &a,
            &r,
            &u,
            &tag,
            (&gadget_parameters.k).try_into().unwrap(),
            &s,
            &base,
        );

        assert_eq!(u, &a * MatZq::from((&sample, &q)));
    }

    /// ensure that the sampled vector is actually a preimage of `f_a`, i.e.
    /// `Ax = u` for `q` not being a power of two
    #[test]
    #[ignore = "WIP/TODO: SampleD not yet implemented"]
    fn working_with_tag_identity_q_not_power_two() {
        let base = Z::from(2);
        let s = Q::from(1);
        let q = Modulus::try_from_z(&Z::from(16777215)).unwrap();
        let gadget_parameters = GadgetParameters::init_default(3, &q);

        let (a, r) = gen_trapdoor_default(3, &q);

        let tag = MatZq::identity(&gadget_parameters.n, &gadget_parameters.n, &q).unwrap();

        let u = MatZq::sample_uniform(&gadget_parameters.n, 1, &q).unwrap();

        let sample = sample_d_efficient(
            3,
            &a,
            &r,
            &u,
            &tag,
            (&gadget_parameters.k).try_into().unwrap(),
            &s,
            &base,
        );

        assert_eq!(u, a * MatZq::from((&sample, &q)));
    }

    /// ensure that the sampled vector is actually a preimage of `f_a`, i.e.
    /// `Ax = u` with an invertible tag and `q` being a power of two
    #[test]
    #[ignore = "WIP/TODO: tag_invert not yet implemented"]
    fn working_with_tag_q_power_two() {
        let base = Z::from(2);
        let s = Q::from(1);
        let q = Modulus::try_from_z(&Z::from(2048)).unwrap();

        let gadget_parameters = GadgetParameters::init_default(7, &q);
        let tag = generate_invertible_tag(7, &q);
        let a_bar = MatZq::sample_uniform(
            &gadget_parameters.n,
            &gadget_parameters.m_bar,
            &gadget_parameters.q,
        )
        .unwrap();

        let (a, r) = gen_trapdoor(&gadget_parameters, &a_bar, &tag).unwrap();

        let u = MatZq::sample_uniform(&gadget_parameters.n, 1, &q).unwrap();

        let sample = sample_d_efficient(
            3,
            &a,
            &r,
            &u,
            &tag,
            (&gadget_parameters.k).try_into().unwrap(),
            &s,
            &base,
        );

        assert_eq!(u, a * MatZq::from((&sample, &q)));
    }

    /// ensure that the sampled vector is actually a preimage of `f_a`, i.e.
    /// `Ax = u` with an invertible tag and `q` not being a power of two
    #[test]
    #[ignore = "WIP/TODO: tag_invert not yet implemented, SampleD not yet implemented"]
    fn working_with_tag_q_not_power_two() {
        let base = Z::from(2);
        let s = Q::from(1);
        let q = Modulus::try_from_z(&Z::from(12345)).unwrap();

        let gadget_parameters = GadgetParameters::init_default(7, &q);
        let tag = generate_invertible_tag(7, &q);
        let a_bar = MatZq::sample_uniform(
            &gadget_parameters.n,
            &gadget_parameters.m_bar,
            &gadget_parameters.q,
        )
        .unwrap();

        let (a, r) = gen_trapdoor(&gadget_parameters, &a_bar, &tag).unwrap();

        let u = MatZq::sample_uniform(&gadget_parameters.n, 1, &q).unwrap();

        let sample = sample_d_efficient(
            3,
            &a,
            &r,
            &u,
            &tag,
            (&gadget_parameters.k).try_into().unwrap(),
            &s,
            &base,
        );

        assert_eq!(u, a * MatZq::from((&sample, &q)));
    }
}

/// `compute_out_vec` is implicitly tested when testing `bucket_sampling`
#[cfg(test)]
mod test_bucket_sampling {
    use super::bucket_sampling;
    use crate::sample::g_trapdoor::gadget_classical::gen_gadget_vec;
    use qfall_math::{integer::Z, rational::Q, traits::GetEntry};

    /// ensure that the sampled values actually map to `u` if multiplied with `g^t`
    #[test]
    fn working() {
        let base = Z::from(4);
        let n = Z::from(5);
        let k = 20;
        let s = Q::from(1);

        for u_ in [0, 1, 2, 10, 1141, 1241531, 15419513] {
            let u = Z::from(u_);

            let sample = bucket_sampling(&u, &n, &s, k, &base);
            let gadget_vec = gen_gadget_vec(k, &base).unwrap();

            assert_eq!(
                u,
                (gadget_vec.transpose() * sample).get_entry(0, 0).unwrap()
            );
        }
    }
}

#[cfg(test)]
mod test_fill_buckets {
    use super::fill_buckets;
    use qfall_math::{integer::Z, integer_mod_q::Modulus, rational::Q};

    /// ensure that the buckets for which we call the function actually contain an
    /// element after calling the function
    #[test]
    fn working() {
        let base = Z::from(5);
        let base_modulus = Modulus::try_from_z(&base).unwrap();
        let n = Z::from(6);
        let s = Q::from(3);

        for i in 0..5 {
            let mut buckets = vec![Vec::new(); i64::try_from(&base).unwrap() as usize];
            fill_buckets(&n, &s, i, &base_modulus, &mut buckets);

            assert!(!buckets[i].is_empty())
        }
    }
}
