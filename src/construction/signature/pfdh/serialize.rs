// Copyright © 2023 Phil Milewski
//
// This file is part of qFALL-crypto.
//
// qFALL-crypto is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

//! Allows to Deserialize an arbitrary [`PFDH`] instantiation

use super::PFDH;
use crate::{construction::hash::HashInto, primitive::psf::PSF};
use serde::{
    de::{Error, MapAccess, Visitor},
    Deserialize, Serialize,
};
use std::{fmt, marker::PhantomData};

impl<'de, A, Trapdoor, Domain, Range, T, Hash> Deserialize<'de>
    for PFDH<A, Trapdoor, Domain, Range, T, Hash>
where
    Domain: Serialize + for<'a> Deserialize<'a>,
    T: PSF<A, Trapdoor, Domain, Range> + Serialize + for<'a> Deserialize<'a>,
    Hash: HashInto<Range> + Serialize + for<'a> Deserialize<'a>,
{
    #[allow(non_camel_case_types)]
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        /// This enum defines the content of the struct to be generated using [`Deserialize`]
        const FIELDS: &[&str] = &["psf", "hash", "randomness_length"];
        #[derive(Deserialize)]
        #[serde(field_identifier, rename_all = "lowercase")]
        enum Field {
            Psf,
            Hash,
            Randomness_Length,
        }

        /// This visitor iterates over the strings content and collects all possible fields.
        /// It sets the corresponding values of the struct based on the values found.
        struct StructVisitor<A, Trapdoor, Domain, Range, T, Hash> {
            a: PhantomData<A>,
            trapdoor: PhantomData<Trapdoor>,
            domain: PhantomData<Domain>,
            range: PhantomData<Range>,
            t: PhantomData<T>,
            hash: PhantomData<Hash>,
        }
        impl<'de, A, Trapdoor, Domain, Range, T, Hash> Visitor<'de>
            for StructVisitor<A, Trapdoor, Domain, Range, T, Hash>
        where
            Domain: Serialize + for<'a> Deserialize<'a>,
            T: PSF<A, Trapdoor, Domain, Range> + Serialize + for<'a> Deserialize<'a>,
            Hash: HashInto<Range> + Serialize + for<'a> Deserialize<'a>,
        {
            type Value = PFDH<A, Trapdoor, Domain, Range, T, Hash>;

            fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
                formatter.write_str("struct $type")
            }

            fn visit_map<V>(self, mut map: V) -> Result<Self::Value, V::Error>
            where
                V: MapAccess<'de>,
            {
                let mut psf = None;
                let mut hash = None;
                let mut randomness_length = None;
                while let Some(key) = map.next_key()? {
                    match key {
                        Field::Psf => {
                            if psf.is_some() {
                                return Err(Error::duplicate_field("psf"));
                            }
                            psf = Some(map.next_value()?);
                        }
                        Field::Hash => {
                            if hash.is_some() {
                                return Err(Error::duplicate_field("hash"));
                            }
                            hash = Some(map.next_value()?);
                        }
                        Field::Randomness_Length => {
                            if randomness_length.is_some() {
                                return Err(Error::duplicate_field("randomness_length"));
                            }
                            randomness_length = Some(map.next_value()?);
                        }
                    }
                }

                Ok(PFDH {
                    psf: Box::new(psf.unwrap()),
                    hash: Box::new(hash.unwrap()),
                    randomness_length: randomness_length.unwrap(),
                    _a_type: PhantomData,
                    _trapdoor_type: PhantomData,
                    _range_type: PhantomData,
                    _domain_type: PhantomData,
                })
            }
        }

        let struct_visitor: StructVisitor<A, Trapdoor, Domain, Range, T, Hash> = StructVisitor {
            a: PhantomData,
            trapdoor: PhantomData,
            domain: PhantomData,
            range: PhantomData,
            t: PhantomData,
            hash: PhantomData,
        };
        deserializer.deserialize_struct("PFDH", FIELDS, struct_visitor)
    }
}

#[cfg(test)]
mod test_deserialization {
    use crate::{
        construction::{
            hash::sha256::HashMatZq,
            signature::{SignatureScheme, PFDH},
        },
        primitive::psf::PSFGPV,
    };
    use qfall_math::{integer::MatZ, integer_mod_q::MatZq, rational::MatQ};

    /// Ensure that deserialization works.
    #[allow(clippy::type_complexity)]
    #[test]
    fn deserialize_gpv() {
        let mut pfdh = PFDH::init_gpv(2, 127, 20, 1233);

        let m = "Hello World!";
        let (pk, sk) = pfdh.gen();
        let signature = pfdh.sign(m.to_owned(), &sk, &pk);

        let pfdh_string = serde_json::to_string(&pfdh).expect("Unable to create a json object");
        let pfdh_2: Result<PFDH<MatZq, (MatZ, MatQ), MatZ, MatZq, PSFGPV, HashMatZq>, _> =
            serde_json::from_str(&pfdh_string);

        assert!(pfdh_2.is_ok());

        //ensure signing still works
        let mut pfdh_2 = pfdh_2.unwrap();
        let signature_2 = pfdh_2.sign(m.to_owned(), &sk, &pk);

        //ensure verification still works
        assert!(pfdh_2.vfy(m.to_string(), &signature, &pk));
        assert!(pfdh_2.vfy(m.to_string(), &signature_2, &pk));
    }
}
