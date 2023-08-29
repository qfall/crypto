// Copyright © 2023 Marvin Beckmann
//
// This file is part of qFALL-crypto.
//
// qFALL-crypto is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

//! Allows to Deserialize an arbitrary [`Fdh`] instantiation

use crate::{primitive::hash::HashInto, sample::distribution::psf::PSF};
use serde::{
    de::{Error, MapAccess, Visitor},
    Deserialize, Serialize,
};
use std::{fmt, marker::PhantomData};

use super::Fdh;
impl<'de, A, Trapdoor, Domain, Range, T, Hash> Deserialize<'de>
    for Fdh<A, Trapdoor, Domain, Range, T, Hash>
where
    Domain: Serialize + for<'a> Deserialize<'a>,
    T: PSF<A, Trapdoor, Domain, Range> + Serialize + for<'a> Deserialize<'a>,
    Hash: HashInto<Range> + Serialize + for<'a> Deserialize<'a>,
{
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: serde::Deserializer<'de>,
    {
        /// This enum defines the content of the struct to be generated using [`Deserialize`]
        const FIELDS: &[&str] = &["psf", "storage", "hash"];
        #[derive(Deserialize)]
        #[serde(field_identifier, rename_all = "lowercase")]
        enum Field {
            Psf,
            Storage,
            Hash,
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
            type Value = Fdh<A, Trapdoor, Domain, Range, T, Hash>;

            fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
                formatter.write_str("struct $type")
            }

            fn visit_map<V>(self, mut map: V) -> Result<Self::Value, V::Error>
            where
                V: MapAccess<'de>,
            {
                let mut psf = None;
                let mut storage = None;
                let mut hash = None;
                while let Some(key) = map.next_key()? {
                    match key {
                        Field::Psf => {
                            if psf.is_some() {
                                return Err(Error::duplicate_field("psf"));
                            }
                            psf = Some(map.next_value()?);
                        }
                        Field::Storage => {
                            if storage.is_some() {
                                return Err(Error::duplicate_field("storage"));
                            }
                            storage = Some(map.next_value()?);
                        }
                        Field::Hash => {
                            if hash.is_some() {
                                return Err(Error::duplicate_field("hash"));
                            }
                            hash = Some(map.next_value()?);
                        }
                    }
                }

                Ok(Fdh {
                    psf: Box::new(psf.unwrap()),
                    storage: storage.unwrap(),
                    hash: Box::new(hash.unwrap()),
                    _a_type: PhantomData,
                    _trapdoor_type: PhantomData,
                    _range_type: PhantomData,
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
        deserializer.deserialize_struct("Fdh", FIELDS, struct_visitor)
    }
}

#[cfg(test)]
mod test_deserialization {
    use crate::{
        construction::signature::{fdh::Fdh, SignatureScheme},
        primitive::hash::HashMatZq,
        sample::distribution::psf::gpv::PSFGPV,
    };
    use qfall_math::{integer::MatZ, integer_mod_q::MatZq, rational::MatQ};

    /// Ensure that deserialization works.
    #[allow(clippy::type_complexity)]
    #[test]
    fn deserialize_gpv() {
        let mut fdh = Fdh::init_gpv(2, 127, 20);

        // fill one entry in the HashMap
        let m = "Hello World!";
        let (pk, sk) = fdh.gen();
        let _ = fdh.sign(m.to_owned(), &sk, &pk);

        let fdh_string = serde_json::to_string(&fdh).expect("Unable to create a json object");
        let fdh_2: Result<Fdh<MatZq, (MatZ, MatQ), MatZ, MatZq, PSFGPV, HashMatZq>, _> =
            serde_json::from_str(&fdh_string);

        assert!(fdh_2.is_ok());
    }
}
