// Copyright © 2023 Niklas Siemer
//
// This file is part of qFALL-crypto.
//
// qFALL-crypto is free software: you can redistribute it and/or modify it under
// the terms of the Mozilla Public License Version 2.0 as published by the
// Mozilla Foundation. See <https://mozilla.org/en-US/MPL/2.0/>.

//! This module provides the trait a struct should implement if it is an
//! instance of a public key encryption scheme. Furthermore, it contains
//! cryptographic schemes implementing the `PKEncryption` trait.

mod dual_regev;

// This trait is just an example implementation. Please think it further through
// and adjust it before the first implementation of a pk encryption scheme
pub trait PKEncryption<PublicParameters, SecurityParameter, PublicKey, SecretKey, Message, Cipher> {
    fn gen(&self, security_parameters: SecurityParameter) -> (PublicKey, SecretKey);
    fn enc(&self, pk: PublicKey, message: Message) -> Cipher;
    fn dec(&self, sk: SecretKey, cipher: Cipher) -> Message;
}
