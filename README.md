# qFALL-crypto
[![made-with-rust](https://img.shields.io/badge/Made%20with-Rust-1f425f.svg)](https://www.rust-lang.org/)
[![CI](https://github.com/qfall/crypto/actions/workflows/push.yml/badge.svg?branch=dev)](https://github.com/qfall/crypto/actions/workflows/pull_request.yml)
[![License: MPL 2.0](https://img.shields.io/badge/License-MPL_2.0-brightgreen.svg)](https://opensource.org/licenses/MPL-2.0)

This repository is currently being developed by the project group [qFALL - quantum resistant fast lattice library](https://cs.uni-paderborn.de/cuk/lehre/veranstaltungen/ws-2022-23/project-group-qfall) in the winter term 2022 and summer term 2023 by the Codes and Cryptography research group in Paderborn.

The main objective of this project is to provide researchers and students with the possibility to easily and quickly prototype (lattice-based) cryptography.

## Disclaimer
Currently, we are in the development phase and interfaces might change.
Feel free to check out the current progress, but be aware, that the content will
change in the upcoming weeks and months. An official release will most likely be published in the second half of 2023.

## Quick-Start

Please refer to [our website](https://qfall.github.io/) as central information point.

To install and add our library to your project, please refer to [our tutorial](https://qfall.github.io/book/index.html).
It provides a step-by-step guide to install the required libraries and gives further insights in the usage of our crates.

## What does qFALL-crypto offer?

qFALL-crypto offers a variety of implementations of cryptographic schemes, constructions, and primitives.
We provide a brief overview in the following list.
For a more detailed description, please refer to [our tutorial section](https://qfall.github.io/book/crypto/features.html).

Full-fledged Cryptographic Features
- [Public Key Encryption](https://github.com/qfall/crypto/blob/dev/src/construction/pk_encryption.rs)
    - [LWE Encryption](https://github.com/qfall/crypto/blob/dev/src/construction/pk_encryption/regev.rs)
    - [Dual LWE Encryption](https://github.com/qfall/crypto/blob/dev/src/construction/pk_encryption/dual_regev.rs)
    - [LPR Encryption](https://github.com/qfall/crypto/blob/dev/src/construction/pk_encryption/lpr.rs)
    - [Ring-based LPR Encryption](https://github.com/qfall/crypto/blob/dev/src/construction/pk_encryption/ring_lpr.rs)
    - [CCA-secure Encryption](https://github.com/qfall/crypto/blob/dev/src/construction/pk_encryption/ccs_from_ibe.rs)
- [Signatures](https://github.com/qfall/crypto/blob/dev/src/construction/signature.rs)
    - [Full-Domain Hash (FDH)](https://github.com/qfall/crypto/blob/dev/src/construction/signature/fdh.rs)
    - [Probabilistic FDH (PFDH)](https://github.com/qfall/crypto/blob/dev/src/construction/signature/pfdh.rs)
    - [Ring-based FDH](https://github.com/qfall/crypto/blob/dev/src/construction/signature/fdh/gpv_ring.rs)
- [Identity Based Encryption](https://github.com/qfall/crypto/blob/dev/src/construction/identity_based_encryption.rs)
    - [From Dual LWE Encryption](https://github.com/qfall/crypto/blob/dev/src/construction/identity_based_encryption/dual_regev_ibe.rs)
- [Hash Functions](https://github.com/qfall/crypto/blob/dev/src/construction/hash.rs)
    - [SIS-Hash Function](https://github.com/qfall/crypto/blob/dev/src/construction/hash/sis.rs)
    - [SHA-256-based Hash](https://github.com/qfall/crypto/blob/dev/src/construction/hash/sha256.rs)

Building Blocks and Primitives
- [Preimage Samplable Functions (PSF)](https://github.com/qfall/crypto/blob/dev/src/primitive/psf.rs)
- [Trapdoors](https://github.com/qfall/crypto/blob/dev/src/sample/g_trapdoor.rs)
    - [G-trapdoor incl. short basis](https://github.com/qfall/crypto/blob/dev/src/sample/g_trapdoor/gadget_classical.rs)
    - [Ring-based G-trapdoor incl. short basis](https://github.com/qfall/crypto/blob/dev/src/sample/g_trapdoor/gadget_ring.rs)

## License
This library is distributed under the **Mozilla Public License Version 2.0** which can be found here [License](https://github.com/qfall/crypto/blob/dev/LICENSE).
Permissions of this weak copyleft license are conditioned on making available source code of licensed files and modifications of those files under the same license (or in certain cases, one of the GNU licenses). Copyright and license notices must be preserved. Contributors provide an express grant of patent rights. However, a larger work using the licensed work may be distributed under different terms and without source code for files added in the larger work.

## Citing

Please use the following bibtex entry to cite [qFALL-crypto](https://github.com/qfall/crypto):

```text
@misc{qFALL-crypto,
    author = {Porzenheim, Laurens and Beckmann, Marvin and Kramer, Paul and Milewski, Phil and Moog, Sven and Schmidt, Marcel and Siemer, Niklas}
    title = {qFALL-crypto v0.0},
    howpublished = {Online: \url{https://github.com/qfall/crypto}},
    month = Mar,
    year = 2023,
    note = {University Paderborn,  Codes and Cryptography}
}
```

## Get in Touch
One can contact the members of the project group with our mailing list `pg-qfall(at)lists.upb.de`.
