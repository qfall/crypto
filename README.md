# qFALL-crypto
[![made-with-rust](https://img.shields.io/badge/Made%20with-Rust-1f425f.svg)](https://www.rust-lang.org/)
[![CI](https://github.com/qfall/crypto/actions/workflows/push.yml/badge.svg?branch=dev)](https://github.com/qfall/crypto/actions/workflows/pull_request.yml)
[![License: MPL 2.0](https://img.shields.io/badge/License-MPL_2.0-brightgreen.svg)](https://opensource.org/licenses/MPL-2.0)
<TODO Badge for Documentation>
<TODO Badge for Code Coverage>
<TODO Badge for Website>

This repository is currently being developed by the project group [qFALL - quantum resistant fast lattice library](https://cs.uni-paderborn.de/cuk/lehre/veranstaltungen/ws-2022-23/project-group-qfall) in the winter term 2022 and summer term 2023 by the Codes and Cryptography research group in Paderborn.

The main objective of this project is to provide researchers and students with the possibility to easily and quickly prototype (lattice-based) cryptography.

## Disclaimer
Currently, we are in the development phase and interfaces might change.
Feel free to check out the current progress, but be aware, that the content will
change in the upcoming weeks and months. An official release will be published in the second half of 2023.

## Installation
In order to use this project one needs to have an [installation of Rust](https://www.rust-lang.org/tools/install). Since the underlying [math](https://github.com/qfall/math) crate internally uses [flint-sys](https://crates.io/crates/flint-sys)
which itself uses [gmp](https://gmplib.org/manual/), we are currently restricted to usage on Mac, Linux and Linux subsystems under Windows. For a subsystem under Windows, one additionally is required to have installed m4 and a C-compiler.

Since our project isn't yet published there is no option to find it on Rust's library collection on [crates.io](https://crates.io/).
If you want to include this project in your own Rust project, you can 
include a link to our version on the `dev` branch in your `Cargo.toml`. 

```text
qfall_crypto = { git = "https://github.com/qfall/crypto", branch="dev" }
```

Be aware that the external libraries in our project have to be compiled at the first installation,
which may take about 30 minutes. After the first installation it should be working fine.

<!-- TODO ## What does qFALL-crypto offer? -->

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
