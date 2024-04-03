# nf-core/metatdenovo: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.1.0 - [2024-02-15]

### `Added`
- [#271](<[https://github.com/nf-core/metatdenovo/issues/271](https://github.com/nf-core/metatdenovo/issues/271)>) - Added flavor to SPADES modules

### `Changed`

- [#268](<[https://github.com/nf-core/ampliseq/pull/681](https://github.com/nf-core/metatdenovo/pull/268)>) - Don't save so many intermediate Megahit files by default

### `Fixed`

- [#269](<[https://github.com/nf-core/ampliseq/pull/681](https://github.com/nf-core/metatdenovo/pull/269)>) - Make Transdecoder work better with `-resume`

### `Dependencies`

### `Deprecated`

## v1.0.1 - [2024-04-02]

### `Fixed`

- [#277](https://github.com/nf-core/metatdenovo/pull/277) - Fix mistake in how `--eukulele_db` parameter is handled. Remove possibility to use a list of dbs in the same run.
- [#277](https://github.com/nf-core/metatdenovo/pull/277) - Gzip user provided assembly files to avoid overwriting by assuming they're already zipped.

## v1.0.0 - [2024-02-15]

Initial release of nf-core/metatdenovo, created with the [nf-core](https://nf-co.re/) template.
