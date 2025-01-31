# nf-core/metatdenovo: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.1.0dev - [date]

Initial release of nf-core/metatdenovo, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- [#320](<[https://github.com/nf-core/metatdenovo/pull/320](https://github.com/nf-core/metatdenovo/pull/320)>) improvments to Diamond taxonomy plus documentation
- [#312](<[https://github.com/nf-core/metatdenovo/pull/312](https://github.com/nf-core/metatdenovo/pull/312)>) added taxonomy directly with Diamond, see `--diamond_dbs`
- [#286](<[https://github.com/nf-core/metatdenovo/pull/286](https://github.com/nf-core/metatdenovo/pull/286)>) added an option to save the fasta file output from formatspades.nf module
- [#285](<[https://github.com/nf-core/metatdenovo/pull/285](https://github.com/nf-core/metatdenovo/pull/285)>) added nf-test for default settings.
- [#280](<[https://github.com/nf-core/metatdenovo/issues/280](https://github.com/nf-core/metatdenovo/issues/280)>) - Added minid option to bbmap_align module. Now the threshold for mapping a read to a contig is an identity of 0.9. The previous version of nf-core/metatdenovo used the default for BBMap, 0.76. This version might hence give slightly different results than the previous.
- [#271](<[https://github.com/nf-core/metatdenovo/issues/271](https://github.com/nf-core/metatdenovo/issues/271)>) - Added flavor to SPADES modules

### `Changed`

- [#311](<[https://github.com/nf-core/metatdenovo/pull/311](https://github.com/nf-core/metatdenovo/pull/311)>) - update modules and subworkflows
- [#295](<[https://github.com/nf-core/metatdenovo/pull/295](https://github.com/nf-core/metatdenovo/pull/295)>) - Update documentation
- [#292](<[https://github.com/nf-core/metatdenovo/pull/292](https://github.com/nf-core/metatdenovo/pull/292)>) - Specify memory to Megahit process
- [#290](<[https://github.com/nf-core/metatdenovo/pull/290](https://github.com/nf-core/metatdenovo/pull/290)>) - Template update to v2.14.1
- [#283](<[https://github.com/nf-core/metatdenovo/pull/283](https://github.com/nf-core/metatdenovo/pull/283)>) - Updated documentation about download databases manually
- [#268](<[https://github.com/nf-core/metatdenovo/pull/268](https://github.com/nf-core/metatdenovo/pull/268)>) - Don't save so many intermediate Megahit files by default

### `Fixed`

- [#305](<[https://github.com/nf-core/ampliseq/pull/681](https://github.com/nf-core/metatdenovo/pull/305)>) - Make EUKulele counts output optional as it's not always created
- [#269](<[https://github.com/nf-core/ampliseq/pull/681](https://github.com/nf-core/metatdenovo/pull/269)>) - Make Transdecoder work better with `-resume`

### `Dependencies`

### `Deprecated`

## v1.0.1 - [2024-04-02]

### `Fixed`

- [#277](https://github.com/nf-core/metatdenovo/pull/277) - Fix mistake in how `--eukulele_db` parameter is handled. Remove possibility to use a list of dbs in the same run.
- [#277](https://github.com/nf-core/metatdenovo/pull/277) - Gzip user provided assembly files to avoid overwriting by assuming they're already zipped.

## v1.0.0 - [2024-02-15]

Initial release of nf-core/metatdenovo, created with the [nf-core](https://nf-co.re/) template.
