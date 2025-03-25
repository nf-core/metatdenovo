# nf-core/metatdenovo: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.1.3 - [date]

### `Added`

### `Changed`

-  - Updated some objects in nextflow_schema.json (@herich0).

### `Fixed`

### `Dependencies`

### `Deprecated`

## v1.1.1 - [2025-03-13]

## v1.1.2 - [date]

### `Added`

### `Changed`

- [#352](https://github.com/nf-core/metatdenovo/pull/352) - removed "removed":"" from objects in nextflow_schema.json as it is not necessary (@erikrikarddaniel).

### `Fixed`

- [#352](https://github.com/nf-core/metatdenovo/pull/352) - Assign less memory to BBNorm to avoid getting killed (@erikrikarddaniel).

### `Dependencies`

### `Deprecated`

## v1.1.1 - [2025-03-13]

### `Added`

### `Changed`

- [#352](https://github.com/nf-core/metatdenovo/pull/352) - Updated some modules (@erikrikarddaniel).

### `Fixed`

- [#352](https://github.com/nf-core/metatdenovo/pull/352) - Assign less memory to BBNorm to avoid getting killed (@erikrikarddaniel).

### `Dependencies`

### `Deprecated`

## v1.1.0 - [2025-02-25]

### `Added`

- [#331](https://github.com/nf-core/metatdenovo/pull/331) - Added nf-tests.
- [#320](https://github.com/nf-core/metatdenovo/pull/320) - added taxonomy directly with Diamond, part 2
- [#312](https://github.com/nf-core/metatdenovo/pull/312) - added taxonomy directly with Diamond, see `--diamond_dbs`.
- [#286](https://github.com/nf-core/metatdenovo/pull/286) - added an option to save the fasta file output from formatspades.nf module.
- [#285](https://github.com/nf-core/metatdenovo/pull/285) - added nf-test for default settings.
- [#280](https://github.com/nf-core/metatdenovo/issues/280) - Added minid option to bbmap_align module. Now the threshold for mapping a read to a contig is an identity of 0.9. The previous version of nf-core/metatdenovo used the default for BBMap, 0.76. This version might hence give slightly different results than the previous.
- [#271](https://github.com/nf-core/metatdenovo/issues/271) - Added flavor to SPADES modules.

### `Changed`

- [#332](https://github.com/nf-core/metatdenovo/pull/332) - Rearranged tree structure for local modules and local subworkflows.
- [#330](https://github.com/nf-core/metatdenovo/pull/330) - Update Usage.md.
- [#326](https://github.com/nf-core/metatdenovo/pull/326) - Clean up overall stats table.
- [#323](https://github.com/nf-core/metatdenovo/pull/323) - Modified param names for input of assembly and ORFs; added name params for output file naming.
- [#323](https://github.com/nf-core/metatdenovo/pull/323) - Removed default for `assembler` and `orf_caller` parameters.
- [#318](https://github.com/nf-core/metatdenovo/pull/318) - Template 3.2.0 update.
- [#311](https://github.com/nf-core/metatdenovo/pull/311) - Update modules and subworkflows.
- [#295](https://github.com/nf-core/metatdenovo/pull/295) - Update documentation.
- [#292](https://github.com/nf-core/metatdenovo/pull/292) - Specify memory to Megahit process.
- [#290](https://github.com/nf-core/metatdenovo/pull/290) - Template update to v2.14.1.
- [#283](https://github.com/nf-core/metatdenovo/pull/283) - Updated documentation about download databases manually.
- [#268](https://github.com/nf-core/metatdenovo/pull/268) - Don't save so many intermediate Megahit files by default.

### `Fixed`

- [#328](https://github.com/nf-core/metatdenovo/pull/328) - Fix BBDuk was passing only one sample.
- [#326](https://github.com/nf-core/metatdenovo/pull/326) - Fix resources for test cases.
- [#326](https://github.com/nf-core/metatdenovo/pull/326) - Fix output file names for Eukulele and Kofamscan.
- [#321](https://github.com/nf-core/metatdenovo/pull/321) - Fix how params.sequence_filter was called in BBDuk module.
- [#305](https://github.com/nf-core/metatdenovo/pull/305) - Make EUKulele counts output optional as it's not always created.
- [#269](https://github.com/nf-core/metatdenovo/pull/269) - Make Transdecoder work better with `-resume`.

### `Dependencies`

### `Deprecated`

## v1.0.1 - [2024-04-02]

### `Fixed`

- [#277](https://github.com/nf-core/metatdenovo/pull/277) - Fix mistake in how `--eukulele_db` parameter is handled. Remove possibility to use a list of dbs in the same run.
- [#277](https://github.com/nf-core/metatdenovo/pull/277) - Gzip user provided assembly files to avoid overwriting by assuming they're already zipped.

## v1.0.0 - [2024-02-15]

Initial release of nf-core/metatdenovo, created with the [nf-core](https://nf-co.re/) template.
