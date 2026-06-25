# nf-core/metatdenovo: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.3.1 - [YYYY-mm-dd]

### `Added`

- [#439](https://github.com/nf-core/metatdenovo/pull/439) - Add paper citation(@erikrikarddaniel)

### `Changed`

- [#442](https://github.com/nf-core/metatdenovo/pull/442) - Increase default memory for Megahit process (@erikrikarddaniel)
- [#440](https://github.com/nf-core/metatdenovo/pull/440) - Remove "uncl." from taxon names in EUKulele `summary_tables` output (@erikrikarddaniel)
- [#440](https://github.com/nf-core/metatdenovo/pull/440) - Make R-package versions specific and move containers to Seqera-hosted (@erikrikarddaniel)
- [#439](https://github.com/nf-core/metatdenovo/pull/439) - Update more software versions (@erikrikarddaniel)
- [#439](https://github.com/nf-core/metatdenovo/pull/439) - Move pipeline to topic channels for versions and better syntax compliance (almost "strict") (@erikrikarddaniel)
- [#435](https://github.com/nf-core/metatdenovo/pull/435) - Template update 4.0.2 (@danilodileo)
- [#430](https://github.com/nf-core/metatdenovo/pull/430) - Nextflow lint (@danilodileo)
- [#429](https://github.com/nf-core/metatdenovo/pull/429) - Module update to nf-core tools 3.5.2 (@danilodileo)
- [#428](https://github.com/nf-core/metatdenovo/pull/428) - Template update to nf-core tools 3.5.2 (@danilodileo)
- [#416](https://github.com/nf-core/metatdenovo/pull/416) - Better content pipeline integration tests (@danilodileo)

### `Fixed`

- [#440](https://github.com/nf-core/metatdenovo/pull/440) - Add database name to eukulele process labels, closes issue [#417](https://github.com/nf-core/metatdenovo/issues/417) (@erikrikarddaniel)
- [#440](https://github.com/nf-core/metatdenovo/pull/440) - Remove "cds." from transdecoder orf names in counts summary table, closes issue [#418](https://github.com/nf-core/metatdenovo/issues/418) (@erikrikarddaniel)
- [#440](https://github.com/nf-core/metatdenovo/pull/440) - Make sure FastQC output is included in the MultiQC report, closes issue [#422](https://github.com/nf-core/metatdenovo/issues/422) (@erikrikarddaniel)
- [#440](https://github.com/nf-core/metatdenovo/pull/440) - Improve documentation of input samplesheet fields (@erikrikarddaniel)
- [#439](https://github.com/nf-core/metatdenovo/pull/439) - Fix download of eggnog database as mentioned in [#423](https://github.com/nf-core/metatdenovo/issues/423) (@erikrikarddaniel)
- [#439](https://github.com/nf-core/metatdenovo/pull/439) - Remove dependency of `versions.yml` presence for eggnog and kofamscan databases (@erikrikarddaniel)

### `Dependencies`

| Tool         | Previous version | New version |
| ------------ | ---------------- | ----------- |
| cat          | 2.3.4            | 2.8         |
| samtools     | 1.21             | 1.23        |
| subread      | 2.0.6            | 2.1.1       |
| trim-galore  | 0.6.10           | 2.1.0       |
| r-base       |                  | 4.5.3       |
| r-dplyr      |                  | 1.2.1       |
| r-readr      |                  | 2.2.0       |
| r-purrr      |                  | 1.2.2       |
| r-tidyr      |                  | 1.3.2       |
| r-stringi    |                  | 1.8.7       |
| r-stringr    |                  | 1.6.0       |
| r-data.table | 1.14.8           | 1.17.8      |
| r-dtplyr     | 1.3.1            | 1.3.3       |

(R packages without previous versions above were used but did not have specified versions as they were used as dependencies of r-tidyverse 2.0.0 which led to drifts in versions.)

### `Deprecated`

## v1.3.0 - [2025-08-29]

### `Added`

### `Changed`

- [#406](https://github.com/nf-core/metatdenovo/pull/406) - Updating modules and removing warnings before release 1.3.0 (@danilodileo)
- [#405](https://github.com/nf-core/metatdenovo/pull/405) - Upgrade EUKulele to 2.1.2. This appears to fix problems with downloads of certain databases (@erikrikarddaniel)
- [#404](https://github.com/nf-core/metatdenovo/pull/404) - Added new reference for SPAdes in CITATIONS.md (@danilodileo)
- [#394](https://github.com/nf-core/metatdenovo/pull/394) - allow unzipped input files (@erikrikarddaniel)
- [#389](https://github.com/nf-core/metatdenovo/pull/389) - template update to nf-core tools 3.3.2 plus module updates (@erikrikarddaniel)

### `Fixed`

- [#402](https://github.com/nf-core/metatdenovo/pull/402) - improve documentation for download of FigShare Diamond files (@erikrikarddaniel)
- [#402](https://github.com/nf-core/metatdenovo/pull/402) - allow the BBNorm process to only use 0.8 of the allocated memory not to fail on oversubscription of memory (@erikrikarddaniel)
- [#400](https://github.com/nf-core/metatdenovo/pull/400) - fix problems with `COLLECT_STATS` when single end reads are used; closes [#396](https://github.com/nf-core/metatdenovo/issues/396) (@erikrikarddaniel)
- [#398](https://github.com/nf-core/metatdenovo/pull/398) - make sure the EUKulele database directory is created if it doesn't exist (@erikrikarddaniel)
- [#391](https://github.com/nf-core/metatdenovo/pull/391),[#392](https://github.com/nf-core/metatdenovo/pull/392) - update the documentation and fix some inconsistencies in which output files are saved (@erikrikarddaniel)
- [#390](https://github.com/nf-core/metatdenovo/pull/390) - remove resource limits on full scale AWS tests to make it work (@erikrikarddaniel)

### `Dependencies`

### `Deprecated`

## v1.2.0 - [2025-06-18]

### `Added`

- [#373](https://github.com/nf-core/metatdenovo/pull/373) - Add module to save tsv with unique ORF Kofamscan hits to `<outdir>/summary_tables` (@erikrikarddaniel)
- [#366](https://github.com/nf-core/metatdenovo/pull/366) - Save amino acid sequences for HMMER hits (@erikrikarddaniel)

### `Changed`

- [#368](https://github.com/nf-core/metatdenovo/pull/368) - Added eukulele database name in filenames (@m3hdad)
- [#367](https://github.com/nf-core/metatdenovo/pull/367) - Gzip Transdecoder output (@erikrikarddaniel)
- [#359](https://github.com/nf-core/metatdenovo/pull/359) - Updated some descriptions and error messages in the json schema for better readability. Also made the input validation stricter in the hopes of preventing more errors during the pipeline run. (@herich0)
- [#352](https://github.com/nf-core/metatdenovo/pull/352) - Updated some modules (@erikrikarddaniel)

### `Fixed`

- [#380](https://github.com/nf-core/metatdenovo/pull/380) - Fix malformatted versions in two modules (@erikrikarddaniel)
- [#378](https://github.com/nf-core/metatdenovo/pull/378) - Add more nf-test tests (@erikrikarddaniel)
- [#377](https://github.com/nf-core/metatdenovo/pull/377) - Updated default nf-test (@erikrikarddaniel)
- [#376](https://github.com/nf-core/metatdenovo/pull/376) - Template update to nf-core tools 3.3.1 (@erikrikarddaniel)
- [#372](https://github.com/nf-core/metatdenovo/pull/372) - Fix bug in overall stats table creation for certain sample names (@erikrikarddaniel)
- [#371](https://github.com/nf-core/metatdenovo/pull/371) - Template update to nf-core tools 3.2.1 (@erikrikarddaniel)
- [#363](https://github.com/nf-core/metatdenovo/pull/363) - Handle duplicate names in taxonomies better (@erikrikarddaniel)
- [#362](https://github.com/nf-core/metatdenovo/pull/362) - Ensure correct Transdecoder publishing and test assertions (@m3hdad)
- [#361](https://github.com/nf-core/metatdenovo/pull/361) - Ensure `COLLECT_STATS` executes properly when trimming is skipped (@m3hdad).

### `Dependencies`

### `Deprecated`

## v1.1.1 - [2025-03-13]

### `Added`

### `Changed`

- [#364](https://github.com/nf-core/metatdenovo/pull/364) - Use `wget` not `gnu-wget` to fetch KofamScan database to improve arm64 support (@dslarm)
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
