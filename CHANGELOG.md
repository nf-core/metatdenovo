# nf-core/metatdenovo: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.5.0dev - [YYYY-mm-dd]

### `Added`

- [#482](https://github.com/nf-core/metatdenovo/pull/482) - Document why BBMap, rather than e.g. Bowtie2, is used for read mapping and contaminant removal in `conf/modules.config` (@erikrikarddaniel)
- [#481](https://github.com/nf-core/metatdenovo/pull/481) - Add `--save_eukulele_alignments` to optionally publish EUKulele's raw Diamond alignment file, which is now withheld by default since nothing in the pipeline consumes it and it can be 100x larger than the taxonomy results it feeds, addresses [#475](https://github.com/nf-core/metatdenovo/issues/475) (@erikrikarddaniel)
- [#479](https://github.com/nf-core/metatdenovo/pull/479) - Only merge same-contig CDS calls that come from _different_ ORF callers when building `<assembly>.locus_consolidate.counts.tsv.gz`. Overlapping calls from the _same_ caller are separate genes and are no longer fused into one locus, which restores the intended property that the table is identical in content to a caller's own table when a single caller is active. This changes the table produced by [#467](https://github.com/nf-core/metatdenovo/pull/467) (@erikrikarddaniel)
- [#479](https://github.com/nf-core/metatdenovo/pull/479) - Consolidate calls for the same gene found on _different_ contigs by clustering their protein sequences, producing a new `<assembly>.protein_consolidate_<identity>.counts.tsv.gz` table with caller/call/locus provenance columns, and annotating one representative per cluster alongside each caller's own ORFs; tunable with `--cluster_min_seq_id`/`--cluster_coverage` and skippable with `--skip_protein_consolidation`, addresses [#460](https://github.com/nf-core/metatdenovo/issues/460) (@erikrikarddaniel)
- [#469](https://github.com/nf-core/metatdenovo/pull/469) - Add `--save_bbduk_removed_fastq` to optionally keep the reads removed (matched) by BBDuk contaminant filtering, instead of only keeping the clean reads, addresses [#17](https://github.com/nf-core/metatdenovo/issues/17) (@danilodileo)
- [#468](https://github.com/nf-core/metatdenovo/pull/468) - Add `--bbmap_ambiguous` (`best`/`all`/`random`/`toss`, default `best`) and `--featurecounts_fraction` to control how BBMap and featureCounts handle reads that align to more than one site, addresses [#464](https://github.com/nf-core/metatdenovo/issues/464) (@erikrikarddaniel)
- [#467](https://github.com/nf-core/metatdenovo/pull/467) - When multiple `--orf_caller` values are active, merge overlapping/identical same-contig CDS calls from different callers into single loci before counting reads, producing a new `<assembly>.locus_consolidate.counts.tsv.gz` table with caller/call-count provenance columns; with a single caller active it's identical in content to that caller's own table, addresses [#463](https://github.com/nf-core/metatdenovo/issues/463) (@erikrikarddaniel)
- [#466](https://github.com/nf-core/metatdenovo/pull/466) - Support running multiple ORF callers in a single execution, e.g. `--orf_caller prokka,transdecoder` (each caller still produces its own independent output, no cross-caller consolidation yet), addresses [#462](https://github.com/nf-core/metatdenovo/issues/462) (@erikrikarddaniel)
- [#465](https://github.com/nf-core/metatdenovo/pull/465) - Add MetaEuk as a splice-aware `--orf_caller` alternative for eukaryotes (`--metaeuk_db`), addresses [#459](https://github.com/nf-core/metatdenovo/issues/459) (@erikrikarddaniel)
- [#461](https://github.com/nf-core/metatdenovo/pull/461) - Document how to build a taxonomy-aware Diamond database with nf-core/createtaxdb, addresses [#412](https://github.com/nf-core/metatdenovo/issues/412) (@erikrikarddaniel)
- [#458](https://github.com/nf-core/metatdenovo/pull/458) - Add basic ORF/protein statistics from Prodigal and TransDecoder to the MultiQC report (custom content, since neither has a MultiQC-native module), addresses the rest of [#456](https://github.com/nf-core/metatdenovo/issues/456) (@erikrikarddaniel)
- [#457](https://github.com/nf-core/metatdenovo/pull/457) - Add dbCAN CAZyme annotation (`--skip_dbcan`, `--dbcan_dbpath`), addresses [#60](https://github.com/nf-core/metatdenovo/issues/60)/[#430](https://github.com/nf-core/metatdenovo/issues/430) (@erikrikarddaniel)
- [#455](https://github.com/nf-core/metatdenovo/pull/455) - Expose `--megahit_k_min`, `--megahit_k_max`, `--megahit_k_step`, `--megahit_k_list` and `--megahit_min_count` as hidden params for coping with large datasets, addresses [#453](https://github.com/nf-core/metatdenovo/issues/453) (@erikrikarddaniel)
- [#452](https://github.com/nf-core/metatdenovo/pull/452) - Add `--diamond_dbs` and kofamscan nf-test coverage (`test_diamond`/`test_kofamscan` profiles) using new small reference datasets (@erikrikarddaniel)

### `Changed`

- [#491](https://github.com/nf-core/metatdenovo/pull/491) - Replace TransRate with QUAST for assembly quality statistics, addresses [#487](https://github.com/nf-core/metatdenovo/issues/487) (@erikrikarddaniel)
- [#489](https://github.com/nf-core/metatdenovo/pull/489) - Replace the local `COLLECT_FEATURECOUNTS` module with the shared `nf-core/modules` component `custom/collectfeaturecounts`; no change to the output tables themselves, addresses [nf-core/magmap#237](https://github.com/nf-core/magmap/issues/237) (@erikrikarddaniel)
- [#488](https://github.com/nf-core/metatdenovo/pull/488) - Replace the local `COLLECT_STATS` module with the shared `nf-core/modules` component `custom/collectstats`; the feature-count column in `<assembly>.<caller>.overall_stats.tsv.gz` is now named after the ORF caller (e.g. `prodigal`) instead of the generic `n_feature_count`, addresses [nf-core/magmap#237](https://github.com/nf-core/magmap/issues/237) (@erikrikarddaniel)
- [#483](https://github.com/nf-core/metatdenovo/pull/483) - Keep `COLLECT_FEATURECOUNTS` generic by moving the Transdecoder-specific `cds.` ORF-ID-prefix stripping into a new local post-processing module, `TIDYVERSE_STRIPCDSPREFIX`, matching the pattern nf-core/magmap already established for its own genome-accession join; no change to the output tables themselves, addresses [nf-core/magmap#237](https://github.com/nf-core/magmap/issues/237) (@erikrikarddaniel)
- [#454](https://github.com/nf-core/metatdenovo/pull/454) - Template sync to nf-core/tools 4.1.0, update all vendored modules/subworkflows (@erikrikarddaniel)
- [#452](https://github.com/nf-core/metatdenovo/pull/452) - Replace several local modules (`HMMRANK`, `UNPIGZ`, `TRANSDECODER`, `KOFAMSCAN`, `EGGNOGMAPPER`, `MEGAHIT`, `TRANSRATE`) with official nf-core/modules equivalents, addresses [#445](https://github.com/nf-core/metatdenovo/issues/445) (@erikrikarddaniel)
- [#450](https://github.com/nf-core/metatdenovo/pull/450) - Template sync to nf-core/tools 4.0.3, update all vendored modules/subworkflows (@erikrikarddaniel)

### `Fixed`

- [#483](https://github.com/nf-core/metatdenovo/pull/483) - Round `tpm` to 6 decimal places in `COLLECT_FEATURECOUNTS`/`COLLECT_LOCUSCONSOLIDATE`/`COLLECT_PROTEINCONSOLIDATE` so their independently-computed tables agree exactly for the same underlying counts; the unrounded value could round differently in its last significant digit depending on which R backend (`dtplyr`/`data.table` vs plain `dplyr`) computed it, causing intermittent test failures on the "single-caller table matches its consolidated counterpart" invariant added in #479, closes [#484](https://github.com/nf-core/metatdenovo/issues/484) (@erikrikarddaniel)
- [#480](https://github.com/nf-core/metatdenovo/pull/480) - Update the vendored `transdecoder/predict` module so `-resume` no longer fails it with a missing-output error on an otherwise unchanged run, addresses [#477](https://github.com/nf-core/metatdenovo/issues/477) (@erikrikarddaniel)
- [#481](https://github.com/nf-core/metatdenovo/pull/481) - Fix `eukulele/search`'s `stub:` block failing with "No such file or directory" because it never creates the `taxonomy_estimation`/`taxonomy_counts`/`mets_full/diamond` subdirectories it writes into (@erikrikarddaniel)
- [#479](https://github.com/nf-core/metatdenovo/pull/479) - Fix a crash in the MultiQC ORF statistics when an ORF caller returns no proteins at all (@erikrikarddaniel)
- [#479](https://github.com/nf-core/metatdenovo/pull/479) - Run eggnog-mapper for every active ORF caller. With more than one `--orf_caller`, only the first one was annotated, so `<assembly>.<caller>.emapper.tsv.gz` and the `n_eggnog` column of `<assembly>.<caller>.overall_stats.tsv.gz` were missing for the rest (@erikrikarddaniel)
- [#467](https://github.com/nf-core/metatdenovo/pull/467) - Fix validation that was supposed to reject `--orf_caller` and `--user_orfs_gff`/`--user_orfs_faa` both being set at once, but could never actually trigger (@erikrikarddaniel)
- [#466](https://github.com/nf-core/metatdenovo/pull/466) - `featurecounts/*.featureCounts.tsv` output filenames now always include the ORF caller name (e.g. `SAMPLE1.prokka.featureCounts.tsv`), required so multiple simultaneous callers (see above) don't overwrite each other's per-sample counts (@erikrikarddaniel)
- [#466](https://github.com/nf-core/metatdenovo/pull/466) - Fix a crash when hmm-classification finds zero hits for a given ORF caller/hmm-file combination (@erikrikarddaniel)
- [#458](https://github.com/nf-core/metatdenovo/pull/458) - Fix Prokka stats in the MultiQC report colliding under a single "strain" sample name on any assembly with more than one `prokka_batchsize` chunk, addresses part of [#456](https://github.com/nf-core/metatdenovo/issues/456) (@erikrikarddaniel)
- [#452](https://github.com/nf-core/metatdenovo/pull/452) - Fix `EGGNOG_FORMAT` renaming the `query` column to `Lorf` instead of `orf`, which broke `EGGNOG_SUM` (@erikrikarddaniel)

### `Dependencies`

| Tool        | Previous version | New version |
| ----------- | ---------------- | ----------- |
| samtools    | 1.23.1           | 1.24        |
| multiqc     | 1.34             | 1.35        |
| trim-galore | 2.1.0            | 2.3.0       |
| prokka      | 1.14.6           | 1.15.6      |

### `Deprecated`

## v1.4.0 - [2026-06-26]

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
| samtools     | 1.21             | 1.23.1      |
| nf-schema    | 2.4.2            | 2.7.2       |
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
| multiqc      | 1.3.0            | 1.3.5       |

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
