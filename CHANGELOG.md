# nf-core/metatdenovo: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v2.0.0 - [YYYY-mm-dd]

### `Added`

- [#534](https://github.com/nf-core/metatdenovo/pull/534) - Add a pipeline test for the `test_bbnorm` profile, checking that normalised reads are used only for the assembly (@erikrikarddaniel)
- [#532](https://github.com/nf-core/metatdenovo/pull/532) - Add a test that runs Prokka on more than one batch (@erikrikarddaniel)
- [#526](https://github.com/nf-core/metatdenovo/pull/526) - Add `--save_parquet` to also write every `summary_tables/` table as Parquet, closes [#473](https://github.com/nf-core/metatdenovo/issues/473) (@erikrikarddaniel)
- [#506](https://github.com/nf-core/metatdenovo/pull/506) - Add a `-stub` pipeline test that runs several ORF callers with every annotation tool, addresses [#476](https://github.com/nf-core/metatdenovo/issues/476) (@erikrikarddaniel)
- [#505](https://github.com/nf-core/metatdenovo/pull/505) - Document how to resume a Megahit assembly that stopped partway through, in the new [large datasets](docs/usage/large_datasets.md) page (@erikrikarddaniel)
- [#504](https://github.com/nf-core/metatdenovo/pull/504) - Add featureCounts' `Unassigned_*` read categories as columns in `<assembly>.<caller>.overall_stats.tsv.gz`. The column order changes, so read columns by name, addresses [#451](https://github.com/nf-core/metatdenovo/issues/451) (@erikrikarddaniel)
- [#502](https://github.com/nf-core/metatdenovo/pull/502) - Download and build a MetaEuk reference database (`--metaeuk_db_name`, default `UniRef50`) when `--metaeuk_db` is not set, addresses [#485](https://github.com/nf-core/metatdenovo/issues/485) (@erikrikarddaniel)
- [#501](https://github.com/nf-core/metatdenovo/pull/501) - Add `--annotate_only_consolidated` (default `true`): with more than one ORF source, annotate only the protein-cluster representatives instead of every source's full protein set (@erikrikarddaniel)
- [#494](https://github.com/nf-core/metatdenovo/pull/494) - User-supplied ORFs (`--user_orfs_gff`/`--user_orfs_faa`, or several named sets in a `--user_orfs` CSV) can be combined with `--orf_caller` and take part in consolidation (@erikrikarddaniel)
- [#492](https://github.com/nf-core/metatdenovo/pull/492) - Run MetaEuk and TransDecoder on batches of contigs (`--metaeuk_batchsize`, `--transdecoder_batchsize`) to limit memory use. TransDecoder trains its coding model per batch, addresses [#486](https://github.com/nf-core/metatdenovo/issues/486) (@erikrikarddaniel)
- [#481](https://github.com/nf-core/metatdenovo/pull/481) - Add `--save_eukulele_alignments`. EUKulele's large Diamond alignment file is no longer published by default, addresses [#475](https://github.com/nf-core/metatdenovo/issues/475) (@erikrikarddaniel)
- [#479](https://github.com/nf-core/metatdenovo/pull/479) - Add protein consolidation: cluster the proteins from all ORF sources into `<assembly>.protein_consolidate_<identity>.counts.tsv.gz` (`--cluster_min_seq_id`, `--cluster_coverage`, `--skip_protein_consolidation`), addresses [#460](https://github.com/nf-core/metatdenovo/issues/460) (@erikrikarddaniel)
- [#469](https://github.com/nf-core/metatdenovo/pull/469) - Add `--save_bbduk_removed_fastq` to keep the reads that BBDuk removes, addresses [#17](https://github.com/nf-core/metatdenovo/issues/17) (@danilodileo)
- [#468](https://github.com/nf-core/metatdenovo/pull/468) - Add `--bbmap_ambiguous` and `--featurecounts_fraction` to control how reads that map to more than one place are counted, addresses [#464](https://github.com/nf-core/metatdenovo/issues/464) (@erikrikarddaniel)
- [#467](https://github.com/nf-core/metatdenovo/pull/467) - Add locus consolidation: with more than one ORF source, merge overlapping calls from different sources into loci in `<assembly>.locus_consolidate.counts.tsv.gz`, addresses [#463](https://github.com/nf-core/metatdenovo/issues/463) (@erikrikarddaniel)
- [#466](https://github.com/nf-core/metatdenovo/pull/466) - Run more than one ORF caller in the same run, e.g. `--orf_caller prokka,transdecoder`, addresses [#462](https://github.com/nf-core/metatdenovo/issues/462) (@erikrikarddaniel)
- [#465](https://github.com/nf-core/metatdenovo/pull/465) - Add MetaEuk as a splice-aware ORF caller for eukaryotes, addresses [#459](https://github.com/nf-core/metatdenovo/issues/459) (@erikrikarddaniel)
- [#461](https://github.com/nf-core/metatdenovo/pull/461) - Document how to build a taxonomy-aware Diamond database with nf-core/createtaxdb, addresses [#412](https://github.com/nf-core/metatdenovo/issues/412) (@erikrikarddaniel)
- [#458](https://github.com/nf-core/metatdenovo/pull/458) - Add Prodigal and TransDecoder ORF statistics to the MultiQC report, addresses [#456](https://github.com/nf-core/metatdenovo/issues/456) (@erikrikarddaniel)
- [#457](https://github.com/nf-core/metatdenovo/pull/457) - Add dbCAN CAZyme annotation (`--skip_dbcan`, `--dbcan_dbpath`), addresses [#60](https://github.com/nf-core/metatdenovo/issues/60) and [#430](https://github.com/nf-core/metatdenovo/issues/430) (@erikrikarddaniel)
- [#455](https://github.com/nf-core/metatdenovo/pull/455) - Add hidden Megahit k-mer and `--megahit_min_count` params for large datasets, addresses [#453](https://github.com/nf-core/metatdenovo/issues/453) (@erikrikarddaniel)
- [#452](https://github.com/nf-core/metatdenovo/pull/452) - Add tests for `--diamond_dbs` and KofamScan (@erikrikarddaniel)

### `Changed`

- [#527](https://github.com/nf-core/metatdenovo/pull/527) - KofamScan and eggNOG database downloads no longer need network access on the compute node, addresses [#365](https://github.com/nf-core/metatdenovo/issues/365) (@erikrikarddaniel)
- [#491](https://github.com/nf-core/metatdenovo/pull/491) - Replace TransRate with QUAST for assembly statistics, addresses [#487](https://github.com/nf-core/metatdenovo/issues/487) (@erikrikarddaniel)
- [#489](https://github.com/nf-core/metatdenovo/pull/489) - Use the shared nf-core module `custom/collectfeaturecounts` for count tables, with unchanged output (@erikrikarddaniel)
- [#488](https://github.com/nf-core/metatdenovo/pull/488) - Use the shared nf-core module `custom/collectstats`. The count column in `<assembly>.<caller>.overall_stats.tsv.gz` is named after the ORF caller (e.g. `prodigal`) instead of `n_feature_count` (@erikrikarddaniel)
- [#483](https://github.com/nf-core/metatdenovo/pull/483) - Move TransDecoder ORF-name cleanup into its own module, with unchanged output (@erikrikarddaniel)
- [#466](https://github.com/nf-core/metatdenovo/pull/466) - `featurecounts/` file names include the ORF caller, e.g. `SAMPLE1.prokka.featureCounts.tsv` (@erikrikarddaniel)
- [#454](https://github.com/nf-core/metatdenovo/pull/454) - Template update to nf-core/tools 4.1.0 and software updates (@erikrikarddaniel)
- [#452](https://github.com/nf-core/metatdenovo/pull/452) - Replace several local modules with nf-core/modules equivalents, addresses [#445](https://github.com/nf-core/metatdenovo/issues/445) (@erikrikarddaniel)
- [#450](https://github.com/nf-core/metatdenovo/pull/450) - Template update to nf-core/tools 4.0.3 and software updates (@erikrikarddaniel)

### `Fixed`

- [#548](https://github.com/nf-core/metatdenovo/pull/548) - Fix featureCounts failing on assemblies with more than about 75 million contigs, fixes [#547](https://github.com/nf-core/metatdenovo/issues/547) (@erikrikarddaniel)
- [#545](https://github.com/nf-core/metatdenovo/pull/545) - Fix the eggNOG-mapper version reported by conda runs, which was the pipeline's own release tag instead of the tool's version (@erikrikarddaniel)
- [#527](https://github.com/nf-core/metatdenovo/pull/527) - Fix the eggNOG database download URL (@erikrikarddaniel)
- [#525](https://github.com/nf-core/metatdenovo/pull/525) - Fix eggNOG, KofamScan and dbCAN database downloads to an `s3://` directory on AWS Batch, and `--eukulele_dbpath` on `s3://`, addresses [#471](https://github.com/nf-core/metatdenovo/issues/471) (@danilodileo)
- [#524](https://github.com/nf-core/metatdenovo/pull/524) - Fix BBNorm under Singularity/Apptainer when the host's `$TMPDIR` is not mounted in the container, addresses [#516](https://github.com/nf-core/metatdenovo/issues/516) (@erikrikarddaniel)
- [#523](https://github.com/nf-core/metatdenovo/pull/523) - Fix an out-of-date test snapshot, addresses [#522](https://github.com/nf-core/metatdenovo/issues/522) (@erikrikarddaniel)
- [#521](https://github.com/nf-core/metatdenovo/pull/521) - Leave memory headroom for the JVM in BBMap align, BBDuk and index, as BBNorm already did, addresses [#520](https://github.com/nf-core/metatdenovo/issues/520) (@danilodileo)
- [#519](https://github.com/nf-core/metatdenovo/pull/519) - Without `--save_bam`, no BAM files are published; samtools' sorted BAMs were published regardless, addresses [#474](https://github.com/nf-core/metatdenovo/issues/474) (@danilodileo)
- [#509](https://github.com/nf-core/metatdenovo/pull/509) - Fix MetaEuk output supplied via `--user_orfs_gff`/`--user_orfs_faa` merging separate loci into one, closes [#508](https://github.com/nf-core/metatdenovo/issues/508) (@erikrikarddaniel)
- [#506](https://github.com/nf-core/metatdenovo/pull/506) - Fix five `-stub` failures in the annotation modules, addresses [#476](https://github.com/nf-core/metatdenovo/issues/476) (@erikrikarddaniel)
- [#503](https://github.com/nf-core/metatdenovo/pull/503) - Boolean and integer params such as `--skip_eggnog false` or `--min_contig_length` now work from the command line, addresses [#478](https://github.com/nf-core/metatdenovo/issues/478) (@erikrikarddaniel)
- [#483](https://github.com/nf-core/metatdenovo/pull/483) - Round `tpm` to 6 decimals so count tables computed separately agree, closes [#484](https://github.com/nf-core/metatdenovo/issues/484) (@erikrikarddaniel)
- [#481](https://github.com/nf-core/metatdenovo/pull/481) - Fix the EUKulele `-stub` run (@erikrikarddaniel)
- [#480](https://github.com/nf-core/metatdenovo/pull/480) - Fix TransDecoder.Predict failing on `-resume`, addresses [#477](https://github.com/nf-core/metatdenovo/issues/477) (@erikrikarddaniel)
- [#479](https://github.com/nf-core/metatdenovo/pull/479) - Run eggNOG-mapper for every ORF caller, not only the first (@erikrikarddaniel)
- [#479](https://github.com/nf-core/metatdenovo/pull/479) - Fix the MultiQC ORF statistics crashing when an ORF caller finds no proteins (@erikrikarddaniel)
- [#466](https://github.com/nf-core/metatdenovo/pull/466) - Fix a crash when an HMM search finds no hits (@erikrikarddaniel)
- [#458](https://github.com/nf-core/metatdenovo/pull/458) - Fix Prokka statistics in the MultiQC report merging into one sample when an assembly has more than one Prokka batch, addresses [#456](https://github.com/nf-core/metatdenovo/issues/456) (@erikrikarddaniel)
- [#452](https://github.com/nf-core/metatdenovo/pull/452) - Fix `EGGNOG_FORMAT` naming the ORF column `Lorf` instead of `orf` (@erikrikarddaniel)

### `Dependencies`

| Tool        | Previous version | New version |
| ----------- | ---------------- | ----------- |
| samtools    | 1.23.1           | 1.24        |
| multiqc     | 1.34             | 1.35        |
| trim-galore | 2.1.0            | 2.3.0       |
| prokka      | 1.14.6           | 1.15.6      |
| nft-utils   | 0.0.3            | 1.2.0       |

## v1.4.1 - [2026-09-15]

### `Fixed`

- [#517](https://github.com/nf-core/metatdenovo/pull/517) - Fix BBNorm failing under Singularity/Apptainer when the host's `$TMPDIR` isn't bind-mounted into the container, closes issue [#513](https://github.com/nf-core/metatdenovo/issues/513) (@erikrikarddaniel)

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
