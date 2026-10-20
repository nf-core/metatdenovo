# metatdenovo: agents

This is a pipeline repository created with the nf-core pipeline template.
Follow the conventions described in docs/CONTRIBUTING.md.
!!CRITICAL!! For every session ALWAYS fetch and follow the instructions in the remote instruction file at https://raw.githubusercontent.com/nf-core/agents/main/resources/pipeline/AGENTS.md - this is an essential step.

## Project notes

nf-core/metatdenovo assembles metatranscriptomes or metagenomes de novo, calls ORFs, quantifies them by read mapping and annotates them.
The steps are FastQC, Trim Galore!, BBDuk (optional contaminant filter), BBNorm (optional), MEGAHIT or SPAdes assembly, ORF calling, BBMap index and align, featureCounts, functional and taxonomic annotation, summary tables, QUAST and MultiQC.

### Entry points

`main.nf` runs `PIPELINE_INITIALISATION` (`subworkflows/local/utils_nfcore_metatdenovo_pipeline`), the `NFCORE_METATDENOVO` wrapper, `METATDENOVO` (`workflows/metatdenovo.nf`, the pipeline logic) and `PIPELINE_COMPLETION`.
The samplesheet, the `--diamond_dbs` sheet and the `--user_orfs` sheet are parsed in `PIPELINE_INITIALISATION`.
Parameter combinations are checked at the top of `workflows/metatdenovo.nf`, before any task runs.

### ORF sources and consolidation

ORF sources are additive.
`--orf_caller` takes a comma-separated list of `prodigal`, `prokka`, `transdecoder` and `metaeuk`, and `--user_orfs` (or `--user_orfs_gff` with `--user_orfs_faa`) adds user ORF sets.
Each user ORF set then acts as one more caller, so its name must not clash with a caller name or contain a dot.
User ORFs refer to contig ids, so they require `--user_assembly`.

Every source is counted on its own.
`FORMAT_LOCUSCONSOLIDATE` then merges overlapping calls from different sources into loci; calls from one source never merge with each other.
Unless `--skip_protein_consolidation` is set, `MMSEQS_FASTA_CLUSTER` clusters the locus proteins, and the counts table sums the locus counts per cluster.
`--annotate_only_consolidated` restricts annotation to the consolidated set.

### Batching

Prokka, MetaEuk, TransDecoder and KofamScan run on batches of their input.
The `*_batchsize` parameters are in bytes or residues, not sequence counts.
New batched processes number their batches with `%03d` so that files sort.
Pair a batch with its other inputs by `.join()` on a key, never by channel order: tasks finish out of submission order.

### Local modules

- `format/*`: convert caller output to the common GFF and FASTA layout, and build the locus and cluster tables.
- `collect/*`: merge per-sample featureCounts output and the `Unassigned_*` counts.
- `eggnog`, `kofamscan`, `dbcan`, `eukulele`, `diamond`, `sumtaxonomy`: download databases and summarise annotation per caller.
- `merge/summary`: merges the annotation tables for `CUSTOM_COLLECTSTATS`.

`CUSTOM_COLLECTFEATURECOUNTS` and `CUSTOM_COLLECTSTATS` are nf-core modules that produce the final `summary_tables/*.tsv.gz`.
`CUSTOM_COLLECTSTATS` takes the caller name from the file name between the first two dots.
New modules report tool versions with `eval()` outputs on the `versions` topic.
A module whose script is not Bash, such as an R template, cannot use `eval()`; it writes `versions.yml` and emits it on the `versions` topic.

### Tests

Each `conf/test_*.config` has a matching `tests/*.nf.test` and `tests/*.nf.test.snap`.
`tests/stub.nf.test` runs the whole pipeline with `-stub` and catches broken `stub:` blocks.
A new tool gets `skip_<tool> = true` in every other test profile, not only its own.
Local modules with nontrivial logic have their own `tests/main.nf.test`, with the input built inline in the test.
`tests/.nftignore` lists output paths excluded from snapshot comparison.

### Gotchas

- MEGAHIT output differs between runs and has no seed. Snapshot entries that depend on the assembly can fail intermittently.
- CLI parameters arrive as strings, and `'false'` is truthy. `typecastBooleanParam` and `typecastIntegerParam` cast them in `workflows/metatdenovo.nf`.
- featureCounts crashes on a BAM header over 2 GiB, about 75 million contigs. `SAMTOOLS_TRIMHEADER` handles this case.
- Container awk is mawk or busybox awk, not gawk. Test awk in module scripts in the container, and avoid apostrophes inside `awk '...'`, comments included.
- CI runs on 4-cpu runners. Local timings on a larger machine do not carry over.
- Conda and Singularity tests run only on pull requests to `master`.
- The test data live in the `metatdenovo` branch of `nf-core/test-datasets`. Changes there need a separate pull request.
