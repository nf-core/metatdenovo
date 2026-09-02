# nf-core/metatdenovo: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/metatdenovo/usage](https://nf-co.re/metatdenovo/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

Metatdenovo is a workflow primarily designed for annotation of metatranscriptomes and metagenomics for which reference genomes are not available.
The approach is to first create an assembly, then call genes and finally quantify and annotate the genes.
Since the workflow includes gene callers and annotation tools and databases for prokaryotes, eukaryotes and viruses, the workflow should be suitable for all
organism groups and mixed communities can be handled by trying different gene callers and comparing the results.

While the rationale for writing the workflow was metatranscriptomes, there is nothing in the workflow that precludes use for single organisms rather than
communities nor genomes rather than transcriptomes.
Instead, the workflow should be usable for any project in which a de novo assembly followed by quantification and annotation is suitable.

If you're working with a large project -- many samples, deep sequencing, or both -- and expect (or hit) memory problems during assembly, see [Coping with large datasets](large_datasets.md) for concrete params to try, and in which order.

## Running the workflow

### Quickstart

A typical command for running the workflow is:

```bash
nextflow run nf-core/metatdenovo -profile docker --outdir results/ --input samples.csv
```

### Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline.
Use this parameter to specify its location.
It must be a comma-separated file with 3 columns, and a header row as shown in the examples below.
The input fastq files need to have names ending with `.fastq` or `.fq`, followed by `.gz` if gzipped.

```csv title="samplesheet.csv"
sample,fastq_1,fastq_2
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
CONTROL_REP1,AEG588A1_S1_L003_R1_001.fastq.gz,AEG588A1_S1_L003_R2_001.fastq.gz
CONTROL_REP1,AEG588A1_S1_L004_R1_001.fastq.gz,AEG588A1_S1_L004_R2_001.fastq.gz
```

#### Full samplesheet

<!-- I commented out text about single-end samples as we don't know whether this works yet. -->
<!-- The pipeline will auto-detect whether a sample is single- or paired-end using the information provided in the samplesheet. The samplesheet can have as many columns as you desire, however, there is a strict requirement for the first 3 columns to match those defined in the table below. -->

<!-- A final samplesheet file consisting of both single- and paired-end data may look something like the one below. This is for 6 samples, where `TREATMENT_REP3` has been sequenced twice. -->

A final samplesheet file consisting of samples taken at time 0 and 24 in triplicate may look like the one below.

```csv title="samplesheet.csv"
sample,fastq_1,fastq_2
T0a,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
T0b,AEG588A2_S2_L002_R1_001.fastq.gz,AEG588A2_S2_L002_R2_001.fastq.gz
T0c,AEG588A3_S3_L002_R1_001.fastq.gz,AEG588A3_S3_L002_R2_001.fastq.gz
T24a,AEG588A4_S1_L002_R1_001.fastq.gz,AEG588A4_S1_L002_R2_001.fastq.gz
T24b,AEG588A5_S2_L002_R1_001.fastq.gz,AEG588A5_S2_L002_R2_001.fastq.gz
T24c,AEG588A6_S3_L002_R1_001.fastq.gz,AEG588A6_S3_L002_R2_001.fastq.gz
```

| Column    | Description                                                                                                                                                                            |
| --------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`  | Custom sample name. This entry will be identical for multiple sequencing libraries/runs from the same sample. Spaces in sample names are automatically converted to underscores (`_`). |
| `fastq_1` | Full path to FastQ file for Illumina short reads 1. The file has to have the extension ".fastq" or ".fq", followed by ".gz" if gzipped.                                                |
| `fastq_2` | Full path to FastQ file for Illumina short reads 2. The file has to have the extension ".fastq" or ".fq", followed by ".gz" if gzipped.                                                |

An [example samplesheet](../assets/samplesheet.csv) has been provided with the pipeline.

#### Multiple runs of the same sample

The `sample` identifiers have to be the same when you have re-sequenced the same sample more than once e.g. to increase sequencing depth. The pipeline will concatenate the raw reads before performing any downstream analysis. Below is an example for the same sample sequenced across 3 lanes:

```console
sample,fastq_1,fastq_2
T0a,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
T0a,AEG588A1_S1_L003_R1_001.fastq.gz,AEG588A1_S1_L003_R2_001.fastq.gz
T0a,AEG588A1_S1_L004_R1_001.fastq.gz,AEG588A1_S1_L004_R2_001.fastq.gz
```

#### Paired end, single end or a mix is allowed

The samplesheet may contain either paired end reads (as in the examples above), single end reads or a mixture.
For single end reads, supply the read file in the `fastq_1` field and leave `fastq_2` empty.
_Do not mix single end and pairs for the same sample!_

### Filter/remove sequences from the samples (e.g. rRNA sequences with SILVA database)

The pipeline can remove potential contaminants using the BBduk program.
Specify a fasta file, gzipped or not, with `--sequence_filter sequences.fasta`.
Use `--save_bbduk_removed_fastq` to also keep the reads that were removed (i.e. matched the filter reference) instead of only keeping the filtered/clean reads.
For further documentation, see the [BBduk official website](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/).

### Digital normalization

Metatdenovo can perform "digital normalization" of the reads before the assembly.
This will reduce coverage of highly abundant sequences and remove sequences that are below a threshold.
This is normally used if the data set is too large to assemble but can potentially improve an assembly by reducing the coverage of very abundant sequences.
N.B. digitally normalized reads are used only for the assembly and the full set of sequences will be used for quantification.
To turn on digital normalization, use the `--bbnorm` parameter and, if required, adjust the `--bbnorm_target` and `--bbnorm_min` parameters.

Digital normalization directly and predictably shrinks the assembler's memory and time footprint, since that scales with the size of the k-mer graph, which in turn scales with the (normalized) k-mer diversity and depth.
It's a coverage-flattening transform, though, so it comes with a real trade-off: k-mers from low-abundance organisms can fall below `--bbnorm_min` and get discarded outright, which means those organisms can be missing from the assembly entirely, not just under-represented in it.
It can also distort strain-level variation for the same reason.
Since normalized reads are only used for the assembly (see the note above), this risk is contained to "what gets assembled" -- it does not bias the abundance counts of whatever does make it into the assembly.

> Please, check the [bbnorm](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbnorm-guide/) documentation for further information about these programs and how digital normalization works. Remember to check [Parameters](https://nf-co.re/metatdenovo/parameters) page for the all options that can be used for this step.

See [Coping with large datasets](large_datasets.md) for concrete `--bbnorm_target`/`--bbnorm_min` starting values, how this combines with the [Assembler options](#assembler-options) below, and in which order to try them.

### Assembler options

By default, the pipeline uses Megahit (`--assembler megahit`) to assemble the cleaned and trimmed reads to create the reference contigs.
Megahit is fast and it does not require a lot of memory to run, making it ideal for large sets of samples.

If Megahit still runs out of memory on very large datasets, a few hidden parameters let you tune its k-mer graph construction: `--megahit_min_count`, `--megahit_k_min`, `--megahit_k_max` and `--megahit_k_step` (or `--megahit_k_list` instead of the min/max/step triplet).
Raising `--megahit_min_count` prunes low-frequency, often erroneous k-mers before the graph is built, and is a more surgical adjustment than the k-mer options.
Raising `--megahit_k_min` skips the smallest, most memory-hungry k-mer iterations, but at the cost of sensitivity to low-coverage or short reads -- treat it more as a last resort.
None of these parameters have a pipeline default; when left unset, Megahit uses its own built-in defaults, which are recorded in its own log file (under `megahit/`) for each run.
See the [Megahit documentation](https://github.com/voutcn/megahit) for the full meaning of these options, and [Coping with large datasets](large_datasets.md) for concrete starting values and where these fit relative to [digital normalization](#digital-normalization) above.

The workflow also supports Spades (`--assembler spades` ) as an alternative.
The default "flavour" of Spades is set to RNA, but this can be changed using the `--spades_flavor` parameter (see [parameter documentation](/metatdenovo/parameters/#spades_flavor))

You can also choose to input contigs from an assembly that you made outside the pipeline using the `--user_assembly file.fna.gz` (where `file.fna.gz` is the name of a fasta file with contigs) parameter.
When you use your own assembly, the name of this -- used in output file names -- can be set using the `--user_assembly_name` parameter.

### ORF caller options

By default, the pipeline uses prodigal (`--orf_caller prodigal` ) to call genes/ORFs from the assembly.
This is suitable for prokaryotes, as is the Prokka alternative (`--orf_caller prokka`).
The latter uses Prodigal to identify ORFs, but also identifies other types of features.
It also performs functional annotation of ORFs.

For eukaryotic genes, we recommend users to use Transdecoder (`--orf_caller transdecoder`) to call ORFs.

MetaEuk (`--orf_caller metaeuk`) is a second eukaryote-targeted alternative.
Unlike Transdecoder, which calls ORFs directly on the assembled contigs/transcripts, MetaEuk is splice-aware: it aligns contigs against a reference protein database and can call a single gene model spanning an intron, which matters for assemblies that include intron-containing genomic sequence alongside spliced transcripts.
It requires a reference protein database, set with `--metaeuk_db` -- either a protein fasta file or a directory containing an mmseqs2-formatted database.
There's no default, so `--metaeuk_db` must be set whenever `--orf_caller metaeuk` is used.

MetaEuk's own `metaeuk databases` command (outside the pipeline) can download and format several standard amino acid reference databases for this purpose, for example:

- `UniRef50`/`UniRef90`/`UniRef100`
- `UniProtKB/Swiss-Prot`
- `GTDB`
- `NR`

```bash
metaeuk databases UniRef50 metaeuk_uniref50 tmp
```

The resulting `metaeuk_uniref50` can then be passed directly as `--metaeuk_db`.
Run `metaeuk databases -h` for the full list.

#### Running more than one ORF caller

`--orf_caller` also accepts a comma-separated list, e.g. `--orf_caller prokka,transdecoder`, to run more than one caller in the same execution -- useful for mixed prokaryote/eukaryote communities.
Each active caller still runs independently and produces its own complete set of output (GFF, protein FASTA, `summary_tables/<assembly>.<caller>.*`) under its own name, exactly as if it had been run alone.
In addition, the pipeline consolidates calls that different callers made for the same gene, at two levels; see [Consolidating calls for the same gene](#consolidating-calls-for-the-same-gene) below and the [Output docs](output.md).

##### Comparing callers: mind the minimum ORF length

Do not read the per-caller tables, or the `callers` column of the consolidated tables, as a like-for-like comparison of what each caller found.
The callers do not use comparable minimum ORF lengths, and the difference is large enough to dominate a naive count.

TransDecoder's `LongOrfs` step defaults to a 100 aa minimum, so it cannot report anything shorter.
Prodigal, which Prokka wraps and which this pipeline runs in metagenome mode, has no equivalent floor and readily calls short ORFs, including partial ones running off the ends of fragmentary contigs.
A caller-specific count therefore mixes "this caller found a gene the other missed" with "the other caller was never eligible to call it".

The pipeline's own full-size test makes the size of the effect concrete (its output is published with the release, so the numbers below can be inspected).
Running `--orf_caller prokka,transdecoder` on a metatranscriptome gave 68829 consolidated loci:

| Loci              | Count | Mean protein length |
| ----------------- | ----- | ------------------- |
| Called by both    | 23034 | ~245 aa             |
| TransDecoder only | 29844 | ~168 aa             |
| Prokka only       | 15951 | ~80 aa              |

83% of the Prokka-only loci are shorter than 100 aa, i.e. below the length TransDecoder would ever report.
Comparing only what both callers could have found leaves about 29800 TransDecoder-only against about 2700 Prokka-only -- an order of magnitude apart, where the raw counts suggest a factor of two.

If you want to compare callers, apply a common length threshold first.
The short Prodigal calls are not necessarily wrong -- small proteins are real, and so are partial genes at contig edges -- they are simply not something the other caller was in a position to confirm.

#### Provide your own ORFs

You can add one or more sets of pre-called ORFs to the pipeline with `--user_orfs orfs.csv`, a comma-separated
file with a header row and three columns: `name`, `gff`, `faa`. Each row supplies one named set of ORFs (a gff
file and its matching amino acid fasta), and is treated exactly like another `--orf_caller` value from that
point on -- it takes part in locus consolidation, protein consolidation and feature counting alongside
whichever built-in callers are also active.

`--user_orfs` is additive with `--orf_caller`, not mutually exclusive with it: provide `--orf_caller`,
`--user_orfs`, or both, but at least one of them is required. Every row's `name` must be unique, and cannot
collide with an active `--orf_caller` value.

```csv title="orfs.csv"
name,gff,faa
my_caller,orfs.gff.gz,orfs.faa.gz
```

### Read mapping and quantification options

After assembly and ORF calling, reads are mapped back to the assembly with BBMap and quantified per ORF/CDS with featureCounts, producing the `summary_tables/<assembly>.<orfcaller>.counts.tsv.gz` tables (and, when multiple ORF callers are active, the locus-consolidated table described above).

A read can align equally well to more than one site in the assembly.
`--bbmap_ambiguous` controls how BBMap handles that:

| bbmap_ambiguous | Behaviour                           |
| --------------- | ----------------------------------- |
| best (default)  | Use the first best-scoring site     |
| random          | Pick one top-scoring site at random |
| toss            | Treat the read as unmapped          |
| all             | Retain every top-scoring site       |

`best`, `random` and `toss` all guarantee at most one reported alignment per read, which is exactly what the counts tables above -- and the locus consolidation described above ([issue #463](https://github.com/nf-core/metatdenovo/issues/463)), and any future cross-contig consolidation -- assume when summing counts per ORF or locus: a read is never counted more than once.
`all` deliberately breaks that guarantee: a read can legitimately count toward more than one ORF.

`--featurecounts_fraction` controls how featureCounts treats reads with more than one reported alignment.
By default (`false`), a read that maps to N sites contributes a full count to each of the N ORFs it hits, so per-sample totals inflate under `--bbmap_ambiguous all`.
Setting `--featurecounts_fraction` makes each such read contribute 1/N to each site instead, keeping totals close to the number of reads actually sequenced.
It has no effect unless `--bbmap_ambiguous all` is also set, since `best`/`random`/`toss` never produce more than one reported alignment per read to begin with.

A typical command line enabling both:

```bash
nextflow run nf-core/metatdenovo -profile docker --outdir results/ --input samplesheet.csv --assembler megahit --orf_caller prodigal --bbmap_ambiguous all --featurecounts_fraction
```

### Consolidating calls for the same gene

The same gene can be called more than once: by two ORF callers on one contig, or -- in a mixed prokaryote/eukaryote assembly -- on two different contigs, when a gene is assembled both from genomic DNA and from its transcript.
Counting each call separately would count one gene's reads several times, so the pipeline reports counts at three levels, from most raw to most consolidated.

| Table                                                     | Consolidation                                                                               |
| --------------------------------------------------------- | ------------------------------------------------------------------------------------------- |
| `<assembly>.<orfcaller>.counts.tsv.gz`                    | None: exactly what that caller found                                                        |
| `<assembly>.locus_consolidate.counts.tsv.gz`              | Overlapping CDS calls from different callers **on the same contig**, merged before counting |
| `<assembly>.protein_consolidate_<identity>.counts.tsv.gz` | Calls **on different contigs** whose proteins cluster together                              |

The three are produced side by side, so you can choose how much consolidation to trust; the more consolidated levels also carry provenance columns (`callers`, `n_calls`, and at the protein level `n_loci` and `loci`) so you can see what was merged into each row.

Locus consolidation is coordinate-based and close to unambiguous.
Protein consolidation is inference-based: it links a splice-aware genomic call to a transcript-derived call for the same gene, which have no coordinate relationship at all but converge on nearly the same protein.
It can also merge genuinely distinct paralogs, which is why the default identity is deliberately high.

Reads are never remapped onto a cluster representative.
Each read keeps aligning to the contig it actually came from and the per-locus counts are summed afterwards, so a cluster member that differs from the representative in nucleotide sequence -- through synonymous variation, say -- is not under-counted.

`--cluster_min_seq_id` sets the minimum protein identity for two calls to be treated as one gene, and `--cluster_coverage` the minimum fraction of both proteins the alignment must cover.
Because these are separate, `--cluster_min_seq_id 1.0` is less strict than it sounds: calls that agree exactly where they align but are trimmed differently at the ends still cluster.
Conversely, a splice-aware caller whose exon boundary is off by a few residues in the middle of a gene falls below 1.0, which is why the default is `0.99` rather than exact identity.

The identity is part of the output name, rounded to whole percent, so runs at different settings don't overwrite each other: `--cluster_min_seq_id 0.99` gives `protein_consolidate_99` and `0.9` gives `protein_consolidate_90`.

Unlike locus consolidation, protein consolidation is **not** a no-op when only one ORF caller is active: near-identical proteins on different contigs are merged whichever caller found them, which also collapses assembly redundancy.
Use `--skip_protein_consolidation` to turn it off.

How much it changes depends mostly on how many ORF callers are active, and secondarily on the community.

With a **single caller** there is no cross-caller redundancy to find, so the step only collapses cases where one caller called near-identical proteins on two different contigs.
On the pipeline's own small prokaryotic test dataset that is nothing at all: 4397 loci produce 4397 clusters, and the protein-consolidated table is identical in content to the locus-consolidated one.

With **two callers** it does what it is for.
On a full-size metatranscriptome run with `--orf_caller prokka,transdecoder` at the default `--cluster_min_seq_id 0.99`, 68829 loci produced 68300 clusters.
476 of those clusters held more than one locus, covering 1005 loci between them (mean 2.1 members), so 1.5% of loci were grouped with at least one other and the feature count fell by 529.
Of the 476, 315 linked calls from different callers -- the genomic-versus-transcript case this exists for -- while 91 merged two Prokka calls and 70 two TransDecoder calls, i.e. redundancy within a single caller.

So the effect is modest in absolute terms, a few percent of loci at most, but two thirds of it is the intended cross-caller kind rather than redundancy collapse.
If you run a single caller and are watching runtime, `--skip_protein_consolidation` costs you little.
Lowering `--cluster_min_seq_id` is what makes the step start merging paralogs and strain variants, so change it deliberately rather than to "make something happen".

```bash
nextflow run nf-core/metatdenovo -profile docker --outdir results/ --input samplesheet.csv --assembler megahit --orf_caller metaeuk,transdecoder --metaeuk_db /path/to/db --cluster_min_seq_id 0.95
```

### Taxonomic annotation options

Metatdenovo uses two different programs for taxonomy annotation: EUKulele and Diamond.

#### Taxonomic annotation with EUKulele

EUKulele can be run with different reference datasets.
The default dataset is PhyloDB (`--eukulele_db phylodb` ) which works for mixed communities of prokaryotes and eukaryotes.
Other database options for running the pipeline are MMETSP (`--eukulele_db mmetsp`; for marine protists) and GTDB (`--eukulele_db gtdb`; for prokarytes
[under development]).

Options:

- PhyloDB: default, covers both prokaryotes and eukaryotes
- MMETSP: marine protists
- GTDB: prokaryotes, both bacteria and archaea

You can also provide your own database, see the [EUKulele documentation](https://eukulele.readthedocs.io/en/latest/#) documentation.

Databases are automatically downloaded by the workflow, but if you already have them available you can use the `--eukulele_dbpath path/to/db` pointing
to the root directory of the EUKulele databases.
(The default for this parameter is `eukulele`.)

> Please, check the [EUKulele documentation](https://eukulele.readthedocs.io/en/latest/#) for more information about the databases.

#### Taxonomic annotation with Diamond

The Diamond taxonomy-annotation process uses Diamond database files (`.dmnd` files) that have been prepared with taxonomy information.
Currently we are only supplying a single standard databases, for GTDB release R09-RS220.
This is provided in collaboration with SciLifeLab Data Center and can be downloaded from here: [GTDB (R09RS220) taxonomy database](https://figshare.scilifelab.se/articles/dataset/nf-core_metatdenovo_taxonomy/28211678), DOI: https://doi.org/10.17044/scilifelab.28211678.
We hope to add more later.

To make your own database, you will need to collect four files: a protein fasta file, the `names.dmp` and `nodes.dmp` files from an
NCBI-style taxon dump plus a mapping file in which protein accessions are translated into taxon ids.

##### Building a database with nf-core/createtaxdb

The recommended way to build your own taxonomy-aware Diamond database is with [nf-core/createtaxdb](https://nf-co.re/createtaxdb), which wraps `diamond makedb` and its taxonomy inputs into a reproducible pipeline of its own.
Below is a worked example that builds and validates a database from [MarFERReT](https://zenodo.org/records/10170983) v1.1, a curated marine microbial eukaryote protein reference -- useful if your community has a substantial eukaryotic fraction not well represented in NCBI RefSeq/GTDB (see also [Coping with large datasets](large_datasets.md) and [issue #459](https://github.com/nf-core/metatdenovo/issues/459) for related eukaryote-focused work).

A minimal `samplesheet.csv`:

```csv title="samplesheet.csv"
id,taxid,fasta_dna,fasta_aa
MarFERReT_v1_1,1,,https://zenodo.org/records/10170983/files/MarFERReT.v1.proteins.faa.gz
```

`taxid` is a required samplesheet column but doesn't matter for Diamond -- the actual per-protein taxonomy comes from `--prot2taxid` below, not this column.
`fasta_dna` can be left empty; the schema only requires one of `fasta_dna`/`fasta_aa`.

A minimal `params.yml`:

```yaml title="params.yml"
input: samplesheet.csv
outdir: results
dbname: marferret_v1.1
build_diamond: true
prot2taxid: https://zenodo.org/records/10170983/files/MarFERReT.v1.taxonomies.tab.gz
nodesdmp: /path/to/nodes.dmp
namesdmp: /path/to/names.dmp
```

Then run:

```bash
nextflow run nf-core/createtaxdb -r dev -profile docker -params-file params.yml
```

A few gotchas worth knowing before you try this:

1. Zenodo's `/api/records/<id>/files/<name>/content` URLs fail `schema_input.json`'s filename regex (it must end in `.fasta`/`.fa`/`.faa`, optionally `.gz`, with nothing after).
   Use the plain `zenodo.org/records/<id>/files/<name>` form instead -- same content, correct extension at the end.
2. `nodesdmp`/`namesdmp` need pre-extracted `.dmp` files -- there's no tar.gz auto-extraction for these two, unlike the FASTA/`prot2taxid` inputs, which Nextflow happily streams straight from a URL.
   Download and extract an [NCBI taxonomy dump](ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz) first.
3. The 4-column accession2taxid-format file (accession, accession.version, taxid, gi) goes in via `--prot2taxid`, which is what actually reaches `DIAMOND_MAKEDB`'s `taxonmap` input.
4. `--taxonmap`/`--taxonnodes`/`--taxonnames` are build-time only.
   Querying the resulting `.dmnd` with `diamond blastp --outfmt 102` doesn't need them again, but also only emits `qseqid taxid evalue`, not a human-readable lineage string.
   The taxid-to-full-lineage expansion is a metatdenovo-side step (`parse_with_taxdump` in `diamond_dbs.csv`, see below), not something createtaxdb/Diamond itself produces.

Once you have a `.dmnd` plus the `names.dmp`/`nodes.dmp` you used to build it, wire them into `diamond_dbs.csv` exactly as described below -- createtaxdb's job ends at producing the database, not at configuring metatdenovo to use it.

##### Building an NCBI NR/RefSeq database with nf-core/createtaxdb

> [!WARNING]
> Unlike the MarFERReT example above, this one has **not** been run end-to-end -- it's the createtaxdb equivalent of the manual NCBI NR procedure in the next section below, built by direct analogy, not validated in practice.
> NR is a much larger database than MarFERReT (the manual procedure's own comment says the FASTA download alone "takes a loooong time"), so expect a correspondingly long build and high memory use if you try this.
> One thing this example can't confirm without actually running it: the manual procedure below pipes the FASTA through `sed '/^>/s/ .*//'` to strip descriptive text after each accession before `diamond makedb` sees it, because the taxonmap lookup needs to match on the bare accession.
> It's not verified here whether createtaxdb's own `DIAMOND_MAKEDB` wrapping does that same stripping internally or expects pre-cleaned input -- if the taxonmap join comes out empty or wrong, that header-cleaning step is the first thing to check.

Using the same NCBI sources as the manual procedure below, a `samplesheet.csv`:

```csv title="samplesheet.csv"
id,taxid,fasta_dna,fasta_aa
ncbi_nr,1,,ftp://ftp.ncbi.nih.gov/blast/db/FASTA/nr.gz
```

and `params.yml`:

```yaml title="params.yml"
input: samplesheet.csv
outdir: results
dbname: refseq
build_diamond: true
prot2taxid: ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.FULL.gz
nodesdmp: /path/to/nodes.dmp
namesdmp: /path/to/names.dmp
```

then, same as above:

```bash
nextflow run nf-core/createtaxdb -r dev -profile docker -params-file params.yml
```

`nodesdmp`/`namesdmp` still need pre-extracting from the [taxonomy dump](ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz) first, same as gotcha 2 above.
`dbname: refseq` matches the naming already used for this database in the `diamond_dbs.csv` example below, though note that NCBI NR and NCBI RefSeq are not literally the same underlying data -- this follows the manual procedure's own choice of source, not a claim that NR is a perfect stand-in for RefSeq specifically.

##### Building a database manually

Alternatively, you can run `diamond makedb` yourself without createtaxdb.
As an example, you can download the [NCBI NR database in fasta format](ftp://ftp.ncbi.nih.gov/blast/db/FASTA/nr.gz), the
[taxonomy dump](ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz) and the [mapping file](ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/).
(_Note_ that the mapping file is actually several files that need to be first downloaded, the concatenated into one.)
After this, you can run `diamond makedb` with the proper parameters to create your database.
Here's the whole procedure for the NCBI NR example (worked for us in January 2025; things might change):

```bash
# Make sure you have Diamond and wget installed

# Download the protein NR fasta file (takes a looong time)
wget ftp://ftp.ncbi.nih.gov/blast/db/FASTA/nr.gz

# Download and untar the taxonomy dump
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
tar xfz taxdump.tar.gz

# Download all the individual mapping files and concatenate them (also takes long)
wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.FULL.gz

# Create the Diamond database
gunzip -c nr.gz | sed '/^>/s/ .*//' | diamond makedb --taxonmap prot.accession2taxid.FULL.gz --taxonnames names.dmp --taxonnodes nodes.dmp --db ncbi-nr.taxonomy.dmnd
```

We are also, in collaboration with SciLifeLab Data Center, providing a [GTDB (R09RS220) taxonomy database](https://figshare.scilifelab.se/articles/dataset/nf-core*metatdenovo_taxonomy/28211678), DOI: https://doi.org/10.17044/scilifelab.28211678.

_Note_: If you can't download the files from FigShare from the command line with `wget` or `curl`, try URLs looking like
`https://ndownloader.figshare.com/files/52095806` instead of the ones you get from the web page.
Also note, that the files do not get the correct file names when you download from the command line.
Use a command like `wget -O gtdb-r220.names.dmp https://ndownloader.figshare.com/files/52095806`.

After creating one or more databases, you can provide them to the pipeline by filling out a file looking like the below and providing that to the pipeline
with `--diamond_dbs diamond_dbs.csv` (see the parameter documentation).
(The example is a `.csv`, but `.tsv`, `.json` and `.yml` are also supported.)

```csv title="diamond_dbs.csv"
db,dmnd_path,taxdump_names,taxdump_nodes,ranks,parse_with_taxdump
gtdb,diamond-taxonomy/gtdb_r220_repr.dmnd,diamond-taxonomy/gtdb_taxdump/names.dmp,diamond-taxonomy/gtdb_taxdump/nodes.dmp,domain;phylum;class;order;genus;species;strain,
refseq,diamond-taxonomy/refseq_protein.taxonomy.dmnd,diamond-taxonomy/ncbi_taxdump/names.dmp,diamond-taxonomy/ncbi_taxdump/nodes.dmp,,true
```

The `db`, `dmnd_path`, `taxdump_names` and `taxdump_nodes` fields are all required, while the remaining two are optional.
We strongly recommend that you download files for nf-core/metatdenovo to your own machine rather than specifying remote urls when running the pipeline since the files are large.
If you specify urls, the files will be downloaded by Nextflow on each run of the pipeline.

The pipeline will run three or four processes:

1. Align protein sequences of ORFs to the database with `diamond blastp`.
2. Add taxonomic lineage information with `taxonkit lineage`.
3. Format the output nicely and, if given a list of ranks in the `ranks` field of the db file, parse the taxonomy string into individual taxa.
4. If the `parse_with_taxdump` field is set to `true`, try to parse the taxonomy string using the `taxdump_names` and `taxdump_nodes` files.

Output from the first two steps will be saved in `results/diamond_taxonomy` whereas the second two goes to `results/summary_tables`.

_Note_ that some databases (e.g. from NCBI) output different numbers of taxa in their taxonomy strings.
These cannot be parsed by enumerating the ranks in `ranks`.
The only option is to set `parse_with_taxdump` to true.
This process is a bit unstable though.

### Functional annotation options

Besides the functional annotation that the gene caller Prokka gives (see above) there are two general purpose functional annotation
programs available in the workflow: the [eggNOG-mapper](http://eggnog-mapper.embl.de/) and
[KofamScan](https://github.com/takaram/kofam_scan).
Both are suitable for both prokaryotic and eukaryotic genes and both are run by default, but can be skipped using the `--skip_eggnog` and
`--skip_kofamscan` options respectivelly.
The tools use large databases which are downloaded automatically but paths can be provided by the user through the `--eggnog_dbpath directory`
and `--kofam_dir dir` parameters respectively.
It is practical to let the pipeline download the files on the first run, and then reuse the data by setting the parameters.

:::note
Currently, the standard download procedure for the eggNOG database using the `download_eggnog_data.py` tool (v.2.1.9) doesn't work because the domain it tries to download from doesn't exist.
Since release 1.4.0, this pipeline therefore uses `wget` to fetch files from [the current download site](http://eggnog6.embl.de/download/emapperdb-5.0.2).
:::

A third functional annotation option is CAZyme annotation using [dbCAN](https://bcb.unl.edu/dbCAN2/) (`run_dbcan`), which is also run by
default and can be skipped with `--skip_dbcan`.
Its database is downloaded automatically, with the path settable through `--dbcan_dbpath directory` (the same "let the pipeline download it
once, then reuse the data" advice above applies here too).
Since the pipeline only runs dbCAN's protein-mode CAZyme annotation (not its gene-cluster/CGC analysis), the download skips the CGC-related
database assets to save space and time.

A more targeted annotation option offered by the workflow is the possibility for the user to provide a set of
[HMMER HMM profiles](http://eddylab.org/software/hmmer/Userguide.pdf) through the `--hmmdir dir` or `hmmfiles file0.hmm,file1.hmm,...,filen.hmm`
parameters.
Each HMM file will be used to search the amino acid sequences of the ORF set and the results will be summarized in a tab separated file in
which each ORF-HMM combination will be ranked according to score and E-value.

#### How to manually download the databases for functional annotation

There are some cases (e.g. offline run) where you prefer to download the databases before running the pipeline. Currently, `eggnog-mapper`, `kofamscan` and `dbcan` use databases that can be downloaded.

##### Eggnog databases

For `eggnog-mapper` the easiest way is to use `download_eggnog_data.py` provided when you install eggnog-mapper locally (documentation
[here](https://github.com/eggnogdb/eggnog-mapper/wiki/eggNOG-mapper-v2.1.5-to-v2.1.12#user-content-Installation)).

First, install eggnog-mapper:

```bash
conda install -c bioconda -c conda-forge eggnog-mapper
```

Then, you can download all databases available

```bash
 download_eggnog_data.py
```

You can select which database you want to download (read eggnog-mapper docs) but you need to be sure you will store them in a directory that
will be called with the option `--eggnog_dbpath`

##### Kofamscan databases

You can use `wget` to download the file in a new directory that will be used with `--kofamscan_dbpath`

```bash
wget https://www.genome.jp/ftp/db/kofam/ko_list.gz
gunzip ko_list.gz

wget https://www.genome.jp/ftp/db/kofam/profiles.tar.gz
tar -zxf profiles.tar.gz
```

##### dbCAN database

You can use `run_dbcan` itself (installed alongside the pipeline's dbCAN container, or via `pip install run-dbcan`) to download the database
into a directory that will be used with `--dbcan_dbpath`.
The `--no-cgc` flag matches what the pipeline itself uses, since only protein-mode CAZyme annotation is run, not gene-cluster analysis:

```bash
run_dbcan database --db_dir dbcan --aws_s3 --no-cgc
```

## Example pipeline command with some common features

```bash
nextflow run nf-core/metatdenovo -profile docker --input samplesheet.csv --assembler spades --orf_caller prokka --eggnog --eukulele_db gtdb
```

In this example, we are running metatdenovo with `spades` as assembler, `prokka` as ORF caller, `eggnog` for functional annotation and EUKulele with the GTDB database for taxonomic annotation.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these
in a parameter `yml` or `json` file and specify this to the pipeline with `-params-file params.yml`.

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/running/run-pipelines#configuring-pipelines), other infrastructural tweaks (such as output directories), or module arguments (args).

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run nf-core/metatdenovo -profile docker -params-file params.yaml
```

with:

```yaml title="params.yaml"
input: 'samplesheet.csv'
assembler: 'spades'
orf_caller: 'prokka'
eggnog: true
eukulele_db: 'gtdb'
<...>
```

You can also generate such `YAML` or `JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/metatdenovo
```

### Reproducibility

It is a good idea to specify the pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/metatdenovo releases page](https://github.com/nf-core/metatdenovo/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducibility, you can use share and reuse [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> [!TIP]
> If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen)

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> [!IMPORTANT]
> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to check if your system is supported, please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer environment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://charliecloud.io/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow ` 24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the pipeline steps, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher resources request (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#set-max-resources) and [customise process resources](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#customize-process-resources) section of the nf-core website.

### Custom Containers

In some cases, you may wish to change the container or conda environment used by a pipeline steps for a particular tool. By default, nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However, in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#update-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/running/configuration/nextflow-for-your-system#modifying-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
