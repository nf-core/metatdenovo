# ![nf-core/metatdenovo](docs/images/nf-core-metatdenovo_logo_light.png#gh-light-mode-only) ![nf-core/metatdenovo](docs/images/nf-core-metatdenovo_logo_dark.png#gh-dark-mode-only)

[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/metatdenovo/results)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A522.10.1-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Nextflow Tower](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Nextflow%20Tower-%234256e7)](https://tower.nf/launch?pipeline=https://github.com/nf-core/metatdenovo)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23metatdenovo-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/metatdenovo)[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/nf_core)[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

**nf-core/metatdenovo** is a bioinformatics best-practice analysis pipeline for Assembly and annotation of metatranscriptomic data, both prokaryotic and eukaryotic.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

On release, automated continuous integration tests run the pipeline on a full-sized dataset on the AWS cloud infrastructure. This ensures that the pipeline runs on AWS, has sensible resource allocation defaults set to run on real-world datasets, and permits the persistent storage of results to benchmark between pipeline releases and other analysis sources.The results obtained from the full-sized test can be viewed on the [nf-core website](https://nf-co.re/metatdenovo/results).

## Pipeline summary

![nf-core/rnaseq metro map](docs/images/nf-core-metatdenovo_metro_map.png)

1. Read QC ([`FastQC`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/))
2. Present QC for raw reads ([`MultiQC`](http://multiqc.info/))
3. Quality trimming and adapters removal for raw reads [`Trimm Galore!`](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)
4. Filter sequences with [`BBduk`](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/)
5. Normalize the sequences with KHMER or [`BBnorm`](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbnorm-guide/)
6. Merge trimmed, pair-end reads ( [`Seqtk`](https://github.com/lh3/seqtk))
7. Choice of de novo assembly:
   1. [`RNAspade`](https://cab.spbu.ru/software/rnaspades/) suggested for Eukaryotes de novo assembly
   2. [`Megahit`](https://github.com/voutcn/megahit) suggested for Prokaryotes de novo assembly
8. Choice of orf caller:
   1. [`TransDecoder`](https://github.com/TransDecoder/TransDecoder) suggested for Eukaryotes
   2. [`Prokka`](https://github.com/tseemann/prokka) suggested for Prokaryotes
   3. [`Prodigal`](https://github.com/hyattpd/Prodigal) suggested for Prokaryotes
9. Quantification of genes identified in assemblies:
   1. generate index of assembly [`BBmap index`](https://sourceforge.net/projects/bbmap/)
   2. Mapping cleaned reads to the assembly for quantification [`BBmap`](https://sourceforge.net/projects/bbmap/)
   3. Get raw counts per each gene present in the assembly [`Featurecounts`](http://subread.sourceforge.net) -> TSV table with collected featurecounts output
10. Functional annotation:
   1. [`Eggnog`](https://github.com/eggnogdb/eggnog-mapper) -> Reformat TSV output "eggnog table"
   2. [`KOfamscan`](https://github.com/takaram/kofam_scan)
   3. [`HMMERsearch`](https://www.ebi.ac.uk/Tools/hmmer/search/hmmsearch) -> Ranking orfs based on HMMprofile with [`Hmmrank`](https://github.com/erikrikarddaniel/hmmrank)
11. Taxonomical annotation:
   1. [`EUKulele`](https://github.com/AlexanderLabWHOI/EUKulele) -> Reformat TSV output "Reformat_tax.R"
   2. [`CAT`](https://github.com/dutilh/CAT)
12. Choice of functional annotation: 1. [`Eggnog-mapper`](http://eggnog-mapper.embl.de) 2. [`Run-DBcan`](https://github.com/linnabrown/run_dbcan) 3. [`Hmmsearch`](https://www.ebi.ac.uk/Tools/hmmer/search/hmmsearch). Besides searching the ORFs, each ORF's hits will be ranked.
   10 Summary statistics table. Collect_stats.R

## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=22.10.1`)

2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) (you can follow [this tutorial](https://singularity-tutorial.github.io/01-installation/)), [`Podman`](https://podman.io/), [`Shifter`](https://nersc.gitlab.io/development/shifter/how-to-use/) or [`Charliecloud`](https://hpc.github.io/charliecloud/) for full pipeline reproducibility _(you can use [`Conda`](https://conda.io/miniconda.html) both to install Nextflow itself and also to manage software within pipelines. Please only use it within pipelines as a last resort; see [docs](https://nf-co.re/usage/configuration#basic-configuration-profiles))_.

3. Download the pipeline and test it on a minimal dataset with a single command:

   ```bash
   nextflow run nf-core/metatdenovo -profile test,YOURPROFILE --outdir <OUTDIR>
   ```

   Note that some form of configuration will be needed so that Nextflow knows how to fetch the required software. This is usually done in the form of a config profile (`YOURPROFILE` in the example command above). You can chain multiple config profiles in a comma-separated string.

   > - The pipeline comes with config profiles called `docker`, `singularity`, `podman`, `shifter`, `charliecloud` and `conda` which instruct the pipeline to use the named tool for software management. For example, `-profile test,docker`.
   > - Please check [nf-core/configs](https://github.com/nf-core/configs#documentation) to see if a custom config file to run nf-core pipelines already exists for your Institute. If so, you can simply use `-profile <institute>` in your command. This will enable either `docker` or `singularity` and set the appropriate execution settings for your local compute environment.
   > - If you are using `singularity`, please use the [`nf-core download`](https://nf-co.re/tools/#downloading-pipelines-for-offline-use) command to download images first, before running the pipeline. Setting the [`NXF_SINGULARITY_CACHEDIR` or `singularity.cacheDir`](https://www.nextflow.io/docs/latest/singularity.html?#singularity-docker-hub) Nextflow options enables you to store and re-use the images from a central location for future pipeline runs.
   > - If you are using `conda`, it is highly recommended to use the [`NXF_CONDA_CACHEDIR` or `conda.cacheDir`](https://www.nextflow.io/docs/latest/conda.html) settings to store the environments in a central location for future pipeline runs.

4. Start running your own analysis!

   <!-- TODO nf-core: Update the example "typical command" below used to run the pipeline -->

   ```bash
   nextflow run nf-core/metatdenovo --input samplesheet.csv --outdir <OUTDIR> --genome GRCh37 -profile <docker/singularity/podman/shifter/charliecloud/conda/institute>
   ```

## Documentation

The nf-core/metatdenovo pipeline comes with documentation about the pipeline [usage](https://github.com/LNUc-EEMiS/metatdenovo/blob/dev/docs/usage.md), [parameters](https://nf-co.re/metatdenovo/parameters) and [output](https://nf-co.re/metatdenovo/output).

## Credits

nf-core/metatdenovo was originally written by Danilo Di Leo (@danilodileo), Emelie Nilsson (@emnilsson) & Daniel Lundin (@erikrikarddaniel).

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#metatdenovo` channel](https://nfcore.slack.com/channels/metatdenovo) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

<!-- If you use  nf-core/metatdenovo for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
