# nf-core/metatdenovo: Usage

## :warning: Please read this documentation on the nf-core website: [https://nf-co.re/metatdenovo/usage] (the link is not working) (https://nf-co.re/metatdenovo/usage)

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

## Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. It must be a comma-separated file with 3 columns, and a header row as shown in the examples below

```bash
--input '[path to samplesheet file]'
```

### Multiple runs of the same sample

The `sample` identifiers have to be the same when you have re-sequenced the same sample more than once e.g. to increase sequencing depth. The pipeline will concatenate the raw reads before performing any downstream analysis. Below is an example for the same sample sequenced across 3 lanes:

```console
sample,fastq_1,fastq_2
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
CONTROL_REP1,AEG588A1_S1_L003_R1_001.fastq.gz,AEG588A1_S1_L003_R2_001.fastq.gz
CONTROL_REP1,AEG588A1_S1_L004_R1_001.fastq.gz,AEG588A1_S1_L004_R2_001.fastq.gz
```

### Full samplesheet

The pipeline will auto-detect whether a sample is single- or paired-end using the information provided in the samplesheet. The samplesheet can have as many columns as you desire, however, there is a strict requirement for the first 3 columns to match those defined in the table below.

A final samplesheet file consisting of both single- and paired-end data may look something like the one below. This is for 6 samples, where `TREATMENT_REP3` has been sequenced twice.

```console
sample,fastq_1,fastq_2
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
CONTROL_REP2,AEG588A2_S2_L002_R1_001.fastq.gz,AEG588A2_S2_L002_R2_001.fastq.gz
CONTROL_REP3,AEG588A3_S3_L002_R1_001.fastq.gz,AEG588A3_S3_L002_R2_001.fastq.gz
TREATMENT_REP1,AEG588A4_S4_L003_R1_001.fastq.gz,
TREATMENT_REP2,AEG588A5_S5_L003_R1_001.fastq.gz,
TREATMENT_REP3,AEG588A6_S6_L003_R1_001.fastq.gz,
TREATMENT_REP3,AEG588A6_S6_L004_R1_001.fastq.gz,
```

| Column    | Description                                                                                                                                                                            |
| --------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`  | Custom sample name. This entry will be identical for multiple sequencing libraries/runs from the same sample. Spaces in sample names are automatically converted to underscores (`_`). |
| `fastq_1` | Full path to FastQ file for Illumina short reads 1. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                             |
| `fastq_2` | Full path to FastQ file for Illumina short reads 2. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                             |

An [example samplesheet](../assets/samplesheet.csv) has been provided with the pipeline.

## Filter/remove sequences from the samples (e.g. rRNA sequences with SILVA database)

The pipeline can remove potential contaminants using the BBduk program.
Specify a fasta file, gzipped or not, with the --sequence_filter sequences.fasta parameter.
For further documentation, see the [BBduk official website](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/).

## Digital normalization

Metatdenovo can perform "digital normalization" on the reads BEFORE the assembly.
This will reduce coverage of highly abundant sequences and remove sequences that are below a threshold, and can be useful if the data set is too large to assemble but also potentially improve an assembly.
N.B. the digital normalization is done only for the assembly and the non-normalized sequences will be used for quantification.
There is one option for digital normalization in the pipeline:

- bbnorm (`--bbnorm`)

> Please, check the [bbnorm](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbnorm-guide/) documentation for further information about these programs and how digital normalization works. Remember to check [Parameters](https://nf-co.re/metatdenovo/parameters) page for the full option that can be used for this step

## Assembler options

By default, the pipeline uses Megahit (i.e. `--assembler megahit`) to assemble the cleaned and trimmed FastQ reads to create the reference genome.
Megahit is fast and it requires a not a lot of memory to run, typically is suggested to be used with prokaryotic samples.
The pipeline allows you to choose another assembler RNAspades, (i.e. `--assembler rnaspades` ), that is usually suggested to use for eukaryotic samples.
You can also choose to input contigs from an assembly that you made outside the pipeline using the `--assembly file.fna` (where `file.fna` is the name of a fasta file with contigs) option.

> N.B. you can use `Megahit` for eukaryotic samples too, we just suggest what is the best option according to our experience (literature?).

## Orf caller options

By default, the pipeline uses prodigal (i.e. `--orf_caller prodigal` ) to generate the genome feature file (.gff) and to generate gene structure from the assembly.

Other orf caller options for running the pipeline are:

- Prokka (`--orf_caller prokka`)

- Transdecoder (`--orf_caller transdecoder`)

> N.B. Prokka and prodigal are suggested to run with prokaryotes while transdecoder is specific for eukaryotes.

## Taxonomical annotation options

Metatdenovo uses `EUKulele` as the main program for taxonomy annotation. `EUKulele` can be run with different reference datasets. The default dataset is PhyloDB (i.e. `--eukulele_db phylodb` ) which works for mixed communities of prokaryotes and eukaryotes.

Other databases options for running the pipeline are:

- MMETSP (`--eukulele_db mmetsp`)

- GTDB (`--eukulele_db gtdb`) [under development]

PhyloDB and GTDB are recommended for prokaryotic datasets and MMETSP for eukaryotes, although PhyoDB can be also recognize eukaryotes and can be used for this purpose.

If you already have these databases ready in your working directory, you can point to the folder so the pipeline will not download the database (e.g. `--eukulele_dbpath your/path/database/`). N.B. When you are using a custom database, don't specify the `--eukulele_db` option. The pipeline will provide a default name for the database to avoid that EUKulele will try to download a new database.

> Please, check the [EUKulele documentation](https://eukulele.readthedocs.io/en/latest/#) for more information about the databases.

An alternative to EUKulele is the CAT program. In contrast to EUKulele that annotates open reading frames (ORFs), CAT annotates the contigs from the assembly.

CAT is uses Prodigal to call ORFs and DIAMOND for the alignment to a reference database. Subsequently, DIAMOND hits for individual ORFs are translated by CAT into contig annotations.

The database can be generated with the option `--cat_db_generate` or you can provide a prepared database that you downloaded from [CAT website](https://tbb.bio.uu.nl/bastiaan/CAT_prepare/).
Check the also the [options]() documentation to learn how to configure CATproperly.

> Please, check the `CAT` documentation for more information about the database cited [HERE](https://github.com/dutilh/CAT)

## Functional annotation options

By default, the metatdenovo pipeline will perform a functional annotation with the [eggNOG-mapper](http://eggnog-mapper.embl.de/) program. In order to run the other programs, you will need to specify them as additional options.

These options are:

- [Eggnog](https://github.com/eggnogdb/eggnog-mapper/wiki) (`--eggnog_dbpath`)

- [hmmsearch](http://eddylab.org/software/hmmer/Userguide.pdf) (`--hmmdir` or `--hmmfiles`)

- [kofamscan](https://github.com/takaram/kofam_scan) (`--kofam_dir`)

All the options can run in the same time (e.g. `nextflow run main.nf -profile test,docker --eggnog --hmmdir hmms/ `) but each program has its own options that you will need to read carefully before running the pipeline.
You can find more information about the different options in the [parameters page](https://nf-co.re/metatdenovo/parameters).
For details about individual programs used, see their respective home pages.

If an Eggnog, kofam, or EUKulele database is already available, they can be specified with the above commands to skip the automatic download that the pipeline performs.

If you don't want run eggNOG-mapper, you will need to add the flag `--skip_eggnog`, otherwise metatdenovo will run the program automatically.

## Example pipeline command with some common features

`nextflow run lnuc-eemis/metatdenovo -profile docker --input samplesheet.csv --assembler rnaspades --orf_caller transdecoder --eggnog`

In this example, we are running metatdenovo with `rnaspades` as assembler, `transdecoder` as ORF caller and `eggnog` for functional annotation.

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run nf-core/metatdenovo --input ./samplesheet.csv --outdir ./results -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

:::warning
Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).
:::

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run nf-core/metatdenovo -profile docker -params-file params.yaml
```

with `params.yaml` containing:

```yaml
input: './samplesheet.csv'
outdir: './results/'
genome: 'GRCh37'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull nf-core/metatdenovo
```

### Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/metatdenovo releases page](https://github.com/nf-core/metatdenovo/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. For example, at the bottom of the MultiQC reports.

To further assist in reproducbility, you can use share and re-use [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

:::tip
If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.
:::

## Core Nextflow arguments

:::note
These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).
:::

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

:::info
We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.
:::

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer enviroment.

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
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher requests (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Containers

In some cases you may wish to change which container or conda environment a step of the pipeline uses for a particular tool. By default nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Azure Resource Requests

To be used with the `azurebatch` profile by specifying the `-profile azurebatch`.
We recommend providing a compute `params.vm_type` of `Standard_D16_v3` VMs by default but these options can be changed if required.

Note that the choice of VM size depends on your quota and the overall workload during the analysis.
For a thorough list, please refer the [Azure Sizes for virtual machines in Azure](https://docs.microsoft.com/en-us/azure/virtual-machines/sizes).

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
