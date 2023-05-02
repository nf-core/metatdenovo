/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowMetatdenovo.initialise(params, log)

// Validate parameters for orf_caller:
ORF_CALLER_PRODIGAL     = 'prodigal'
ORF_CALLER_PROKKA       = 'prokka'

// Validate parameters for assembler:
MEGAHIT   = 'megahit'

def valid_params = [
    orf_caller      : [ ORF_CALLER_PRODIGAL, ORF_CALLER_PROKKA ],
    assembler       : [ MEGAHIT ]
]

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

// set an empty multiqc channel
ch_multiqc_files = Channel.empty()

// If the user supplied hmm files, we will run hmmsearch and then rank the results.
// Create a channel for hmm files.
ch_hmmrs = Channel.empty()
if ( params.hmmdir ) {
    Channel
        .fromPath(params.hmmdir + params.hmmpattern)
        .set { ch_hmmrs }
} else if ( params.hmmfiles ) {
    Channel
        .of( params.hmmfiles.split(',') )
        .map { [ file(it) ] }
        .set { ch_hmmrs }
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: local
//
include { MEGAHIT_INTERLEAVED              } from '../modules/local/megahit/interleaved'
include { COLLECT_FEATURECOUNTS            } from '../modules/local/collect_featurecounts'
include { COLLECT_STATS                    } from '../modules/local/collect_stats'
include { UNPIGZ as UNPIGZ_CONTIGS         } from '../modules/local/unpigz'
include { UNPIGZ as UNPIGZ_GFF             } from '../modules/local/unpigz'
include { MERGE_TABLES                     } from '../modules/local/merge_summary_tables'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'

//
// SUBWORKFLOW: Consisting of local modules
//
include { EGGNOG            } from '../subworkflows/local/eggnog'
include { HMMCLASSIFY       } from '../subworkflows/local/hmmclassify'
include { PROKKA_SUBSETS    } from '../subworkflows/local/prokka_subsets'
include { DIGINORM          } from '../subworkflows/local/diginorm'
include { FASTQC_TRIMGALORE } from '../subworkflows/local/fastqc_trimgalore'
include { PRODIGAL          } from '../subworkflows/local/prodigal'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { BBMAP_INDEX                                } from '../modules/nf-core/bbmap/index/main'
include { BBMAP_ALIGN                                } from '../modules/nf-core/bbmap/align/main'
include { BBMAP_BBNORM                               } from '../modules/nf-core/bbmap/bbnorm/main'
include { SUBREAD_FEATURECOUNTS as FEATURECOUNTS_CDS } from '../modules/nf-core/subread/featurecounts/main'
include { CAT_FASTQ 	          	                 } from '../modules/nf-core/cat/fastq/main'
include { FASTQC                                     } from '../modules/nf-core/fastqc/main'
include { MULTIQC                                    } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS                } from '../modules/nf-core/custom/dumpsoftwareversions/main'

//
// SUBWORKFLOWS: Installed directly from nf-core/modules
//
include { BAM_SORT_STATS_SAMTOOLS                    } from '../subworkflows/nf-core/bam_sort_stats_samtools/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow METATDENOVO {

    ch_versions = Channel.empty()

    // STEP 0: Validate input
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    .reads
    .map {
        meta, fastq ->
            new_id = meta.id - ~/_T\d+/
            [ meta + [id: new_id], fastq ]
    }
    .groupTuple()
    .branch {
        meta, fastq ->
            single  : fastq.size() == 1
                return [ meta, fastq.flatten() ]
            multiple: fastq.size() > 1
                return [ meta, fastq.flatten() ]
    }
    .set { ch_fastq }
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    // 
    // MODULE: Concatenate FastQ files from same sample if required
    //
    CAT_FASTQ (
        ch_fastq.multiple
    )
    .reads
    .mix(ch_fastq.single)
    .set { ch_cat_fastq }

    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first().ifEmpty(null))

    // Step 1 + Step 3: FastQC and Trim Galore!
    // SUBWORKFLOW: Read QC and trim adapters
    //
    FASTQC_TRIMGALORE (
        ch_cat_fastq,
        params.skip_fastqc || params.skip_qc,
        params.skip_trimming
    )
    ch_versions = ch_versions.mix(FASTQC_TRIMGALORE.out.versions)

    ch_collect_stats = ch_cat_fastq.collect { it[0].id }.map { [ [ id:"${params.assembler}.${params.orf_caller}" ], it ] }
    if ( params.skip_trimming ) {
        ch_collect_stats
            .map { [ it[0], it[1], [] ] }
            .set { ch_collect_stats }
    } else {
        ch_collect_stats
            .combine(FASTQC_TRIMGALORE.out.trim_log.collect { it[1][0] }.map { [ it ] })
            .set { ch_collect_stats }
    }

    // Step 4
    // Remove host sequences, bowtie2 align to Bos taurus
    //

    // Step 5
    // rRNA remove (sortmerna)
    //

    // Step 6
    // Deduplication with BBdup
    //

    // Step 7
    // Filter by taxa with Kraken2
    //

    // Step 8
    // SUBWORKFLOW: Perform digital normalization. There are two options: khmer or BBnorm. The latter is faster.
    //
    ch_pe_reads_to_assembly = Channel.empty()
    ch_se_reads_to_assembly = Channel.empty()

    if ( ! params.assembly ) {
        if ( params.diginorm ) {
            DIGINORM(ch_interleaved.collect { meta, fastq -> fastq }, [], 'all_samples')
            ch_versions = ch_versions.mix(DIGINORM.out.versions)
            ch_pe_reads_to_assembly = DIGINORM.out.pairs
            ch_se_reads_to_assembly = DIGINORM.out.singles
        } else if ( params.bbnorm) {
                BBMAP_BBNORM(ch_interleaved.collect { meta, fastq -> fastq }.map {[ [id:'all_samples', single_end:true], it ] } )
                ch_pe_reads_to_assembly = BBMAP_BBNORM.out.fastq.map { it[1] }
                ch_se_reads_to_assembly = Channel.empty()
        } else {
            ch_pe_reads_to_assembly = ch_interleaved.map { meta, fastq -> fastq }
            ch_se_reads_to_assembly = Channel.empty()
        }
    }

    // Step 9
    // MODULE: Run Megahit on all interleaved fastq files
    //
    if ( params.assembly ) {
        Channel
            .value ( [ [ id: 'user_assembly' ], file(params.assembly) ] )
            .set { ch_assembly_contigs }
    } else if ( params.assembler == MEGAHIT ) {
        MEGAHIT_INTERLEAVED(
            ch_pe_reads_to_assembly.collect().ifEmpty([]),
            ch_se_reads_to_assembly.collect().ifEmpty([]),
            'megahit_assembly'
        )
        MEGAHIT_INTERLEAVED.out.contigs
            .map { [ [ id: 'megahit' ], it ] }
            .set { ch_assembly_contigs }
        ch_versions = ch_versions.mix(MEGAHIT_INTERLEAVED.out.versions)
    }

    // Step 10
    // Clustering with CD-HIT-EST
    // (for concatenating multiple assemblies)
    //

    // Step 11
    // Call ORFs
    //
    ch_gff = Channel.empty()
    ch_aa  = Channel.empty()

    //
    // SUBWORKFLOW: Run PROKKA_SUBSETS on assmebly output, but split the fasta file in chunks of 10 MB, then concatenate and compress output.
    //
    if ( params.orf_caller == ORF_CALLER_PROKKA ) {
        PROKKA_SUBSETS(ch_assembly_contigs)
        UNPIGZ_GFF(PROKKA_SUBSETS.out.gff.map { [ [id: "${params.orf_caller}.${it[0].id}"], it[1] ] })
        ch_versions      = ch_versions.mix(PROKKA_SUBSETS.out.versions)
        ch_gff           = UNPIGZ_GFF.out.unzipped
        ch_aa            = PROKKA_SUBSETS.out.faa
        ch_multiqc_files = ch_multiqc_files.mix(PROKKA_SUBSETS.out.prokka_log.collect{it[1]}.ifEmpty([]))
    }

    //
    // MODULE: Run PRODIGAL on assembly output.
    //
    if ( params.orf_caller == ORF_CALLER_PRODIGAL ) {
        PRODIGAL( ch_assembly_contigs.map { [ [id: 'prodigal.megahit' ], it[1] ] } )
        UNPIGZ_GFF(PRODIGAL.out.gff.map { [ [id: "${params.orf_caller}.${it[0].id}"], it[1] ] })
        ch_gff          = UNPIGZ_GFF.out.unzipped
        ch_aa           = PRODIGAL.out.faa
        ch_versions     = ch_versions.mix(PRODIGAL.out.versions)
    }

    //
    // MODULE: Create a BBMap index
    //
    BBMAP_INDEX(ch_assembly_contigs.map { it[1] })
    ch_versions   = ch_versions.mix(BBMAP_INDEX.out.versions)

    //
    // MODULE: Call BBMap with the index once per sample
    //
    BBMAP_ALIGN ( ch_clean_reads, BBMAP_INDEX.out.index )
    ch_versions = ch_versions.mix(BBMAP_ALIGN.out.versions)

    // Step 12
    // Quantification (salmon, rsem, or bbmap)
    //

    // Step 13
    // SUBWORKFLOW: classify ORFs with a set of hmm files
    //
    ch_hmmrs
        .combine(ch_aa)
        .map { [ [id: it[0].baseName ], it[0], it[2] ] }
        .set { ch_hmmclassify }
    HMMCLASSIFY ( ch_hmmclassify )
    ch_versions = ch_versions.mix(HMMCLASSIFY.out.versions)

    //
    // MODULE: FeatureCounts. Create a table for each samples that provides raw counts as result of the alignment.
    //

    BAM_SORT_STATS_SAMTOOLS ( BBMAP_ALIGN.out.bam, ch_assembly_contigs.map { it[1] } )
    ch_versions = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS.out.versions)

    // if ( orf_caller == 
    BAM_SORT_STATS_SAMTOOLS.out.bam
        .combine(ch_gff.map { it[1] } )
        .set { ch_featurecounts }

    ch_collect_stats
        .combine(BAM_SORT_STATS_SAMTOOLS.out.idxstats.collect { it[1]}.map { [ it ] })
        .set { ch_collect_stats }

    FEATURECOUNTS_CDS ( ch_featurecounts)
    ch_versions       = ch_versions.mix(FEATURECOUNTS_CDS.out.versions)

    //
    // MODULE: Collect featurecounts output counts in one table
    //
    FEATURECOUNTS_CDS.out.counts
        .collect() { it[1] }
        .map { [ [ id:'all_samples'], it ] }
        .set { ch_collect_feature }

    COLLECT_FEATURECOUNTS ( ch_collect_feature )
    ch_versions           = ch_versions.mix(COLLECT_FEATURECOUNTS.out.versions)
    ch_fcs_for_stats      = COLLECT_FEATURECOUNTS.out.counts.collect { it[1]}.map { [ it ] }
    ch_fcs_for_summary    = COLLECT_FEATURECOUNTS.out.counts.map { it[1]}
    ch_collect_stats
        .combine(ch_fcs_for_stats)
        .set { ch_collect_stats }

    //
    // SUBWORKFLOW: run eggnog_mapper on the ORF-called amino acid sequences
    //
    if ( ! params.skip_eggnog ) {
        EGGNOG(ch_aa, ch_fcs_for_summary )
        ch_versions = ch_versions.mix(EGGNOG.out.versions)
        ch_merge_tables = EGGNOG.out.sumtable
    } else {
        ch_aa
            .map { [ it[0], [] ] }
            .set { ch_merge_tables }
    }

    // Step 14
    // Kraken2 taxonomical annotation
    //

    // Step 15
    // MODULE: Collect statistics from mapping analysis
    //
    if( !params.skip_eggnog  || !params.skip_eukulele ) {
        MERGE_TABLES ( ch_merge_tables )
        ch_collect_stats
            .combine(MERGE_TABLES.out.merged_table.collect{ it[1]}.map { [ it ] })
            .set { ch_collect_stats }
    } else {
        ch_collect_stats
            .map { [ it[0], it[1], it[2], it[3], it[4], it[5], [] ] }
            .set { ch_collect_stats }
    }

    COLLECT_STATS(ch_collect_stats)
    ch_versions     = ch_versions.mix(COLLECT_STATS.out.versions)

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    // Step 2
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowMetatdenovo.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowMetatdenovo.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_TRIMGALORE.out.trim_zip.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(BAM_SORT_STATS_SAMTOOLS.out.idxstats.collect{it[1]}.ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(FEATURECOUNTS_CDS.out.summary.collect{it[1]}.ifEmpty([]))


    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
