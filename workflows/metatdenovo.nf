/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowMetatdenovo.initialise(params, log)

// Validate parameters Orfcaller:
def valid_params = [
    orf_caller  : ['prodigal']
]

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
========================================================================================
    CONFIG FILES
========================================================================================
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

include { MEGAHIT_INTERLEAVED } from '../modules/local/megahit/interleaved.nf' addParams( options: modules['megahit'] )
include { UNPIGZ as UNPIGZ_MEGAHIT_CONTIGS } from '../modules/local/unpigz.nf' addParams( options: modules['unpigz_megahit_contigs'])

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check' addParams( options: [:] )

//
// SUBWORKFLOW: Adapted from rnaseq!
//
def trimgalore_options    = modules['trimgalore']
trimgalore_options.args  += params.trim_nextseq > 0 ? Utils.joinModuleArgs(["--nextseq ${params.trim_nextseq}"]) : ''

include { FASTQC_TRIMGALORE } from '../subworkflows/local/fastqc_trimgalore' addParams( fastqc_options: modules['fastqc'], trimgalore_options: trimgalore_options )

//
// SUBWORKFLOW: Perform digital normalization
//
def diginorm_normalizebymedian_options    = modules['diginorm_normalizebymedian']
diginorm_normalizebymedian_options.args  += Utils.joinModuleArgs(["-C ${params.diginorm_C} -k ${params.diginorm_k}"])
def diginorm_filterabund_options          = modules['diginorm_filterabund']
diginorm_filterabund_options.args        += Utils.joinModuleArgs(["-C ${params.diginorm_C} -Z ${params.diginorm_C}"])
def diginorm_extractpairedreads_options   = modules['diginorm_extractpairedreads']

include { DIGINORM } from '../subworkflows/local/diginorm' addParams(
    diginorm_normalizebymedian_options:  diginorm_normalizebymedian_options, 
    diginorm_filterabund_options:        diginorm_filterabund_options,
    diginorm_extractpairedreads_options: diginorm_extractpairedreads_options
)

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

def multiqc_options   = modules['multiqc']
multiqc_options.args += params.multiqc_title ? Utils.joinModuleArgs(["--title \"$params.multiqc_title\""]) : ''

def prodigal_options   = modules['prodigal']
prodigal_options.args += params.prodigal_trainingfile ? Utils.joinModuleArgs("-t $params.prodigal_trainingfile") : ""

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC        } from '../modules/nf-core/modules/fastqc/main'        addParams( options: modules['fastqc'] )
include { BBMAP_BBDUK   } from '../modules/nf-core/modules/bbmap/bbduk/main'   addParams( options: modules['bbduk'] )
include { BBMAP_INDEX   } from '../modules/nf-core/modules/bbmap/index/main'   addParams( options: modules['bbmap_index'] )
include { BBMAP_ALIGN   } from '../modules/nf-core/modules/bbmap/align/main'   addParams( options: modules['bbmap_align'] )
include { SEQTK_MERGEPE } from '../modules/nf-core/modules/seqtk/mergepe/main' addParams( options: modules['seqtk_mergepe'] )
include { PROKKA        } from '../modules/nf-core/modules/prokka/main' addParams( options: modules['prokka'] )
include { PRODIGAL      } from '../modules/nf-core/modules/prodigal/main' addParams( options: prodigal_options )
include { MULTIQC       } from '../modules/nf-core/modules/multiqc/main' addParams( options: multiqc_options   )
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/modules/custom/dumpsoftwareversions/main'  addParams( options: [publish_files : ['_versions.yml':'']] )

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

// Info required for completion email and summary
def multiqc_report = []

workflow METATDENOVO {

    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        ch_input
    )
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    //
    // SUBWORKFLOW: Read QC and trim adapters
    //
    FASTQC_TRIMGALORE (
        INPUT_CHECK.out.reads,
        params.skip_fastqc || params.skip_qc,
        params.skip_trimming
    )
    ch_versions = ch_versions.mix(FASTQC_TRIMGALORE.out.versions)

    //
    // MODULE: Run BBDuk to clean out whatever sequences the user supplied via params.sequence_filter
    //
    if ( params.sequence_filter ) {
        BBMAP_BBDUK ( FASTQC_TRIMGALORE.out.reads, params.sequence_filter )
        ch_clean_reads  = BBMAP_BBDUK.out.reads
        ch_bbduk_logs = BBMAP_BBDUK.out.log.map { it[1] }
        ch_versions   = ch_versions.mix(BBMAP_BBDUK.out.versions)
    } else {
        ch_clean_reads  = FASTQC_TRIMGALORE.out.reads
        ch_bbduk_logs = []
    }

    //
    // MODULE: Interleave sequences
    //
    SEQTK_MERGEPE(ch_clean_reads)
    ch_versions = ch_versions.mix(SEQTK_MERGEPE.out.versions)

    //
    // SUBWORKFLOW: Perform digital normalization
    //
    ch_reads_to_assembly = Channel.empty()
    if ( params.diginorm ) {
        DIGINORM(SEQTK_MERGEPE.out.reads.collect { meta, fastq -> fastq }, [], 'all_samples')
        ch_versions = ch_versions.mix(SEQTK_MERGEPE.out.versions)
        ch_pe_reads_to_assembly = DIGINORM.out.pairs
        ch_se_reads_to_assembly = DIGINORM.out.singles
    } else {
        ch_pe_reads_to_assembly = SEQTK_MERGEPE.out.reads.map { meta, fastq -> fastq }
        ch_se_reads_to_assembly = []
    }

    //
    // MODULE: Run Megahit on all interleaved fastq files
    //
    MEGAHIT_INTERLEAVED(
        ch_pe_reads_to_assembly.collect(), 
        ch_se_reads_to_assembly.collect(), 
        'all_samples'
    )
    ch_assembly_contigs = MEGAHIT_INTERLEAVED.out.contigs
    ch_versions = ch_versions.mix(MEGAHIT_INTERLEAVED.out.versions)

    //
    // MODULE: Create a BBMap index
    //
    BBMAP_INDEX(ch_assembly_contigs)
    ch_versions   = ch_versions.mix(BBMAP_INDEX.out.versions)

    //
    // MODULE: Call BBMap with the index once per sample
    //
    BBMAP_ALIGN(ch_clean_reads, BBMAP_INDEX.out.index)
    ch_versions   = ch_versions.mix(BBMAP_ALIGN.out.versions)

    //
    // MODULE: Run PROKKA on Megahit output, but split the fasta file in chunks of 10 MiB
    //
    PROKKA(MEGAHIT_INTERLEAVED.out.contigs.splitFasta( size: 10.MB, file: true).map { [[id: 'part_of_all_samples'], it] }, [], [] )
    ch_versions = ch_versions.mix(PROKKA.out.versions)

    //
    // MODULE: Call Prodigal
    //
    UNPIGZ_MEGAHIT_CONTIGS(ch_assembly_contigs)
    ch_versions = ch_versions.mix(UNPIGZ_MEGAHIT_CONTIGS.out.versions)
    
    ch_prodigal = Channel.empty()
    if( params.orf_caller == 'prodigal' ) {
        PRODIGAL(
            UNPIGZ_MEGAHIT_CONTIGS.out.unzipped.collect { [ [ id: 'all_samples' ], it ] },
            'gff'
        )
        ch_prodigal_gff = PRODIGAL.out.gene_annotations
        ch_prodigal_aa  = PRODIGAL.out.amino_acid_fasta
        ch_prodigal_fna = PRODIGAL.out.nucleotide_fasta
        ch_versions = ch_versions.mix(PRODIGAL.out.versions)
    }
    
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowMetatdenovo.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(Channel.from(ch_multiqc_config))
    ch_multiqc_files = ch_multiqc_files.mix(ch_multiqc_custom_config.collect().ifEmpty([]))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    // Make sure we integrate FASTQC output from FASTQC_TRIMGALORE here!!!
    //ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

    MULTIQC (
        ch_multiqc_files.collect()
    )
    multiqc_report = MULTIQC.out.report.toList()
    ch_versions    = ch_versions.mix(MULTIQC.out.versions)
}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
========================================================================================
    THE END
========================================================================================
*/
