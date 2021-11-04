/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowMetatdenovo.initialise(params, log)

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
def diginorm_normalizebymedian_options   = modules['diginorm_normalizebymedian']
diginorm_normalizebymedian_options.args += Utils.joinModuleArgs(["-C ${params.diginorm_C} -k ${params.diginorm_k}"])
def diginorm_filterabund_options         = modules['diginorm_filterabund']
diginorm_filterabund_options.args       += Utils.joinModuleArgs(["-C ${params.diginorm_C} -Z ${params.diginorm_C}"])

include { DIGINORM } from '../subworkflows/local/diginorm' addParams(
    diginorm_normalizebymedian_options: diginorm_normalizebymedian_options, 
    diginorm_filterabund_options: diginorm_filterabund_options
)

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

def multiqc_options   = modules['multiqc']
multiqc_options.args += params.multiqc_title ? Utils.joinModuleArgs(["--title \"$params.multiqc_title\""]) : ''

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC  } from '../modules/nf-core/modules/fastqc/main'  addParams( options: modules['fastqc'] )
include { SEQTK_MERGEPE } from '../modules/nf-core/modules/seqtk/mergepe/main' addParams( options: modules['seqtk_mergepe'] )
//include { PRODIGAL } from '../modules/nf-core/modules/prodigal/main' addParams( options: modules['prodigal'] )
include { MULTIQC } from '../modules/nf-core/modules/multiqc/main' addParams( options: multiqc_options   )
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
    // MODULE: Run FastQC
    //
    FASTQC (
        INPUT_CHECK.out.reads
    )
    ch_versions = ch_versions.mix(FASTQC.out.versions.first())

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
    // MODULE: Interleave sequences
    //
    SEQTK_MERGEPE(FASTQC_TRIMGALORE.out.reads)
    ch_versions = ch_versions.mix(SEQTK_MERGEPE.out.versions)

    //
    // SUBWORKFLOW: Perform digital normalization
    //
    ch_reads_to_assembly = Channel.empty()
    if ( params.diginorm ) {
        DIGINORM(SEQTK_MERGEPE.out.reads.collect { meta, fastq -> fastq }, [], 'all_samples')
        ch_versions = ch_versions.mix(SEQTK_MERGEPE.out.versions)
        ch_reads_to_assembly = DIGINORM.out.reads
    } else {
        ch_reads_to_assembly = SEQTK_MERGEPE.out.reads.map { meta, fastq -> fastq }
    }

    //
    // MODULE: Run Megahit on all interleaved fastq files
    //
    MEGAHIT_INTERLEAVED(ch_reads_to_assembly.collect(), 'all_samples')
    ch_versions = ch_versions.mix(MEGAHIT_INTERLEAVED.out.versions)

//    //
//    // MODULE: Call Prodigal
//    //
//    PRODIGAL([ [id: 'full_assembly' ], MEGAHIT_INTERLEAVED.out.contigs)
//    ch_versions = ch_versions.mix(PRODIGAL.out.versions)

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
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]}.ifEmpty([]))

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
