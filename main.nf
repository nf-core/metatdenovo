/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/metatdenovo
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/metatdenovo
    Website: https://nf-co.re/metatdenovo
    Slack  : https://nfcore.slack.com/channels/metatdenovo
----------------------------------------------------------------------------------------
*/

<<<<<<< HEAD
=======
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { METATDENOVO  } from './workflows/metatdenovo'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_metatdenovo_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_metatdenovo_pipeline'
include { getGenomeAttribute      } from './subworkflows/local/utils_nfcore_metatdenovo_pipeline'

>>>>>>> TEMPLATE
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

<<<<<<< HEAD
include { METATDENOVO             } from './workflows/metatdenovo'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_metatdenovo_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_metatdenovo_pipeline'
=======
// TODO nf-core: Remove this line if you don't need a FASTA file
//   This is an example of how to use getGenomeAttribute() to fetch parameters
//   from igenomes.config using `--genome`
params.fasta = getGenomeAttribute('fasta')
>>>>>>> TEMPLATE

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow NFCORE_METATDENOVO {

    take:
    samplesheet // channel: samplesheet read in from --input
<<<<<<< HEAD
    diamond_dbs // channel: paths to Diamond taxonomy databases, read from --diamond_dbs
=======
>>>>>>> TEMPLATE

    main:

    //
    // WORKFLOW: Run pipeline
    //
    METATDENOVO (
<<<<<<< HEAD
        samplesheet,
        diamond_dbs
=======
        samplesheet
>>>>>>> TEMPLATE
    )
    emit:
    multiqc_report = METATDENOVO.out.multiqc_report // channel: /path/to/multiqc_report.html
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
<<<<<<< HEAD
        params.diamond_dbs
=======
        params.help,
        params.help_full,
        params.show_hidden
>>>>>>> TEMPLATE
    )

    //
    // WORKFLOW: Run main workflow
    //
    NFCORE_METATDENOVO (
<<<<<<< HEAD
        PIPELINE_INITIALISATION.out.samplesheet,
        PIPELINE_INITIALISATION.out.diamond_paths
    )

=======
        PIPELINE_INITIALISATION.out.samplesheet
    )
>>>>>>> TEMPLATE
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        NFCORE_METATDENOVO.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
