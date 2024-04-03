/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input, params.multiqc_config ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }
// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

WorkflowMetatdenovo.initialise(params, log)

// Deal with user-supplied assembly to make sure output names are correct
if ( params.assembly ) {
    assembler = 'user_assembly'
} else {
    assembler = params.assembler
}

// Deal with params from user-supplied ORFs, and set orf_caller correctly
if ( params.gff && params.protein_fasta ) {
    orf_caller = 'user_orfs'
} else if ( params.gff && ! params.protein_fasta ) {
    error 'When supplying ORFs, both --gff and --protein_fasta must be specified, --protein_fasta file is missing!'
} else if ( params.protein_fasta && ! params.gff ) {
    error 'When supplying ORFs, both --gff and --protein_fasta must be specified, --gff file is missing!'
} else {
    orf_caller = params.orf_caller
}

// set an empty multiqc channel
ch_multiqc_files = Channel.empty()

// If the user supplied hmm files, we will run hmmsearch and then rank the results.
// Create a channel for hmm files.
ch_hmmrs = Channel.empty()
if ( params.hmmdir ) {
    Channel
        .fromPath(params.hmmdir + params.hmmpattern, checkIfExists: true)
        .set { ch_hmmrs }
} else if ( params.hmmfiles ) {
    Channel
        .fromList( params.hmmfiles.tokenize(',') )
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
include { WRITESPADESYAML                  } from '../modules/local/writespadesyaml'
include { MEGAHIT_INTERLEAVED              } from '../modules/local/megahit/interleaved'
include { COLLECT_FEATURECOUNTS            } from '../modules/local/collect_featurecounts'
include { COLLECT_STATS                    } from '../modules/local/collect_stats'
include { FORMATSPADES                     } from '../modules/local/formatspades'
include { UNPIGZ as UNPIGZ_CONTIGS         } from '../modules/local/unpigz'
include { UNPIGZ as UNPIGZ_GFF             } from '../modules/local/unpigz'
include { MERGE_TABLES                     } from '../modules/local/merge_summary_tables'
include { TRANSRATE                        } from '../modules/local/transrate'

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'

//
// SUBWORKFLOW: Consisting of local modules
//
include { EGGNOG            } from '../subworkflows/local/eggnog'
include { SUB_EUKULELE      } from '../subworkflows/local/eukulele'
include { HMMCLASSIFY       } from '../subworkflows/local/hmmclassify'
include { PROKKA_SUBSETS    } from '../subworkflows/local/prokka_subsets'
include { TRANSDECODER      } from '../subworkflows/local/transdecoder'
include { FASTQC_TRIMGALORE } from '../subworkflows/local/fastqc_trimgalore'
include { PRODIGAL          } from '../subworkflows/local/prodigal'
include { KOFAMSCAN         } from '../subworkflows/local/kofamscan'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { BBMAP_BBDUK                                } from '../modules/nf-core/bbmap/bbduk/main'
include { BBMAP_INDEX                                } from '../modules/nf-core/bbmap/index/main'
include { BBMAP_ALIGN                                } from '../modules/nf-core/bbmap/align/main'
include { BBMAP_BBNORM                               } from '../modules/nf-core/bbmap/bbnorm/main'
include { SEQTK_MERGEPE                              } from '../modules/nf-core/seqtk/mergepe/main'
include { SUBREAD_FEATURECOUNTS as FEATURECOUNTS_CDS } from '../modules/nf-core/subread/featurecounts/main'
include { SPADES                                     } from '../modules/nf-core/spades/main'
include { SEQTK_SEQ as SEQTK_SEQ_CONTIG_FILTER       } from '../modules/nf-core/seqtk/seq/main'
include { CAT_FASTQ            	                     } from '../modules/nf-core/cat/fastq/main'
include { FASTQC                                     } from '../modules/nf-core/fastqc/main'
include { MULTIQC                                    } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS                } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { PIGZ_COMPRESS as PIGZ_ASSEMBLY             } from '../modules/nf-core/pigz/compress/main'

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

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        file(params.input)
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

    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first())

    //
    // SUBWORKFLOW: Read QC and trim adapters
    //
    FASTQC_TRIMGALORE (
        ch_cat_fastq,
        params.skip_fastqc || params.skip_qc,
        params.skip_trimming
    )
    ch_versions = ch_versions.mix(FASTQC_TRIMGALORE.out.versions)
    ch_collect_stats = ch_cat_fastq.collect { meta, fasta -> meta.id }.map { [ [ id:"${assembler}.${orf_caller}" ], it ] }
    if ( params.skip_trimming ) {
        ch_collect_stats
            .map { meta, samples -> [ meta, samples, [] ] }
            .set { ch_collect_stats }

    } else {
        if ( params.se_reads ) {
            ch_collect_stats
                .combine(FASTQC_TRIMGALORE.out.trim_log.collect { meta, report -> report }.map { [ it ] })
                .set { ch_collect_stats }
        } else {
            ch_collect_stats
                .combine(FASTQC_TRIMGALORE.out.trim_log.collect { meta, report -> report[0] }.map { [ it ] })
                .set { ch_collect_stats }
        }
    }

    //
    // MODULE: Run BBDuk to clean out whatever sequences the user supplied via params.sequence_filter
    //
    if ( params.sequence_filter ) {
        BBMAP_BBDUK ( FASTQC_TRIMGALORE.out.reads, params.sequence_filter )
        ch_clean_reads  = BBMAP_BBDUK.out.reads
        ch_bbduk_logs = BBMAP_BBDUK.out.log.collect { meta, log ->  log }.map { [ it ] }
        ch_versions   = ch_versions.mix(BBMAP_BBDUK.out.versions.first())
        ch_collect_stats
            .combine(ch_bbduk_logs)
            .set {ch_collect_stats}
        ch_multiqc_files = ch_multiqc_files.mix(BBMAP_BBDUK.out.log.collect{ meta, log -> log })
    } else {
        ch_clean_reads  = FASTQC_TRIMGALORE.out.reads
        ch_bbduk_logs = Channel.empty()
        ch_collect_stats
            .map { meta, samples, report -> [ meta, samples, report, [] ] }
            .set { ch_collect_stats }
    }

    //
    // MODULE: Interleave sequences for assembly
    //
    // DL & DDL: We can probably not deal with single end input
    ch_interleaved = Channel.empty()
    if ( ! params.assembly ) {
        if ( params.se_reads) {
            ch_single_end = ch_clean_reads
        } else {
            SEQTK_MERGEPE(ch_clean_reads)
            ch_interleaved = SEQTK_MERGEPE.out.reads
            ch_versions    = ch_versions.mix(SEQTK_MERGEPE.out.versions)
        }
    }

    //
    // SUBWORKFLOW: Perform digital normalization. There are two options: khmer or BBnorm. The latter is faster.
    //
    if ( ! params.assembly ) {
        if ( params.se_reads ) {
            if ( params.bbnorm ) {
                BBMAP_BBNORM(ch_single_end.collect { meta, fastq -> fastq }.map {[ [id:'all_samples', single_end:true], it ] } )
                ch_se_reads_to_assembly = BBMAP_BBNORM.out.fastq.map { meta, fasta -> fasta }
                ch_pe_reads_to_assembly = Channel.empty()
                ch_versions    = ch_versions.mix(BBMAP_BBNORM.out.versions)
            } else {
                ch_se_reads_to_assembly = ch_single_end.map { meta, fastq -> fastq }
                ch_pe_reads_to_assembly = Channel.empty()
            }
        }
        else if ( params.bbnorm ) {
            BBMAP_BBNORM(ch_interleaved.collect { meta, fastq -> fastq }.map {[ [id:'all_samples', single_end:true], it ] } )
            ch_pe_reads_to_assembly = BBMAP_BBNORM.out.fastq.map { meta, fasta -> fasta }
            ch_se_reads_to_assembly = Channel.empty()
            ch_versions    = ch_versions.mix(BBMAP_BBNORM.out.versions)
        } else {
            ch_pe_reads_to_assembly = ch_interleaved.map { meta, fastq -> fastq }
            ch_se_reads_to_assembly = Channel.empty()
        }
    }

    //
    // MODULE: Run Megahit or RNAspades on all interleaved fastq files
    //
    if ( params.assembly ) {
        // If the input assembly is not gzipped, do that since all downstream calls assume this
        if ( ! params.assembly.endsWith('.gz') ) {
            PIGZ_ASSEMBLY(Channel.fromPath(params.assembly).map { [ [ id:params.assembly ], it ] } )
            PIGZ_ASSEMBLY.out.archive.first().set { ch_assembly_contigs }
        } else {
            Channel
                .value ( [ [ id: 'user_assembly' ], file(params.assembly) ] )
                .set { ch_assembly_contigs }
        }
    } else if ( assembler == 'rnaspades' ) {
        // 1. Write a yaml file for Spades
        WRITESPADESYAML (
            ch_pe_reads_to_assembly.toList(),
            ch_se_reads_to_assembly.toList()
        )
        ch_versions    = ch_versions.mix(WRITESPADESYAML.out.versions)
        // 2. Call the module with a channel with all fastq files plus the yaml
        ch_pe_reads_to_assembly
            .mix(ch_se_reads_to_assembly)
            .collect()
            .map { [ [ id:'rnaspades' ], it, [], [] ] }
            .set { ch_spades }
        SPADES (
            ch_spades,
            WRITESPADESYAML.out.yaml,
            []
        )
        ch_assembly = SPADES.out.transcripts
        ch_versions = ch_versions.mix(SPADES.out.versions)
        FORMATSPADES( ch_assembly )
        ch_assembly_contigs = FORMATSPADES.out.assembly
        ch_versions    = ch_versions.mix(FORMATSPADES.out.versions)
    } else if ( assembler == 'megahit' ) {
        MEGAHIT_INTERLEAVED(
            ch_pe_reads_to_assembly.toList(),
            ch_se_reads_to_assembly.toList(),
            'megahit_assembly'
        )
        MEGAHIT_INTERLEAVED.out.contigs
            .map { [ [ id: 'megahit' ], it ] }
            .set { ch_assembly_contigs }
        ch_versions = ch_versions.mix(MEGAHIT_INTERLEAVED.out.versions)
    } else { 
        error 'Assembler not specified!' 
    }

    // If the user asked for length filtering, perform that with SEQTK_SEQ (the actual length parameter is used in modules.config)
    if ( params.min_contig_length > 0 ) {
        SEQTK_SEQ_CONTIG_FILTER ( ch_assembly_contigs )
        ch_assembly_contigs = SEQTK_SEQ_CONTIG_FILTER.out.fastx
        ch_versions = ch_versions.mix(SEQTK_SEQ_CONTIG_FILTER.out.versions)
    }

    //
    // Call ORFs
    //
    ch_gff      = Channel.empty()
    ch_protein  = Channel.empty()

    //
    // SUBWORKFLOW: Run PROKKA_SUBSETS on assmebly output, but split the fasta file in chunks of 10 MB, then concatenate and compress output.
    //
    if ( params.orf_caller == 'prokka' ) {
        PROKKA_SUBSETS(ch_assembly_contigs, params.prokka_batchsize)
        ch_versions      = ch_versions.mix(PROKKA_SUBSETS.out.versions)
        ch_protein       = PROKKA_SUBSETS.out.faa
        ch_multiqc_files = ch_multiqc_files.mix(PROKKA_SUBSETS.out.prokka_log)

        UNPIGZ_GFF(PROKKA_SUBSETS.out.gff.map { meta, gff -> [ [id: "${params.orf_caller}.${meta}"], gff ] })
        ch_gff           = UNPIGZ_GFF.out.unzipped
        ch_versions      = ch_versions.mix(UNPIGZ_GFF.out.versions)
    }

    //
    // MODULE: Run PRODIGAL on assembly output.
    //
    if ( orf_caller == 'prodigal' ) {
        PRODIGAL( ch_assembly_contigs.map { meta, contigs -> [ [id: "${assembler}.${orf_caller}"], contigs  ] } )
        ch_protein      = PRODIGAL.out.faa
        ch_versions     = ch_versions.mix(PRODIGAL.out.versions)

        UNPIGZ_GFF(PRODIGAL.out.gff.map { meta, gff -> [ [id: "${meta.id}.${orf_caller}"], gff ] })
        ch_gff          = UNPIGZ_GFF.out.unzipped
        ch_versions     = ch_versions.mix(UNPIGZ_GFF.out.versions)
    }

    //
    // SUBWORKFLOW: run TRANSDECODER. Orf caller alternative for eukaryotes.
    //
    if ( orf_caller == 'transdecoder' ) {
        TRANSDECODER ( ch_assembly_contigs.map { meta, contigs -> [ [id: "transdecoder.${meta.id}" ], contigs ] } )
        ch_gff      = TRANSDECODER.out.gff
        ch_protein  = TRANSDECODER.out.pep
        ch_versions = ch_versions.mix(TRANSDECODER.out.versions)
    }

    // Populate channels if the user provided the orfs
    if ( orf_caller == 'user_orfs' ) {
        Channel
            .value ( [ [ id: "${assembler}.${orf_caller}" ], file(params.gff) ] )
            .set { ch_gff }
        Channel
            .value ( [ [ id: "${assembler}.${orf_caller}" ], file(params.protein_fasta) ] )
            .set { ch_protein }
    }

    //
    // MODULE: Create a BBMap index
    //
    BBMAP_INDEX(ch_assembly_contigs.map { meta, contigs -> contigs })
    ch_versions   = ch_versions.mix(BBMAP_INDEX.out.versions)

    //
    // MODULE: Call BBMap with the index once per sample
    //
    BBMAP_ALIGN ( ch_clean_reads, BBMAP_INDEX.out.index )
    ch_versions = ch_versions.mix(BBMAP_ALIGN.out.versions)

    //
    // SUBWORKFLOW: classify ORFs with a set of hmm files
    //
    ch_hmmrs
        .combine(ch_protein)
        .map { hmm, meta, protein ->[ [ id: "${assembler}.${orf_caller}" ], hmm, protein ] }
        .set { ch_hmmclassify }
    HMMCLASSIFY ( ch_hmmclassify )
    ch_versions = ch_versions.mix(HMMCLASSIFY.out.versions)

    //
    // MODULE: FeatureCounts. Create a table for each samples that provides raw counts as result of the alignment.
    //

    BAM_SORT_STATS_SAMTOOLS ( BBMAP_ALIGN.out.bam, ch_assembly_contigs )
    ch_versions = ch_versions.mix(BAM_SORT_STATS_SAMTOOLS.out.versions)

    BAM_SORT_STATS_SAMTOOLS.out.bam
        .combine(ch_gff.map { meta, bam -> bam } )
        .set { ch_featurecounts }

    ch_collect_stats
        .combine(BAM_SORT_STATS_SAMTOOLS.out.idxstats.collect { meta, idxstats -> idxstats }.map { [ it ] } )
        .set { ch_collect_stats }

    FEATURECOUNTS_CDS ( ch_featurecounts)
    ch_versions       = ch_versions.mix(FEATURECOUNTS_CDS.out.versions)

    //
    // MODULE: Collect featurecounts output counts in one table
    //
    FEATURECOUNTS_CDS.out.counts
        .collect() { meta, featurecounts -> featurecounts }
        .map { featurecounts -> [ [ id:"${assembler}.${orf_caller}" ], featurecounts ] }
        .set { ch_collect_feature }

    COLLECT_FEATURECOUNTS ( ch_collect_feature )
    ch_versions           = ch_versions.mix(COLLECT_FEATURECOUNTS.out.versions)
    ch_fcs_for_stats      = COLLECT_FEATURECOUNTS.out.counts.collect { meta, tsv -> tsv }.map { [ it ] }
    ch_fcs_for_summary    = COLLECT_FEATURECOUNTS.out.counts.map { meta, tsv -> tsv }
    ch_collect_stats
        .combine(ch_fcs_for_stats)
        .set { ch_collect_stats }

    //
    // SUBWORKFLOW: run eggnog_mapper on the ORF-called amino acid sequences
    //
    if ( ! params.skip_eggnog ) {
        EGGNOG(ch_protein, ch_fcs_for_summary)
        ch_versions = ch_versions.mix(EGGNOG.out.versions)
        ch_merge_tables = EGGNOG.out.sumtable
    } else {
        ch_protein
            .map { meta, protein -> [ meta, [] ] }
            .set { ch_merge_tables }
    }


    //
    // SUBWORKFLOW: run kofamscan on the ORF-called amino acid sequences
    //
    if( !params.skip_kofamscan ) {
        ch_protein
            .map { meta, protein -> [ meta, protein ] }
            .set { ch_kofamscan }
        KOFAMSCAN( ch_kofamscan, ch_fcs_for_summary)
        ch_versions = ch_versions.mix(KOFAMSCAN.out.versions)
        ch_kofamscan_summary = KOFAMSCAN.out.kofamscan_summary.collect{ meta, tsv -> tsv }
        ch_merge_tables
            .combine( ch_kofamscan_summary )
            .set { ch_merge_tables }
    } else {
        ch_merge_tables
            .map { meta, tsv -> [ meta, tsv, [] ] }
            .set { ch_merge_tables }
    }

    // set up contig channel to use in CAT and TransRate
    UNPIGZ_CONTIGS(ch_assembly_contigs)
    ch_unzipped_contigs = UNPIGZ_CONTIGS.out.unzipped
    ch_versions = ch_versions.mix(UNPIGZ_CONTIGS.out.versions)

    //
    // MODULE: Use TransRate to judge assembly quality, piped into MultiQC
    //
    TRANSRATE(ch_unzipped_contigs)
    ch_versions = ch_versions.mix(TRANSRATE.out.versions)

    //
    // SUBWORKFLOW: Eukulele
    //
    ch_eukulele_db = Channel.empty()
    if( ! params.skip_eukulele ) {
        // Create a channel for EUKulele either with a named database or not. The latter means a user-provided database in a directory.
        if ( params.eukulele_db ) {
            Channel
                .of ( params.eukulele_db )
                .map { [ it, file(params.eukulele_dbpath) ] }
                .set { ch_eukulele_db }
        } else {
            Channel.fromPath(params.eukulele_dbpath, checkIfExists: true)
                .map { [ [], it ] }
                .set { ch_eukulele_db }
        }
        ch_protein
            .map { meta, protein -> [ [ id:"${meta.id}" ], protein ] }
            .combine( ch_eukulele_db )
            .set { ch_eukulele }
        SUB_EUKULELE( ch_eukulele, ch_fcs_for_summary )
        ch_taxonomy_summary = SUB_EUKULELE.out.taxonomy_summary.collect{ meta, tsv -> tsv }
        ch_versions = ch_versions.mix(SUB_EUKULELE.out.versions)
        ch_merge_tables
            .combine( ch_taxonomy_summary )
            .set { ch_merge_tables }
    } else {
        ch_merge_tables
            .map { meta, tsv1, tsv2 -> [ meta, tsv1, tsv2, [] ] }
            .set { ch_merge_tables }
    }

    //
    // MODULE: Collect statistics from mapping analysis
    //
    if( !params.skip_eggnog  || !params.skip_eukulele || !params.skip_kofamscan) {
        MERGE_TABLES ( ch_merge_tables )
        ch_collect_stats
            .combine(MERGE_TABLES.out.merged_table.collect{ meta, tblout -> tblout }.map { [ it ] })
            .set { ch_collect_stats }
        ch_versions       = ch_versions.mix(MERGE_TABLES.out.versions)
    } else {
        ch_collect_stats
            .map { meta, samples, report, tsv, idxstats, counts -> [ meta, samples, report, tsv, idxstats, counts, [] ] }
            .set { ch_collect_stats }
    }

    COLLECT_STATS(ch_collect_stats)
    ch_versions     = ch_versions.mix(COLLECT_STATS.out.versions)

    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowMetatdenovo.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowMetatdenovo.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC_TRIMGALORE.out.trim_zip.collect{ meta, zip -> zip })
    ch_multiqc_files = ch_multiqc_files.mix(TRANSRATE.out.assembly_qc.collect{ meta, tbl -> tbl })
    ch_multiqc_files = ch_multiqc_files.mix(BAM_SORT_STATS_SAMTOOLS.out.idxstats.collect{ meta, idxstats -> idxstats })
    ch_multiqc_files = ch_multiqc_files.mix(FEATURECOUNTS_CDS.out.summary.collect{ meta, summary -> summary })


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
    NfcoreTemplate.dump_parameters(workflow, params)
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

workflow.onError {
    if (workflow.errorReport.contains("Process requirement exceeds available memory")) {
        println("🛑 Default resources exceed availability 🛑 ")
        println("💡 See here on how to configure pipeline: https://nf-co.re/docs/usage/configuration#tuning-workflow-resources 💡")
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
