// Exit if the user provides both --assembler and --user_assembly, or both --orf_caller and either of the two --user_orfs,
// or none of them
if ( ( params.assembler && params.user_assembly ) || ( ! params.assembler && ! params.user_assembly ) ) {
    error "Provide either `--assembler` or `--user_assembly`!"
}
if ( ( params.orf_caller && ( params.user_orfs_gff ) ) && ( ! params.orf_caller && ! params.user_orfs_gff ) ) {
    error "Provide either `--orf_caller` or `--user_orfs_gff`/`--user_orfs_faa`!"
}

// Exit if the user forgot one of the two `--user_orfs_*`
if ( params.user_orfs_gff && ! params.user_orfs_faa ) {
    error 'When supplying ORFs, both --user_orfs_gff and --user_orfs_faa must be specified, --user_orfs_faa file is missing!'
} else if ( params.user_orfs_faa && ! params.user_orfs_gff ) {
    error 'When supplying ORFs, both --user_orfs_gff and --user_orfs_faa must be specified, --user_orfs_gff file is missing!'
}

// Exit if the user set params.assembler plus any of params.user_orfs_*
if ( params.assembler && ( params.user_orfs_gff || params.user_orfs_faa ) ) {
    error "You can't input your own ORFs (`--user_orfs_*`) if you call for assembly with `--assembler`."
}

// Deal with user-supplied assembly to make sure output names are correct
assembler     = params.assembler
assembly_name = params.assembler ?: params.user_assembly_name

// Deal with params from user-supplied ORFs, and set orf_caller correctly
orf_caller = params.orf_caller
orfs_name  = params.orf_caller ?: params.user_orfs_name

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
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: local
//
include { COLLECT_FEATURECOUNTS              } from '../modules/local/collect/featurecounts/main'
include { COLLECT_STATS                      } from '../modules/local/collect/stats/main'
include { DIAMOND_BLASTP as DIAMOND_TAXONOMY } from '../modules/local/diamond/blastp/main'
include { FORMATSPADES                       } from '../modules/local/format/spades/main'
include { MEGAHIT_INTERLEAVED                } from '../modules/local/megahit/interleaved/main'
include { MERGE_TABLES                       } from '../modules/local/merge/summary/main'
include { FORMAT_DIAMOND_TAX_RANKLIST        } from '../modules/local/diamond/format_tax/ranklist/main'
include { FORMAT_DIAMOND_TAX_TAXDUMP         } from '../modules/local/diamond/format_tax/taxdump/main'
include { SUMTAXONOMY as SUM_DIAMONDTAX      } from '../modules/local/sumtaxonomy/main'
include { TRANSDECODER                       } from '../modules/local/transdecoder/main'
include { TRANSRATE                          } from '../modules/local/transrate/main'
include { UNPIGZ as UNPIGZ_CONTIGS           } from '../modules/local/unpigz/main'
include { UNPIGZ as UNPIGZ_GFF               } from '../modules/local/unpigz/main'
include { WRITESPADESYAML                    } from '../modules/local/spades/writeyaml/main'


//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { validateInputSamplesheet       } from '../subworkflows/local/utils_nfcore_metatdenovo_pipeline'

//
// SUBWORKFLOW: Consisting of local modules
//
include { EGGNOG                  } from '../subworkflows/local/eggnog/main'
include { SUB_EUKULELE            } from '../subworkflows/local/eukulele/main'
include { HMMCLASSIFY             } from '../subworkflows/local/hmmclassify/main'
include { PROKKA_SUBSETS          } from '../subworkflows/local/prokka/subsets/main'
include { FASTQC_TRIMGALORE       } from '../subworkflows/local/fastqc/trimgalore/main'
include { PRODIGAL                } from '../subworkflows/local/prodigal/main'
include { KOFAMSCAN               } from '../subworkflows/local/kofamscan/main'
include { PIPELINE_INITIALISATION } from '../subworkflows/local/utils_nfcore_metatdenovo_pipeline'
include { PIPELINE_COMPLETION     } from '../subworkflows/local/utils_nfcore_metatdenovo_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { BBMAP_ALIGN                                } from '../modules/nf-core/bbmap/align/main'
include { BBMAP_BBDUK                                } from '../modules/nf-core/bbmap/bbduk/main'
include { BBMAP_BBNORM                               } from '../modules/nf-core/bbmap/bbnorm/main'
include { BBMAP_INDEX                                } from '../modules/nf-core/bbmap/index/main'
include { CAT_FASTQ            	                     } from '../modules/nf-core/cat/fastq/main'
include { FASTQC                                     } from '../modules/nf-core/fastqc/main'
include { MULTIQC                                    } from '../modules/nf-core/multiqc/main'
include { PIGZ_COMPRESS as PIGZ_ASSEMBLY             } from '../modules/nf-core/pigz/compress/main'
include { PIGZ_COMPRESS as PIGZ_DIAMOND_LINEAGE      } from '../modules/nf-core/pigz/compress/main'
include { PIGZ_COMPRESS as PIGZ_TRANSDECODER_BED     } from '../modules/nf-core/pigz/compress/main'
include { PIGZ_COMPRESS as PIGZ_TRANSDECODER_CDS     } from '../modules/nf-core/pigz/compress/main'
include { PIGZ_COMPRESS as PIGZ_TRANSDECODER_GFF     } from '../modules/nf-core/pigz/compress/main'
include { PIGZ_COMPRESS as PIGZ_TRANSDECODER_PEP     } from '../modules/nf-core/pigz/compress/main'
include { SEQTK_MERGEPE                              } from '../modules/nf-core/seqtk/mergepe/main'
include { SEQTK_SEQ as SEQTK_SEQ_CONTIG_FILTER       } from '../modules/nf-core/seqtk/seq/main'
include { SPADES                                     } from '../modules/nf-core/spades/main'
include { SUBREAD_FEATURECOUNTS as FEATURECOUNTS_CDS } from '../modules/nf-core/subread/featurecounts/main'
include { TAXONKIT_LINEAGE                           } from '../modules/nf-core/taxonkit/lineage/main'

//
// SUBWORKFLOWS: Installed directly from nf-core/modules
//
include { paramsSummaryMap                           } from 'plugin/nf-schema'
include { paramsSummaryMultiqc                       } from '../subworkflows/nf-core/utils_nfcore_pipeline/'
include { softwareVersionsToYAML                     } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { BAM_SORT_STATS_SAMTOOLS                    } from '../subworkflows/nf-core/bam_sort_stats_samtools/main'
include { UTILS_NEXTFLOW_PIPELINE                    } from '../subworkflows/nf-core/utils_nextflow_pipeline/main'
include { UTILS_NFCORE_PIPELINE                      } from '../subworkflows/nf-core/utils_nfcore_pipeline/main'
include { methodsDescriptionText                     } from '../subworkflows/local/utils_nfcore_metatdenovo_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow METATDENOVO {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    ch_diamond_dbs // channel: paths to Diamond taxonomy databases, read from --diamond_dbs

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // Form a fastq channel from the samplesheet channel
    // DL: I'm not sure which parts are still required after nf-schema. The branch { } certainly is needed.
    ch_fastq = ch_samplesheet
        .flatMap { meta, fastq_files ->
            if (fastq_files.size() <= 2) {
                return [[ meta.id, [meta], fastq_files ]]
            } else {
                def pairs = fastq_files.collate(2)
                return [[ meta.id, pairs.collect { meta + [id: "${meta.id}_${pairs.indexOf(it) + 1}"] }, fastq_files ]]
            }
        }
        .map { id, metas, fastq_files ->
            // Ensure single_end is set correctly in meta
            def updatedMetas = metas.collect { it + [single_end: (fastq_files.size() / metas.size() == 1)] }
            return [id, updatedMetas, fastq_files]
        }
        .map { validateInputSamplesheet(it) }
        .branch {
            meta, fastqs ->
                single  : fastqs.size() == 1
                    return [ meta, fastqs ]
                multiple: fastqs.size() > 1
                    return [ meta, fastqs ]
        }

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
    ch_collect_stats = ch_cat_fastq.collect { meta, fasta -> meta.id }.map { [ [ id:"${assembly_name}.${orfs_name}" ], it ] }
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
        BBMAP_BBDUK ( FASTQC_TRIMGALORE.out.reads, Channel.fromPath(params.sequence_filter).first() )
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
    if ( ! params.user_assembly ) {
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
    if ( ! params.user_assembly ) {
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
    // MODULE: Run Megahit or Spades on all interleaved fastq files
    //
    if ( params.user_assembly ) {
        // If the input assembly is not gzipped, do that since all downstream calls assume this
        if ( ! params.user_assembly.endsWith('.gz') ) {
            PIGZ_ASSEMBLY(Channel.fromPath(params.user_assembly).map { [ [ id:params.user_assembly ], it ] } )
            PIGZ_ASSEMBLY.out.archive.first().set { ch_assembly_contigs }
        } else {
            Channel
                .value ( [ [ id: assembly_name ], file(params.user_assembly) ] )
                .set { ch_assembly_contigs }
        }
    } else if ( assembler == 'spades' ) {
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
            .map { it -> [ [ id: assembly_name ], it, [], [] ] }
            .set { ch_spades }
        SPADES (
            ch_spades,
            WRITESPADESYAML.out.yaml,
            []
        )

        SPADES.out.transcripts
            .ifEmpty { [] }
            .combine(SPADES.out.contigs.ifEmpty { [] } )
            .set { ch_assembly }
        ch_versions = ch_versions.mix(SPADES.out.versions)

        FORMATSPADES( ch_assembly.first() )
        ch_assembly_contigs = FORMATSPADES.out.assembly
        ch_versions    = ch_versions.mix(FORMATSPADES.out.versions)
    } else if ( assembler == 'megahit' ) {
        MEGAHIT_INTERLEAVED(
            ch_pe_reads_to_assembly.toList(),
            ch_se_reads_to_assembly.toList(),
            'megahit_assembly'
        )
        MEGAHIT_INTERLEAVED.out.contigs
            .map { it -> [ [ id: assembly_name ], it ] }
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

        //UNPIGZ_GFF(PROKKA_SUBSETS.out.gff.map { meta, gff -> [ [id: "${orfs_name}.${meta.id}"], gff ] })
        UNPIGZ_GFF(PROKKA_SUBSETS.out.gff)
        ch_gff           = UNPIGZ_GFF.out.unzipped
        ch_versions      = ch_versions.mix(UNPIGZ_GFF.out.versions)
    }

    //
    // MODULE: Run PRODIGAL on assembly output.
    //
    if ( orf_caller == 'prodigal' ) {
        PRODIGAL( ch_assembly_contigs.map { meta, contigs -> [ [id: "${assembly_name}.${orfs_name}"], contigs  ] } )
        ch_protein      = PRODIGAL.out.faa
        ch_versions     = ch_versions.mix(PRODIGAL.out.versions)

        //UNPIGZ_GFF(PRODIGAL.out.gff.map { meta, gff -> [ [id: "${meta.id}"], gff ] })
        UNPIGZ_GFF(PRODIGAL.out.gff)
        ch_gff          = UNPIGZ_GFF.out.unzipped
        ch_versions     = ch_versions.mix(UNPIGZ_GFF.out.versions)
    }

    //
    // SUBWORKFLOW: run TRANSDECODER. Orf caller alternative for eukaryotes.
    //
    if ( orf_caller == 'transdecoder' && ! ( params.user_orfs_gff ) ) {
        TRANSDECODER ( ch_assembly_contigs.map { meta, contigs -> [ [id: "${assembly_name}.${orfs_name}" ], contigs ] } )
        ch_gff      = TRANSDECODER.out.gff
        ch_protein  = TRANSDECODER.out.pep
        ch_versions = ch_versions.mix(TRANSDECODER.out.versions)

        PIGZ_TRANSDECODER_BED(TRANSDECODER.out.bed)
        ch_versions = ch_versions.mix(PIGZ_TRANSDECODER_BED.out.versions)
        PIGZ_TRANSDECODER_CDS(TRANSDECODER.out.cds)
        ch_versions = ch_versions.mix(PIGZ_TRANSDECODER_CDS.out.versions)
        PIGZ_TRANSDECODER_GFF(TRANSDECODER.out.gff)
        ch_versions = ch_versions.mix(PIGZ_TRANSDECODER_GFF.out.versions)
        PIGZ_TRANSDECODER_PEP(TRANSDECODER.out.pep)
        ch_versions = ch_versions.mix(PIGZ_TRANSDECODER_PEP.out.versions)
    }

    // Populate channels if the user provided the orfs
    if ( params.user_orfs_faa && params.user_orfs_gff ) {
        Channel
            .value ( [ [ id: "${assembly_name}.${orfs_name}" ], file(params.user_orfs_gff) ] )
            .set { ch_gff }
        Channel
            .value ( [ [ id: "${assembly_name}.${orfs_name}" ], file(params.user_orfs_faa) ] )
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
        .map { hmm, meta, protein ->[ [ id: "${assembly_name}.${orfs_name}" ], hmm, protein ] }
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
        .map { featurecounts -> [ [ id:"${assembly_name}.${orfs_name}" ], featurecounts ] }
        .set { ch_collect_feature }

    COLLECT_FEATURECOUNTS ( ch_collect_feature )
    ch_versions           = ch_versions.mix(COLLECT_FEATURECOUNTS.out.versions)
    ch_fcs_for_stats      = COLLECT_FEATURECOUNTS.out.counts.collect { meta, tsv -> tsv }.map { [ it ] }
    ch_fcs_for_summary    = COLLECT_FEATURECOUNTS.out.counts.map { meta, tsv -> tsv }
    ch_collect_stats
        .combine(ch_fcs_for_stats)
        .set { ch_collect_stats }

    ch_merge_tables = Channel.empty()
    //
    // SUBWORKFLOW: run eggnog_mapper on the ORF-called amino acid sequences
    //
    if ( ! params.skip_eggnog ) {
        EGGNOG(ch_protein, ch_fcs_for_summary)
        ch_versions = ch_versions.mix(EGGNOG.out.versions)
        ch_merge_tables = ch_merge_tables.mix ( EGGNOG.out.sumtable.map { meta, tsv -> tsv } )
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
        ch_merge_tables = ch_merge_tables.mix ( KOFAMSCAN.out.kofamscan_summary.map { meta, tsv -> tsv } )
    }

    // set up contig channel to use in TransRate
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
        ch_versions = ch_versions.mix(SUB_EUKULELE.out.versions)
        ch_merge_tables = ch_merge_tables.mix ( SUB_EUKULELE.out.taxonomy_summary.map { meta, tsv -> tsv } )
    }

    //
    // Call Diamond for taxonomy with amino acid sequences
    //
    DIAMOND_TAXONOMY(
        ch_protein,
        ch_diamond_dbs.map { [ it[0], it[1] ] },
        102,
        []
    )
    ch_versions     = ch_versions.mix(DIAMOND_TAXONOMY.out.versions)

    // Create a unified channel of the output from Diamond together with the diamond db info to
    // make sure the channels are synchronized before calling TAXONKIT_LINEAGE
    ch_taxonkit_lineage = DIAMOND_TAXONOMY.out.tsv
        .map { it -> [ [ id: it[0].db ], [ id: "${it[0].id}.${it[0].db}.lineage", db: it[0].db ], it[1] ] }
        .join(ch_diamond_dbs)

    TAXONKIT_LINEAGE(
        ch_taxonkit_lineage.map { it -> [ it[1], [], it[2] ] },
        [],
        ch_taxonkit_lineage.map { it -> it[4] },
        ch_taxonkit_lineage.map { it -> it[5] }
    )
    ch_versions     = ch_versions.mix(TAXONKIT_LINEAGE.out.versions)

    PIGZ_DIAMOND_LINEAGE(
        TAXONKIT_LINEAGE.out.tsv
    )
    ch_versions     = ch_versions.mix(PIGZ_DIAMOND_LINEAGE.out.versions)

    FORMAT_DIAMOND_TAX_RANKLIST(
        PIGZ_DIAMOND_LINEAGE.out.archive
            .map { it -> [ [ id: it[0].db ], it[0], it[1] ] }
            .join(ch_diamond_dbs)
            .map { it -> [ [ id: it[1].id - ".lineage" + ".diamond", db: it[1].db ], it[2], it[6] ] }
    )
    ch_versions     = ch_versions.mix(FORMAT_DIAMOND_TAX_RANKLIST.out.versions)

    FORMAT_DIAMOND_TAX_TAXDUMP(
        PIGZ_DIAMOND_LINEAGE.out.archive
            .map { it -> [ [ id: it[0].db ], it[0], it[1] ] }
            .join(ch_diamond_dbs.filter { it -> it[5] })
            .map { it -> [ [ id: it[1].id - ".lineage" + ".diamond", db: it[1].db ], it[2], it[4], it[5], it[6] ] }
    )
    ch_versions     = ch_versions.mix(FORMAT_DIAMOND_TAX_TAXDUMP.out.versions)

    SUM_DIAMONDTAX(
        FORMAT_DIAMOND_TAX_RANKLIST.out.taxonomy
            .map { it -> [ it[0], it[0].db, it[1] ] },
        ch_fcs_for_summary,
        'diamondtax'
    )
    ch_versions     = ch_versions.mix(SUM_DIAMONDTAX.out.versions)

    ch_merge_tables = ch_merge_tables.mix ( SUM_DIAMONDTAX.out.taxonomy_summary.map { meta, tsv -> tsv } )

    //
    // MODULE: Collect statistics from mapping analysis
    //
    MERGE_TABLES (
        ch_merge_tables
            .collect()
            .map { it -> [ [ id: "${assembly_name}.${orfs_name}" ], it ] }
    )
    MERGE_TABLES.out.merged_table
        //.view { "merged0: ${it}" }
        //.collect { meta, tblout -> tblout }
        //.view { "merged1: ${it}" }
        //.map { meta, tblout -> [ tblout ] }
        //.view { "merged2: ${it}" }
    ch_collect_stats = ch_collect_stats
        .combine(
            Channel.empty()
                .mix ( MERGE_TABLES.out.merged_table.map { meta, tblout -> [ tblout ] } )
                .ifEmpty { [ [] ] }
                //.map { [ it ] }
        )
    ch_versions     = ch_versions.mix(MERGE_TABLES.out.versions)

    COLLECT_STATS(ch_collect_stats)
    ch_versions     = ch_versions.mix(COLLECT_STATS.out.versions)

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'metatdenovo_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
