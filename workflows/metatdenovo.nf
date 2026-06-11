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
include { PIGZ_COMPRESS as PIGZ_PE_READS_FWD         } from '../modules/nf-core/pigz/compress/main'
include { PIGZ_COMPRESS as PIGZ_PE_READS_REV         } from '../modules/nf-core/pigz/compress/main'
include { PIGZ_COMPRESS as PIGZ_SE_READS             } from '../modules/nf-core/pigz/compress/main'
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
    multiqc_config
    multiqc_logo
    multiqc_methods_description
    outdir

    main:

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


    // If the user supplied hmm files, we will run hmmsearch and then rank the results.
    // Create a channel for hmm files.
    ch_hmmrs = channel.empty()
    if ( params.hmmdir ) {
        channel
            .fromPath(params.hmmdir + params.hmmpattern, checkIfExists: true)
            .set { ch_hmmrs }
    } else if ( params.hmmfiles ) {
        channel
            .fromList( params.hmmfiles.tokenize(',') )
            .map { hmmfile -> [ file(hmmfile) ] }
            .set { ch_hmmrs }
    }

    def ch_versions = channel.empty()
    def ch_multiqc_files = channel.empty()

    // Form a fastq channel from the samplesheet channel
    // DL: I'm not sure which parts are still required after nf-schema. The branch { } certainly is needed.
    ch_fastq = ch_samplesheet
        .flatMap { meta, fastq_files ->
            if (fastq_files.size() <= 2) {
                return [[ meta.id, [meta], fastq_files ]]
            } else {
                def pairs = fastq_files.collate(2)
                return [[ meta.id, pairs.collect { pair -> meta + [id: "${meta.id}_${pairs.indexOf(pair) + 1}"] }, fastq_files ]]
            }
        }
        /** DL: In my testing, this fails as entries come in with single meta, multiple read files when appear for single ends
        .map { id, metas, fastq_files ->
            // Ensure single_end is set correctly in meta
            def updatedMetas = metas
                .collect { it + [single_end: (fastq_files.size() / metas.size() == 1)] }
            return [id, updatedMetas, fastq_files]
        }
        **/
        .map { row -> validateInputSamplesheet(row) }
        .branch {
            meta, fastqs ->
                single  : ( meta.single_end && fastqs.size() == 1 ) || ( ! meta.single_end && fastqs.size == 2 )
                    return [ meta, fastqs ]
                multiple: true
                    return [ meta, fastqs ]
        }

    //
    // MODULE: Concatenate FastQ files from the same sample if required
    //
    CAT_FASTQ (
        ch_fastq.multiple
    )

    ch_versions = channel.empty()
    ch_multiqc_files = channel.empty()
    //
    // Gzip unzipped read files
    //
    // We're only doing this for samples having a single row in the sample sheet since those with more than one
    // were gzipped by CAT_FASTQ above.
    //

    // Paired end, forward
    fwd = ch_fastq.single
        .filter { meta, _f -> ! meta.single_end }
        .map { meta, fastqs -> [ meta, fastqs[0] ] }
        .branch {
            meta, fastqs ->
                zipped  : fastqs.name.endsWith('.gz')
                    return [ meta, fastqs ]
                unzipped: true
                    return [ meta, fastqs ]
        }
    PIGZ_PE_READS_FWD(fwd.unzipped)

    // Paired end, reverse
    rev = ch_fastq.single
        .filter { meta, _f -> ! meta.single_end }
        .map { meta, fastqs -> [ meta, fastqs[1] ] }
        .branch {
            meta, fastqs ->
                zipped  : fastqs.name.endsWith('.gz')
                    return [ meta, fastqs ]
                unzipped: true
                    return [ meta, fastqs ]
        }
    PIGZ_PE_READS_REV(rev.unzipped)

    // Single end
    se = ch_fastq.single
        .filter { meta, _f -> meta.single_end }
        .map { meta, fastqs -> [ meta, fastqs[0] ] }
        .branch {
            meta, fastqs ->
                zipped  : fastqs.name.endsWith('.gz')
                    return [ meta, fastqs ]
                unzipped: true
                    return [ meta, fastqs ]
        }
    PIGZ_SE_READS(se.unzipped)

    // Join the three channels with the originally zipped to form a new ch_fastq of the same structure as the original
    ch_fastq = fwd.zipped.concat(PIGZ_PE_READS_FWD.out.archive)
        .join(rev.zipped.concat(PIGZ_PE_READS_REV.out.archive))
        .map { meta, fwd_read, rev_read -> [ meta, [ fwd_read, rev_read ] ] }
        .concat(
            se.zipped
                .concat(PIGZ_SE_READS.out.archive)
                .map { meta, fastq -> [ meta, [ fastq ] ] }
        )
        .concat(CAT_FASTQ.out.reads)

    //
    // SUBWORKFLOW: Read QC and trim adapters
    //
    FASTQC_TRIMGALORE (
        ch_fastq,
        params.skip_fastqc || params.skip_qc,
        params.skip_trimming
    )

    ch_collect_stats = ch_fastq
        .collect { meta, _fasta -> meta }
        .map { metas -> [ [ id:"${assembly_name}.${orfs_name}" ], metas ] }

    if ( params.skip_trimming ) {
        ch_collect_stats = ch_collect_stats
            .map { meta, samples -> [ meta, samples, [] ] }

    } else {
        ch_collect_stats = ch_collect_stats
            .combine(
                FASTQC_TRIMGALORE.out.trim_log
                    .collect { _meta, report ->
                        if ( report in List ) {
                            report[0]
                        } else {
                            report
                        }
                    }
                    .map { report -> [ report ] }
            )
    }

    //
    // MODULE: Run BBDuk to clean out whatever sequences the user supplied via params.sequence_filter
    //
    if ( params.sequence_filter ) {
        BBMAP_BBDUK ( FASTQC_TRIMGALORE.out.reads, channel.fromPath(params.sequence_filter).first() )
        ch_clean_reads  = BBMAP_BBDUK.out.reads
        ch_bbduk_logs = BBMAP_BBDUK.out.log.collect { _meta, log ->  log }.map { log -> [ log ] }
        ch_collect_stats = ch_collect_stats.combine(ch_bbduk_logs)
        ch_multiqc_files = ch_multiqc_files.mix(BBMAP_BBDUK.out.log.collect{ _meta, log -> log })
    } else {
        ch_clean_reads  = FASTQC_TRIMGALORE.out.reads
        ch_bbduk_logs = channel.empty()
        ch_collect_stats = ch_collect_stats
            .map { meta, samples, report -> [ meta, samples, report, [] ] }
    }

    //
    // MODULE: Interleave sequences for assembly
    //
    ch_interleaved = channel.empty()
    if ( ! params.user_assembly ) {
        SEQTK_MERGEPE(ch_clean_reads)
        ch_interleaved = SEQTK_MERGEPE.out.reads
    }

    //
    // SUBWORKFLOW: Perform digital normalization.
    //
    if ( ! params.user_assembly ) {
        if ( params.bbnorm ) {
            BBMAP_BBNORM(
                ch_interleaved
                    .collect { _meta, fastq -> fastq }
                    .map { fastq -> [ [id:'all_samples', single_end:true], fastq ] }
            )
            ch_pe_reads_to_assembly = BBMAP_BBNORM.out.fastq.map { _meta, fasta -> fasta }
            ch_se_reads_to_assembly = channel.empty()
        } else {
            ch_pe_reads_to_assembly = ch_interleaved
                .filter { meta, _fastq -> ! meta.single_end }
                .map { _meta, fastq -> fastq }
            ch_se_reads_to_assembly = ch_interleaved
                .filter { meta, _fastq -> meta.single_end }
                .map { _meta, fastq -> fastq }
        }
    }

    //
    // MODULE: Run Megahit or Spades on all interleaved fastq files
    //
    if ( params.user_assembly ) {
        // If the input assembly is not gzipped, do that since all downstream calls assume this
        if ( ! params.user_assembly.endsWith('.gz') ) {
            PIGZ_ASSEMBLY(
                channel
                    .fromPath(params.user_assembly)
                    .map { path -> [ [ id:params.user_assembly ], path ] }
            )
            ch_assembly_contigs = PIGZ_ASSEMBLY.out.archive.first()
        } else {
            ch_assembly_contigs = channel
                .value ( [ [ id: assembly_name ], file(params.user_assembly) ] )
        }
    } else if ( assembler == 'spades' ) {
        // 1. Write a yaml file for Spades
        WRITESPADESYAML (
            ch_pe_reads_to_assembly.toList(),
            ch_se_reads_to_assembly.toList()
        )

        // 2. Call the module with a channel with all fastq files plus the yaml
        ch_spades = ch_pe_reads_to_assembly
            .mix(ch_se_reads_to_assembly)
            .collect()
            .map { it -> [ [ id: assembly_name ], it, [], [] ] }
        SPADES (
            ch_spades,
            WRITESPADESYAML.out.yaml,
            []
        )

        ch_spades_assembly = SPADES.out.transcripts
            .ifEmpty { [] }
            .combine(SPADES.out.contigs.ifEmpty { [] } )

        FORMATSPADES( ch_spades_assembly.first() )
        ch_assembly_contigs = FORMATSPADES.out.assembly
    } else if ( assembler == 'megahit' ) {
        MEGAHIT_INTERLEAVED(
            ch_pe_reads_to_assembly.toList(),
            ch_se_reads_to_assembly.toList(),
            'megahit_assembly'
        )
        ch_assembly_contigs = MEGAHIT_INTERLEAVED.out.contigs
            .map { it -> [ [ id: assembly_name ], it ] }
    } else {
        error 'Assembler not specified!'
    }

    // If the user asked for length filtering, perform that with SEQTK_SEQ (the actual length parameter is used in modules.config)
    if ( params.min_contig_length > 0 ) {
        SEQTK_SEQ_CONTIG_FILTER ( ch_assembly_contigs )
        ch_assembly_contigs = SEQTK_SEQ_CONTIG_FILTER.out.fastx
    }

    //
    // Call ORFs
    //
    ch_gff      = channel.empty()
    ch_protein  = channel.empty()

    //
    // SUBWORKFLOW: Run PROKKA_SUBSETS on assmebly output, but split the fasta file in chunks of 10 MB, then concatenate and compress output.
    //
    if ( params.orf_caller == 'prokka' ) {
        PROKKA_SUBSETS(ch_assembly_contigs, params.prokka_batchsize)
        ch_protein       = PROKKA_SUBSETS.out.faa
        ch_multiqc_files = ch_multiqc_files.mix(PROKKA_SUBSETS.out.prokka_log)

        //UNPIGZ_GFF(PROKKA_SUBSETS.out.gff.map { meta, gff -> [ [id: "${orfs_name}.${meta.id}"], gff ] })
        UNPIGZ_GFF(PROKKA_SUBSETS.out.gff)
        ch_gff           = UNPIGZ_GFF.out.unzipped
    }

    //
    // MODULE: Run PRODIGAL on assembly output.
    //
    if ( orf_caller == 'prodigal' ) {
        PRODIGAL( ch_assembly_contigs.map { _meta, contigs -> [ [id: "${assembly_name}.${orfs_name}"], contigs  ] } )
        ch_protein      = PRODIGAL.out.faa
        UNPIGZ_GFF(PRODIGAL.out.gff)
        ch_gff          = UNPIGZ_GFF.out.unzipped
    }

    //
    // SUBWORKFLOW: run TRANSDECODER. Orf caller alternative for eukaryotes.
    //
    if ( orf_caller == 'transdecoder' ) {
        TRANSDECODER ( ch_assembly_contigs.map { _meta, contigs -> [ [id: "${assembly_name}.${orfs_name}" ], contigs ] } )
        ch_gff      = TRANSDECODER.out.gff
        ch_protein  = TRANSDECODER.out.pep

        PIGZ_TRANSDECODER_BED(TRANSDECODER.out.bed)
        PIGZ_TRANSDECODER_CDS(TRANSDECODER.out.cds)
        PIGZ_TRANSDECODER_GFF(TRANSDECODER.out.gff)
        PIGZ_TRANSDECODER_PEP(TRANSDECODER.out.pep)
    }

    // Populate channels if the user provided the orfs
    if ( params.user_orfs_faa && params.user_orfs_gff ) {
        ch_gff = channel.value ( [ [ id: "${assembly_name}.${orfs_name}" ], file(params.user_orfs_gff) ] )
        ch_protein = channel.value ( [ [ id: "${assembly_name}.${orfs_name}" ], file(params.user_orfs_faa) ] )
    }

    //
    // MODULE: Create a BBMap index
    //
    BBMAP_INDEX(ch_assembly_contigs.map { _meta, contigs -> contigs })

    //
    // MODULE: Call BBMap with the index once per sample
    //
    BBMAP_ALIGN ( ch_clean_reads, BBMAP_INDEX.out.index )

    //
    // SUBWORKFLOW: classify ORFs with a set of hmm files
    //
    ch_hmmclassify = ch_hmmrs
        .combine(ch_protein)
        .map { hmm, _meta, protein ->[ [ id: "${assembly_name}.${orfs_name}" ], hmm, protein ] }
    HMMCLASSIFY ( ch_hmmclassify )

    //
    // MODULE: FeatureCounts. Create a table for each samples that provides raw counts as result of the alignment.
    //
    BAM_SORT_STATS_SAMTOOLS (
        BBMAP_ALIGN.out.bam,
        ch_assembly_contigs.map { meta, fasta -> [meta, fasta, []] }
    )

    ch_featurecounts = BAM_SORT_STATS_SAMTOOLS.out.bam
        .combine(ch_gff.map { _meta, bam -> bam } )

    ch_collect_stats = ch_collect_stats
        .combine(BAM_SORT_STATS_SAMTOOLS.out.idxstats.collect { _meta, idxstats -> idxstats }.map { idxstats -> [ idxstats ] } )

    FEATURECOUNTS_CDS ( ch_featurecounts)

    //
    // Collate and save software versions
    //
    def topic_versions = channel.topic("versions")
        .distinct()
        .branch { entry ->
            versions_file: entry instanceof Path
            versions_tuple: true
        }

    def topic_versions_string = topic_versions.versions_tuple
        .map { process, tool, version ->
            [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
        }
        .groupTuple(by:0)
        .map { process, tool_versions ->
            tool_versions.unique().sort()
            "${process}:\n${tool_versions.join('\n')}"
        }

    def ch_collated_versions = softwareVersionsToYAML(ch_versions.mix(topic_versions.versions_file))
        .mix(topic_versions_string)
        .collectFile(
            storeDir: "${outdir}/pipeline_info",
            name: 'nf_core_'  +  'metatdenovo_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    //
    // MODULE: Collect featurecounts output counts in one table
    //
    ch_collect_feature = FEATURECOUNTS_CDS.out.counts
        .collect() { _meta, featurecounts -> featurecounts }
        .map { featurecounts -> [ [ id:"${assembly_name}.${orfs_name}" ], featurecounts ] }

    COLLECT_FEATURECOUNTS ( ch_collect_feature )
    ch_versions           = ch_versions.mix(COLLECT_FEATURECOUNTS.out.versions)
    ch_fcs_for_stats      = COLLECT_FEATURECOUNTS.out.counts.collect { _meta, tsv -> tsv }.map { tsv -> [ tsv ] }
    ch_fcs_for_summary    = COLLECT_FEATURECOUNTS.out.counts.map { _meta, tsv -> tsv }
    ch_collect_stats = ch_collect_stats.combine(ch_fcs_for_stats)

    // Initialize ch_merge_tables that will be populated with tables from annotation tools and used by the MERGE_TABLES module which output will then be passed to the COLLECT_STATS module
    ch_merge_tables = channel.empty()

    //
    // SUBWORKFLOW: run eggnog_mapper on the ORF-called amino acid sequences
    //
    if ( ! params.skip_eggnog ) {
        EGGNOG(ch_protein, ch_fcs_for_summary)
        ch_merge_tables = ch_merge_tables.mix ( EGGNOG.out.sumtable.map { _meta, tsv -> tsv } )
    }

    //
    // SUBWORKFLOW: run kofamscan on the ORF-called amino acid sequences
    //
    if( !params.skip_kofamscan ) {
        ch_kofamscan = ch_protein.map { meta, protein -> [ meta, protein ] }
        KOFAMSCAN( ch_kofamscan, ch_fcs_for_summary)
        ch_merge_tables = ch_merge_tables.mix ( KOFAMSCAN.out.kofamscan_summary.map { _meta, tsv -> tsv } )
    }

    // set up contig channel to use in TransRate
    UNPIGZ_CONTIGS(ch_assembly_contigs)
    ch_unzipped_contigs = UNPIGZ_CONTIGS.out.unzipped

    //
    // MODULE: Use TransRate to judge assembly quality, piped into MultiQC
    //
    TRANSRATE(ch_unzipped_contigs)

    //
    // SUBWORKFLOW: Eukulele
    //
    ch_eukulele_db = channel.empty()
    if( ! params.skip_eukulele ) {
        // Make sure the eukulele_dbpath exists
        d = new File("${params.eukulele_dbpath}")
        if ( ! d.exists() ) {
            d.mkdirs()
        }

        // Create a channel for EUKulele either with a named database or not. The latter means a user-provided database in a directory.
        if ( params.eukulele_db ) {
            ch_eukulele_db = channel
                .of ( params.eukulele_db )
                .map { db -> [ db, file(params.eukulele_dbpath) ] }
        } else {
            ch_eukulele_db = channel.fromPath(params.eukulele_dbpath, checkIfExists: true)
                .map { path -> [ [], path ] }
        }
        ch_eukulele = ch_protein
            .map { meta, protein -> [ [ id:"${meta.id}" ], protein ] }
            .combine( ch_eukulele_db )
        SUB_EUKULELE( ch_eukulele, ch_fcs_for_summary )

        ch_merge_tables = ch_merge_tables.mix ( SUB_EUKULELE.out.taxonomy_summary.map { _meta, tsv -> tsv } )
    }

    //
    // Call Diamond for taxonomy with amino acid sequences
    //
    DIAMOND_TAXONOMY(
        ch_protein,
        ch_diamond_dbs.map { db -> [ db[0], db[1] ] },
        102,
        []
    )

    // Create a unified channel of the output from Diamond together with the diamond db info to
    // make sure the channels are synchronized before calling TAXONKIT_LINEAGE
    ch_taxonkit_lineage = DIAMOND_TAXONOMY.out.tsv
        .map { it -> [ [ id: it[0].db ], [ id: "${it[0].id}.${it[0].db}.lineage", db: it[0].db ], it[1] ] }
        .join(ch_diamond_dbs)

    TAXONKIT_LINEAGE(
        ch_taxonkit_lineage.map { it -> [ it[1], [], it[2] ] },
        ch_taxonkit_lineage.map { it -> it[4] },
        ch_taxonkit_lineage.map { it -> it[5] }
    )

    PIGZ_DIAMOND_LINEAGE(
        TAXONKIT_LINEAGE.out.tsv
    )

    FORMAT_DIAMOND_TAX_RANKLIST(
        PIGZ_DIAMOND_LINEAGE.out.archive
            .map { it -> [ [ id: it[0].db ], it[0], it[1] ] }
            .join(ch_diamond_dbs)
            .map { it -> [ [ id: it[1].id - ".lineage" + ".diamond", db: it[1].db ], it[2], it[6] ] }
    )

    FORMAT_DIAMOND_TAX_TAXDUMP(
        PIGZ_DIAMOND_LINEAGE.out.archive
            .map { it -> [ [ id: it[0].db ], it[0], it[1] ] }
            .join(ch_diamond_dbs.filter { it -> it[5] })
            .map { it -> [ [ id: it[1].id - ".lineage" + ".diamond", db: it[1].db ], it[2], it[4], it[5], it[6] ] }
    )

    SUM_DIAMONDTAX(
        FORMAT_DIAMOND_TAX_RANKLIST.out.taxonomy
            .map { it -> [ it[0], it[0].db, it[1] ] },
        ch_fcs_for_summary,
        'diamondtax'
    )

    ch_merge_tables = ch_merge_tables.mix ( SUM_DIAMONDTAX.out.taxonomy_summary.map { _meta, tsv -> tsv } )

    //
    // MODULE: Collect statistics from mapping analysis
    //
    MERGE_TABLES (
        ch_merge_tables
            .collect()
            .map { it -> [ [ id: "${assembly_name}.${orfs_name}" ], it ] }
    )
    MERGE_TABLES.out.merged_table

    ch_collect_stats = ch_collect_stats
        .combine(
            channel.empty()
                .mix ( MERGE_TABLES.out.merged_table.map { _meta, tblout -> [ tblout ] } )
                .ifEmpty { [ [] ] }
        )

    COLLECT_STATS(ch_collect_stats)

    //
    // MODULE: MultiQC
    //
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    def ch_summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def ch_workflow_summary = channel.value(paramsSummaryMultiqc(ch_summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    def ch_multiqc_custom_methods_description = multiqc_methods_description
        ? file(multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    def ch_methods_description = channel.value(methodsDescriptionText(ch_multiqc_custom_methods_description))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml', sort: true))
    MULTIQC(
        ch_multiqc_files.flatten().collect().map { files ->
            [
                [id: 'metatdenovo'],
                files,
                multiqc_config
                    ? file(multiqc_config, checkIfExists: true)
                    : file("${projectDir}/assets/multiqc_config.yml", checkIfExists: true),
                multiqc_logo ? file(multiqc_logo, checkIfExists: true) : [],
                [],
                [],
            ]
        }
    )

    emit:
    multiqc_report = MULTIQC.out.report.map { _meta, report -> [report] }.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
