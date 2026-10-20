//
// Assign taxonomy to proteins with Diamond, against every database, in batches of proteins
//

include { DIAMOND_BLASTP as DIAMOND_TAXONOMY         } from '../../../modules/nf-core/diamond/blastp/main'
include { TAXONKIT_LINEAGE                           } from '../../../modules/nf-core/taxonkit/lineage/main'
include { PIGZ_COMPRESS as PIGZ_DIAMOND_LINEAGE      } from '../../../modules/nf-core/pigz/compress/main'
include { FORMAT_DIAMOND_TAX_RANKLIST                } from '../../../modules/local/diamond/format_tax/ranklist/main'
include { FORMAT_DIAMOND_TAX_TAXDUMP                 } from '../../../modules/local/diamond/format_tax/taxdump/main'
include { SUMTAXONOMY as SUM_DIAMONDTAX              } from '../../../modules/local/sumtaxonomy/main'

workflow DIAMONDTAX {

    take:
    proteins    // channel: [ val(meta), path(faa) ]; meta.caller names the ORF source
    dbs         // channel: [ val(meta), path(dmnd), path(names), path(nodes), val(ranks), val(parse_with_taxdump) ]
    fcs         // channel: [ val(meta), path(fcs) ]; meta.caller must match the proteins'
    batchsize   // integer: residues per Diamond batch (splitFasta size counts sequence characters only)

    main:

    // splitFasta cannot read a zero-byte gzip, which stub runs produce
    def ch_proteins_split = proteins
        .branch { _meta, faa ->
            empty: faa.size() == 0
            other: true
        }

    def ch_batches = ch_proteins_split.other
        .splitFasta(size: batchsize, file: true, elem: 1)
        .map { meta, faa -> [ meta, faa, faa.baseName.tokenize('.').last() as Integer ] }
        .mix( ch_proteins_split.empty.map { meta, faa -> [ meta, faa, 1 ] } )
        .map { meta, faa, n -> [ meta + [ id: String.format('%s.%03d', meta.id, n), batch: n, parent: meta ], faa ] }

    // E-values depend on the database size, not the number of queries
    def ch_diamond_input = ch_batches.combine( dbs.map { db -> [ db[0], db[1] ] } )

    DIAMOND_TAXONOMY(
        ch_diamond_input.map { pm, protein, _dm, _db -> [ pm, protein ] },
        ch_diamond_input.map { _pm, _protein, dm, db -> [ dm, db ] },
        102,
        []
    )

    // Sorted by batch number, so the lineage file lists proteins in input order
    def ch_diamond_tsv = DIAMOND_TAXONOMY.out.tsv
        .map { meta, tsv -> [ meta.parent + [ db: meta.db ], [ meta.batch, tsv ] ] }
        .groupTuple()
        .map { meta, batches -> [ meta, batches.sort { b -> b[0] }.collect { b -> b[1] } ] }

    // Not .join(): callers share db names, and join mishandles duplicate keys
    def ch_taxonkit_lineage = ch_diamond_tsv
        .map { meta, tsvs -> [ [ id: "${meta.id}.${meta.db}.lineage", db: meta.db, caller: meta.caller ], tsvs ] }
        .combine(dbs)
        .filter { meta, _tsvs, dbMeta, _dmnd, _names, _nodes, _ranks, _parse -> meta.db == dbMeta.id }

    TAXONKIT_LINEAGE(
        ch_taxonkit_lineage.map { it -> [ it[0], [], it[1] ] },
        ch_taxonkit_lineage.map { it -> it[4] },
        ch_taxonkit_lineage.map { it -> it[5] }
    )

    PIGZ_DIAMOND_LINEAGE(
        TAXONKIT_LINEAGE.out.tsv
    )

    FORMAT_DIAMOND_TAX_RANKLIST(
        PIGZ_DIAMOND_LINEAGE.out.archive
            .combine(dbs)
            .filter { archiveMeta, _archiveFile, dbMeta, _dmnd, _names, _nodes, _ranks, _parse -> archiveMeta.db == dbMeta.id }
            .map { archiveMeta, archiveFile, _dbMeta, _dmnd, _names, _nodes, ranks, _parse -> [ [ id: archiveMeta.id - ".lineage" + ".diamond", db: archiveMeta.db, caller: archiveMeta.caller ], archiveFile, ranks ] }
    )

    FORMAT_DIAMOND_TAX_TAXDUMP(
        PIGZ_DIAMOND_LINEAGE.out.archive
            .combine(dbs.filter { it -> it[5] })
            .filter { archiveMeta, _archiveFile, dbMeta, _dmnd, _names, _nodes, _ranks, _parse -> archiveMeta.db == dbMeta.id }
            .map { archiveMeta, archiveFile, _dbMeta, _dmnd, names, nodes, ranks, _parse -> [ [ id: archiveMeta.id - ".lineage" + ".diamond", db: archiveMeta.db, caller: archiveMeta.caller ], archiveFile, names, nodes, ranks ] }
    )

    // Not .join(): dbs share caller names
    def ch_diamondtax_sum_input = FORMAT_DIAMOND_TAX_RANKLIST.out.taxonomy
        .combine( fcs )
        .filter { meta, _taxonomy, fcsMeta, _fcs -> meta.caller == fcsMeta.caller }
        .map { meta, taxonomy, _fcsMeta, fcs_file -> [ meta, meta.db, taxonomy, fcs_file ] }

    SUM_DIAMONDTAX(ch_diamondtax_sum_input, 'diamondtax')

    emit:
    tsv               = ch_diamond_tsv                            // channel: [ val(meta), [ path(tsv) ] ], one file per batch
    lineage           = PIGZ_DIAMOND_LINEAGE.out.archive          // channel: [ val(meta), path(tsv.gz) ]
    taxonomy          = FORMAT_DIAMOND_TAX_RANKLIST.out.taxonomy  // channel: [ val(meta), path(tsv.gz) ]
    taxonomy_taxdump  = FORMAT_DIAMOND_TAX_TAXDUMP.out.taxonomy   // channel: [ val(meta), path(tsv.gz) ]
    taxonomy_summary  = SUM_DIAMONDTAX.out.taxonomy_summary       // channel: [ val(meta), path(tsv.gz) ]
}
