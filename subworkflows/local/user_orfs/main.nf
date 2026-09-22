//
// Normalise user-supplied ORFs. A re-supplied raw MetaEuk gff/faa pair gets the same rewrite as an internal MetaEuk run.
//

include { FORMAT_METAEUKFAA  } from '../../../modules/local/format/metaeukfaa/main'
include { FORMAT_METAEUK_GFF } from '../../../modules/local/format/metaeuk/main'

// Test ID=, not TCS_ID=: normalised MetaEuk output keeps TCS_ID=.
// Comments are skipped before the 20-line limit; a GFF3 header has one ##sequence-region per contig.
def hasIdAttribute(gff) {
    def reader = gff.name.endsWith('.gz') ?
        new java.util.zip.GZIPInputStream(gff.newInputStream()).newReader() :
        gff.newReader()
    reader.withCloseable { r ->
        r.lines()
            .filter { line -> ! line.startsWith('#') }
            .limit(20)
            .anyMatch { line -> line.contains('\tID=') || line.contains(';ID=') }
    }
}

// Raw MetaEuk headers have >= 7 pipe-delimited fields, rewritten ones 4. Cross-checks the gff.
def hasRawMetaeukHeader(faa) {
    def reader = faa.name.endsWith('.gz') ?
        new java.util.zip.GZIPInputStream(faa.newInputStream()).newReader() :
        faa.newReader()
    reader.withCloseable { r ->
        def header = r.lines().filter { line -> line.startsWith('>') }.findFirst()
        header.isPresent() && header.get().split(/\|/).length >= 7
    }
}

workflow USER_ORFS {

    take:
    user_orfs // channel: [ val(meta), path(gff), path(faa) ]

    main:

    ch_branched = user_orfs
        .map { meta, gff, faa ->
            def raw_gff = ! hasIdAttribute(gff)
            def raw_faa = hasRawMetaeukHeader(faa)
            if (raw_gff != raw_faa) {
                error "USER_ORFS: ${meta.id}'s gff and faa disagree on whether they're raw MetaEuk output or already normalised (gff looks ${raw_gff ? 'raw' : 'normalised'}, faa looks ${raw_faa ? 'raw' : 'normalised'}) -- supply a matched gff/faa pair from the same source."
            }
            [ meta, gff, faa, raw_gff ]
        }
        .branch { _meta, _gff, _faa, raw_gff ->
            needs_format: raw_gff
            generic: true
        }

    FORMAT_METAEUK_GFF ( ch_branched.needs_format.map { meta, gff, _faa, _raw -> [ meta, gff ] } )
    FORMAT_METAEUKFAA  ( ch_branched.needs_format.map { meta, _gff, faa, _raw -> [ meta, faa ] } )

    ch_formatted = FORMAT_METAEUK_GFF.out.format_gff
        .join( FORMAT_METAEUKFAA.out.format_faa )
        .mix( ch_branched.generic.map { meta, gff, faa, _raw -> [ meta, gff, faa ] } )

    emit:
    gff = ch_formatted.map { meta, gff, _faa -> [ meta, gff ] } // channel: [ val(meta), path(gff) ]
    faa = ch_formatted.map { meta, _gff, faa -> [ meta, faa ] } // channel: [ val(meta), path(faa) ]
}
