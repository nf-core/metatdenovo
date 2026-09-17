//
// Normalise user-supplied ORFs (--user_orfs_gff/--user_orfs_faa) before they join the pipeline's
// own callers. A GFF/FASTA pair re-supplied from a MetaEuk run (e.g. re-fed after MetaEuk's own
// resume/cache broke on a large run) is raw MetaEuk output, not yet through the same
// FORMAT_METAEUK_GFF/FORMAT_METAEUKFAA rewrite a pipeline-internal `--orf_caller metaeuk` run
// always gets first -- so it carries `TCS_ID=`/`Target_ID=` attributes instead of a plain `ID=`.
//

include { FORMAT_METAEUKFAA  } from '../../../modules/local/format/metaeukfaa/main'
include { FORMAT_METAEUK_GFF } from '../../../modules/local/format/metaeuk/main'

// A GFF already carrying a proper ID= attribute needs no further work, whether that's from a
// pipeline-internal MetaEuk run (which prepends one in front of its own untouched
// Target_ID=/TCS_ID= attributes -- so those remain present even after formatting) or any other
// caller with its own GFF3 ID=. Checking for TCS_ID='s mere presence instead would misfire on
// that same already-normalised MetaEuk output and reformat it a second time. Reads at most 20
// lines, decompressing if needed, regardless of how large the underlying file is -- these can be
// tens of GB for a real assembly.
def hasIdAttribute(gff) {
    def reader = gff.name.endsWith('.gz') ?
        new java.util.zip.GZIPInputStream(gff.newInputStream()).newReader() :
        gff.newReader()
    reader.withCloseable { r ->
        r.lines().limit(20).anyMatch { line -> line.contains('\tID=') || line.contains(';ID=') }
    }
}

// A raw MetaEuk header has >= 7 pipe-delimited fields (see FORMAT_METAEUKFAA); an already-rewritten
// one has exactly 4. Used only to cross-check against the gff's own raw/normalised call below --
// gff and faa are supposed to be a matched pair from the same source, and disagreeing on which
// state they're in means one of them is not.
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
