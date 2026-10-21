// Shared nf-test helpers for the software_mqc_versions.yml snapshot (see nf-test's lib folder convention).
class VersionsYaml {

    // Database provenance strings (kofamscan_db, dbcan_db, eukulele_db) embed the real retrieval
    // mtime, which is "today" whenever the pipeline runs, so it can't be snapshotted literally.
    // Mask it the same way nft-utils' removeNextflowVersion masks the Nextflow version.
    static maskRetrievedDate(obj) {
        if (obj instanceof Map) {
            return obj.collectEntries { k, v -> [ k, maskRetrievedDate(v) ] }
        } else if (obj instanceof String) {
            return obj.replaceAll(/retrieved=\d{4}-\d{2}-\d{2}/, 'retrieved=[DATE]')
        }
        return obj
    }
}
