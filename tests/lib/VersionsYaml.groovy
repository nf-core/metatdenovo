// Shared nf-test helpers for the software_mqc_versions.yml snapshot (see nf-test's lib folder convention).
class VersionsYaml {

    // Retrieval dates change between test runs; mask them as removeNextflowVersion masks the Nextflow version.
    static maskRetrievedDate(obj) {
        if (obj instanceof Map) {
            return obj.collectEntries { k, v -> [ k, maskRetrievedDate(v) ] }
        } else if (obj instanceof String) {
            return obj.replaceAll(/retrieved=\d{4}-\d{2}-\d{2}/, 'retrieved=[DATE]')
        }
        return obj
    }
}
