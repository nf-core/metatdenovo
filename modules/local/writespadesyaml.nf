process WRITESPADESYAML {
    tag "spades.yaml"
    label 'process_single'

    conda "conda-forge::pigz=2.3.4=2.6"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pigz:2.3.4' :
        'quay.io/biocontainers/pigz:2.3.4' }"

    input:
    //tuple val(meta), path(pe), path(se)
    path(pe)
    path(se)

    output:
    //tuple val(meta), path("*.yaml"), emit: yaml
    path("*.yaml")     , emit: yaml
    path "versions.yml"            , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    //def pe_split = pe.join('",\n            "')
    def pe_split = pe.join('", "')
    """
    cat <<-YAML > spades.yaml
    [
      {
        orientation: "fr",
            type: "paired-end",
            interlaced reads: [ "${pe_split}" ]
      }
    ]
    YAML

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        writespadesyaml: \$(echo \$(bash --version | grep 'GNU bash' | sed 's/.*version //' | sed 's/ .*//'))
    END_VERSIONS
    """
}
