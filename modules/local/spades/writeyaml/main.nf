process WRITESPADESYAML {
    tag "spades.yaml"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/pigz:2.3.4' :
        'biocontainers/pigz:2.3.4' }"

    input:
    path(pe)
    path(se)

    output:
    path("*.yaml")     , emit: yaml

    when:
    task.ext.when == null || task.ext.when

    script:
    def read_list = []
    if ( pe ) read_list.add('{ orientation: "fr", type: "paired-end", interlaced reads: [ "' + pe.join('", "') + '" ] }')
    if ( se ) read_list.add('{ type: "single", single reads: [ "' + se.join('", "') + '" ] }')
    def reads = read_list.join(", ")
    """
    cat <<-YAML > spades.yaml
    [ $reads ]
    YAML
    """

    stub:
    """
    touch spades.yaml
    """
}
