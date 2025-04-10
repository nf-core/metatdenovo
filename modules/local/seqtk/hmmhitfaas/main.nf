process SEQTK_HMMHITFAAS {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/seqtk:1.4--h577a1d6_3':
        'biocontainers/seqtk:1.4--h577a1d6_3' }"

    input:
    tuple val(meta), path(hmmrank), path(faa)

    output:
    tuple val(meta), path("hits/*.faa.gz"), emit: faas
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir hits/

    # Loop over all unique profiles found in hmmrank, and call seqtk subseq for all orfs matching that profile
    for profile in \$(gunzip -c $hmmrank | grep -v '^profile' | cut -f 1 | sort -u); do
        seqtk subseq $faa <(gunzip -c $hmmrank | grep "^\${profile}" | cut -f 2) | gzip -c > hits/${prefix}.\${profile}.faa.gz
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqtk: \$(echo \$(seqtk 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    mkdir hits/
    echo "" | gzip -c > hits/${prefix}.faa.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        seqtk: \$(echo \$(seqtk 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
    END_VERSIONS
    """
}
