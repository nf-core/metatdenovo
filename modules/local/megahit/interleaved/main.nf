process MEGAHIT_INTERLEAVED {
    tag "$assembly"
    label 'process_high'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-0f92c152b180c7cd39d9b0e6822f8c89ccb59c99:8ec213d21e5d03f9db54898a2baeaf8ec729b447-0' :
        'biocontainers/mulled-v2-0f92c152b180c7cd39d9b0e6822f8c89ccb59c99:8ec213d21e5d03f9db54898a2baeaf8ec729b447-0' }"

    input:
    path intl_pe_reads
    path se_reads
    val  assembly

    output:
    path("megahit_out/*.contigs.fa.gz")                            , emit: contigs
    path("megahit_out/*.log")                                      , emit: log
    path("megahit_out/intermediate_contigs/k*.contigs.fa.gz")      , emit: k_contigs
    path("megahit_out/intermediate_contigs/k*.addi.fa.gz")         , emit: addi_contigs
    path("megahit_out/intermediate_contigs/k*.local.fa.gz")        , emit: local_contigs
    path("megahit_out/intermediate_contigs/k*.final.contigs.fa.gz"), emit: kfinal_contigs
    path "versions.yml"                                            , emit: versions

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: ''
    single_ends = se_reads ? "--read ${se_reads.join(',')}" : ""
    pair_ends = intl_pe_reads ? "--12 ${intl_pe_reads.join(',')}" : ""

    """
    megahit \\
        $pair_ends \\
        ${single_ends} \\
        -t $task.cpus \\
        -m ${task.memory.toBytes()} \\
        $args \\
        --out-prefix $assembly

    pigz \\
        --no-name \\
        -p $task.cpus \\
        $args2 \\
        megahit_out/*.fa \\
        megahit_out/intermediate_contigs/*.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        megahit_interleaved: \$(echo \$(megahit -v 2>&1) | sed 's/MEGAHIT v//')
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: ''
    """
    mkdir -p megahit_out/intermediate_contigs
    echo "" | gzip > megahit_out/${assembly}.contigs.fa.gz
    touch megahit_out/${assembly}.log
    echo "" | gzip > megahit_out/intermediate_contigs/k21.contigs.fa.gz
    echo "" | gzip > megahit_out/intermediate_contigs/k21.addi.fa.gz
    echo "" | gzip > megahit_out/intermediate_contigs/k21.local.fa.gz
    echo "" | gzip > megahit_out/intermediate_contigs/k21.final.contigs.fa.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        megahit_interleaved: 1.2.9
    END_VERSIONS
    """
}
