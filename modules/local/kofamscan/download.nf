process KOFAMSCAN_DOWNLOAD {
    tag "KEGG data"
    label 'process_long'

    conda "bioconda::gnu-wget=1.18"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/gnu-wget:1.18--hed695b0_4':
        'quay.io/biocontainers/gnu-wget:1.18--hed695b0_4' }"

    input:
    path kofam_dir

    output:
    path "${kofam_dir}/ko_list" , emit: ko_list
    path "${kofam_dir}/profiles", emit: koprofiles
    path "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    if [ ! -e ${kofam_dir}/ko_list ]; then
        wget https://www.genome.jp/ftp/db/kofam/ko_list.gz
        gunzip -c ko_list.gz > ${kofam_dir}/ko_list
    fi

    if [ ! -e ${kofam_dir}/koprofiles ]; then
        wget -P ${kofam_dir} -c https://www.genome.jp/ftp/db/kofam/profiles.tar.gz
        gunzip -c kofam/profiles.tar.gz | tar vxf - -C kofam/
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        wget: \$(wget --version | grep 'GNU Wget' | sed 's/GNU Wget \\([0-9.]\\+\\) .*/\\1/')
    END_VERSIONS
    """
}
