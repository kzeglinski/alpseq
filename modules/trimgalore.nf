// this module runs trimgalore for read trimming and fastqc
// adapted from the nf-core module here: https://github.com/kzeglinski/nabseq_nf/blob/main/modules/local/cutadapt.nf
process trimgalore {
    tag "$sequence_id"
    label 'process_medium'
    //publishDir "${params.out_dir}/trimmed_reads", mode: 'copy', pattern: "*.fq.gz"

    conda "bioconda::trim-galore=0.6.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/trim-galore:0deab4fc8772f146' :
        'community.wave.seqera.io/library/trim-galore:3a7047006bbc4999' }"

    input:
    tuple val(sequence_id), path(reads)
    val(adapter_r1)
    val(adapter_r2)

    output:
    tuple val(sequence_id), path("*{3prime,5prime,trimmed,val}*.fq.gz"), emit: reads
    path("*.zip")                                                      , emit: fastqc_zip , optional: true
    path "versions.yml"                                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // Calculate number of --cores for TrimGalore based on value of task.cpus
    // See: https://github.com/FelixKrueger/TrimGalore/blob/master/Changelog.md#version-060-release-on-1-mar-2019
    // See: https://github.com/nf-core/atacseq/pull/65
    def cores = 1
    if (task.cpus) {
        cores = (task.cpus as int) - 4
        if (cores < 1) cores = 1
        if (cores > 8) cores = 8
    }

    // Added soft-links to original fastqs for consistent naming in MultiQC
    def prefix = task.ext.prefix ?: "${sequence_id}"
    """
    [ ! -f  ${prefix}_1.fastq.gz ] && ln -s ${reads[0]} ${prefix}_1.fastq.gz
    [ ! -f  ${prefix}_2.fastq.gz ] && ln -s ${reads[1]} ${prefix}_2.fastq.gz
    trim_galore \\
        $args \\
        --cores $cores \\
        --paired \\
        --gzip \\
        --quality 20 \\
        --fastqc \\
        -a $adapter_r1 \\
        -a2 $adapter_r2 \\
        ${prefix}_1.fastq.gz \\
        ${prefix}_2.fastq.gz
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trimgalore: \$(echo \$(trim_galore --version 2>&1) | sed 's/^.*version //; s/Last.*\$//')
        cutadapt: \$(cutadapt --version)
    END_VERSIONS
    """
}

process trimgalore_sanger {
    tag "$sequence_id"
    label 'process_medium'
    //publishDir "${params.out_dir}/trimmed_reads", mode: 'copy', pattern: "*.fq.gz"

    conda "bioconda::trim-galore=0.6.7"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/trim-galore:0deab4fc8772f146' :
        'community.wave.seqera.io/library/trim-galore:3a7047006bbc4999' }"

    input:
    tuple val(sequence_id), path(reads)
    val(adapter_r1)

    output:
    tuple val(sequence_id), path("*{3prime,5prime,trimmed,val}*.fq.gz"), emit: reads
    path("*.zip")                                                      , emit: fastqc_zip , optional: true
    path "versions.yml"                                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    // Calculate number of --cores for TrimGalore based on value of task.cpus
    // See: https://github.com/FelixKrueger/TrimGalore/blob/master/Changelog.md#version-060-release-on-1-mar-2019
    // See: https://github.com/nf-core/atacseq/pull/65
    def cores = 1
    if (task.cpus) {
        cores = (task.cpus as int) - 4
        if (cores < 1) cores = 1
        if (cores > 8) cores = 8
    }

    // Added soft-links to original fastqs for consistent naming in MultiQC
    def prefix = task.ext.prefix ?: "${sequence_id}"
    """
    [ ! -f  ${prefix}.fastq.gz ] && ln -s ${reads[0]} ${prefix}.fastq.gz
    trim_galore \\
        $args \\
        --cores $cores \\
        --gzip \\
        --quality 20 \\
        --fastqc \\
        -a $adapter_r1 \\
        ${prefix}.fastq.gz
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        trimgalore: \$(echo \$(trim_galore --version 2>&1) | sed 's/^.*version //; s/Last.*\$//')
        cutadapt: \$(cutadapt --version)
    END_VERSIONS
    """
}