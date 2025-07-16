// flash to merge R1 and R2
// adapted from nf-core module: https://github.com/nf-core/modules/blob/master/modules/nf-core/flash/main.nf

process flash {
    tag "$sequence_id"
    label 'process_medium'
    //publishDir "${params.out_dir}/merged_reads", mode: 'copy', pattern: "*.fastq.gz"

    conda "bioconda::flash=1.2.11"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/flash:89b360363f36dab3' :
        'community.wave.seqera.io/library/flash:9791f2226797de01' }"

    input:
    tuple val(sequence_id), path(reads)
    val(maximum_overlap)

    output:
    tuple val(sequence_id), path("*.extendedFrags.fastq.gz"), emit: combined_reads
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${sequence_id}"
    """
    flash \\
        $args \\
        -o ${prefix} \\
        -M $maximum_overlap \\
        -z \\
        ${reads[0]} \\
        ${reads[1]}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        flash: \$(echo \$(flash --version 2>&1) | sed 's/^.*FLASH v//; s/ .*\$//')
    END_VERSIONS
    """
}