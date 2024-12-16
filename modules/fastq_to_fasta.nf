// convert fastq to fasta (needed for IgBLAST)
process fastq_to_fasta {
    tag "$sequence_id"
    label 'process_medium'

    conda (params.enable_conda ? 'bioconda::seqkit=2.3.1' : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/seqkit:1c7b907eba99f587' :
        'community.wave.seqera.io/library/seqkit:825acf14813a21d5' }"

    input:
    tuple val(sequence_id), path(reads)

    output:
    tuple val(sequence_id), path("*.fasta")

    script:
    """
    read_base_name=\$(basename "$reads" .fastq)
    seqkit fq2fa --threads $task.cpus $reads -o "\${read_base_name}.fasta"
    """
}