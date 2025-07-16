// module for running Jakob's find-cdr3
// general layout is based on the nf-core modules
process find_cdr3 {
    tag "$sample_id"
    label 'process_high'
    publishDir "${params.out_dir}/original_annotation_tsv", mode: 'copy', pattern: "*_annotated.tsv"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'library://kzeglinski/nanologix/find-cdr3:v0.0.2' :
        'MAKE A DOCKER CONTAINER!!!!' }"

    input:
    tuple val(sample_id), path(reads)
    val igblast_databases

    output:
    tuple val(sample_id), path('*_annotated.tsv'), emit: annotation

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    find-cdr3 \
    --input-reads $reads \
    --reference-fasta ${igblast_databases}/databases/imgt_alpaca_ighv \
    --threads $task.cpus \
    --output-tsv ${sample_id}_annotated.tsv
    """

}