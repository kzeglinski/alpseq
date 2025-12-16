// module for running Jakob's find-cdr3
// general layout is based on the nf-core modules
process matchbox_annotate {
    tag "$sample_id"
    label 'process_low'
    publishDir "${params.out_dir}/original_annotation", mode: 'copy', pattern: "*_annotation.csv"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'library://kzeglinski/kzeglinski/nanologix-matchbox:v0.0.4' :
        'ghcr.io/jakob-schuster/matchbox:sha256-774786ff07c5d9d16d1fb64d8329c9c2cf9fd0fe3d89856e2a2672133e0c3fae' }"

    input:
    tuple val(sample_id), path(reads), path(read_sample)
    val igblast_databases
    val mb_scripts
    val fwr4_seq
    val mb_error_rate
    val num_v_genes

    output:
    tuple val(sample_id), path('*_mb_annotation.csv'), emit: annotation

    when:
    task.ext.when == null || task.ext.when

    script:    
    """
    read_base_name=\$(basename "$reads" .fasta)
    # make codon csv file
    echo "codon\nTGT\nTGC\nGGT\nGGC\nGGA\nGGG" > codons.csv

    # headers for the output files
    echo "id,v_gene,cdr3_nt,cdr3_aa,full_seq_nt,full_seq_aa" > \${read_base_name}_mb_annotation.csv

    # then run matchbox
    matchbox --script-file ${mb_scripts}/generate_reference_counts.mb $read_sample -e 0.2 --threads ${task.cpus} --args "references = '${igblast_databases}/databases/imgt_alpaca_ighv', fwr4 = '$fwr4_seq'" -o "."
    sort -nk2 names.csv -t, | tail -n $num_v_genes | cut -f1 -d, > names_sorted.csv
    sed -i '1iname' names_sorted.csv

    # first generate the trimmed reference file
    cp ${igblast_databases}/databases/imgt_alpaca_ighv reference.fasta
    matchbox --script-file ${mb_scripts}/preprocess_references.mb reference.fasta --threads ${task.cpus} -e 0.2 > trimmed_reference.fasta

    # then process the data
    # (ideally, matchbox should also be able to add the header to this csv output)
    matchbox --script-file ${mb_scripts}/find_cdr3_rescue.mb $reads -e $mb_error_rate --threads ${task.cpus} --args "fwr4 = '$fwr4_seq'" >> \${read_base_name}_mb_annotation.csv
    """
}