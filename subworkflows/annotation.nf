include { fastq_to_fasta } from '../modules/fastq_to_fasta'
include { igblast } from '../modules/igblast'
include { find_cdr3 } from '../modules/find_cdr3'

process merge_chunked_tsvs {
    tag "$sequence_id"
    label 'process_low'
    publishDir "${params.out_dir}/original_annotation_tsv", mode: 'copy', pattern: "*_merged.tsv"

    input:
    tuple val(sequence_id), path(tsvs)

    output:
    tuple val(sequence_id), path("*_merged.tsv")

    script:
    """
    cat *.tsv > ${sequence_id}_merged.tsv
    """
}

workflow annotation {
    take:
        trimmed_and_merged_reads
        igblast_databases
        use_igblast

    main:

        if (use_igblast == false) {
            // use jakob's cdr3 finder
            // cdr3_finder()
            final_tsvs = find_cdr3(trimmed_and_merged_reads, igblast_databases).annotation

        } else {
            // use igblast
            // split reads into chunks bc igblast is slow
            trimmed_and_merged_reads
                .splitFastq(by: params.chunk_size, file: true)
                .set{chunked_merged_reads}


            // convert fastq to fasta (needed to run igblast)
            chunked_merged_fastas = fastq_to_fasta(chunked_merged_reads)

            // set up environment variables
            igdata_dir="${igblast_databases}/igdata/"
            igblastdb_dir="${igblast_databases}/databases/"

            // annotate reads using igblast
            annotated_tsv = igblast(
                chunked_merged_fastas,
                igblast_databases,
                igdata_dir,
                igblastdb_dir).airr_table

            // merge the tsvs
            grouped_tsvs = annotated_tsv.groupTuple(by: 0) // group by first element (the sample ID)
            final_tsvs = merge_chunked_tsvs(grouped_tsvs) // cat all of them
        }

    emit:
        final_tsvs
 }