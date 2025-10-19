include { fastq_to_fasta } from '../modules/helpers'
include { sample_1000 } from '../modules/helpers'
include { igblast } from '../modules/igblast'
include { matchbox_annotate } from '../modules/matchbox_annotate'

process merge_chunked_tsvs {
    tag "$sequence_id"
    label 'process_low'
    publishDir "${params.out_dir}/original_annotation", mode: 'copy', pattern: "*_merged.tsv"

    input:
    tuple val(sequence_id), path(tsvs)

    output:
    tuple val(sequence_id), path("*_merged.tsv")

    script:
    """
    cat *.tsv > ${sequence_id}_merged.tsv
    """
}

process merge_chunked_csvs {
    tag "$sequence_id"
    label 'process_low'
    publishDir "${params.out_dir}/original_annotation", mode: 'copy', pattern: "*_merged.csv"

    input:
    tuple val(sequence_id), path(tsvs)

    output:
    tuple val(sequence_id), path("*_merged.csv")

    script:
    """
    cat *.csv > ${sequence_id}_merged.csv
    """
}

workflow annotation {
    take:
        trimmed_and_merged_reads
        igblast_databases
        use_igblast
        mb_scripts

    main:

        if (use_igblast == false) {
            // use jakob's cdr3 finder
            // split reads into chunks
            trimmed_and_merged_reads
                .splitFastq(by: 1000000, file: true)
                .set{chunked_merged_reads}
            // sample 1000 reads
            reads_with_sample = sample_1000(chunked_merged_reads)

            // use matchbox
            annotated_csv = matchbox_annotate(reads_with_sample, igblast_databases, mb_scripts).annotation

            // merge the tsvs
            grouped_csvs = annotated_csv.groupTuple(by: 0) // group by first element (the sample ID)
            final_annotation = merge_chunked_csvs(grouped_csvs) // cat all of them

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
            final_annotation = merge_chunked_tsvs(grouped_tsvs) // cat all of them
        }

    emit:
        final_annotation
 }