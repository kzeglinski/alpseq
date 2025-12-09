include { multiqc } from '../modules/multiqc'

process percentage_passing_trim_merge {
    tag "$sequence_id"
    label 'process_medium'

    input:
    tuple val(sequence_id), path(before_trimming), path(after_trimming)

    output:
    path("*.tsv"), emit: num_reads

    script:
    """
    # get number of reads before
    num_reads_before=\$(zcat ${before_trimming} | echo \$((`wc -l`/4)))

    # get number of reads after
    num_reads_after=\$(zcat ${after_trimming} | echo \$((`wc -l`/4)))

    # get proportion
    proportion=\$(echo "scale=3 ; \$num_reads_after / \$num_reads_before" | bc)

    # get percentage
    percentage=\$(echo "scale=3 ; \$proportion * 100" | bc)

    # write to text file
    printf "$sequence_id\t\$num_reads_before\t\$num_reads_after\t\$percentage\n" >> ${sequence_id}_pass_trim_merge.tsv
    """
}

process cat_all_percentage_passing_trim_merge {
    label 'process_low'
    publishDir "${params.out_dir}", mode: 'copy', pattern: "percentage_passing_trim_merge.tsv"

    input:
    path(individual_passing_trim_merge)

    output:
    path("*.tsv"), emit: percentage_passing_trim_merge

    script:
    """
    # write header
    printf "sequence_id\tnum_reads_before\tnum_reads_after\tpercentage_passing_trim_merge\n" >> header.tsv

    # cat them
    cat header.tsv ${individual_passing_trim_merge} > percentage_passing_trim_merge.tsv
    """
}

workflow quality_control {
    take:
        sample_read_pairs
        trimmed_and_merged_reads
        fastqc_reports

    main:
        // make multiqc report
        multiqc(fastqc_reports)
        multiqc_report = multiqc.out.report
        multiqc_plots = multiqc.out.plots

        // determine the percentage of reads that pass trimming/merging
        // we need a tuple like [sample name, r1 file, trimmed_merged_file]
        sample_read_pairs
            .map(it -> [it[0], it[1][0]])
            .combine(trimmed_and_merged_reads, by: 0)
            .set{before_after_trimming}

        percentage_passing_trim_merge(before_after_trimming)

        // collect the results so our cat process only runs once
        individual_passing_trim_merge = percentage_passing_trim_merge.out.num_reads.collect()
        cat_all_percentage_passing_trim_merge(individual_passing_trim_merge)
        percentage_passing_trim_merge = cat_all_percentage_passing_trim_merge.out.percentage_passing_trim_merge

    emit:
        multiqc_report
        percentage_passing_trim_merge
        multiqc_plots
}