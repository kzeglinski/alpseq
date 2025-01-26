process run_r_script {
    tag "$sequence_id"
    label 'process_high'
    publishDir "${params.out_dir}/processed_tsv", mode: 'copy', pattern: "*.tsv"
    container "library://kzeglinski/nanologix/nanologix-report:v0.4.0"

    input:
    tuple val(sequence_id), path(merged_tsv)

    output:
    tuple val(sequence_id), path('*.tsv'), emit: processed_tsv

    script:
    """
    #!/usr/bin/env Rscript
    library(stringr)
    library(dplyr)
    library(tidyr)
    library(vroom)

    working_dir <- getwd()

    # CURRENT TSV FILE
    sample_id <- '$sequence_id'

    # reading in the data
    entire_data <- vroom(
        '$merged_tsv',
        col_select = c("productive", "complete_vdj","stop_codon", "vj_in_frame", "v_frameshift", "v_call", "d_call", "j_call", "cdr3_aa", "sequence_alignment_aa", "sequence_alignment",
        "v_cigar", "j_cigar")
    ) %>%
    mutate(across(c(productive, complete_vdj, stop_codon, vj_in_frame, v_frameshift),
                as.logical))

    # categories of productivity
    entire_data %>%
        select(c(productive, stop_codon, vj_in_frame, v_frameshift)) %>%
        drop_na() %>%
        summarise(count = n(), .by = c(productive, stop_codon, vj_in_frame, v_frameshift)) %>%
        mutate(percentage = count / sum(count) * 100) %>%
        mutate(sample_id = sample_id) -> productivity_categories

    # germline V genes
    entire_data %>%
        filter(productive == TRUE) %>%
        select(v_call) %>%
        table() %>%
        as_tibble() %>%
        mutate(percentage = n / sum(n) * 100) %>%
        select(-n) %>%
        mutate(sample_id = sample_id) -> v_calls

    # germline D genes
    entire_data %>%
        filter(productive == TRUE) %>%
        select(d_call) %>%
        table() %>%
        as_tibble() %>%
        mutate(percentage = n / sum(n) * 100) %>%
        select(-n) %>%
        mutate(sample_id = sample_id) -> d_calls

    # germline J genes
    entire_data %>%
        filter(productive == TRUE) %>%
        select("j_call") %>%
        table() %>%
        as_tibble() %>%
        mutate(percentage = n / sum(n) * 100) %>%
        select(-n) %>%
        mutate(sample_id = sample_id) -> j_calls

    # clone information
    entire_data %>%
        filter(productive == TRUE) %>%
        select(c(cdr3_aa, sequence_alignment_aa, sequence_alignment)) %>%
        filter(!str_detect(sequence_alignment_aa, "X")) %>% # remove sequences with Xs
        drop_na() %>%
        summarise(count = n(),
            .by = c(cdr3_aa, sequence_alignment_aa, sequence_alignment)) %>%
        group_by(cdr3_aa) %>%
        summarise(
            cdr3_aa = cdr3_aa[1],
            sequence_alignment_aa = sequence_alignment_aa[which.max(count)],
            sequence_alignment = sequence_alignment[which.max(count)],
            cdr3_count = sum(count)) %>%
        filter(cdr3_count > 1) %>% # remove CDR3 with count of 1
        filter(nchar(sequence_alignment_aa) > 100) %>% # remove sequences with less than 100 amino acids
        mutate(proportion = cdr3_count / sum(cdr3_count)) %>%
        mutate(cdr3_cpm = proportion * 1000000) %>%
        mutate(sample_id = sample_id) %>%
        arrange(desc(cdr3_cpm)) -> clone_information

    # CDR3 samples for saturation plot data
    # make a vector of cdr3s to sample from
    all_cdr3_vector <- rep(clone_information[["cdr3_aa"]], times = clone_information[["cdr3_count"]])

    # do the sampling
    cdr3_samples <- tibble(
        sample_size = round(seq(0, length(all_cdr3_vector), length.out = 100)),
        sample_num = rep(sample_id, 100))
    cdr3_samples[["number_unique"]] <- sapply(
        cdr3_samples[["sample_size"]],
        function(x) length(unique(sample(all_cdr3_vector, size = x, replace = FALSE))))

    # calculate the derivative
    cdr3_samples[["derivative"]] <- c(0,
    diff(cdr3_samples[["number_unique"]]) / diff(cdr3_samples[["sample_size"]]))

    # nucleotide sequences
    nucleotide_sequences <- entire_data %>%
        filter(productive == TRUE) %>%
        mutate(sample_id = sample_id) %>%
        select(c(cdr3_aa, sequence_alignment, v_cigar, j_cigar))

    # write out
    vroom_write(productivity_categories, "${sequence_id}_productivity_categories.tsv")
    vroom_write(v_calls, "${sequence_id}_v_calls.tsv")
    vroom_write(d_calls, "${sequence_id}_d_calls.tsv")
    vroom_write(j_calls, "${sequence_id}_j_calls.tsv")
    vroom_write(clone_information, "${sequence_id}_clone_information.tsv")
    vroom_write(cdr3_samples, "${sequence_id}_saturation_plot_data.tsv")
    vroom_write(nucleotide_sequences, "${sequence_id}_nucleotide_sequences.tsv")
    """
}

workflow r_processing {
    take:
        igblast_tsvs

    main:
        // convert fastq to fasta (needed to run igblast)
        processed_tsv = run_r_script(igblast_tsvs)

    emit:
        processed_tsv
 }
