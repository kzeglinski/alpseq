// adapted from https://github.com/stevekm/nextflow-demos/blob/master/parse-samplesheet/main.nf
workflow parse_sample_sheet {
    take:
        fastq_dir
        sample_sheet

    main:
        sample_sheet_channel = Channel.fromPath(file(sample_sheet))
        check_sample_sheet(sample_sheet_channel)

        sample_sheet_checked = check_sample_sheet.out.sample_sheet_checked
        sample_sheet_deduplicated = check_sample_sheet.out.sample_sheet_deduplicated

        sample_sheet_deduplicated
        .splitCsv(header: true, sep: ',')
        .map{row ->
            def sample_ID = row['sample_num']
            def reads1 = "${fastq_dir}/" + row['r1_file_name']
            def reads2 = "${fastq_dir}/" + row['r2_file_name']
            return [ sample_ID, [reads1, reads2]]
        }
        .set { samples_R1_R2 } // set of all fastq R1 R2 per sample

        sample_sheet_checked
        .splitCsv(header: true, sep: ',')
        .map{row ->
            def report_ID = row['report_id']
            def sample_ID = row['sample_num']
            return [ sample_ID, report_ID ]
        }
        .set { report_sample }

    emit:
        samples_R1_R2
        report_sample
        sample_sheet_checked

}

process check_sample_sheet {
    tag "checking sample sheet"
    label 'process_low'
    container "library://kzeglinski/nanologix/nanologix-report:v0.4.0"

    input:
    path(sample_sheet)

    output:
    path('sample_sheet_checked.csv'), emit: sample_sheet_checked
    path('sample_sheet_no_duplicates.csv'), emit: sample_sheet_deduplicated

    script:
    """
    #!/usr/bin/env Rscript
    sample_sheet_file <- '$sample_sheet'

    # sample sheet checker
    library(vroom)

    # check that sample sheet exists
    if(!file.exists(sample_sheet_file)){
        stop("Sample sheet does not exist at ", sample_sheet_file, "\n Please check the path and try again.")
    }

    # read it in
    sample_sheet <- vroom(sample_sheet_file)

    # check all columns present (need either panning_id OR report_id not both)
    columns_ok <- all(c("sample_num", "round", "r1_file_name", "r2_file_name") %in% colnames(sample_sheet))

    if(!columns_ok){
        stop("Sample sheet is missing one or more of the required columns: sample_num, round, r1_file_name, r2_file_name")
    }

    columns_ok <- "panning_id" %in% colnames(sample_sheet) || "report_id" %in% colnames(sample_sheet)

    if(!columns_ok){
        stop("Sample sheet needs either a panning_id column (to generate a single report containing multiple different pans) or a report_id column (to generate multiple reports for each pan)")
    }

    # if no panning_id column, make it all ones
    if(!"panning_id" %in% colnames(sample_sheet)){
        sample_sheet[["panning_id"]] <- 1
    }

    # if no report_id column, make it all ones
    if(!"report_id" %in% colnames(sample_sheet)){
        sample_sheet[["report_id"]] <- 1
    }

    # if no replicate column, make it all NAs
    if(!"replicate" %in% colnames(sample_sheet)){
        sample_sheet[["replicate"]] <- NA
    }

    # check that round is numeric
    if(!all(is.numeric(sample_sheet[["round"]]))){
        stop("Round column must be numeric")
    }

    # check same sample numbers have the same library, antigen, round
    for(i in seq_along(unique(sample_sheet[["report_id"]]))){
            report_id <- unique(sample_sheet[["report_id"]])[i]
            report_subset <- sample_sheet[sample_sheet[["report_id"]] == report_id,]
        for(j in seq_along(unique(sample_sheet[["sample_num"]]))){
            sample_num <- unique(report_subset[["sample_num"]])[j]
            sample_subset <- report_subset[report_subset[["sample_num"]] == sample_num,]
            if(length(unique(sample_subset[["library"]])) > 1){
                stop("Sample number ", sample_num, " has multiple libraries. Each sample number should be unique to a particular library/antigen/round combination")
            }
            if(length(unique(sample_subset[["antigen"]])) > 1){
                stop("Sample number ", sample_num, " has multiple antigens. Each sample number should be unique to a particular library/antigen/round combination")
            }
            if(length(unique(sample_subset[["round"]])) > 1){
                stop("Sample number ", sample_num, " has multiple rounds. Each sample number should be unique to a particular library/antigen/round combination")
            }
        }
    }

    # make sample sheet with no duplicate sample IDs (for pre-processing)
    sample_sheet_no_duplicates <- sample_sheet[!duplicated(sample_sheet[["sample_num"]]),]
    vroom_write(sample_sheet_no_duplicates, "sample_sheet_no_duplicates.csv", delim = ",")
    vroom_write(sample_sheet, "sample_sheet_checked.csv", delim = ",")
    """
}
