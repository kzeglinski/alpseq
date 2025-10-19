process render_report {
    tag "creating report"
    label 'process_high'
    publishDir "${params.out_dir}/report/", mode: 'copy', pattern: "*.html"
   // stageInMode 'copy'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://ghcr.io/kzeglinski/alpseq:0.4.0' :
        'ghcr.io/kzeglinski/alpseq:0.4.0' }"

    input:
    path(sample_sheet)
    path(template_dir)
    path(extensions_dir)
    path(qmd_templates)
    tuple val(report_id), path(report_data)

    output:
    path('*.html'), emit: report

    script:
    """
    #!/usr/bin/env bash

    # need the PWD so the caches don't clash between jobs
    export DENO_DIR="\$PWD"
    export XDG_CACHE_HOME="\$PWD/tmp/quarto_cache_home"
    export XDG_DATA_HOME="\$PWD/tmp/quarto_data_home"
    export XDG_RUNTIME_DIR="\$PWD/tmp/quarto_runtime_dir"

    tar -xvf _extensions.tar
    tar -xvf _template.tar

    quarto render "${report_id}_report.qmd" --log wnp_report.log --output "${report_id}_report.html"

    """
}

process prepare_report_templates {
    tag "preparing report templates"
    label 'process_low'
    stageInMode 'copy'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://ghcr.io/kzeglinski/alpseq:0.4.0' :
        'ghcr.io/kzeglinski/alpseq:0.4.0' }"


    input:
    path(sample_sheet)
    path(qmd_templates)
    tuple val(report_id), path(report_data)
    val(use_igblast)

    output:
    path('*.qmd'), emit: report_templates

    script:
    """
    #!/usr/bin/env Rscript
    library(stringr)
    library(dplyr)
    library(tidyr)
    library(vroom)

    # read in the sample sheet
    metadata <- vroom(fs::dir_ls(glob = "*.csv"),
        col_select = c("sample_num", "library", "antigen", "round", "replicate", "panning_id", "report_id")) %>%
        filter(report_id == "${report_id}")

    if(length(unique(metadata[["panning_id"]])) > 1){
        stop("More than one panning ID found in the sample sheet. Comparisons are not yet supported.")
    } else{
        panning_id <- unique(metadata[["panning_id"]])
    }
    # create replicate ID if required
    if(!all(is.na(metadata[["replicate"]]))){
        metadata <- metadata %>%
            filter(round != 0) %>% # temporarily remove R0 samples
            group_by(library, antigen, round) %>%
            mutate(replicate_id = cur_group_id()) %>%
            # and give our replicates an informative name, not just a number
            mutate(replicate_id_informative =
                paste0(library, "_", antigen, "_round_", round)) %>%
            ungroup() %>%
            bind_rows(filter(metadata, round == 0))
    }

    # wnp_report
    if (as.logical('$use_igblast')) {
        template = readLines("template_wnp_report.qmd")
    } else {
        template = readLines("template_wnp_report_mb.qmd")
    }
    template %>%
        stringr::str_replace(pattern = "param_analysis_name", replace = "${report_id}") %>%
        stringr::str_replace(pattern = "param_panning_id", replace = as.character(panning_id)) %>%
        writeLines(con = "${report_id}_report.qmd")

    """
}

process render_qc_report {
    tag "creating QC report"
    label 'process_low'
    publishDir "${params.out_dir}/report/", mode: 'copy', pattern: "*.html"
   // stageInMode 'copy'
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://ghcr.io/kzeglinski/alpseq:0.4.0' :
        'ghcr.io/kzeglinski/alpseq:0.4.0' }"

    input:
    path(processed_tsv_for_qc_report)
    path(sample_sheet)
    path(multiqc_plots)
    path(percentage_passing_trim_merge)
    path(template_dir)
    path(extensions_dir)
    path(qmd_templates)
    val(use_igblast)

    output:
    path('*.html'), emit: report

    script:
    """
    #!/usr/bin/env bash

    export DENO_DIR="\$PWD"
    export XDG_CACHE_HOME="\$PWD/tmp/quarto_cache_home"
    export XDG_DATA_HOME="\$PWD/tmp/quarto_data_home"
    export XDG_RUNTIME_DIR="\$PWD/tmp/quarto_runtime_dir"

    tar -xvf _extensions.tar
    tar -xvf _template.tar

    if [ $use_igblast ]; then
        quarto render qc_report.qmd --log qc_report.log
    else
        quarto render qc_report_mb.qmd --log qc_report.log
    fi
    """
}