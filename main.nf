#!/usr/bin/env nextflow

/*
 * A nextflow pipeline for the pre-processing of
 * 2x300 illumina sequencing data from nanobodies
 */

version = "v0.3.0"

if(params.help == true){
log.info """
꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙
 ∩~~~∩
ξ ･×･ ξ
ξ ~  ξ   		★ NanoLogix ★
ξ　  ξ              v0.3.0
ξ　  “~~~~~~_
ξ　          ξ
 ξ ξ ξ~~~ξ ξ
 ξ_ξξ_ξ  ξ_ξξ_ξ
꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙

Usage: nextflow run ./nanobody_preprocessing/main.nf --fastq_dir [input path] --sample_sheet [sample sheet]
--help                : prints this help message
Required arguments:
--out_dir             : where the output files will be written to (default: "$projectDir/results)
--fastq_dir           : where the input fastq files are located
--sample_sheet        : location of the .csv sample sheet (format: sample_num,library,antigen,round,replicate)
Optional (only needed for advanced users)
--igblast_databases   : location of the igblast databases (default: "$projectDir/igblast_refs/")
--adapter_r1          : pattern for trimgalore to trim off R1 (default: "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA")
--adapter_r2          : pattern for trimgalore to trim off R2 (default: "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT")
--maximum_overlap     : maximum overlap (in bp) of the two reads (default: 200)
--chunk_size          : size of chunks to use for IgBLAST processing (default: 100000)
"""
System.exit(0)
}

// TODO: parameter validation
// validate that these files/directories exist
fastq_dir = params.fastq_dir
sample_sheet = params.sample_sheet
igblast_databases = params.igblast_databases
template_dir = params.template_dir
extensions_dir = params.extensions_dir
quarto_base_yaml = params.quarto_base_yaml

// validate that these are strings
adapter_r1 = params.adapter_r1
adapter_r2 = params.adapter_r2
sequence_trim_5p = params.sequence_trim_5p
sequence_trim_3p = params.sequence_trim_3p

// validate these are numbers with appropriate limits
maximum_overlap = params.maximum_overlap
chunk_size = params.chunk_size

log.info """
꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙
 ∩~~~∩
ξ ･×･ ξ
ξ ~  ξ   		★ NanoLogix ★
ξ　  ξ              v0.1.0
ξ　  “~~~~~~_
ξ　          ξ
 ξ ξ ξ~~~ξ ξ
 ξ_ξξ_ξ  ξ_ξξ_ξ
꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙
★ read directory           : ${params.fastq_dir}
★ sample sheet             : ${params.sample_sheet}
★ output directory         : ${params.out_dir}
"""

/*
 * Bring in modules
 */
include { parse_sample_sheet } from './subworkflows/file_import'
include { trim_merge } from './subworkflows/trim_merge'
include { quality_control } from './subworkflows/quality_control'
include { annotation } from './subworkflows/annotation'
include { r_processing } from './subworkflows/r_processing'
include { render_report } from './subworkflows/reporting'
include { prepare_report_templates } from './subworkflows/reporting'
include { render_qc_report } from './subworkflows/reporting'

/*
 * Run the workflow
 */

workflow{
    // read the sample sheet to associate sample names and fastq files
    // output is a tuple with [sample_num, [R1, R2]]
    parse_sample_sheet(fastq_dir, sample_sheet)
    sample_read_pairs = parse_sample_sheet.out.samples_R1_R2
    report_sample = parse_sample_sheet.out.report_sample //association bw sample ID and report ID

    // TODO: avoid re-running shared R0s but now cbf

    // trim and merge the data
    trim_merge(sample_read_pairs, adapter_r1, adapter_r2, maximum_overlap)
    trimmed_and_merged_reads = trim_merge.out.trimmed_and_merged_reads
    fastqc_reports = trim_merge.out.fastqc_reports.collect()

    // perform quality control (multiQC on the fastQC output),
    // as well as simple checking of how many reads made it
    // through the trim/merge process
    quality_control(sample_read_pairs, trimmed_and_merged_reads, fastqc_reports)

    multiqc_report = quality_control.out.multiqc_report
    percentage_passing_trim_merge = quality_control.out.percentage_passing_trim_merge
    multiqc_plots = quality_control.out.multiqc_plots

    // annotation using IgBLAST
    // first let's split each file into n chunks (default: 100,000)
    trimmed_and_merged_reads
        .splitFastq(by: chunk_size, file: true)
        .set{chunked_merged_reads}

    igblast_tsvs = annotation(chunked_merged_reads, igblast_databases)

    // R processing of the IgBLAST output
    r_processing(igblast_tsvs)
    processed_tsv = r_processing.out.processed_tsv.collect(flat: false).flatMap{ it -> tuple(it[0], *it[1..-1]) } // why does this work? who knows
    processed_tsv
    .flatMap{ it -> it[1..-1] }
    .collect()
    .set{processed_tsv_for_qc_report}

    report_sample
        .combine(processed_tsv, by: 0) // associate TSVs with reports
        .map{ it -> tuple(it[1], *it[2])} //don't need sample ID anymore, use *it[2] to unlist
        .groupTuple(by: 0)
        // now it's like [report_id, [v_gene_files], [nucleotide files], etc...]
        // we want all TSVs in one list so they're easy to stage
        // also this is more resilient to changes in the number/type of TSVs
        .map{it -> tuple(it[0], it[1..-1].flatten())}
        .set{report_data}

    // reporting
    // edit the templates to include the parameters
    original_qmd_templates = Channel.fromPath("$projectDir/modules/report/*.qmd").collect()

    prepare_report_templates(
        sample_sheet,
        original_qmd_templates,
        report_data)

    edited_qmd_templates = prepare_report_templates.out.report_templates.collect()
edited_qmd_templates.view()

    render_report(
        sample_sheet,
        template_dir,
        extensions_dir,
        edited_qmd_templates,
        report_data)

    render_qc_report(
        processed_tsv_for_qc_report,
        sample_sheet,
        multiqc_plots,
        percentage_passing_trim_merge,
        template_dir,
        extensions_dir,
        original_qmd_templates)

}
/* workflow.onComplete{
log.info """
꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙
 ∩~~~∩
ξ ･×･ ξ         ~~ THANK-YOU FOR RUNNING ~~
ξ ~  ξ   		★ NanoLogix ★
ξ　  ξ              v0.3.0
ξ　  “~~~~~~_
ξ　          ξ
 ξ ξ ξ~~~ξ ξ
 ξ_ξξ_ξ  ξ_ξξ_ξ
꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙
You might like to look at the following:
★ output directory                             : ${params.out_dir}
★ multi QC report                              : ${params.out_dir}/multiqc_report.html
★ how many reads passed the trimming & merging : ${params.out_dir}/percentage_passing_trim_merge.tsv
"""
}
 */
