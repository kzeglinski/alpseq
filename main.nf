#!/usr/bin/env nextflow

/*
 * A nextflow pipeline for the pre-processing of
 * 2x300 illumina sequencing data from nanobodies
 */

version = "v0.4.0"

if(params.help == true){
log.info """
꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙
 ∩~~~∩
ξ ･×･ ξ
ξ ~  ξ   		  ★ alpseq v0.4.0 ★
ξ　  ξ              
ξ　  “~~~~~~_
ξ　          ξ
 ξ ξ ξ~~~ξ ξ
 ξ_ξξ_ξ  ξ_ξξ_ξ
꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙

Usage: nextflow run ./nanobody_preprocessing/main.nf --read_dir [input path] --sample_sheet [sample sheet]
--help                : prints this help message
Required arguments:
--out_dir             : where the output files will be written to (default: "$projectDir/results")
--read_dir            : where the input fastq files are located
--sample_sheet        : location of the .csv sample sheet (format: sample_num,library,antigen,round,replicate)
Other helpful settings:
--use_igblast         : whether or not to use IgBLAST for annotation (default: true)
--qc_only             : whether to perfom panning analysis (enrichment etc) or just quality control (default: false)
Optional (only needed for advanced users)
--igblast_databases   : location of the igblast databases (default: "$projectDir/igblast_refs/")
--adapter_r1          : pattern for trimgalore to trim off R1 (default: "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA")
--adapter_r2          : pattern for trimgalore to trim off R2 (default: "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT")
--maximum_overlap     : maximum overlap (in bp) of the two reads (default: 200)
--chunk_size          : size of chunks to use for annotation (default: 25000)
"""
System.exit(0)
}

// TODO: parameter validation
// validate that these files/directories exist
read_dir = params.read_dir
sample_sheet = params.sample_sheet
igblast_databases = params.igblast_databases
mb_scripts = params.mb_scripts
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

// boolean
use_igblast = params.use_igblast

log.info """
꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙
 ∩~~~∩
ξ ･×･ ξ
ξ ~  ξ   		★ alpseq v0.3.0 ★
ξ　  ξ              
ξ　  “~~~~~~_
ξ　          ξ
 ξ ξ ξ~~~ξ ξ
 ξ_ξξ_ξ  ξ_ξξ_ξ
꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙
★ read directory           : ${params.read_dir}
★ sample sheet             : ${params.sample_sheet}
★ output directory         : ${params.out_dir}
★ only perform qc?         : ${params.qc_only}
★ use igblast              : ${params.use_igblast}
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

// to do: allow for non gzipped input (trimgalore expects it rn)

// make a process that uses tracy to basecall all ab1 files in a directory
// then cat the resulting fastq files into a single fastq file and gzip
process basecall_sanger {
    tag "basecalling sanger"
    label 'process_low'
    container "geargenomics/tracy:latest"

    input:
        tuple val(sample_num), path(read_dir)

    output:
        tuple val(sample_num), path ('*.fastq.gz'), emit: basecalled_fastq

    script:
    """
    #!/usr/bin/env bash

    # basecall all ab1 files in the directory
    # need to loop over the files in the directory
    # and run tracy on each one
    # this is because tracy doesn't support wildcards
    for file in ${read_dir}/*.ab1; do
        # get the sample name from the file name
        file_name=\$(basename "\$file" .ab1)
        # run tracy on the file
        tracy basecall \$file --format fastq --otype consensus --output \${file_name}.fastq
    done

    # cat all the fastq files into a single file and gzip it
    cat *.fastq | gzip > ${sample_num}_basecalled.fastq.gz
    """
}

workflow{
    // read the sample sheet to associate sample names and fastq files
    // output is a tuple with [sample_num, [R1, R2]]
    parse_sample_sheet(read_dir, sample_sheet)
    sample_read_pairs = parse_sample_sheet.out.samples_with_reads
    report_sample = parse_sample_sheet.out.report_sample //association bw sample ID and report ID
    sample_sheet_checked = parse_sample_sheet.out.sample_sheet_checked.first() //sample sheet with optional columns filled out. have to use first() to make it a value channel
    sample_sheet_for_pan_reports = parse_sample_sheet.out.sample_sheet_for_pan_reports.first() // sample sheet with only the samples needed for the pan reports

    // if sanger, need to basecall the ab1 files
    if (params.sanger_mode) {
        basecall_sanger(sample_read_pairs)
        basecalled_reads = basecall_sanger.out.basecalled_fastq
        trim_merge(basecalled_reads, adapter_r1, adapter_r2, maximum_overlap)
    } else {
        trim_merge(sample_read_pairs, adapter_r1, adapter_r2, maximum_overlap)
    }

    // trim and merge the data
    trimmed_and_merged_reads = trim_merge.out.trimmed_and_merged_reads
    fastqc_reports = trim_merge.out.fastqc_reports.collect()

    // perform quality control (multiQC on the fastQC output),
    // as well as simple checking of how many reads made it
    // through the trim/merge process
    quality_control(sample_read_pairs, trimmed_and_merged_reads, fastqc_reports)

    multiqc_report = quality_control.out.multiqc_report
    percentage_passing_trim_merge = quality_control.out.percentage_passing_trim_merge
    multiqc_plots = quality_control.out.multiqc_plots

    // annotation using
    annotated_tsvs = annotation(trimmed_and_merged_reads, igblast_databases, use_igblast, mb_scripts)

    // R processing of the IgBLAST output
    r_processing(annotated_tsvs, use_igblast)


    processed_tsv = r_processing.out.processed_tsv.collect(flat: false).flatMap{it -> it}
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
        .collect(flat: false)
        .flatMap()
        .map{it -> tuple(it[0], it[1..-1].flatten().unique())}
        .set{report_data}

    // reporting
    // edit the templates to include the parameters
    original_qmd_templates = Channel.fromPath("$projectDir/modules/report/*.qmd").collect()

    // only render the pan report if qc_only is false
    if(!params.qc_only){
        prepare_report_templates(
        sample_sheet_for_pan_reports,
        original_qmd_templates,
        report_data,
        use_igblast)

        edited_qmd_templates = prepare_report_templates.out.report_templates.collect()

        render_report(
        sample_sheet_for_pan_reports,
        template_dir,
        extensions_dir,
        edited_qmd_templates,
        report_data)
    }

    render_qc_report(
        processed_tsv_for_qc_report,
        sample_sheet_checked,
        multiqc_plots,
        percentage_passing_trim_merge,
        template_dir,
        extensions_dir,
        original_qmd_templates,
        use_igblast)

}
