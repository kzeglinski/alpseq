#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * A nextflow pipeline for the pre-processing of
 * 2x300 illumina sequencing data from nanobodies
 */

version = "v0.5.0"

/*
 * help message
 */

if(params.help == true){
log.info """
꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙
 ∩~~~∩
ξ ･×･ ξ
ξ ~  ξ   		  ★ alpseq v0.5.0 ★
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
--fwr4_seq            : nt sequence of fwr4 (only used for annotation with matchbox)
--maximum_overlap     : maximum overlap (in bp) of the two reads (default: 200)
--chunk_size          : size of chunks to use for annotation (default: 25000)
"""
System.exit(0)
}

/*
 * initial log message
 */

log.info """
꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙꧙
 ∩~~~∩
ξ ･×･ ξ
ξ ~  ξ   		★ alpseq v0.5.0 ★
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
 * parameter validation
 */

// can't work out how to have these functions in a separate file so they just have to be here for now >:(
/*
 * Validate a filesystem path:
 *   - type 'file'  → must exist and be a file
 *   - type 'dir'   → must exist and be a directory
 */
def validate_path = { String name, def path, String type ->
    if (!path) {
        exit 1, "ERROR: Missing required parameter --${name}"
    }

    def obj = file(path)

    if (type == "file") {
        if (!obj.exists() || !obj.isFile()) {
            exit 1, "ERROR: Required file does not exist or is not a file: --${name} (${path})"
        }
    } else if (type == "dir") {
        if (!obj.exists() || !obj.isDirectory()) {
            exit 1, "ERROR: Required directory does not exist or is not a directory: --${name} (${path})"
        }
    } else {
        exit 1, "ERROR: Unknown validation type '${type}' for --${name}; expected 'file' or 'dir'"
    }
}

/*
 * Validate that a parameter is a DNA sequence:
 *   - must be a string
 *   - only A/T/G/C allowed (case-insensitive)
 */
def validate_dna_seq = { String name, def seq ->
    if (!seq) {
        exit 1, "ERROR: Missing required parameter --${name}"
    }

    if (!(seq instanceof String)) {
        exit 1, "ERROR: Parameter --${name} must be a string, but got: ${seq.getClass().getSimpleName()}"
    }

    // Regex ensures only A/T/G/C (case insensitive)
    if (!(seq ==~ /(?i)^[ATGC]+$/)) {
        exit 1, """|
        ERROR: Invalid DNA sequence for --${name}
        Provided: '${seq}'
        Allowed characters: A, T, G, C (case-insensitive)
        """.stripIndent()
    }
}

/*
 * Validate that a parameter is a DNA sequence:
 *   - parameter must exist
 *   - empty string ("") is allowed and considered valid
 *   - if non-empty, only A/T/G/C allowed (case-insensitive)
 */
def validate_dna_seq_allow_empty = { String name, def seq ->
    if (seq == null) {
        exit 1, "ERROR: Missing required parameter --${name}"
    }

    if (!(seq instanceof String)) {
        exit 1, "ERROR: Parameter --${name} must be a string, but got: ${seq.getClass().getSimpleName()}"
    }

    // Allow empty string
    if (seq.trim() == '') {
        return  // empty DNA sequence is considered valid
    }

    // Validate non-empty sequences
    if (!(seq ==~ /(?i)^[ATGC]+$/)) {
        exit 1, """|
        ERROR: Invalid DNA sequence for --${name}
        Provided: '${seq}'
        Allowed characters: A, T, G, C (case-insensitive), or empty string
        """.stripIndent()
    }
}

/*
 * Validate that a parameter is boolean (true/false)
 */
def validate_bool = { String name, def value ->
    if (value == null) {
        exit 1, "ERROR: Missing required boolean parameter --${name}"
    }

    if (!(value instanceof Boolean)) {
        exit 1, "ERROR: Parameter --${name} must be boolean (true/false), but got: ${value} (${value.getClass().getSimpleName()})"
    }
}

/*
 * Validate that a parameter is a number within a specified inclusive range
 * min and max may be ints or floats
 */
def validate_number = { String name, def value, def min, def max ->
    if (value == null) {
        exit 1, "ERROR: Missing required numeric parameter --${name}"
    }

    // Accept Integer or Float/Double
    if (!(value instanceof Number)) {
        exit 1, "ERROR: Parameter --${name} must be a number, but got: ${value} (${value.getClass().getSimpleName()})"
    }

    if (value < min || value > max) {
        exit 1, "ERROR: Parameter --${name} must be between ${min} and ${max} (inclusive), but got: ${value}"
    }
}

// files or directories that must exist    
validate_path('sample_sheet', params.sample_sheet, 'file')
validate_path('template_dir', params.template_dir, 'file') // kind of dumb the param is called template_dir but its a tar, should change at some point
validate_path('extensions_dir', params.extensions_dir, 'file')
validate_path('mb_scripts', params.mb_scripts, 'dir')
validate_path('igblast_databases', params.igblast_databases, 'dir')
validate_path('read_dir', params.read_dir, 'dir')

// validate DNA sequence parameters
validate_dna_seq('adapter_r1', params.adapter_r1)
validate_dna_seq('adapter_r2', params.adapter_r2)
validate_dna_seq('fwr4_seq', params.fwr4_seq)

// DNA sequence or empty string allowed
validate_dna_seq_allow_empty('sequence_trim_5p', params.sequence_trim_5p)
validate_dna_seq_allow_empty('sequence_trim_3p', params.sequence_trim_3p)

// numbers, within expected range
validate_number('mb_error_rate', params.mb_error_rate, 0, 0.4)
validate_number('num_v_genes', params.num_v_genes, 1, 999)
validate_number('maximum_overlap', params.maximum_overlap, 1, 999)
validate_number('chunk_size', params.chunk_size, 100, 100000)

// boolean
validate_bool('use_igblast', params.use_igblast)
validate_bool('sanger_mode', params.sanger_mode)
validate_bool('enable_conda', params.enable_conda)
validate_bool('help', params.help)
validate_bool('qc_only', params.qc_only)

log.info "✓ Parameter validation successful!"

/*
 * bring in modules
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
 * run the workflow
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
    parse_sample_sheet(params.read_dir, params.sample_sheet)
    sample_read_pairs = parse_sample_sheet.out.samples_with_reads
    report_sample = parse_sample_sheet.out.report_sample //association bw sample ID and report ID
    sample_sheet_checked = parse_sample_sheet.out.sample_sheet_checked.first() //sample sheet with optional columns filled out. have to use first() to make it a value channel
    sample_sheet_for_pan_reports = parse_sample_sheet.out.sample_sheet_for_pan_reports.first() // sample sheet with only the samples needed for the pan reports

    // if sanger, need to basecall the ab1 files
    if (params.sanger_mode) {
        basecall_sanger(sample_read_pairs)
        basecalled_reads = basecall_sanger.out.basecalled_fastq
        trim_merge(basecalled_reads, params.adapter_r1, params.adapter_r2, params.maximum_overlap)
    } else {
        trim_merge(sample_read_pairs, params.adapter_r1, params.adapter_r2, params.maximum_overlap)
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
    annotated_tsvs = annotation(trimmed_and_merged_reads, params.igblast_databases, params.use_igblast, params.mb_scripts, params.fwr4_seq, params.mb_error_rate, params.num_v_genes, params.igblast_db_name)

    // R processing of the IgBLAST output
    r_processing(annotated_tsvs, params.use_igblast)


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
        params.use_igblast)

        edited_qmd_templates = prepare_report_templates.out.report_templates.collect()

        render_report(
        sample_sheet_for_pan_reports,
        params.template_dir,
        params.extensions_dir,
        edited_qmd_templates,
        report_data)
    }

    render_qc_report(
        processed_tsv_for_qc_report,
        sample_sheet_checked,
        multiqc_plots,
        percentage_passing_trim_merge,
        params.template_dir,
        params.extensions_dir,
        original_qmd_templates,
        params.use_igblast)

}
