include { trimgalore } from '../modules/trimgalore'
include { trimgalore_sanger } from '../modules/trimgalore'
include { flash } from '../modules/flash'

workflow trim_merge {
    take:
        sample_read_pairs
        adapter_r1
        adapter_r2
        maximum_overlap

    main:
        if (params.sanger_mode) {
            // trim + qc using trimgalore
            trimgalore_sanger(sample_read_pairs, adapter_r1)
            // no need to merge reads
            trimmed_and_merged_reads = trimgalore_sanger.out.reads
            fastqc_reports = trimgalore_sanger.out.fastqc_zip
        } else {
            // trim + qc using trimgalore
            trimgalore(sample_read_pairs, adapter_r1, adapter_r2)
            trimmed_reads = trimgalore.out.reads
            fastqc_reports = trimgalore.out.fastqc_zip
            // merging using FLASH
            flash(trimmed_reads, maximum_overlap)
            trimmed_and_merged_reads = flash.out.combined_reads
        }

    emit:
        trimmed_and_merged_reads
        fastqc_reports
}