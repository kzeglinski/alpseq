// run multiqc
// adapted from nf-core module: https://github.com/nf-core/modules/blob/master/modules/nf-core/multiqc/main.nf
process multiqc {
    label 'process_medium'
    publishDir "${params.out_dir}", mode: 'copy', pattern: "multiqc_report.html"

    conda "bioconda::multiqc=1.14"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/pip_multiqc:3f9d6693c55189c9' :
        'community.wave.seqera.io/library/pip_multiqc:97bcaba1a9b5034c' }"

    input:
    path  multiqc_files, stageAs: "?/*"

    output:
    path "multiqc_report.html"      , emit: report
    path "multiqc_plots/png/*.png"  , emit: plots
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    """
    multiqc \\
        --force \\
        --export \\
        $args \\
        .

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$( multiqc --version | sed -e "s/multiqc, version //g" )
    END_VERSIONS
    """
}