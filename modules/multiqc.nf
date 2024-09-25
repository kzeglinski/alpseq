// run multiqc
// adapted from nf-core module: https://github.com/nf-core/modules/blob/master/modules/nf-core/multiqc/main.nf
process multiqc {
    label 'process_medium'
    publishDir "${params.out_dir}", mode: 'copy', pattern: "multiqc_report.html"

    conda "bioconda::multiqc=1.14"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/multiqc%3A1.24.1--pyhdfd78af_0' :
        'quay.io/biocontainers/multiqc:1.24.1--pyhdfd78af_0' }"

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