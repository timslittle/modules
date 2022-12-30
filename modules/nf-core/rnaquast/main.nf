process RNAQUAST {
    tag "$meta.id"
    label 'process_medium'

    conda (params.enable_conda ? "bioconda::rnaquast=2.2.1" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/rnaquast:2.2.1--h9ee0642_0':
        'quay.io/biocontainers/rnaquast:2.2.1--h9ee0642_0' }"

    input:
    tuple val(meta), path(assembly)

    output:
    tuple val(meta), path("comparison_output/*"), emit: rnaquast_dir
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def reference = params.reference ? "--reference ${params.reference}" : ""
    def gtf = params.gtf ? "--gtf ${params.gtf}" : "" 
    def busco_lineage = params.busco_lineage ? "--busco_lineage ${params.busco_lineage}" : ""
    def gene_mark = params.gene_mark ? "--gene_mark" : ""
    def strandedness = (meta.strandedness == "unstranded") ? "" : "--strand_specific"

    // TODO nf-core: Where possible, a command MUST be provided to obtain the version number of the software e.g. 1.10
    //               If the software is unable to output a version number on the command-line then it can be manually specified
    //               e.g. https://github.com/nf-core/modules/blob/master/modules/nf-core/homer/annotatepeaks/main.nf
    //               Each software used MUST provide the software name and version number in the YAML version file (versions.yml)
    """
    python rnaQUAST.py \\
        --transcipts $assembly \\
        --threads $task.cpus \\
        --labels $prefix \\
        $reference \\
        $gtf \\
        $strandedness \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rnaquast: \$(echo \$(python rnaQuast.py -h 2>&1) | sed 's/^.*rnaQUAST //; s/Using.*\$//' ))
    END_VERSIONS
    """
}
