// TODO nf-core: If in doubt look at other nf-core/modules to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/modules/nf-core/
//               You can also ask for help via your pull request or on the #modules channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A module file SHOULD only define input and output files as command-line parameters.
//               All other parameters MUST be provided using the "task.ext" directive, see here:
//               https://www.nextflow.io/docs/latest/process.html#ext
//               where "task.ext" is a string.
//               Any parameters that need to be evaluated in the context of a particular sample
//               e.g. single-end/paired-end data MUST also be defined and evaluated appropriately.
// TODO nf-core: Software that can be piped together SHOULD be added to separate module files
//               unless there is a run-time, storage advantage in implementing in this way
//               e.g. it's ok to have a single module for bwa to output BAM instead of SAM:
//                 bwa mem | samtools view -B -T ref.fasta
// TODO nf-core: Optional inputs are not currently supported by Nextflow. However, using an empty
//               list (`[]`) instead of a file can be used to work around this issue.

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
    tuple val(meta), path("*.bam"), emit: bam
    path "versions.yml"           , emit: versions

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
        rnaquast: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' ))
    END_VERSIONS
    """
}
