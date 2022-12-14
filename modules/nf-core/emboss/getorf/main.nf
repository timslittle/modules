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

process EMBOSS_GETORF {
    tag "$meta.id"
    label 'process_low'

    conda (params.enable_conda ? "bioconda::emboss=6.6.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/emboss:6.6.0--hf657eab_5':
        'quay.io/biocontainers/emboss:6.6.0--h440b012_4' }"

    input:
    tuple val(meta), path(sequence)
    val out_ext

    output:
    // TODO nf-core: Named file extensions MUST be emitted for ALL output channels
    tuple val(meta), path("*.${out_ext}"), emit: bam
    // TODO nf-core: List additional required output channels/values here
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def osformat2 = args.contains('-osformat2') ? '' : "-osformat2 ${out_ext}"
    def table = params.getorf_table ?: 0
    def minsize = params.getorf_minsize ?: 30
    def maxsize = params.getorf_maxsize ?: 1000000
    def getorf_find = params.getorf_find ?: 0
    // TODO nf-core: Where possible, a command MUST be provided to obtain the version number of the software e.g. 1.10
    //               If the software is unable to output a version number on the command-line then it can be manually specified
    //               e.g. https://github.com/nf-core/modules/blob/master/modules/nf-core/homer/annotatepeaks/main.nf
    //               Each software used MUST provide the software name and version number in the YAML version file (versions.yml)
    // TODO nf-core: It MUST be possible to pass additional parameters to the tool as a command-line string via the "task.ext.args" directive
    // TODO nf-core: If the tool supports multi-threading then you MUST provide the appropriate parameter
    //               using the Nextflow "task" variable e.g. "--threads $task.cpus"
    // TODO nf-core: Please replace the example samtools command below with your module's command
    // TODO nf-core: Please indent the command appropriately (4 spaces!!) to help with readability ;)
    """
    getorf \\
    -table ${table} \\
    -minsize ${getorf_minsize} \\
    -maxsize ${getorf_maxsize} \\
    -find ${getorf_find} \\
    -sequence ${sequence} \\
    ${outfmt2} \\
    -outseq ${prefix}.${out_ext} \\
    $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        emboss: \$(echo \$(getorf -version 2>&1) | sed 's/EMBOSS://')
    END_VERSIONS
    """
}
