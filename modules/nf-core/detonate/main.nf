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

process DETONATE {
    tag "$meta_rnaseq.id"
    label 'process_high'

    conda (params.enable_conda ? "bioconda::detonate=1.11" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/detonate:1.11--hae1ec2f_3':
        'quay.io/biocontainers/detonate:1.11--hae1ec2f_3' }"

    input:
    tuple val(meta_rnaseq), path(rnaseq)
    tuple val(meta_fasta), path(fasta)

    output:
    tuple val(meta), path("*.score"), emit: score
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta_rnaseq.id}"
    def read_length = params.read_length ?: 76 
    if(!read_length){ log.warn("detonate: Average read length not specified, defaulting to 76nt") }
    def transcript_length_parameters = params.transcript_length_parameters ? "--transcript-length-parameters ${params.transcript_length_parameters}" : ""

    if (meta_rnaseq.single_end){
        rnaseq_reads = rnaseq
    } else {
        rnaseq_reads = "--paired_end ${rnaseq}"
    }
    // TODO nf-core: Where possible, a command MUST be provided to obtain the version number of the software e.g. 1.10
    //               If the software is unable to output a version number on the command-line then it can be manually specified
    //               e.g. https://github.com/nf-core/modules/blob/master/modules/nf-core/homer/annotatepeaks/main.nf
    //               Each software used MUST provide the software name and version number in the YAML version file (versions.yml)
    // TODO nf-core: It MUST be possible to pass additional parameters to the tool as a command-line string via the "task.ext.args" directive
    """
    ./rsem-eval/rsem-eval-calculate-score \\
        $rnaseq_reads \\
        $fasta \\
        $prefix \\
        $read_length \\
        $transcript_length_parameters \\
        -p $task.cpus

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        detonate: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//' ))
    END_VERSIONS
    """
}
