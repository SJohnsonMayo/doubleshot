process FASTP {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/fastp_env.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastp:0.23.2--h5f740d0_3' :
        'biocontainers/fastp:0.23.2--h5f740d0_3' }"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}_1.fastq.gz"), path("${meta.id}_2.fastq.gz"), emit: processed_reads
    path "${meta.id}_fastp.html", emit: html
    path "${meta.id}_fastp.json", emit: json
    path "${meta.id}_read_counts.txt", emit: read_counts

    script:
    """

    
    fastp \\
        --in1 ${reads[0]} --in2 ${reads[1]} \\
        --out1 ${meta.id}_1.fastq.gz --out2 ${meta.id}_2.fastq.gz \\
        --html ${meta.id}_fastp.html --json ${meta.id}_fastp.json


    initial_reads=\$(grep '"total_reads":' ${meta.id}_fastp.json | awk -F":" 'NR==1 {print \$2/2}' | tr -d ',')
    post_fastp=\$(grep '"total_reads":' ${meta.id}_fastp.json | awk -F":" 'NR==2 {print \$2/2}' | tr -d ',')

    read_counts="${meta.id}\t\${initial_reads}\t\${post_fastp}"
    echo \$read_counts > ${meta.id}_read_counts.txt
    """
}

process BBT {
    tag "$meta.id"
    label 'process_medium'
    label 'process_biobloomtools'

    conda "${moduleDir}/bbt_env.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biobloomtools:2.3.5--h077b44d_6' :
        'biocontainers/biobloomtools:2.3.5--h077b44d_6' }"

    input:
    tuple val(meta), path(reads1), path(reads2)
    path bloom_filter_files

    output:
    tuple val(meta), path("${meta.id}_noMatch_1.fq"), path("${meta.id}_noMatch_2.fq"), emit: filtered_reads
    path("${meta.id}_summary.tsv"), emit: bbt_summary

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def bf_files = bloom_filter_files.findAll { it.toString().endsWith('.bf') }.join(' ')
    
    """
    # Assuming you have a reference bloom filter file
    biobloomcategorizer -t 16 -p ${prefix} -f "${bf_files}"  -e -i --fq ${reads1} ${reads2}

    """
}

process SORTMERNA {
    tag "${meta.id}"
    label 'process_high'

    conda "${moduleDir}/sortmerna_env.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sortmerna:4.3.7--hdbdd923_1' :
        'biocontainers/sortmerna:4.3.7--hdbdd923_1' }"
 
    input:
    tuple val(meta), path(reads1), path(reads2)
    path ref_files

    output:
    tuple val(meta), path("*.non_rRNA_fwd.fq"), path("*.non_rRNA_rev.fq") , emit: filtered_reads
    tuple val(meta), path("*.rRNA_fwd.fq"), path("*.rRNA_rev.fq"), emit: rRNA_reads
    tuple val(meta), path("*.log"), emit: logs
    path "versions.yml", emit: versions

    script:
    def args = task.cpus ? "--threads ${task.cpus}" : ""
    def prefix = meta.id
    
    // Get just the .fasta files for the --ref arguments
   
    def ref_args = ref_files.findAll { it.toString().endsWith('.fasta') }
                           .collect { "--ref ${it}" }
                           .join(" ") 
    """
    # Run SortMeRNA
    sortmerna \\
        ${ref_args} \\
        --reads ${reads1} \\
        --reads ${reads2} \\
        --workdir sorterna_tmp \\
        --aligned ${prefix}.rRNA \\
        --other ${prefix}.non_rRNA \\
        --fastx \\
        --out2 \\
        ${args}

    # Compress outputs
    ## gzip ${prefix}.non_rRNA.fastq
    ## gzip ${prefix}.rRNA.fastq
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        sortmerna: \$(sortmerna --version 2>&1 | sed 's/SortMeRNA version //; s/,.*\$//')
    END_VERSIONS
    """
}
