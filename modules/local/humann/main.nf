process HUMANN {
    tag "$meta.id"
    label 'process_high'

    conda "${moduleDir}/humann3_env.yml"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/humann:3.6.1--pyh7cba7a3_1' :
        'biocontainers/humann:3.6.1--pyh7cba7a3_1' }"

    
    input:
    tuple val(meta), path(reads1), path(reads2)

    path "protein_db/*"     // Stage all protein db files into a directory called protein_db
    path "nucleotide_db/*"  // Stage all nucleotide db files into a directory called nucleotide_db
    path map_files        // Stage mapping files

    output:
    path("*_genefamilies.tsv"), emit: gene_families_rpk
    path("*_pathabundance.tsv"), emit: path_abundance_rpk
    path("*_go.tsv"), emit: go_rpk
    path("*_ko.tsv"), emit: ko_rpk
    path("*_ec.tsv"), emit: ec_rpk
    
    path("*_genefamilies.cpm.tsv"), emit: gene_families_cpm
    path("*_pathabundance.cpm.tsv"), emit: path_abundance_cpm
    path("*_go.cpm.tsv"), emit: go_cpm
    path("*_ko.cpm.tsv"), emit: ko_cpm
    path("*_ec.cpm.tsv"), emit: ec_cpm

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = params.run_mode == 'test' ? params.tools.humann.test.args : params.tools.humann.prod.args
    def nuc_db = params.run_mode == 'prod' ? "--nucleotide-database nucleotide_db" : ""
 

    """
    cat ${reads1} ${reads2} > combined.fastq
    humann3 -vvv \\
        --input combined.fastq \\
        --output . \\
        --output-basename ${prefix} \\
        --threads 16 \\
        --protein-database protein_db \\
        ${nuc_db} \\
        ${args} \\

    rm combined.fastq
  
    humann_regroup_table -i ${prefix}_genefamilies.tsv -c map_go_uniref90.txt.gz --output ${prefix}_go.tsv
    humann_regroup_table -i ${prefix}_genefamilies.tsv -c map_ko_uniref90.txt.gz --output ${prefix}_ko.tsv
    humann_regroup_table -i ${prefix}_genefamilies.tsv -c map_level4ec_uniref90.txt.gz --output ${prefix}_ec.tsv

    humann_renorm_table -i ${prefix}_genefamilies.tsv -p -o ${prefix}_genefamilies.cpm.tsv
    humann_renorm_table -i ${prefix}_pathabundance.tsv -p -o ${prefix}_pathabundance.cpm.tsv
    humann_renorm_table -i ${prefix}_go.tsv -p -o ${prefix}_go.cpm.tsv
    humann_renorm_table -i ${prefix}_ec.tsv -p -o ${prefix}_ec.cpm.tsv
    humann_renorm_table -i ${prefix}_ko.tsv -p -o ${prefix}_ko.cpm.tsv
    
    """
}

process HUMANN_COMBINE_TABLES {
    label 'process_low'

    conda "${moduleDir}/humann3_env.yml"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/humann:3.6.1--pyh7cba7a3_1' :
        'biocontainers/humann:3.6.1--pyh7cba7a3_1' }"
    
    input:
    path(gene_families_rpk)
    path(path_abundance_rpk)
    path(go_rpk)
    path(ko_rpk)
    path(ec_rpk)

    path(gene_families_cpm)
    path(path_abundance_cpm)
    path(go_cpm)
    path(ko_cpm)
    path(ec_cpm)

    output:
    path("combined_*"), emit: humann_combined

    script:
    """
    humann_join_tables \
        --input . \
        --output combined_gene_families.rpk.tsv \
        --file_name genefamilies.tsv

    humann_join_tables \
        --input . \
        --output combined_path_abundance.rpk.tsv \
        --file_name pathabundance.tsv

    humann_join_tables \
        --input . \
        --output combined_go.rpk.tsv \
        --file_name go.tsv

    humann_join_tables \
        --input . \
        --output combined_ko.rpk.tsv \
        --file_name ko.tsv

    humann_join_tables \
        --input . \
        --output combined_ec.rpk.tsv \
        --file_name ec.tsv



    humann_join_tables \
        --input . \
        --output combined_gene_families.cpm.tsv \
        --file_name genefamilies.cpm.tsv

    humann_join_tables \
        --input . \
        --output combined_path_abundance.cpm.tsv \
        --file_name pathabundance.cpm.tsv

    humann_join_tables \
        --input . \
        --output combined_go.cpm.tsv \
        --file_name go.cpm.tsv

    humann_join_tables \
        --input . \
        --output combined_ko.cpm.tsv \
        --file_name ko.cpm.tsv

    humann_join_tables \
        --input . \
        --output combined_ec.cpm.tsv \
        --file_name ec.cpm.tsv

    """
}
