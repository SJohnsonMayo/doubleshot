process SOURMASH {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/sourmash_env.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/sourmash:4.8.11--hdfd78af_0' :
        'biocontainers/sourmash:4.8.11--hdfd78af_0' }"

    input:
    tuple val(meta), path(reads1), path(reads2)
    path db_file    // Stage the zip database file
    path tax_file   // Stage the taxonomy database file

    output:
    tuple val(meta), path("*.sig"), emit: signatures
    path("*.csv"), emit: sourmash_results

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Compute sourmash signature

    sourmash sketch dna -p k=21,k=31,k=51,scaled=1000,abund --merge ${prefix} -o ${prefix}.sig ${reads1} ${reads2}
    sourmash gather ${prefix}.sig ${db_file} -o ${prefix}_gather.csv
    sourmash tax annotate -g ${prefix}_gather.csv -t ${tax_file}
    """
}

process IMPORT_SOURMASH_TO_PHYLOSEQ {
    label 'process_low'

    conda "${moduleDir}/r-sourmash_env.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-aeed62b5d197f017b6a55f1566f3f84e56c9b3f3:b0c847e4fb89c343b04036e33b2daa19c4152cf5-0' :
        'biocontainers/mulled-v2-aeed62b5d197f017b6a55f1566f3f84e56c9b3f3:b0c847e4fb89c343b04036e33b2daa19c4152cf5-0' }"
    
    input:
    path(sourmash_results)

    output:
    path("sourmash_phyloseq.rds")

    script:
    """
    #!/usr/bin/env Rscript
    
    library(tidyverse)
    library(phyloseq)
    
    sourmash_taxonomy_results <- Sys.glob("*with-lineages.csv") %>%
        map_dfr(read_csv, col_types = "ddddddddcccddddcccdc") %>%
        mutate(name = gsub(" .*", "", name))

    tax_table <- sourmash_taxonomy_results %>%
        select(name, lineage) %>%
        distinct() %>%
        separate(lineage, into = c("domain", "phylum", "class", "order", "family", "genus", "species"), sep = ";") %>%
        column_to_rownames("name")

    scaled = 1000

    gather_table <- sourmash_taxonomy_results %>%
        mutate(n_unique_kmers = (unique_intersect_bp / scaled) * average_abund) %>% # calculate the number of uniquely matched k-mers and multiply it by average k-mer abundance
        select(query_name, name, n_unique_kmers) %>% # select only the columns that have information we need
        pivot_wider(id_cols = name, names_from = query_name, values_from = n_unique_kmers) %>% # transform to wide format
        replace(is.na(.), 0) %>% # replace all NAs with 0
        column_to_rownames("name") # move the metagenome sample name to a rowname


    phylo <- phyloseq(otu_table(gather_table, taxa_are_rows = T),
        tax_table(as.matrix(tax_table))
        )
    
    saveRDS(phylo, "sourmash_phyloseq.rds")
    """
}
