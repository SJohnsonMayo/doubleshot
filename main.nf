#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { FASTQC } from './modules/nf-core/fastqc/main'
include { FASTP; BBT; SORTMERNA } from './modules/local/QC/main'
include { SOURMASH; IMPORT_SOURMASH_TO_PHYLOSEQ } from './modules/local/sourmash/main'
include { HUMANN; HUMANN_COMBINE_TABLES } from './modules/local/humann/main'
include { FASTP_READ_COUNTS; BBT_READ_COUNTS } from './modules/local/stats/main'

include { validateParameters; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'

process SAVE_FASTP_SUMMARY {
    input:
    tuple val(meta), val(read_counts)

    output:
    path "fastp_summary.tsv"

    script:
    """
    echo -e "sampleID\tinitial_reads\tpost_fastp" > fastp_summary.tsv
    for sample in ${read_counts}; do
        echo -e "\${sample[0]}\t\${sample[1]}\t\${sample[2]}" >> fastp_summary.tsv
    done
    """
}

workflow {

    if (!params.input) {
        error "The 'input' parameter is not specified. Please provide it in params.yaml."
    }

    samplesheet = Channel.fromList(samplesheetToList(params.input, "assets/schema_input.json"))

    ch_input = samplesheet
        .map { meta, fastq_1, fastq_2, sample_type ->
            def new_meta = meta + [sample_type: sample_type]
            return [ new_meta, [fastq_1, fastq_2]]
        }


    FASTQC ( ch_input )

    ch_fastp = FASTP ( ch_input )

    FASTP_READ_COUNTS(ch_fastp.read_counts.collect())

    ch_bloom_filters = Channel.fromPath("${params.bloom_db_dir}*")  // This will get both .bf and .txt files
    ch_bbt = BBT(ch_fastp.processed_reads, ch_bloom_filters.collect())


    BBT_READ_COUNTS(ch_bbt.bbt_summary.collect())

    ch_bbt.filtered_reads
        .branch {
            meta, reads1, reads2 ->
            metatranscriptomic: meta.sample_type == 'metatranscriptomic'
                return [ meta, reads1, reads2 ]
            metagenomic: meta.sample_type == 'metagenomic'
                return [ meta, reads1, reads2 ]
            unknown: true  // Catch any undefined sample types
                error "Unknown sample type for sample ${meta.id}. Must be either 'metagenomic' or 'metatranscriptomic'"
        }
        .set { ch_typed_reads }


    // Create a channel for reference files
    ch_sortmerna_refs = Channel.fromPath("${params.sortmerna_db_dir}/*")

    ch_sortmerna = SORTMERNA( ch_typed_reads.metatranscriptomic, ch_sortmerna_refs.collect() )

    ch_for_profiling = ch_sortmerna.filtered_reads
        .mix(ch_typed_reads.metagenomic)


    ch_sourmash_db = file(params.sourmash_db)
    ch_sourmash_tax = file(params.sourmash_tax)
    ch_sourmash = SOURMASH ( ch_for_profiling, ch_sourmash_db, ch_sourmash_tax )

    ch_sourmash_phyloseq = IMPORT_SOURMASH_TO_PHYLOSEQ (ch_sourmash.sourmash_results.collect()) 

    ch_map_files = Channel.fromPath("${params.humann_map_dir}/*.txt.gz")
    ch_protein_db = Channel.fromPath("${params.run_mode == 'test' ? params.humann_test_protein_db : params.humann_protein_db}/*")
    ch_nuc_db = params.run_mode == 'prod' ? Channel.fromPath("${params.humann_nuc_db}/*") : Channel.value([])


    ch_humann = HUMANN ( ch_for_profiling, ch_protein_db.collect(), ch_nuc_db.collect(), ch_map_files.collect())
    
    ch_humann_merge = HUMANN_COMBINE_TABLES ( ch_humann.gene_families_rpk.collect(),
                                              ch_humann.path_abundance_rpk.collect(),
                                              ch_humann.go_rpk.collect(),
                                              ch_humann.ko_rpk.collect(),
                                              ch_humann.ec_rpk.collect(),
                                              ch_humann.gene_families_cpm.collect(),
                                              ch_humann.path_abundance_cpm.collect(),
                                              ch_humann.go_cpm.collect(),
                                              ch_humann.ko_cpm.collect(),
                                              ch_humann.ec_cpm.collect() )

}

