process FASTP_READ_COUNTS {
    input:
    path(fastp_read_counts_files)

    output:
    path "fastp_read_counts.txt"

    script:
    """
    # Start with base header
    echo -e "SampleID\tinitial_reads\tpost_fastp" > fastp_read_counts.txt
    
    # Add FASTP read counts
    for file in ${fastp_read_counts_files}; do
        cat \$file >> fastp_read_counts.txt
    done
    
    """
}

process BBT_READ_COUNTS {
    input:
    path(bbt_summaries)

    output:
    path "bbt_read_counts.txt"

    script:
    """
    # Extract sample ID and BBT read count from BBT log
    echo -e "SampleID\tpost_bbt" > bbt_read_counts.txt
    for file in ${bbt_summaries}; do
        misses=\$(grep "noMatch" \$file | awk '{print \$2}')
        id=\$(basename \$file _summary.tsv)
        echo -e "\$id\t\$misses" >> bbt_read_counts.txt
    done
    """
}

