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
    echo -e "SampleID\tGRCm49\tphiX\thuman_hg38\tmultiMatch\tnoMatch" > bbt_read_counts.txt
    for file in ${bbt_summaries}; do
	grcm49=\$(grep "GRCm49" "\$file" | awk '{print \$2}')
	phix=\$(grep "phiX" "\$file" | awk '{print \$2}')
	human=\$(grep "humann_hg38" "\$file" | awk '{print \$2}')
	multiMatch=\$(grep "multiMatch" "\$file" | awk '{print \$2}')
        noMatch=\$(grep "noMatch" \$file | awk '{print \$2}')
        id=\$(basename \$file _summary.tsv)
	echo -e "\$id\t\$grcm49\t\$phix\t\$human\t\$multiMatch\t\$noMatch" >> bbt_read_counts.txt
    done
    """
}

