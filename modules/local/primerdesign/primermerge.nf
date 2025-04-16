process PRIMER_MERGE {
    label "process_single"

    conda "bioconda::csvkit=1.0.6"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/csvkit:2.1.0--2074e44ef35da83c' :
        'community.wave.seqera.io/library/csvkit:2.1.0--89933f91a588cea2' }"

    input:
    path(summaries)

    output:
    path("*.csv")      , emit: combined
    path "versions.yml", emit: versions

    script:
    """
    # Create a directory for the first-line files
    mkdir -p first_lines
    
    # Extract header and first data row from each file
    for csv in ${summaries}; do
        # Get the header
        head -n 1 \$csv > first_lines/header_\$(basename \$csv)
        
        # If the file has more than one line, get the second line (first data row)
        # Otherwise, don't add any data row (just the header)
        if [[ \$(wc -l < \$csv) -gt 1 ]]; then
            head -n 2 \$csv | tail -n 1 >> first_lines/header_\$(basename \$csv)
        fi
    done
    
    # Stack the first-line files
    csvstack first_lines/header_* > all_primers.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        csvkit: \$(csvstack --version | head -n 1)
    END_VERSIONS
    """

    stub:
    """
    touch all_primers.csv
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        csvkit: \$(csvstack --version | head -n 1)
    END_VERSIONS
    """
}