process PRIMER_SELECTION {
    label 'process_single'
    tag "$meta.id"

    conda "bioconda::blast=2.12.0 bioconda::biopython=1.80 conda-forge::pandas=1.3.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/blast_biopython_pandas:88b324a4076bf67d' :
        "community.wave.seqera.io/library/blast_biopython_pandas:deaa74d652752c3a" }" 

    input:
    tuple val(meta), val(seq), path(primer_file), path(quality_file), path(specificity_file)

    output:
    tuple val(meta), path("${meta.id}_selected_primers.csv"), emit: selected_primers
    path "versions.yml", emit: versions

    script:
    """
    primerselection.py \\
        --primer_file "${primer_file}" \\
        --quality_file "${quality_file}" \\
        --specificity_file "${specificity_file}" \\
        --meta_id "${meta.id}" \\
        --output_file "${meta.id}_selected_primers.csv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version | sed 's/Python //')
        biopython: \$(python3 -c "import Bio; print(Bio.__version__)")
        pandas: \$(python3 -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """

    stub:
    """
    # Placeholder for the actual script
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(blastn -version | head -n 1 | awk '{print \$2}')
        python: \$(python --version | sed 's/Python //')
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """
}