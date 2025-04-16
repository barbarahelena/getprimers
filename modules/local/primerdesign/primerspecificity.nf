process PRIMER_SPECIFICITY {
    label 'process_single'
    tag "$meta.id"

    conda "bioconda::blast=2.12.0 bioconda::biopython=1.80 conda-forge::pandas=1.3.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/blast_biopython_pandas:88b324a4076bf67d' :
        "community.wave.seqera.io/library/blast_biopython_pandas:deaa74d652752c3a" }" 

    input:
    tuple val(meta), val(seq), path(primer_file), path(reference_genomes)

    output:
    tuple val(meta), path("*_specificity_results.csv"), emit: specificity_results
    tuple val(meta), path("*_specificity_summary.csv"), emit: specificity_summary
    tuple val(meta), path("*_off_target_amplicons.csv"), optional: true, emit: off_target_amplicons
    path "versions.yml", emit: versions

    script:
    """
    primerspecificity.py --primer_file "${primer_file}" \\
        --reference_genomes "${reference_genomes}" \\
        --meta_id "${meta.id}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(blastn -version | head -n 1 | awk '{print \$2}')
        python: \$(python3 --version | sed 's/Python //')
        biopython: \$(python3 -c "import Bio; print(Bio.__version__)")
        pandas: \$(python3 -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """

    stub:
    """
    echo "Gene,Primer_Set,Forward_Primer,Reverse_Primer,Forward_Tm,Reverse_Tm" > "${meta.id}_primer_quality.csv"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(blastn -version | head -n 1 | awk '{print \$2}')
        python: \$(python --version | sed 's/Python //')
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """
}