process PRIMER_DESIGN {
    label 'process_single'
    tag "$meta.id"
    
    conda "bioconda::primer3-py=2.0.3 bioconda::biopython=1.85 conda-forge::pandas=2.2.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/primer3-py_pandas_pip_biopython:e1f5eb2f1ba237b8' :
        'community.wave.seqera.io/library/primer3-py_pandas_pip_biopython:e797e4d1477d8e75' }"

    input:
    tuple val(meta), val(gene_seq)
    val(tf_fow)
    val(tf_rev)
    val(primer_opt_size)
    val(primer_min_size)
    val(primer_max_size)
    val(primer_min_tm)
    val(primer_max_tm)
    val(primer_pair_max_diff_tm)
    val(primer_product_size_range)
    val(primer_task)
    val(primer_opt_tm)
    val(primer_min_gc)
    val(primer_max_gc)
    val(primer_n_left)
    val(primer_n_right)

    output:
    tuple val(meta), val(gene_seq), path("${meta.id}_primers.csv")   , emit: primers
    path "versions.yml"                                              , emit: versions

    script:
    def sequence = gene_seq.first().toLowerCase()
    def TF_FOW = tf_fow ?: 'gaaattaatacgactcactatagggagaccacaacggtttccctctagaaataattttgtttaagaaggagatatacat'
    def TF_REV_COMP = tf_rev ?: 'aggcttacctaatgccaactttgtacaagaaagctg'
    def PRIMER_OPT_SIZE = primer_opt_size ?: 20
    def PRIMER_MIN_SIZE = primer_min_size ?: PRIMER_OPT_SIZE - 2
    def PRIMER_MAX_SIZE = primer_max_size ?: PRIMER_OPT_SIZE + 7
    def PRIMER_OPT_TM = primer_opt_tm ?: 60
    def PRIMER_MIN_TM = primer_min_tm ?: PRIMER_OPT_TM - 3
    def PRIMER_MAX_TM = primer_max_tm ?: PRIMER_OPT_TM + 3
    def PRIMER_PAIR_MAX_DIFF_TM = primer_pair_max_diff_tm ?: 3.0
    def PRIMER_PRODUCT_SIZE_RANGE = primer_product_size_range ?: '[[75,150]]'
    def PRIMER_TASK = primer_task ?: 'generic'
    def PRIMER_MIN_GC = primer_min_gc ?: 30
    def PRIMER_MAX_GC = primer_max_gc ?: 70
    def PRIMER_PICK_LEFT_PRIMER = primer_n_left ?: 1
    def PRIMER_PICK_RIGHT_PRIMER = primer_n_right ?: 1

    """
    python3 <<EOF
    import os
    import sys
    import primer3
    import pandas as pd
    from Bio.Seq import Seq

    # Define universal primer sequences
    TF_FOW = "${TF_FOW}"
    TF_REV_COMP = "${TF_REV_COMP}"

    # Configure primer3 parameters
    primers = primer3.design_primers(
        {
            'SEQUENCE_ID': "${meta.id}",
            'SEQUENCE_TEMPLATE': "${sequence}"
        },
        {
            'PRIMER_TASK': "${PRIMER_TASK}",
            'PRIMER_PICK_LEFT_PRIMER': ${PRIMER_PICK_LEFT_PRIMER},
            'PRIMER_PICK_RIGHT_PRIMER': ${PRIMER_PICK_RIGHT_PRIMER},
            'PRIMER_OPT_SIZE': ${PRIMER_OPT_SIZE},
            'PRIMER_MIN_SIZE': ${PRIMER_MIN_SIZE},
            'PRIMER_MAX_SIZE': ${PRIMER_MAX_SIZE},
            'PRIMER_OPT_TM': ${PRIMER_OPT_TM},
            'PRIMER_MIN_TM': ${PRIMER_MIN_TM},
            'PRIMER_MAX_TM': ${PRIMER_MAX_TM},
            'PRIMER_MIN_GC': ${PRIMER_MIN_GC},
            'PRIMER_MAX_GC': ${PRIMER_MAX_GC},
            'PRIMER_PAIR_MAX_DIFF_TM': ${PRIMER_PAIR_MAX_DIFF_TM},
            'PRIMER_PRODUCT_SIZE_RANGE': ${PRIMER_PRODUCT_SIZE_RANGE}
        }
    )

    # Extract primer results
    results = [
        {
            "Gene": "${meta.id}",
            "Primer_Set": i,
            "Forward_Primer": primers.get(f"PRIMER_LEFT_{i}_SEQUENCE", ""),
            "Reverse_Primer": primers.get(f"PRIMER_RIGHT_{i}_SEQUENCE", ""),
            "Forward_Primer_Adapter": f"{TF_FOW}{primers.get(f'PRIMER_LEFT_{i}_SEQUENCE', '')}",
            "Reverse_Primer_Adapter": f"{primers.get(f'PRIMER_RIGHT_{i}_SEQUENCE', '')}{TF_REV_COMP}",
            "Forward_Tm": primers.get(f"PRIMER_LEFT_{i}_TM", "N/A"),
            "Reverse_Tm": primers.get(f"PRIMER_RIGHT_{i}_TM", "N/A")
        }
        for i in range(int(primers.get("PRIMER_PAIR_NUM_RETURNED", 0)))
    ]

    # Save results to CSV
    pd.DataFrame(results).to_csv("${meta.id}_primers.csv", index=False)
    EOF
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        primer3: \$(primer3_core --version 2>&1 | head -n 1 | sed 's/^.*release //')
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
    
    stub:
    """
    echo "Gene,Primer_Set,Forward_Primer,Reverse_Primer,Forward_Tm,Reverse_Tm" > "${meta.id}_primers.csv"
    echo "${meta.id},0,AAAAAAAAAAAAAAAAAAAA,TTTTTTTTTTTTTTTTTTTT,60.0,60.0" >> "${meta.id}_primers.csv"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        primer3: \$(primer3_core --version 2>&1 | head -n 1 | sed 's/^.*release //')
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}