process PRIMER_QUALITY {
    tag "$meta.id"
    label 'process_single'

    conda "bioconda::primer3-py=2.0.3 bioconda::biopython=1.85 conda-forge::pandas=2.2.3"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/primer3-py_pandas_pip_biopython:e1f5eb2f1ba237b8' :
        'community.wave.seqera.io/library/primer3-py_pandas_pip_biopython:e797e4d1477d8e75' }"
    
    input:
    tuple val(meta), val(seq), path(primer_file)
    val(hairpin_threshold)
    val(dimer_threshold)

    output:
    tuple val(meta), path("${meta.id}_quality.csv"), emit: primerquality

    script:
    def HAIRPIN_THRESHOLD = hairpin_threshold ?: -3.0
    def DIMER_THRESHOLD = dimer_threshold ?: -6.0
    """
    python3 <<EOF
    import pandas as pd
    import primer3

    # Define thresholds for problematic primers
    HAIRPIN_THRESHOLD = ${HAIRPIN_THRESHOLD}  # kcal/mol
    DIMER_THRESHOLD = ${DIMER_THRESHOLD}    # kcal/mol

    def analyze_primers(csv_path, output_path):
        # Read primer data
        df = pd.read_csv(csv_path)

        # Add columns for thermodynamic values and flags
        df['Forward_Hairpin_dG'] = 0.0
        df['Reverse_Hairpin_dG'] = 0.0
        df['Forward_Homodimer_dG'] = 0.0
        df['Reverse_Homodimer_dG'] = 0.0
        df['Heterodimer_dG'] = 0.0
        df['Hairpin_Flag'] = "No"
        df['Dimer_Flag'] = "No"

        # Calculate thermodynamic values for each primer pair
        for idx, row in df.iterrows():
            forward = row['Forward_Primer']
            reverse = row['Reverse_Primer']

            # Calculate hairpins
            fwd_hairpin = primer3.calcHairpin(forward)
            rev_hairpin = primer3.calcHairpin(reverse)

            # Calculate homodimers
            fwd_homodimer = primer3.calcHomodimer(forward)
            rev_homodimer = primer3.calcHomodimer(reverse)

            # Calculate heterodimers
            heterodimer = primer3.calcHeterodimer(forward, reverse)

            # Store ΔG values in the dataframe (convert from millical to kcal)
            df.at[idx, 'Forward_Hairpin_dG'] = fwd_hairpin.dg / 1000.0
            df.at[idx, 'Reverse_Hairpin_dG'] = rev_hairpin.dg / 1000.0
            df.at[idx, 'Forward_Homodimer_dG'] = fwd_homodimer.dg / 1000.0
            df.at[idx, 'Reverse_Homodimer_dG'] = rev_homodimer.dg / 1000.0
            df.at[idx, 'Heterodimer_dG'] = heterodimer.dg / 1000.0

            # Flag problematic primers
            if fwd_hairpin.dg / 1000.0 < HAIRPIN_THRESHOLD or rev_hairpin.dg / 1000.0 < HAIRPIN_THRESHOLD:
                df.at[idx, 'Hairpin_Flag'] = "Yes"
            if (fwd_homodimer.dg / 1000.0 < DIMER_THRESHOLD or
                rev_homodimer.dg / 1000.0 < DIMER_THRESHOLD or
                heterodimer.dg / 1000.0 < DIMER_THRESHOLD):
                df.at[idx, 'Dimer_Flag'] = "Yes"

        # Save results
        df.to_csv(output_path, index=False)

    # Run the analysis
    input_file = "${primer_file}"
    output_file = "${meta.id}_quality.csv"
    analyze_primers(input_file, output_file)
    print(f"Analysis complete. Results saved to {output_file}")
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
    echo "Gene,Primer_Set,Forward_Primer,Reverse_Primer,Forward_Tm,Reverse_Tm" > "${meta.id}_quality.csv"
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        primer3: \$(primer3_core --version 2>&1 | head -n 1 | sed 's/^.*release //')
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}