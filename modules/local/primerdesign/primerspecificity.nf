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
    tuple val(meta), path("${meta.id}_specificity_results.csv")  , emit: specificity_results
    path "versions.yml"                                          , emit: versions

    script:
    """
    python3 <<EOF
    import os
    import pandas as pd
    from Bio.Blast.Applications import NcbiblastnCommandline

    # Read primer sequences from the input file
    primer_file = "${primer_file}"
    primers = pd.read_csv(primer_file)

    # Debug: Print the columns in the primer file
    print("Columns in the primer file:", primers.columns)

    # Check if required columns exist
    required_columns = ['Gene', 'Primer_Set', 'Forward_Primer', 'Reverse_Primer']
    for col in required_columns:
        if col not in primers.columns:
            raise ValueError(f"Missing required column '{col}' in the primer file.")

    # Write primers to a FASTA file for BLAST
    primer_fasta = "${meta.id}_primers.fasta"
    with open(primer_fasta, "w") as f:
        for index, row in primers.iterrows():
            f.write(f">{row['Gene']}_Set{row['Primer_Set']}_F\\n{row['Forward_Primer']}\\n")
            f.write(f">{row['Gene']}_Set{row['Primer_Set']}_R\\n{row['Reverse_Primer']}\\n")

    # Initialize a DataFrame to store combined results
    combined_results = pd.DataFrame()

    # Iterate over each reference genome and run BLAST
    reference_genomes = "${reference_genomes}".split()
    for ref_genome in reference_genomes:
        genome_name = os.path.basename(ref_genome).replace(".fna", "")
        blast_output = f"{genome_name}_blast_results.tsv"

        # Calculate the word size based on the primer length
        primer_length = primers["Forward_Primer"].str.len().max()  # Get the maximum primer length
        word_size = max(7, min(primer_length, 20))  # Ensure word_size is between 7 and 20

        # Construct the BLAST command
        blastn_cline = NcbiblastnCommandline(
            query=primer_fasta,
            subject=ref_genome,
            outfmt="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore",
            out=blast_output,
            task="blastn-short",
            word_size=word_size,
            evalue=10
        )

        # Debug: Print the BLAST command
        print(f"Running BLAST for {genome_name}: {blastn_cline}")

        # Execute the BLAST command
        stdout, stderr = blastn_cline()

        # Debug: Print BLAST output and errors
        print(stdout)
        print(stderr)

        # Parse BLAST results and add genome name
        columns = [
            "Query_ID", "Subject_ID", "Percent_Identity", "Alignment_Length",
            "Mismatches", "Gap_Openings", "Query_Start", "Query_End",
            "Subject_Start", "Subject_End", "E-value", "Bit_Score"
        ]
        blast_results = pd.read_csv(blast_output, sep="\\t", names=columns)
        blast_results["Genome"] = genome_name

        # Append results to the combined DataFrame
        combined_results = pd.concat([combined_results, blast_results], ignore_index=True)

    # Save combined results to CSV
    combined_results.to_csv("${meta.id}_specificity_results.csv", index=False)
    EOF

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        blast: \$(blastn -version | head -n 1 | awk '{print \$2}')
        python: \$(python --version | sed 's/Python //')
        biopython: \$(python -c "import Bio; print(Bio.__version__)")
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