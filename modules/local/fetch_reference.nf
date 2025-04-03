process FETCH_REFERENCE_GENOME {
    tag "${meta.species}"
    label 'process_low'
    storeDir "reference_genomes/"

    conda "bioconda::ncbi-genome-download=0.3.3 conda-forge::python=3.13.2 conda-forge::gzip=1.13"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/ncbi-genome-download_gzip_python:118c33aac0776166' :
        'community.wave.seqera.io/library/ncbi-genome-download_gzip_python:d9744426f7d491f4'}"

    input:
    tuple val(meta), val(gene_seq)
    path(taxdb)

    output:
    tuple val(meta), path("${meta.species}/*.fna")     , emit: reference
    path "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args           = task.ext.args ?: ''
    def species        = meta.species.replaceAll("_", " ")
    """
    # Set a custom cache directory
    export XDG_CACHE_HOME="cache"
    mkdir -p "\$XDG_CACHE_HOME"

    # Generate the taxid list using Python
    python3 << 'EOF'
    import sys

    # Parameters
    species_name = "${species}"
    names_dmp_path = "${taxdb}"
    output_file = "mytaxa.txt"

    try:
        # Open the names.dmp file and search for the species name
        with open(names_dmp_path, "r") as f:
            for line in f:
                fields = line.strip().split("|")
                taxid = fields[0].strip()
                name_txt = fields[1].strip()
                name_class = fields[3].strip()

                # Check if the name matches and is a scientific name
                if name_txt == species_name and name_class == "scientific name":
                    # Write the taxid to the output file
                    with open(output_file, "w") as out:
                        out.write(f"{taxid}\\n")
                    print(f"Taxid for '{species_name}': {taxid}")
                    break
            else:
                raise ValueError(f"Species '{species_name}' not found in names.dmp.")

    except Exception as e:
        sys.stderr.write(f"Error: {str(e)}\\n")
        sys.exit(1)
    EOF

    mytaxid=\$(head -n1 mytaxa.txt)
    # Use ncbi-genome-download with the generated taxid list
    ncbi-genome-download ${args} --formats fasta --taxids mytaxa.txt --output-folder ./ bacteria

    mkdir -p ${meta.species}
    mv refseq/bacteria/*/*.fna.gz ${meta.species}/
    gunzip ${meta.species}/*.fna.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ncbigenomedownload: \$( ncbi-genome-download --version )
    END_VERSIONS
    """
}