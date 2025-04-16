/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { paramsSummaryMap       } from 'plugin/nf-schema'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { PRIMER_DESIGN          } from '../modules/local/primerdesign/primerdesign'
include { PRIMER_QUALITY         } from '../modules/local/primerdesign/primerquality'
include { PRIMER_SPECIFICITY     } from '../modules/local/primerdesign/primerspecificity'
include { PRIMER_SELECTION       } from '../modules/local/primerdesign/primerselection'
include { PRIMER_MERGE           } from '../modules/local/primerdesign/primermerge'
include { BOWTIE2_BUILD          } from '../modules/nf-core/bowtie2/build'
include { FETCH_REFERENCE_GENOME } from '../modules/local/fetch_reference'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow GETPRIMERS {

    take:
    ch_samplesheet // channel: samplesheet read in from --input

    main:
    ch_versions = Channel.empty()

    PRIMER_DESIGN(ch_samplesheet,
        params.tf_fow,
        params.tf_rev,
        params.primer_opt_size,
        params.primer_min_size,
        params.primer_max_size,
        params.primer_min_tm,
        params.primer_max_tm,
        params.primer_pair_max_diff_tm,
        params.primer_product_size_range,
        params.primer_task,
        params.primer_opt_tm,
        params.primer_min_gc,
        params.primer_max_gc,
        params.primer_n_left,
        params.primer_n_right)
    ch_versions = ch_versions.mix(PRIMER_DESIGN.out.versions)

    // Primer Quality Module
    PRIMER_QUALITY(PRIMER_DESIGN.out.primers, params.hairpin_threshold, params.dimer_threshold)
    ch_versions = ch_versions.mix(PRIMER_QUALITY.out.versions)
    
    if (params.organism_type == 'human') {
        // Primer Specificity Module
        ch_specificity = PRIMER_DESIGN.out.primers.combine(params.genome)
        PRIMER_SPECIFICITY(ch_specificity)
        ch_versions = ch_versions.mix(PRIMER_SPECIFICITY.out.versions)
    }
    else if (params.organism_type == 'bacteria') {
        if (params.genome) {
                ch_specificity = PRIMER_DESIGN.out.primers.combine(params.genome)
                
        } else {
                ch_species = ch_samplesheet.map { it[0].species }.unique()
                FETCH_REFERENCE_GENOME(ch_species, params.taxdb)

                ch_primers_by_species = PRIMER_DESIGN.out.primers
                    .map { meta, seq, primers -> [meta.species, [meta, seq, primers]] }
                    .groupTuple()
                    .map { species, list -> [species, list.collect { it[0] }, list.collect { it[1] }, list.collect { it[2] }] }
                ch_specificity = FETCH_REFERENCE_GENOME.out.reference
                                    .combine(ch_primers_by_species, by: 0)
                                    .flatMap { species, reference_files, meta_list, sequence_list, primer_files ->
                                        meta_list.withIndex().collect { meta, index ->
                                            [meta, sequence_list[index], primer_files[index], reference_files]
                                        }
                                    }
        }
        PRIMER_SPECIFICITY(ch_specificity)
        ch_versions = ch_versions.mix(PRIMER_SPECIFICITY.out.versions)
    }  else {
        error "Unsupported organism type: ${params.organism_type}. Please specify 'human' or 'bacteria'."
    }

    ch_primer_output = PRIMER_DESIGN.out.primers
                            .combine(PRIMER_QUALITY.out.primerquality, by: 0)
                            .combine(PRIMER_SPECIFICITY.out.specificity_summary, by: 0)
    PRIMER_SELECTION( ch_primer_output )
    ch_versions = ch_versions.mix(PRIMER_SELECTION.out.versions)

    ch_allprimers = PRIMER_SELECTION.out.selected_primers.map{ it[1] }.collect()
    PRIMER_MERGE( ch_allprimers )
    ch_versions = ch_versions.mix(PRIMER_MERGE.out.versions)
    
    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'getprimers_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        )

    emit:
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
