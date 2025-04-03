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
        [],[],[],[],[],
        [],[],[],[],[],
        [],[],[],[],[])

    // Primer Quality Module
    PRIMER_QUALITY(PRIMER_DESIGN.out.primers, [], [])
    
    if (params.organism_type == 'human') {
        // Primer Specificity Module
        ch_specificity = PRIMER_DESIGN.out.primers.combine(params.genome)
        PRIMER_SPECIFICITY(ch_specificity)
    }
    else if (params.organism_type == 'bacteria') {
        if (params.genome) {
                ch_specificity = PRIMER_DESIGN.out.primers.combine(params.genome)
                
        } else {
                FETCH_REFERENCE_GENOME(ch_samplesheet, params.taxdb)
                ch_specificity = PRIMER_DESIGN.out.primers
                                    .combine(FETCH_REFERENCE_GENOME.out.reference, by: 0) 
        }
        PRIMER_SPECIFICITY(ch_specificity)
    }  else {
        error "Unsupported organism type: ${params.organism_type}. Please specify 'human' or 'bacteria'."
    }
    
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
