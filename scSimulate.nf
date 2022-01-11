#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


////////////////////////////////////////////////////
/* --          VALIDATE INPUTS                 -- */
////////////////////////////////////////////////////

Channel
    .fromPath(params.count_file, checkIfExists:true)
    .map { it -> [ it.getParent().getParent().baseName, it.baseName.tokenize('.')[0], it ] }  
    .set { ch_count }
Channel
    .fromPath(params.features, checkIfExists:true)
    .map { it -> [ it.getParent().getParent().baseName, it ] }
    .set { ch_features }
Channel
    .fromPath(params.nozero_genes, checkIfExists:true)
    .map { it -> [ it.getParent().getParent().baseName, it ] }
    .set { ch_nozero }
ch_count
    .filter{ (it[1] == "absolute") }
    .combine( ch_features, by:0 )
    .set { ch_input2model }  // cell_type, count_type, count_file, features_file


////////////////////////////////////////////////////
/* --          IMPORT LOCAL MODULES            -- */
////////////////////////////////////////////////////

include { MODEL    } from "${baseDir}/modules/simulate.nf"
include { SIMULATE } from "${baseDir}/modules/simulate.nf"
include { CHECK_DROPOUT } from "${baseDir}/modules/simulate.nf"


////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////

workflow {

    /* fix model */
    MODEL(ch_input2model)
    MODEL.out.ch_model
        .combine( Channel.fromList(params.size_factors) )
        .set { ch_model2simulate }

    /* simulate data */
    SIMULATE(ch_model2simulate)

    /* check the simulated data dropout */
    SIMULATE.out.ch_simulated
        .combine( ch_nozero, by:0 )
        .set { ch_simulated2dropout }
    CHECK_DROPOUT(ch_simulated2dropout)

}