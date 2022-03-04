#!/usr/bin/env nextflow

/* model and simulate single cell data using scDesign2 */


nextflow.enable.dsl = 2

include { MODEL_DATA } from "${launchDir}/modules/model-data.nf"
include { SIMULATE_DATA_BYDEPTH } from "${launchDir}/modules/simulate-data.nf"
include { SIMULATE_DATA_BYSLOPE } from "${launchDir}/modules/simulate-data.nf"
include { SIMULATE_DATA_BYSTEP  } from "${launchDir}/modules/simulate-data.nf"



workflow DATA_SIMULATION {

    take:
    ch_count
    ch_features
    ch_barcodes
    ch_barcodes_sim
    ch_size_factor
    ch_slope
    ch_ndata

    main:

    dataset = ch_count.map{ it -> [it[0]] }.flatten() 
    count   = ch_count.map{ it -> [it[3]] }.flatten() 
    
    /* model data */
    MODEL_DATA(dataset, count, ch_features, ch_barcodes_sim)

    /* simulate data */
    ch_dataset2simulate = MODEL_DATA.out.ch_dataset
    ch_model2simulate   = MODEL_DATA.out.ch_model
    SIMULATE_DATA_BYDEPTH(
        ch_dataset2simulate, 
        ch_model2simulate, 
        ch_size_factor
    )
    SIMULATE_DATA_BYSLOPE(
        ch_dataset2simulate,
        ch_model2simulate,
        ch_slope,
        ch_ndata
    )
    SIMULATE_DATA_BYSTEP(
        ch_dataset2simulate,
        ch_model2simulate,
        ch_slope,
        ch_ndata
    )

    SIMULATE_DATA_BYDEPTH.out.ch_dataset
        .mix( SIMULATE_DATA_BYSLOPE.out.ch_dataset )
        .mix( SIMULATE_DATA_BYSTEP.out.ch_dataset )
        .set{ ch_dataset }
    SIMULATE_DATA_BYDEPTH.out.ch_type1
        .mix( SIMULATE_DATA_BYSLOPE.out.ch_type1 )
        .mix( SIMULATE_DATA_BYSTEP.out.ch_type1 )
        .set{ ch_type1 }
    SIMULATE_DATA_BYDEPTH.out.ch_simulated
        .mix( SIMULATE_DATA_BYSLOPE.out.ch_simulated )
        .mix( SIMULATE_DATA_BYSTEP.out.ch_simulated )
        .set{ ch_simulated }

    emit:
    dataset    = ch_dataset
    type1      = ch_type1
    simulated  = ch_simulated
}