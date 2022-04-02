#!/usr/bin/env nextflow

/* model and simulate single cell data using scDesign2 */


nextflow.enable.dsl = 2

include { MODEL_DATA } from "${launchDir}/modules/model-data.nf"
include { SIMULATE_DATA_BYDEPTH; SIMULATE_DATA_BYSLOPE; SIMULATE_DATA_BYSTEP } from "${launchDir}/modules/simulate-data.nf"


workflow DATA_SIMULATION {

    take:
    ch_input  // dataset, experimental, full, absolute, count, features, barcodes
    slope
    ndata
    size_factor
    cell_factor
    do_simulate_byslope
    do_simulate_bydepth
    do_simulate_bystep

    main:
    
    // TODO check why they compute simulations with cell factor 0
    /* model data */
    MODEL_DATA(ch_input)
    MODEL_DATA.out
        .map{ it -> [ 
            it[1].getName().minus('.csv').tokenize('_')[0], 
            it[1].getName().minus('.csv').tokenize('_')[1], 
            it[1].getName().minus('.csv').tokenize('_')[2], 
            it[1].getName().minus('.csv').tokenize('_')[3], 
            it[0], 
            it[1], 
            it[2] 
        ] }
        .set{ ch_model2simulation }
    
    /* simulate data */
    ch_simulated = Channel.empty()
    if (do_simulate_byslope){
        SIMULATE_DATA_BYSLOPE(ch_model2simulation, slope, ndata, cell_factor)
        ch_simulated
            .mix( SIMULATE_DATA_BYSLOPE.out )
            .set{ ch_simulated }
    }
    if (do_simulate_bydepth){ 
        SIMULATE_DATA_BYDEPTH(ch_model2simulation, size_factor)
        ch_simulated
            .mix( SIMULATE_DATA_BYDEPTH.out )
            .set{ ch_simulated }
    }
    if (do_simulate_bystep){
        SIMULATE_DATA_BYSTEP(ch_model2simulation, slope, ndata, cell_factor)
        ch_simulated
            .mix( SIMULATE_DATA_BYSTEP.out )
            .set{ ch_simulated }
    }

    // organize output channel
    ch_simulated
        .map{ it -> [ 
            it[0].getName().minus('.csv.gz').tokenize('_')[0], 
            it[0].getName().minus('.csv.gz').tokenize('_')[1], 
            it[0].getName().minus('.csv.gz').tokenize('_')[2], 
            it[0].getName().minus('.csv.gz').tokenize('_')[3], 
            it[0], 
            it[1], 
            it[2] 
        ] }
        .set{ ch_simulated }
        

    // TODO Compare simulated vs original
    // be careful here, because they might have different dimensions

    emit:
    ch_simulated = ch_simulated   // dataset, simulated tag, full, absolute, simulated count, features, barcodes
}