#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


////////////////////////////////////////////////////
/* --          VALIDATE INPUTS                 -- */
////////////////////////////////////////////////////

// input count data
Channel
    .fromPath(params.input_count_file, checkIfExists:true)
    .map { it -> [ it.getSimpleName().tokenize('_')[0], it.getSimpleName().tokenize('_')[1], it.getSimpleName().tokenize('_')[2], it ] }  // cell type, experimental, absolute, path to count file
    .set { ch_count }
// features
Channel
    .fromPath(params.input_features_file, checkIfExists:true)
    .map { it -> [ it.getParent().baseName, it ] }   // cell type, path to features file
    .set { ch_features }
// barcodes
Channel
    .fromPath(params.input_barcodes_file, checkIfExists:true)
    .map { it -> [ it.getParent().baseName, it ] }  // cell type, path to barcodes file
    .set { ch_barcodes }
// cells from a homogeneous population (eg. same cell cycle phase)
Channel
    .fromPath(params.select_cells_pos_simulation)
    .map { it -> [ it.getParent().baseName, it ] }   // cell type, path to barcode positions file
    .set { ch_barcodes_sim }
// simulation parameters
Channel
    .fromList(params.simulation_dataset)
    .set { ch_simulation_dataset }
Channel
    .fromList(params.simulation_size_factors)
    .set { ch_size_factor }
Channel
    .fromList(params.simulation_slope)
    .combine(Channel.fromList(params.simulate_ndata))
    .multiMap{ it ->
        slope: it[0]
        ndata: it[1]
    }
    .set{ch_slope_ndata}


///////////////////////////////////////////////////////////////
/* --          IMPORT LOCAL MODULES/SUBWORKFLOWS          -- */
///////////////////////////////////////////////////////////////

include { DATA_EXPLORATION } from "${launchDir}/subworkflows/1-data_exploration.nf"
include { DATA_SIMULATION }  from "${launchDir}/subworkflows/2-data_simulation.nf"


workflow{

    // subworkflow: model and simulate data with scDesign2
    ch_count
        .combine(ch_simulation_dataset, by:0)
        .set{ ch_count2simulation }
    ch_features
        .combine(ch_simulation_dataset, by:0)
        .set{ ch_features2simulation }
    ch_barcodes_sim
        .combine(ch_simulation_dataset, by:0)
        .set{ ch_barcodessim2simulation}

    DATA_SIMULATION(
        ch_count2simulation,
        ch_features2simulation,
        ch_barcodessim2simulation,
        ch_size_factor,
        ch_slope_ndata.slope,
        ch_slope_ndata.ndata
    )
    DATA_SIMULATION.out.dataset
        .join( DATA_SIMULATION.out.type1 )
        .combine( Channel.from('absolute') )
        .join( DATA_SIMULATION.out.simulated )
        .set{ ch_simulated }

    // subworkflow: plot basic figures for experimental and simulated data 
    ch_count
        .mix( ch_simulated )
        .set{ ch_count_simulated }
    DATA_EXPLORATION( ch_count_simulated )
}

