#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


////////////////////////////////////////////////////
/* --          VALIDATE INPUTS                 -- */
////////////////////////////////////////////////////

// input count data
Channel
    .fromPath(params.input_count_file, checkIfExists:true)
    .map { it -> [ it.getSimpleName().tokenize('_')[0], it.getSimpleName().tokenize('_')[1], it.getSimpleName().tokenize('_')[2], it.getSimpleName().tokenize('_')[3], it ] }  
    .set { ch_count }  
// features
Channel
    .fromPath(params.input_features_file, checkIfExists:true)
    .map { it -> [ it.getParent().baseName, it ] }   
    .set { ch_features }
// barcodes
Channel
    .fromPath(params.input_barcodes_file, checkIfExists:true)
    .map { it -> [ it.getParent().baseName, it ] } 
    .set { ch_barcodes }
// main input channel
ch_count
    .combine( ch_features, by:0 )
    .combine( ch_barcodes, by:0 )
    .set{ ch_input }  // dataset, experimental, full, absolute, count, features, barcodes


///////////////////////////////////////////////////////////////
/* --          IMPORT LOCAL MODULES/SUBWORKFLOWS          -- */
///////////////////////////////////////////////////////////////

// modules
include { SELECT_NOZERO_GENES; SELECT_BARCODES_SIM } from "${launchDir}/modules/select-data.nf"
include { GET_RELATIVE; GET_RELATIVE_FOR_PLOT      } from "${launchDir}/modules/get-relative.nf"

// subworkflows
include { DATA_SIMULATION     } from "${launchDir}/subworkflows/1-data_simulation.nf"
include { DATA_EXPLORATION    } from "${launchDir}/subworkflows/2-data_exploration.nf"
include { DATA_PROCESSING     } from "${launchDir}/subworkflows/3-data_processing.nf"
include { CORRELATION         } from "${launchDir}/subworkflows/4-correlation.nf"


workflow {


    /* subworkflow: model and simulate data with scDesign2 */
    if (params.do_simulation_byslope || params.do_simulation_bydepth || params.do_simulation_bystep){

        // select the data for simulation
        if (params.select_cells_simulation){
            Channel
                .fromPath( params.select_cells_simulation, checkIfExists:true )
                .map { it -> [ it.getParent().baseName, it ] }
                .set{ ch_barcodes_sim }     
            ch_input
                .combine( ch_barcodes_sim, by:0 )
                .set{ ch_input2select_simulation }   // dataset, experimental, full or nozero, absolute, count, features, barcodes, barcodes_sim
            ch_input2simulation = SELECT_BARCODES_SIM(ch_input2select_simulation)
            ch_input2plot = ch_input2simulation   // homogeneous population
        }else{
            ch_input2simulation = ch_input
            ch_input2plot = Channel.empty()
        }

        // model and simulate data
        DATA_SIMULATION(
            ch_input2simulation,
            params.simulation_slope,
            params.simulation_ndata,
            params.simulation_size_factor,
            params.simulation_cell_factor,
            params.do_simulation_byslope, 
            params.do_simulation_bydepth,
            params.do_simulation_bystep
        )  
        ch_input
            .mix( DATA_SIMULATION.out.ch_simulated )
            .set{ ch_input }    
    } 


    /* module: get non-zero gene datasets */
    if (params.do_nozero_genes){
        SELECT_NOZERO_GENES(ch_input)
        ch_input
            .mix( SELECT_NOZERO_GENES.out )
            .set{ ch_input }
    }
    if (params.do_test){
        ch_input 
            .filter{ it[2] == 'nozero' }
            .set{ ch_input }
    }


    /* subworkflow: plot basic figures for experimental and simulated data */
    // TODO also compute data exploration plots for transformed data
    ch_input2plot
        .mix( ch_input )
        .set{ ch_input2plot }
    GET_RELATIVE_FOR_PLOT(ch_input2plot)
    ch_input2plot
        .mix( GET_RELATIVE_FOR_PLOT.out )
        .set{ ch_input2plot }
    DATA_EXPLORATION( ch_input2plot )  


    /* module: get relative data */
    GET_RELATIVE(ch_input)
    ch_input
        .mix( GET_RELATIVE.out )
        .set{ ch_input }


    /* subworkflow: data processing */
    DATA_PROCESSING(
        ch_input,
        params.method_replace_zero.tokenize(','),
        params.method_transform_data.tokenize(','),
        params.refgenes.plus('NA')
    )

    /* subworkflow: compute correlation */
    CORRELATION(
        DATA_PROCESSING.out,
        params.method_correlation.tokenize(','),
        params.method_evaluation.tokenize(',')
    )

    // TODO plot scatter vs index, colored by reference gene
}

