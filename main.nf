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
    .map { it -> [ it.getSimpleName().tokenize('_')[0], it.getSimpleName().tokenize('_')[1], it ] }   
    .set { ch_features }
// barcodes
Channel
    .fromPath(params.input_barcodes_file, checkIfExists:true)
    .map { it -> [ it.getSimpleName().tokenize('_')[0], it.getSimpleName().tokenize('_')[1], it ] } 
    .set { ch_barcodes }
// main input channel
ch_count
    .combine( ch_features, by:[0,1] )
    .combine( ch_barcodes, by:[0,1] )
    .set{ ch_input }  // dataset, experimental, full, absolute, count, features, barcodes


///////////////////////////////////////////////////////////////
/* --          IMPORT LOCAL MODULES/SUBWORKFLOWS          -- */
///////////////////////////////////////////////////////////////

// modules
// include { GET_REDUCED_DATASET } from "${launchDir}/modules/select-data.nf"
include { SELECT_NOZERO_GENES; FILTER_GENES } from "${launchDir}/modules/select-data.nf"
include { GET_RELATIVE } from "${launchDir}/modules/get-relative.nf"

// subworkflows
include { DATA_SIMULATION     } from "${launchDir}/subworkflows/1-data_simulation.nf"
include { DATA_EXPLORATION    } from "${launchDir}/subworkflows/2-data_exploration.nf"
include { DATA_PROCESSING     } from "${launchDir}/subworkflows/3-data_processing.nf"
include { CORRELATION         } from "${launchDir}/subworkflows/4-correlation.nf"


workflow {


    /* subworkflow: model and simulate data with scDesign2 */
    if ( params.do_simulation ){

        ch_input
            .filter{ it[1] in params.do_simulation_type }
            .set{ ch_input2simulation }

        // model and simulate data
        DATA_SIMULATION(
            ch_input2simulation,
            params.simulation_slope,
            params.simulation_ndata,
            params.simulation_cell_factor,
            params.simulation_depth_factor,
            params.do_simulation
        )  
        ch_input
            .filter{ it[1] == 'experimental' }   // don't want to compute the benchmark on phaseS original data
            .mix( DATA_SIMULATION.out.ch_simulated )
            .set{ ch_input }    
    } 

    /* module: get relative data */
    GET_RELATIVE(ch_input)
    ch_input
        .mix( GET_RELATIVE.out )
        .set{ ch_input }

    /* filter out genes with high zero rate */
    if ( params.do_filter_genes ){
        FILTER_GENES( ch_input )
        ch_input = FILTER_GENES.out
    }

    /* module: get non-zero gene datasets */
    // if ( params.do_nozero_genes ){
    //     GET_REDUCED_DATASET( ch_input )
    //     ch_input
    //         .mix( GET_REDUCED_DATASET.out )
    //         .set{ ch_input }
    // }
    if ( params.do_nozero_genes ){
        SELECT_NOZERO_GENES( ch_input )
        ch_input
            .mix( SELECT_NOZERO_GENES.out )
            .set{ ch_input }
    }

    /* subworkflow: plot basic figures for experimental and simulated data */
    // TODO also compute data exploration plots for transformed data
    if ( params.do_plot_data_exploration ){
        DATA_EXPLORATION( ch_input )  
    }


    /* subworkflow: data processing */
    ch_ori = ch_input
    if ( params.do_data_processing ){
        DATA_PROCESSING(
            ch_input,
            params.method_replace_zero.plus('NA'),
            params.method_transform_data,
            params.refgenes.plus('NA')
        )
    }

    /* subworkflow: compute correlation */
    if ( params.do_correlation ){
        CORRELATION(
            DATA_PROCESSING.out,
            ch_ori,
            params.method_correlation.tokenize(','),
            params.method_evaluation.tokenize(','),
            params.scatter_colorby
        )
    }

    // TODO plot scatter vs index, colored by reference gene
}

