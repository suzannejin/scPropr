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
include { SELECT_NOZERO_GENES; EXP_NOZERO_GENES_2_SELECT } from "${launchDir}/modules/select-data.nf"
include { GET_RELATIVE; GET_RELATIVE_FOR_PLOT } from "${launchDir}/modules/get-relative.nf"

// subworkflows
include { DATA_SIMULATION     } from "${launchDir}/subworkflows/1-data_simulation.nf"
include { DATA_EXPLORATION    } from "${launchDir}/subworkflows/2-data_exploration.nf"
include { DATA_PROCESSING     } from "${launchDir}/subworkflows/3-data_processing.nf"
include { CORRELATION         } from "${launchDir}/subworkflows/4-correlation.nf"


workflow {


    /* subworkflow: model and simulate data with scDesign2 */
    if ( params.do_simulation || params.do_simulation_byslope || params.do_simulation_bydepth || params.do_simulation_bystep){

        ch_input
            .filter{ it[1] in params.do_simulation_type }
            .set{ ch_input2simulation }

        // model and simulate data
        DATA_SIMULATION(
            ch_input2simulation,
            params.simulation_slope,
            params.simulation_ndata,
            params.simulation_size_factor,
            params.simulation_cell_factor,
            params.do_simulation,
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
        // select non zero genes in full experimental data
        SELECT_NOZERO_GENES(
            ch_input.filter{ it[1] == 'experimental' }
        )
        ch_input
            .mix( SELECT_NOZERO_GENES.out )
            .set{ ch_input }
        // do the same for the other dataset following this non zero gene set
        ch_input
            .filter{ it[1] != 'experimental' }
            .combine(
                SELECT_NOZERO_GENES.out.map{ it -> [ it[0], it[5] ] }.unique(),
                by: 0
            )
            .set{ ch_input2select }
        ch_input2select.take(5).view()
        EXP_NOZERO_GENES_2_SELECT( ch_input2select )
        ch_input
            .mix( EXP_NOZERO_GENES_2_SELECT.out )
            .set{ ch_input }
    }

    /* module: get relative data */
    GET_RELATIVE(ch_input)
    ch_input
        .mix( GET_RELATIVE.out )
        .set{ ch_input }


    /* subworkflow: plot basic figures for experimental and simulated data */
    // TODO also compute data exploration plots for transformed data
    if (params.do_plot_data_exploration){
        DATA_EXPLORATION( ch_input )  
    }


    /* subworkflow: data processing */
    ch_ori = ch_input
    DATA_PROCESSING(
        ch_input,
        params.method_replace_zero.tokenize(','),
        params.method_transform_data.tokenize(','),
        params.refgenes.plus('NA')
    )

    /* subworkflow: compute correlation */
    CORRELATION(
        DATA_PROCESSING.out,
        ch_ori,
        params.method_correlation.tokenize(','),
        params.method_evaluation.tokenize(','),
        params.scatter_colorby
    )

    // TODO plot scatter vs index, colored by reference gene
}

