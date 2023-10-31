#!/usr/bin/env nextflow

/* Data processing
 *   - zero handling
 *   - normalization
 *   - log transformation
 */


nextflow.enable.dsl = 2

include { PROCESS_DATA } from "${launchDir}/modules/process-data.nf"
include { PLOT_ORIGINAL_VS_TRANSFORMED; PLOT_LOG_VS_TRANSFORMED } from "${launchDir}/modules/plot-original-vs-transformed.nf"
include { PLOT_ABS_VS_REL_TRANSF; PLOT_LOG2ABS_VS_REL_TRANSF    } from "${launchDir}/modules/plot-abs-vs-rel-transf.nf"


workflow DATA_PROCESSING {

    take:
    ch_input   // dataset, exp_sim, full, abs_rel, count, features, barcodes
    method_replace_zero
    method_transform_data
    refgenes

    main:

    /* process data: zero handling, normalization and log transformation */
    PROCESS_DATA(
        ch_input,
        method_replace_zero,
        method_transform_data,
        refgenes
    )
    PROCESS_DATA.out.count
        .map{ it -> [
            it[0].getName().minus('.csv.gz').tokenize('_'),
            it[1],
            it[2]
        ].flatten() }
        .set{ ch_processed }   // dataset, exp_sim, full, abs_rel, replace_zero, transf_data, refgene, transf_count, features, barcodes


    /* scatter plot original count data vs processed count data */
    ch_input
        .map{ it -> [it[0..4]].flatten() }
        .combine( ch_processed, by:[0,1,2,3])
        .map{ it -> [
            it[0],   // dataset
            it[1],   // exp_sim
            it[2],   // full
            it[3],   // abs_rel
            it[5],   // replace_zero
            it[6],   // transf_data
            it[7],   // refgene
            it[4],   // original count
            it[8],   // transformed count
            it[9],   // features
            it[10]   // barcodes
        ]}
        .groupTuple( by:[0,1,2,3,7,9,10] )
        .set{ ch2plot_original_vs_transformed }
    PLOT_ORIGINAL_VS_TRANSFORMED(ch2plot_original_vs_transformed)


    /* scatter plot log2 count data vs processed count data */
    ch_processed
        .filter{ it[5] == 'log2' }
        .set{ ch_log2 }
    ch_processed
        .filter{ it[5] != 'log2' }
        .set{ ch_transf }
    ch_log2
        .combine(ch_transf, by:[0,1,2,3,4])
        .map{ it -> [
            it[0],   // dataset
            it[1],   // exp_sim
            it[2],   // full
            it[3],   // abs_rel
            it[4],   // replace_zero
            it[10],  // transf_data
            it[11],  // refgene
            it[7],   // log2 count
            it[12],  // transformed count
            it[8],   // features
            it[9]    // barcodes
        ] }
        .groupTuple(by:[0,1,2,3,4,7,9,10])
        .set{ ch2plot_log_vs_transformed }
    PLOT_LOG_VS_TRANSFORMED(ch2plot_log_vs_transformed)


    /* plot transformed count on relative vs absolute data */
    ch_processed
        .filter{ it[3] == 'absolute' }
        .map{ it -> it.remove(3) }
        .set{ ch_abs }
    ch_processed
        .filter{ it[3] == 'relative' }
        .map{ it -> it.remove(3) }
        .set{ ch_rel }
    ch_abs
        .combine(ch_rel, by:[0,1,2,3,4,5])
        .map{ it -> [
            it[0],   // dataset
            it[1],   // exp_sim
            it[2],   // full
            it[3],   // replace_zero
            it[4],   // transf_data
            it[5],   // refgene
            it[6],   // abs count
            it[9]    // rel count
        ] }
        .groupTuple( by:[0,1,2,3] )
        .set{ ch2plot_abs_vs_rel }  
    PLOT_ABS_VS_REL_TRANSF(ch2plot_abs_vs_rel)


    /* plot absolute log2 count vs absolute transformed count */
    ch_processed
        .filter{ it[3] == 'absolute' && it[5] == 'log2' }
        .map{ it -> it.remove(3) }
        .set{ ch_log2abs }
    ch_log2abs
        .combine( ch_rel, by:[0,1,2,3] )
        .map{ it -> [
            it[0],   // dataset
            it[1],   // exp_sim
            it[2],   // full
            it[3],   // replace_zero
            it[9],   // transf_data
            it[10],  // refgene
            it[6],   // log2 abs count
            it[11],  // transformed rel count
            it[7]    // features
        ] }
        .groupTuple( by:[0,1,2,3,6,8])
        .set{ ch2plot_log2abs_vs_rel }
    PLOT_LOG2ABS_VS_REL_TRANSF(ch2plot_log2abs_vs_rel)


    /* TODO check normalized factors per method */
    /* pairwise plot */
    /* distribution */
    /* pairwise correlation */

    /* TODO plot average count vs total rna per cell on transformed data */
    /* plot everything in one plot, colored by method */


    emit:
    ch_processed 
}