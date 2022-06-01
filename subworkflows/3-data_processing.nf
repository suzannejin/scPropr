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
    PROCESS_DATA.out
        .map { it -> [ 
            it[0].getName().minus('.csv.gz').tokenize('_')[0],
            it[0].getName().minus('.csv.gz').tokenize('_')[1],
            it[0].getName().minus('.csv.gz').tokenize('_')[2],
            it[0].getName().minus('.csv.gz').tokenize('_')[3],
            it[0],
            it[1]
         ] }
        .combine( ch_input, by:[0,1,2,3] )
        .map{ it -> [ 
            it[4].getName().minus('.csv.gz').tokenize('_'),
            it[6],
            it[4..5],
            it[7..8]
         ].flatten() }
        .set{ ch_processed }   // dataset, exp_sim, full, abs_rel, replace_zero, transf_data, refgene, original_count, norm_count, norm_list, features, barcodes
    ch_processed
        .map{ it -> [it[0..6], it[10..11], it[7..8]].flatten() }
        .set{ ch_processed2 }  // dataset, exp_sim, full, abs_rel, replace_zero, transf_data, refgene, features, barcodes, original_count, norm_count

    /* scatter plot original count data vs processed count data */
    ch_processed2
        .groupTuple( by:[0,1,2,3,4,7,8,9] )
        .set{ ch2plot_original_vs_transformed }
    PLOT_ORIGINAL_VS_TRANSFORMED(ch2plot_original_vs_transformed)


    /* scatter plot log2 count data vs processed count data */
    ch_processed2
        .filter{ it[5] == 'log2' }
        .map{ it -> [ it[0..3], it[10] ].flatten() }
        .unique()
        .set{ ch_tmp }
    ch_processed2
        .filter{ it[5] != 'log2' }
        .map{ it -> [ it[0..8], it[10] ].flatten() }
        .unique()
        .combine( ch_tmp, by:[0,1,2,3] )
        .unique()
        .map{ it -> [ 
            it[0..8],
            it[10],
            it[9]
         ].flatten() }
        .groupTuple( by:[0,1,2,3,7,8,9] )
        .set{ ch2plot_log_vs_transformed }    // dataset, exp_sim, full, abs_rel, replace_zero, transf_data, refgene, features, barcodes, log2_count, norm_count
    PLOT_LOG_VS_TRANSFORMED(ch2plot_log_vs_transformed)


    /* plot transformed count on relative vs absolute data */
    ch_processed
        .filter{ it[3] == 'absolute' }
        .map{ it -> [ it[0..2], it[4..6], it[8] ].flatten() }
        .unique()
        .set{ ch_abs }
    ch_processed
        .filter{ it[3] == 'relative' }
        .map{ it -> [ it[0..2], it[4..6], it[8] ].flatten() }
        .unique()
        .set{ ch_rel }
    ch_abs
        .combine( ch_rel, by:[0,1,2,3,4,5] )  
        .groupTuple( by:[0,1,2,3] )
        .set{ ch2plot_abs_vs_rel }  
    PLOT_ABS_VS_REL_TRANSF(ch2plot_abs_vs_rel)


    /* plot absolute log2 count vs absolute transformed count */
    ch_processed
        .filter{ it[3] == 'absolute' && it[5] == 'log2' }
        .map{ it -> [ it[0..2], it[8], it[10] ].flatten() }  
        .unique()
        .set{ ch_log2abs }  
    ch_log2abs
        .combine( ch_rel, by:[0,1,2] )   // dataset, exp_sim, full, abslog_count, features, replace_zero, transf_data, refgene, normrel_count
        .map{ it -> [ it[0..2], it[5..7], it[3], it[8], it[4]].flatten() }
        .groupTuple( by:[0,1,2,6,8])
        .set{ ch2plot_log2abs_vs_rel }
    PLOT_LOG2ABS_VS_REL_TRANSF(ch2plot_log2abs_vs_rel)


    /* TODO check normalized factors per method */
    /* pairwise plot */
    /* distribution */
    /* pairwise correlation */

    /* TODO plot average count vs total rna per cell on transformed data */
    /* plot everything in one plot, colored by method */


    emit:
    ch_processed = ch_processed.map{ it -> [ it[0..6], it[8], it[10..11] ].flatten() } 
                    // dataset, exp_sim, full, abs_rel, replace_zero, transf_data, refgene, processed_count, features, barcodes
}