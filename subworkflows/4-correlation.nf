#!/usr/bin/env nextflow

/* Compute and evaluate correlation coefficients */

nextflow.enable.dsl = 2

include { GET_CORRELATION } from "${launchDir}/modules/get-correlation.nf"
include { PLOT_ABS_VS_REL_COR; PLOT_LOG2ABS_VS_REL_COR_COLORED; PLOT_LOG2ABS_VS_REL_COR_COLORED_INDIV; PLOT_LOG2ABS_VS_REL_COR_COLORED_INDIV_GM } from "${launchDir}/modules/plot-abs-vs-rel-cor.nf"
include { EVALUATE_LOG2ABS_VS_REL_COR; EVALUATE_LOG2ABS_VS_REL_COR_FILTER } from "${launchDir}/modules/evaluate-log2abs-vs-rel-cor.nf"


workflow CORRELATION {

    take:
    ch_input  // dataset, exp_sim, full, abs_rel, replace_zero, transf_data, refgene, processed_count, features, barcodes
    ch_ori    // dataset, experimental, full, absolute, count, features, barcodes
    method_cor
    method_eval
    scatter_colorby

    main:

    /* compute correlation coefficients */
    GET_CORRELATION(
        ch_input,
        method_cor
    )
    GET_CORRELATION.out
        .map{ it -> [ 
            it[0].getName().minus('.csv.gz').tokenize('_'),
            it[0],
            it[1]
        ].flatten() }
        .set{ ch_cor }  // dataset, exp_sim, full, abs_rel, method_replace_zero, method_transf_data, refgene, method_cor, cor matrix, features

    /* plot transformed count on relative vs absolute data */
    ch_cor
        .filter{ it[3] == 'absolute' }
        .map{ it -> [ it[0..2], it[4..8] ].flatten() }
        .unique()
        .set{ ch_abs }
    ch_cor
        .filter{ it[3] == 'relative' }
        .map{ it -> [ it[0..2], it[4..8] ].flatten() }
        .unique()
        .set{ ch_rel }  
    ch_abs
        .combine( ch_rel, by:[0,1,2,3,4,5,6] )
        .groupTuple( by:[0,1,2,3] )
        .set{ ch2plot_abs_vs_rel }   
    PLOT_ABS_VS_REL_COR(ch2plot_abs_vs_rel)


    /* plot absolute log2 count vs absolute transformed count */
    ch_ori
        .filter{ it[3] == 'absolute' }
        .map{ it -> [
            it[0],
            it[1],
            it[2],
            it[4]
        ] }
        .set{ ch_ori }
    ch_cor
        .filter { it[3] == 'absolute' && it[5] == 'log2' }
        .map{ it -> [ it[0..2], it[4], it[8..9], it[7] ].flatten() } 
        .unique()
        .set{ ch_log2abs }  // dataset, exp_sim, full, method_replace_zero, cor matrix, features, method_cor
    ch_log2abs
        .combine( ch_rel, by:[0,1,2,3,6] )  
        .map{ it -> [ it[0..3], it[7..8], it[4], it[5], it[9], it[6] ].flatten() }
        .groupTuple( by:[0,1,2,3] )
        .map{ it -> [
            it[0], 
            it[1],
            it[2],
            it[3],
            it[4],
            it[5],
            it[6],
            it[7].unique(),
            it[8],
            it[9][0]   
        ] }
        .set{ ch2plot_log2abs_vs_rel }  // dataset, exp_sim, full, method_replace_zero, method_transf_data, refgene, method_cor, log2abs matrix, rel matrix, features
    ch2plot_log2abs_vs_rel
        .combine( ch_ori, by:[0,1,2] )
        .set{ ch2plot_log2abs_vs_rel_colored }
    PLOT_LOG2ABS_VS_REL_COR_COLORED(
        ch2plot_log2abs_vs_rel_colored,
        scatter_colorby
    )


    /* TODO compute metric */
    ch_log2abs
        .combine( ch_rel, by:[0,1,2,3,6] )
        .map{ it -> [ it[0..3], it[7..8], it[4], it[5], it[9], it[6] ].flatten() }
        .combine( ch_ori, by:[0,1,2])
        .set{ ch2evaluate_log2abs_vs_rel }   // dataset, exp_sim, full, method_replace_zero, method_transf_zero, refgene, method_cor, log2abs matrix, rel matrix, features
    EVALUATE_LOG2ABS_VS_REL_COR(
        ch2evaluate_log2abs_vs_rel, 
        method_eval
    )
    EVALUATE_LOG2ABS_VS_REL_COR_FILTER(
        ch2evaluate_log2abs_vs_rel,
        method_eval
    )


    /* individual plots */
    ch2evaluate_log2abs_vs_rel
        .combine( ch_ori, by:[0,1,2] )
        .set{ ch2plot_log2abs_vs_rel_colored_indiv }
    PLOT_LOG2ABS_VS_REL_COR_COLORED_INDIV(ch2plot_log2abs_vs_rel_colored_indiv)
    PLOT_LOG2ABS_VS_REL_COR_COLORED_INDIV_GM(ch2plot_log2abs_vs_rel_colored_indiv)


    emit:
    ch_cor = ch_cor  // dataset, exp_sim, full, abs_rel, method_replace_zero, method_transf_data, refgene, method_cor, cor matrix, features
}