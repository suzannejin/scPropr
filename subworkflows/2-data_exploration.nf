#!/usr/bin/env nextflow

/* plot initial figures to help with data exploration */


nextflow.enable.dsl = 2

include { PLOT_SEQDEPTH_VIOLIN_BARPLOT } from "${launchDir}/modules/plot-seqdepth-violin-barplot.nf"
include { PLOT_SCATTER_VS_INDEX        } from "${launchDir}/modules/plot-scatter-vs-index.nf"
include { PLOT_HISTOGRAM               } from "${launchDir}/modules/plot-histogram.nf"
include { PLOT_DROPOUT                 } from "${launchDir}/modules/plot-dropout.nf"
include { PLOT_MEAN_VS_VARIANCE        } from "${launchDir}/modules/plot-mean-vs-variance.nf"


workflow DATA_EXPLORATION {

    take:
    ch_input  // dataset, exp_sim, full, abs_rel, count, features, barcodes

    main:

    // TODO run data exploration for phase S experimental dataset as well

    /* plot 
     * 1. reads per cell
     * 2. genes per cell
     * 3. seqdepth barplot
     */ 
    ch_input
        .map{ it -> it[0..4].flatten() }
        .groupTuple(by:[0,2,3])
        .set{ ch_input2plot_seqdepth }
    PLOT_SEQDEPTH_VIOLIN_BARPLOT(ch_input2plot_seqdepth)

    /* plot scatter vs cell index
     *      a. y = library size
     *      b. y = geometric mean
     * NOTE for relative data, mean vs index has always the same value, but the var vs index looks very similar to absolute data, since by dividing the absolute data by the library size we only standardize by the libsize, and we don't account for the variance. This means that our relative data still keeps the original variance structure. Do normalisation methods account for this and still be able to normalize the data based on the data variability (eg. guess unchanged genes considering variance-like measures)?
     * TODO if normalisations do account for the data variability and get the correct answer, then for our benchmark we need a better way to define relative data
     */
    PLOT_SCATTER_VS_INDEX(ch_input)

    /* plot histogram of
     *      a. library size
     *      b. geometric mean
     */
    PLOT_HISTOGRAM(ch_input)

    /* plot dropout */
    PLOT_DROPOUT(ch_input)

    /* TODO plot dropout vs index, colored by reference gene */

    /* plot mean expression vs variance */
    PLOT_MEAN_VS_VARIANCE(ch_input)

    /* TODO plot mean expression vs seqdepth */

    /* TODO plot umi count vs total rna per cell 
     * black line = mean
     * colored region = hinges or whiskers
     */
    

    /* TODO plot log umi count vs total rna per cell 
     * black line = mean
     * colored region = hinges or whiskers
     */
    

    /* correlation between umi counts and total rna per cell */

}