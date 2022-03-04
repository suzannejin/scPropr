#!/usr/bin/env nextflow

/* plot initial figures to help with data exploration */


nextflow.enable.dsl = 2

include { PLOT_SEQDEPTH_VIOLIN_BARPLOT } from "${launchDir}/modules/plot-seqdepth-violin-barplot.nf"
include { PLOT_SCATTER_VS_INDEX } from "${launchDir}/modules/plot-scatter-vs-index.nf"
include { PLOT_HISTOGRAM } from "${launchDir}/modules/plot-histogram.nf"
include { PLOT_DROPOUT } from "${launchDir}/modules/plot-dropout.nf"
include { PLOT_MEAN_VS_VARIANCE } from "${launchDir}/modules/plot-mean-vs-variance.nf"


workflow DATA_EXPLORATION {

    take:
    ch_count // cell type, data type 1, data type 2, path to count file

    main:

    datasets = ch_count.map{ it -> [it[0]] }.collect()
    counts   = ch_count.map{ it -> [it[3]] }.collect()

    dataset = ch_count.map{ it -> [it[0]] }.flatten()  // cell type dataset name
    type1   = ch_count.map{ it -> [it[1]] }.flatten()  // experimental, simulated, merged, etc
    type2   = ch_count.map{ it -> [it[2]] }.flatten()  // absolute or relative
    count   = ch_count.map{ it -> [it[3]] }.flatten()  // path to count file


    /* plot 
     * 1. reads per cell
     * 2. genes per cell
     * 3. seqdepth barplot
     */ 
    PLOT_SEQDEPTH_VIOLIN_BARPLOT(datasets, counts, type1, type2)

    /* plot scatter vs cell index
     *      a. y = library size
     *      b. y = geometric mean
     */
    PLOT_SCATTER_VS_INDEX(dataset, count, type1, type2)

    /* plot histogram of
     *      a. library size
     *      b. geometric mean
     */
    PLOT_HISTOGRAM(dataset, count, type1, type2)

    /* plot dropout */
    PLOT_DROPOUT(dataset, count, type1, type2)

    /* plot mean expression vs variance */
    PLOT_MEAN_VS_VARIANCE(dataset, count, type1, type2)


    /* plot umi count vs total rna per cell 
     * black line = mean
     * colored region = hinges or whiskers
     */
    

    /* plot log umi count vs total rna per cell 
     * black line = mean
     * colored region = hinges or whiskers
     */
    

    /* correlation between umi counts and total rna per cell */

}