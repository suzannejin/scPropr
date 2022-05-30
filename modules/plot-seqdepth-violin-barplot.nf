
/* Plot rna per cell, gene per cell, and sample sequencing depth */

process PLOT_SEQDEPTH_VIOLIN_BARPLOT {
    label 'process_low_short'
    tag "$dataset-$full-$abs_rel"
    container 'suzannejin/scpropr:plot'
    publishDir "${params.outdir}/${dataset}/plot/seqdepth/${dataset}_${full}_${abs_rel}", mode: params.publish_dir_mode

    input:
    tuple val(dataset),
          val(exp_sim),
          val(full),
          val(abs_rel),
          file(count)
    
    output:
    file "seqdepth-violin-barplot.png"
    file ".command.trace"
    file ".command.sh"
    file ".command.log"

    script:
    def exp_sims = exp_sim.join(' ')
    def counts   = count.join(' ')
    """
    plot-seqdepth-violin-barplot.R \
        -i $counts \
        -n $exp_sims \
        -o seqdepth-violin-barplot.png
    """

    stub:
    def exp_sims = exp_sim.join(' ')
    def counts   = count.join(' ')
    """
    echo plot-seqdepth-violin-barplot.R \
        -i $counts \
        -n $exp_sims \
        -o seqdepth-violin-barplot.png
    touch seqdepth-violin-barplot.png
    """
}