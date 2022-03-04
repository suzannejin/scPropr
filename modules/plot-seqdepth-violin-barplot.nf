
/* Plot rna per cell, gene per cell, and sample sequencing depth */

process PLOT_SEQDEPTH_VIOLIN_BARPLOT {
    tag "${type1}_${type2}"
    container 'suzannejin/scpropr:plot'
    publishDir "${params.outdir}/plot/seqdepth", mode: params.publish_dir_mode

    input:
    val datasets
    path counts
    val type1
    val type2

    output:
    path "seqdepth-violin-barplot-${type1}-${type2}.png"

    script:
    def dat = datasets.join(' ')
    """
    plot-seqdepth-violin-barplot.R \
        -i ${counts} \
        -n ${dat} \
        -o .

    mv seqdepth-violin-barplot.png seqdepth-violin-barplot-${type1}-${type2}.png
    """
}