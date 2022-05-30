process PLOT_MEAN_VS_VARIANCE {
    label 'process_low_short'
    tag "$dataset-$exp_sim-$full-$abs_rel"
    container 'suzannejin/scpropr:plot'
    publishDir "${params.outdir}/${dataset}/plot/mean-vs-variance/${dataset}_${exp_sim}_${full}_${abs_rel}", mode: params.publish_dir_mode

    input:
    tuple val(dataset),
          val(exp_sim),
          val(full),
          val(abs_rel),
          file(count),
          file(features),
          file(barcodes)

    output:
    path "mean-vs-variance-bygene.png"
    path "mean-vs-variance-bycell.png"
    path ".command.trace"
    path ".command.sh"
    path ".command.log"

    script:
    """
    plot-mean-vs-variance.R \
        $count \
        $dataset \
        . \
        --bygene
    plot-mean-vs-variance.R \
        $count \
        $dataset \
        . \
        --bycell
    """

    stub:
    """
    touch mean-vs-variance-bygene.png mean-vs-variance-bycell.png
    """
}