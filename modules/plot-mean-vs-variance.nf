process PLOT_MEAN_VS_VARIANCE {
    tag "${dataset}_${type1}_${type2}"
    container 'suzannejin/scpropr:plot'
    publishDir "${params.outdir}/${dataset}/${type1}/plot/expr/${type2}", mode: params.publish_dir_mode

    input:
    val dataset
    path count
    val type1
    val type2

    output:
    path "mean-vs-variance-bygene.png"
    path "mean-vs-variance-bycell.png"

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
}