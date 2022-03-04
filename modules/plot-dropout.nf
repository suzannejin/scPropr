process PLOT_DROPOUT {
    tag "${dataset}_${type1}_${type2}"
    container 'suzannejin/scpropr:plot'
    publishDir "${params.outdir}/${dataset}/${type1}/plot/dropout/${type2}", mode: params.publish_dir_mode

    input:
    val dataset
    path count
    val type1
    val type2

    output:
    path "dropout"
    path "dropout.png"
    path "nozero_cells"
    path "nozero_genes"

    script:
    """
    plot-dropout.R \
        $count \
        $dataset \
        .
    """
}