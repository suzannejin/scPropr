process PLOT_SCATTER_VS_INDEX {
    tag "${dataset}_${type1}_${type2}"
    container 'suzannejin/scpropr:plot'
    publishDir "${params.outdir}/${dataset}/${type1}/plot/scatter/${type2}", mode: params.publish_dir_mode

    input:
    val dataset
    path count
    val type1
    val type2

    output:
    path "scatter-libsize-index*.png"
    path "scatter-geomean-index*.png"
    path "scatter-cellmean-index*.png"
    path "scatter-cellvar-index*.png"
    path "scatter-genemean-index*.png"
    path "scatter-genevar-index*.png"

    script:
    """
    plot-scatter-vs-index.R \
        $count \
        $dataset \
        . \
        --libsize
    plot-scatter-vs-index.R \
        $count \
        $dataset \
        . \
        --geomean
    plot-scatter-vs-index.R \
        $count \
        $dataset \
        . \
        --mean
    plot-scatter-vs-index.R \
        $count \
        $dataset \
        . \
        --var
    plot-scatter-vs-index.R \
        $count \
        $dataset \
        . \
        --mean \
        --bygene
    plot-scatter-vs-index.R \
        $count \
        $dataset \
        . \
        --var \
        --bygene
    """
}