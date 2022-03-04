process PLOT_HISTOGRAM {
    tag "${dataset}_${type1}_${type2}"
    container 'suzannejin/scpropr:plot'
    publishDir "${params.outdir}/${dataset}/${type1}/plot/histogram/${type2}", mode: params.publish_dir_mode

    input:
    val dataset
    path count
    val type1
    val type2

    output:
    path "histogram-libsize*.png"
    path "histogram-geomean*.png"
    path "histogram-cellmean*.png"
    path "histogram-cellvar*.png"
    path "histogram-genemean*.png"
    path "histogram-genevar*.png"

    script:
    """
    plot-histogram.R \
        $count \
        $dataset \
        . \
        --libsize
    plot-histogram.R \
        $count \
        $dataset \
        . \
        --geomean
    plot-histogram.R \
        $count \
        $dataset \
        . \
        --mean 
    plot-histogram.R \
        $count \
        $dataset \
        . \
        --var
    plot-histogram.R \
        $count \
        $dataset \
        . \
        --mean \
        --bygene
    plot-histogram.R \
        $count \
        $dataset \
        . \
        --var \
        --bygene
    """
}