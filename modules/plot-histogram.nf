process PLOT_HISTOGRAM {
    label 'process_low_short'
    tag "$dataset-$exp_sim-$full-$abs_rel"
    container 'suzannejin/scpropr:plot'
    publishDir "${params.outdir}/${dataset}/plot/histogram/${dataset}_${exp_sim}_${full}_${abs_rel}", mode: params.publish_dir_mode

    input:
    tuple val(dataset),
          val(exp_sim),
          val(full),
          val(abs_rel),
          file(count),
          file(features),
          file(barcodes)

    output:
    path "histogram-libsize-percell.png"
    path "histogram-geomean-percell.png"
    path "histogram-mean-percell.png"
    path "histogram-var-percell.png"
    path "histogram-geomean-pergene.png"
    path "histogram-mean-pergene.png"
    path "histogram-var-pergene.png"
    path ".command.trace"
    path ".command.sh"
    path ".command.log"

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
        --geomean \
        --bygene
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

    stub:
    """
    touch histogram-libsize-percell.png
    touch histogram-geomean-percell.png
    touch histogram-mean-percell.png
    touch histogram-var-percell.png
    touch histogram-geomean-pergene.png
    touch histogram-mean-pergene.png
    touch histogram-var-pergene.png
    """
}