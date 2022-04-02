process PLOT_SCATTER_VS_INDEX {
    
    label 'process_low_short'
    tag "$dataset-$exp_sim-$full-$abs_rel"
    container 'suzannejin/scpropr:plot'
    publishDir "${params.outdir}/${dataset}/plot/scatter-vs-index/${dataset}_${exp_sim}_${full}_${abs_rel}", mode: params.publish_dir_mode

    input:
    tuple val(dataset),
          val(exp_sim),
          val(full),
          val(abs_rel),
          file(count),
          file(features),
          file(barcodes)

    output:
    path "scatter-libsize-cellindex*.png"
    path "scatter-geomean-cellindex*.png"
    path "scatter-mean-cellindex*.png"
    path "scatter-var-cellindex*.png"
    path "scatter-geomean-geneindex*.png"
    path "scatter-mean-geneindex*.png"
    path "scatter-var-geneindex*.png"
    path ".command.trace"
    path ".command.sh"

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
        --geomean \
        --bygene
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
    
    stub:
    """
    touch scatter-libsize-cellindex.png
    touch scatter-geomean-cellindex.png
    touch scatter-mean-cellindex.png
    touch scatter-var-cellindex.png
    touch scatter-geomean-geneindex.png
    touch scatter-mean-geneindex.png
    touch scatter-var-geneindex.png
    """
}