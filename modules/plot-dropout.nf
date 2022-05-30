process PLOT_DROPOUT {
    label 'process_low_short'
    tag "$dataset-$exp_sim-$full-$abs_rel"
    container 'suzannejin/scpropr:plot'
    publishDir "${params.outdir}/${dataset}/plot/dropout/${dataset}_${exp_sim}_${full}_${abs_rel}", mode: params.publish_dir_mode

    input:
    tuple val(dataset),
          val(exp_sim),
          val(full),
          val(abs_rel),
          file(count),
          file(features),
          file(barcodes)

    output:
    file "dropout"
    file "dropout.png"
    file "nozero_cells"
    file "nozero_genes"
    file ".command.trace"
    file ".command.sh"
    file ".command.log"

    script:
    """
    plot-dropout.R \
        $count \
        $dataset \
        .
    """

    stub:
    """
    echo plot-dropout.R \
        $count \
        $dataset \
        .
    touch dropout dropout.png nozero_cells nozero_genes
    """
}