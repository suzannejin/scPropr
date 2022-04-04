/* select genes with dropout < threshold to compute procrustes analysis */

process SELECT_GENE_PROCRUSTES{
    label 'process_low_short'
    container 'suzannejin/scpropr:propr'
    tag "$dataset-$exp_sim-$full-$abs_rel-$threshold"
    
    input:
    tuple val(dataset),
          val(exp_sim),
          val(full),
          val(abs_rel),
          path(count),
          path(features),
          path(barcodes)
    val threshold

    output:
    tuple val(dataset),
          val(exp_sim),
          val(full),
          val(abs_rel),
          file("*.gene")

    script:
    """
    select-data.R \
        -i $count \
        -f $features \
        -b $barcodes \
        -o ${dataset}_${exp_sim}_procrustes_${abs_rel}.csv.gz \
        -o2 ${dataset}_${exp_sim}_procrustes_features.csv \
        -o3 ${dataset}_${exp_sim}_procrustes_barcodes.csv \
        --filter_gene $threshold

    awk '{gene=\$0; print gene > gene".gene"}' ${dataset}_${exp_sim}_procrustes_features.csv
    """
}