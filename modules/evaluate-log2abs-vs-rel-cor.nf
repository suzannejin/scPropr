
/* evaluate coefficients computed on absolute vs relative data */

process EVALUATE_LOG2ABS_VS_REL_COR {
    label 'process_high'
    container 'suzannejin/scpropr:propr'
    tag "${dataset}-${exp_sim}-${full}-${method_replace_zero}-${method_transform_data}-${refgene}-${method_cor}-${method_eval}"
    storeDir "${params.outdir}/${dataset}/evaluate/${method_eval}/${dataset}_${exp_sim}_${full}_${method_replace_zero}_${method_transform_data}_${refgene}_${method_cor}"

    input:
    tuple val(dataset), 
          val(exp_sim), 
          val(full), 
          val(method_replace_zero), 
          val(method_transform_data), 
          val(refgene),
          val(method_cor),
          file(abs), 
          file(rel),
          file(features)
    each method_eval

    output:
    file "${dataset}_${exp_sim}_${full}_${method_replace_zero}_${method_transform_data}_${refgene}_${method_cor}_${method_eval}.csv"
    file ".command.trace"
    file ".command.sh"

    script:
    """
    evaluate-log2abs-vs-rel-cor.R \
        --abs $abs \
        --rel $rel \
        --features $features \
        --out ${dataset}_${exp_sim}_${full}_${method_replace_zero}_${method_transform_data}_${refgene}_${method_cor}_${method_eval}.csv \
        --method $method_eval
    sleep 30
    """

    stub:
    """
    echo evaluate-log2abs-vs-rel-cor.R \
        --abs $abs \
        --rel $rel \
        --features $features \
        --out ${dataset}_${exp_sim}_${full}_${method_replace_zero}_${method_transform_data}_${refgene}_${method_cor}_${method_eval}.csv \
        --method $method_eval
    touch ${dataset}_${exp_sim}_${full}_${method_replace_zero}_${method_transform_data}_${refgene}_${method_cor}_${method_eval}.csv
    """
}