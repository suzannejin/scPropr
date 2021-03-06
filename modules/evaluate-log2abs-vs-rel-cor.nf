
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
          file(features),
          file(ori)
    each method_eval

    output:
    file "${dataset}_${exp_sim}_${full}_${method_replace_zero}_${method_transform_data}_${refgene}_${method_cor}_${method_eval}.csv"
    file ".command.trace"
    file ".command.sh"
    file ".command.log"

    when:
    params.do_evaluate_log2abs_vs_rel_cor

    script:
    """
    evaluate-log2abs-vs-rel-cor.R \
        --ori $ori \
        --abs $abs \
        --rel $rel \
        --features $features \
        --filter 1 \
        --out ${dataset}_${exp_sim}_${full}_${method_replace_zero}_${method_transform_data}_${refgene}_${method_cor}_${method_eval}.csv \
        --method $method_eval
    sleep 30
    """

    stub:
    """
    echo evaluate-log2abs-vs-rel-cor.R \
        --ori $ori \
        --abs $abs \
        --rel $rel \
        --features $features \
        --filter 1 \
        --out ${dataset}_${exp_sim}_${full}_${method_replace_zero}_${method_transform_data}_${refgene}_${method_cor}_${method_eval}.csv \
        --method $method_eval
    touch ${dataset}_${exp_sim}_${full}_${method_replace_zero}_${method_transform_data}_${refgene}_${method_cor}_${method_eval}.csv
    """
}

process EVALUATE_LOG2ABS_VS_REL_COR_FILTER {
    label 'process_high'
    container 'suzannejin/scpropr:propr'
    tag "${dataset}-${exp_sim}-${full}-${method_replace_zero}-${method_transform_data}-${refgene}-${method_cor}-${method_eval}"
    storeDir "${params.outdir}/${dataset}/evaluate-filtered/${method_eval}/${dataset}_${exp_sim}_${full}_${method_replace_zero}_${method_transform_data}_${refgene}_${method_cor}"

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
          file(features),
          file(ori)
    each method_eval

    output:
    file "${dataset}_${exp_sim}_${full}_${method_replace_zero}_${method_transform_data}_${refgene}_${method_cor}_${method_eval}.csv"
    file ".command.trace"
    file ".command.sh"
    file ".command.log"

    when:
    params.do_evaluate_log2abs_vs_rel_cor_filter

    script:
    """
    evaluate-log2abs-vs-rel-cor.R \
        --ori $ori \
        --abs $abs \
        --rel $rel \
        --features $features \
        --filter 0.2 \
        --out ${dataset}_${exp_sim}_${full}_${method_replace_zero}_${method_transform_data}_${refgene}_${method_cor}_${method_eval}.csv \
        --method $method_eval
    sleep 30
    """

    stub:
    """
    echo evaluate-log2abs-vs-rel-cor.R \
        --ori $ori \
        --abs $abs \
        --rel $rel \
        --features $features \
        --filter 0.2 \
        --out ${dataset}_${exp_sim}_${full}_${method_replace_zero}_${method_transform_data}_${refgene}_${method_cor}_${method_eval}.csv \
        --method $method_eval
    touch ${dataset}_${exp_sim}_${full}_${method_replace_zero}_${method_transform_data}_${refgene}_${method_cor}_${method_eval}.csv
    """
}