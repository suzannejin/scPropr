
/* compute correlation coefficients */

process GET_CORRELATION {
    
    label 'process_high_long'
    container 'suzannejin/scpropr:propr'
    tag "$dataset-$exp_sim-$full-$abs_rel-$method_replace_zero-$method_transform_data-$refgene-$method_cor"
    storeDir "${params.outdir}/${dataset}/cor/${dataset}_${exp_sim}_${full}_${abs_rel}_${method_replace_zero}_${method_transform_data}_${refgene}_${method_cor}"
    // publishDir "${params.outdir}/${dataset}/cor/${dataset}_${exp_sim}_${full}_${abs_rel}_${method_replace_zero}_${method_transform_data}_${refgene}_${method_cor}", mode: params.publish_dir_mode

    input:
    tuple val(dataset), 
          val(exp_sim), 
          val(full), 
          val(abs_rel), 
          val(method_replace_zero),
          val(method_transform_data),
          val(refgene),
          file(count), 
          file(features), 
          file(barcodes)
    each method_cor

    output:
    tuple file("${dataset}_${exp_sim}_${full}_${abs_rel}_${method_replace_zero}_${method_transform_data}_${refgene}_${method_cor}.csv.gz"),
          file(features), 
          file("${dataset}_${exp_sim}_${full}_${abs_rel}_${method_replace_zero}_${method_transform_data}_${refgene}_${method_cor}.lambda"),
          file(".command.trace"),
          file(".command.sh"),
          file(".command.log")

    script:
    """
    get-correlation.R \
        -i $count \
        -o ${dataset}_${exp_sim}_${full}_${abs_rel}_${method_replace_zero}_${method_transform_data}_${refgene}_${method_cor}.csv.gz \
        -o2 ${dataset}_${exp_sim}_${full}_${abs_rel}_${method_replace_zero}_${method_transform_data}_${refgene}_${method_cor}.lambda \
        -m $method_cor
    sleep 30
    """

    stub:
    """
    echo get-correlation.R \
        -i $count \
        -o ${dataset}_${exp_sim}_${full}_${abs_rel}_${method_replace_zero}_${method_transform_data}_${refgene}_${method_cor}.csv.gz \
        -o2 ${dataset}_${exp_sim}_${full}_${abs_rel}_${method_replace_zero}_${method_transform_data}_${refgene}_${method_cor}.lambda \
        -m $method_cor
    touch ${dataset}_${exp_sim}_${full}_${abs_rel}_${method_replace_zero}_${method_transform_data}_${refgene}_${method_cor}.csv.gz
    touch ${dataset}_${exp_sim}_${full}_${abs_rel}_${method_replace_zero}_${method_transform_data}_${refgene}_${method_cor}.lambda
    """
}