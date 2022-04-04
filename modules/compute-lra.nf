process COMPUTE_LRA {
    label 'process_high'
    container 'suzannejin/scpropr:propr'
    tag "${dataset}-${exp_sim}-${full}-${abs_rel}-${method_replace_zero}"
    storeDir "${params.outdir}/${dataset}/lra/${dataset}_${exp_sim}_${full}_${abs_rel}_${method_replace_zero}"

    input:
    tuple val(dataset),
          val(exp_sim),
          val(full),
          val(abs_rel),
          val(method_replace_zero),
          file(count),
          file(features),
          file(barcodes)
    
    output:
    tuple file(count),
          file(features),
          file(barcodes),
          file("${dataset}_${exp_sim}_${full}_${abs_rel}_${method_replace_zero}_lra.rds"),
          file(".command.trace"),
          file(".command.sh")
    
    script:
    """
    compute-lra.R \
        $count \
        ${dataset}_${exp_sim}_${full}_${abs_rel}_${method_replace_zero}_lra.rds
    """
}