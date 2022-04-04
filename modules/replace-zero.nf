
process REPLACE_ZERO {
    label 'process_high_long'
    container 'suzannejin/scpropr:propr'
    tag "${dataset}-${exp_sim}-${full}-${abs_rel}-${method_replace_zero}"
    storeDir "${params.outdir}/${dataset}/imputed/${dataset}_${exp_sim}_${full}_${abs_rel}_${method_replace_zero}"

    input:
    tuple val(dataset), 
          val(exp_sim), 
          val(full), 
          val(abs_rel), 
          file(count), 
          file(features), 
          file(barcodes)
    each method_replace_zero

    output:
    tuple file("${dataset}_${exp_sim}_${full}_${abs_rel}_${method_replace_zero}.csv.gz"),
          file(features),
          file(barcodes),
          file(".command.trace"),
          file(".command.sh")

    script:
    """
    replace-zero.R \
        -i $count \
        -o ${dataset}_${exp_sim}_${full}_${abs_rel}_${method_replace_zero}.csv.gz \
        --method $method_replace_zero
    """
}