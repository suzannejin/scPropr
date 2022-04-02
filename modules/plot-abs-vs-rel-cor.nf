
/* plot correlation coefficients computed on absolute vs relative data */

process PLOT_ABS_VS_REL_COR {
    label 'process_high_long'
    container 'suzannejin/scpropr:plot'
    tag "${dataset}-${exp_sim}-${full}-${method_replace_zero}"
    publishDir "${params.outdir}/${dataset}/plot/abs-vs-rel-cor/${dataset}_${exp_sim}_${full}_${method_replace_zero}", mode: params.publish_dir_mode

    input:
    tuple val(dataset), 
          val(exp_sim), 
          val(full), 
          val(method_replace_zero), 
          val(method_transform_data), 
          val(refgene),
          val(method_cor),
          file(abs), 
          file(rel)

    output:
    file "${dataset}_${exp_sim}_${full}_${method_replace_zero}.png"
    file ".command.trace"
    file ".command.sh"

    script:
    def abs2 = abs.join(' ')
    def rel2 = rel.join(' ')
    def methods_transf = method_transform_data.join(' ')
    def refgenes = refgene.join(' ')
    def methods_cor = method_cor.join(' ')
    """
    plot-abs-vs-rel-cor.R \
        --abs $abs2 \
        --rel $rel2 \
        --method_transf $methods_transf \
        --refgene $refgenes \
        --method_cor $methods_cor \
        --output ${dataset}_${exp_sim}_${full}_${method_replace_zero}.png
    """

    stub:
    def abs2 = abs.join(' ')
    def rel2 = rel.join(' ')
    def methods_transf = method_transform_data.join(' ')
    def refgenes = refgene.join(' ')
    def methods_cor = method_cor.join(' ')
    """
    echo plot-abs-vs-rel-cor.R \
        --abs $abs2 \
        --rel $rel2 \
        --method_transf $methods_transf \
        --refgene $refgenes \
        --method_cor $methods_cor \
        --output ${dataset}_${exp_sim}_${full}_${method_replace_zero}.png
    touch ${dataset}_${exp_sim}_${full}_${method_replace_zero}.png
    """
}

process PLOT_LOG2ABS_VS_REL_COR {
    label 'process_high_long'
    container 'suzannejin/scpropr:plot'
    tag "${dataset}-${exp_sim}-${full}-${method_replace_zero}"
    publishDir "${params.outdir}/${dataset}/plot/log2abs-vs-rel-cor/${dataset}_${exp_sim}_${full}_${method_replace_zero}", mode: params.publish_dir_mode 

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

    output:
    file "${dataset}_${exp_sim}_${full}_${method_replace_zero}.png"
    file ".command.trace"
    file ".command.sh"

    script:
    def abs2 = abs.join(' ')
    def rel2 = rel.join(' ')
    """
    plot-log2abs-vs-rel-cor.R \
        --abs $abs2 \
        --rel $rel2 \
        --features $features \
        --output ${dataset}_${exp_sim}_${full}_${method_replace_zero}.png
    """

    stub:
    def abs2 = abs.join(' ')
    def rel2 = rel.join(' ')
    """
    echo plot-log2abs-vs-rel-cor.R \
        --abs $abs2 \
        --rel $rel2 \
        --features $features \
        --output ${dataset}_${exp_sim}_${full}_${method_replace_zero}.png
    touch ${dataset}_${exp_sim}_${full}_${method_replace_zero}.png
    """
}