
/* plot transformed/normalized counts on absolute vs relative data */

process PLOT_ABS_VS_REL_TRANSF {
    
    container 'suzannejin/scpropr:plot'
    tag "${dataset}-${exp_sim}-${full}-${method_replace_zero}"
    publishDir "${params.outdir}/${dataset}/plot/abs-vs-rel-transf/${dataset}_${exp_sim}_${full}_${method_replace_zero}", mode: params.publish_dir_mode

    input:
    tuple val(dataset), 
          val(exp_sim), 
          val(full), 
          val(method_replace_zero), 
          val(method_transform_data), 
          val(refgene), 
          file(abs), 
          file(rel)

    output:
    file "${dataset}_${exp_sim}_${full}_${method_replace_zero}.png"
    file ".command.trace"
    file ".command.sh"

    when:
    params.do_plot_abs_vs_rel_transf && params.do_transform_abs

    script:
    def abs2 = abs.join(' ')
    def rel2 = rel.join(' ')
    def methods = method_transform_data.join(' ')
    def refgenes = refgene.join(' ')
    """
    plot-abs-vs-rel-transf.R \
        --abs $abs2 \
        --rel $rel2 \
        --method $methods \
        --refgene $refgenes \
        --output ${dataset}_${exp_sim}_${full}_${method_replace_zero}.png
    """

    stub:
    def abs2 = abs.join(' ')
    def rel2 = rel.join(' ')
    def methods = method_transform_data.join(' ')
    def refgenes = refgene.join(' ')
    """
    echo plot-abs-vs-rel-transf.R \
        --abs $abs2 \
        --rel $rel2 \
        --method $methods \
        --refgene $refgenes \
        --output ${dataset}_${exp_sim}_${full}_${method_replace_zero}.png
    touch ${dataset}_${exp_sim}_${full}_${method_replace_zero}.png
    """
}

process PLOT_LOG2ABS_VS_REL_TRANSF {
    
    container 'suzannejin/scpropr:plot'
    tag "${dataset}-${exp_sim}-${full}-${method_replace_zero}"
    publishDir "${params.outdir}/${dataset}/plot/log2abs-vs-rel-transf/${dataset}_${exp_sim}_${full}_${method_replace_zero}", mode: params.publish_dir_mode

    input:
    tuple val(dataset), 
          val(exp_sim), 
          val(full), 
          val(method_replace_zero), 
          val(method_transform_data), 
          val(refgene), 
          file(abs), 
          file(rel),
          file(features)

    output:
    file "${dataset}_${exp_sim}_${full}_${method_replace_zero}.png"
    file ".command.trace"
    file ".command.sh"

    when:
    params.do_plot_log2abs_vs_rel_transf

    script:
    def abs2 = abs.join(' ')
    def rel2 = rel.join(' ')
    def methods = method_transform_data.join(' ')
    def refgenes = refgene.join(' ')
    """
    plot-log2abs-vs-rel-transf.R \
        --abs $abs2 \
        --rel $rel2 \
        --method $methods \
        --refgene $refgenes \
        --features $features \
        --output ${dataset}_${exp_sim}_${full}_${method_replace_zero}.png
    """

    stub:
    def abs2 = abs.join(' ')
    def rel2 = rel.join(' ')
    def methods = method_transform_data.join(' ')
    def refgenes = refgene.join(' ')
    """
    echo plot-log2abs-vs-rel-transf.R \
        --abs $abs2 \
        --rel $rel2 \
        --method $methods \
        --refgene $refgenes \
        --features $features \
        --output ${dataset}_${exp_sim}_${full}_${method_replace_zero}.png
    touch ${dataset}_${exp_sim}_${full}_${method_replace_zero}.png
    """
}