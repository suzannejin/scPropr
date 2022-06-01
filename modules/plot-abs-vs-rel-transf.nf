
/* plot transformed/normalized counts on absolute vs relative data */

process PLOT_ABS_VS_REL_TRANSF {
    
    container 'suzannejin/scpropr:plot'
    tag "${dataset}-${exp_sim}-${full}"
    publishDir "${params.outdir}/${dataset}/plot/abs-vs-rel-transf/${dataset}_${exp_sim}_${full}", mode: params.publish_dir_mode

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
    file ".command.log"

    when:
    params.do_plot_abs_vs_rel_transf && params.do_transform_abs

    script:
    def abs2           = abs.join(' ')
    def rel2           = rel.join(' ')
    def methods_zero   = method_replace_zero.join(' ')
    def methods_transf = method_transform_data.join(' ')
    def refgenes       = refgene.join(' ')
    """
    plot-abs-vs-rel-transf.R \
        --abs $abs2 \
        --rel $rel2 \
        --method_zero $methods_zero \
        --method_transf $methods_transf \
        --refgene $refgenes \
        --output ${dataset}_${exp_sim}_${full}.png
    """

    stub:
    def abs2           = abs.join(' ')
    def rel2           = rel.join(' ')
    def methods_zero   = method_replace_zero.join(' ')
    def methods_transf = method_transform_data.join(' ')
    def refgenes       = refgene.join(' ')
    """
    echo plot-abs-vs-rel-transf.R \
        --abs $abs2 \
        --rel $rel2 \
        --method_zero $methods_zero \
        --method_transf $methods_transf \
        --refgene $refgenes \
        --output ${dataset}_${exp_sim}_${full}.png
    touch ${dataset}_${exp_sim}_${full}.png
    """
}

process PLOT_LOG2ABS_VS_REL_TRANSF {
    
    container 'suzannejin/scpropr:plot'
    tag "${dataset}-${exp_sim}-${full}"
    publishDir "${params.outdir}/${dataset}/plot/log2abs-vs-rel-transf/${dataset}_${exp_sim}_${full}", mode: params.publish_dir_mode

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
    file "${dataset}_${exp_sim}_${full}.png"
    file ".command.trace"
    file ".command.sh"
    file ".command.log"

    when:
    params.do_plot_log2abs_vs_rel_transf

    script:
    def abs2           = abs.join(' ')
    def rel2           = rel.join(' ')
    def methods_zero   = method_replace_zero.join(' ')
    def methods_transf = method_transform_data.join(' ')
    def refgenes       = refgene.join(' ')
    """
    plot-log2abs-vs-rel-transf.R \
        --abs $abs2 \
        --rel $rel2 \
        --method_zero $methods_zero \
        --method_transf $methods_transf \
        --refgene $refgenes \
        --features $features \
        --output ${dataset}_${exp_sim}_${full}.png
    """

    stub:
    def abs2           = abs.join(' ')
    def rel2           = rel.join(' ')
    def methods_zero   = method_replace_zero.join(' ')
    def methods_transf = method_transform_data.join(' ')
    def refgenes       = refgene.join(' ')
    """
    echo plot-log2abs-vs-rel-transf.R \
        --abs $abs2 \
        --rel $rel2 \
        --method_zero $methods_zero \
        --method_transf $methods_transf \
        --refgene $refgenes \
        --features $features \
        --output ${dataset}_${exp_sim}_${full}.png
    touch ${dataset}_${exp_sim}_${full}.png
    """
}