/* plot original vs transformed data*/

process PLOT_ORIGINAL_VS_TRANSFORMED {

    container 'suzannejin/scpropr:plot'
    tag "${dataset}-${exp_sim}-${full}-${abs_rel}"
    publishDir "${params.outdir}/${dataset}/plot/original-vs-transformed/${dataset}_${exp_sim}_${full}_${abs_rel}", mode: params.publish_dir_mode

    input:
    tuple val(dataset), 
          val(exp_sim), 
          val(full), 
          val(abs_rel), 
          val(method_replace_zero), 
          val(method_transform_data), 
          val(refgene), 
          file(features),
          file(barcodes),
          file(original), 
          file(transformed)

    output:
    file "${dataset}_${exp_sim}_${full}_${abs_rel}.png"
    file "${dataset}_${exp_sim}_${full}_${abs_rel}-logx.png"
    file "${dataset}_${exp_sim}_${full}_${abs_rel}-logy.png"
    file "${dataset}_${exp_sim}_${full}_${abs_rel}-logxy.png"
    file ".command.trace"
    file ".command.sh"
    file ".command.log"

    when:
    params.do_plot_original_vs_transf

    script:
    def methods_zero   = method_replace_zero.join(' ')
    def methods_transf = method_transform_data.join(' ')
    def refgenes       = refgene.join(' ')
    def transformed2   = transformed.join(' ')
    """
    plot-original-vs-transformed.R \
        --original $original \
        --transformed $transformed2 \
        --method_zero $methods_zero \
        --method_transf $methods_transf \
        --refgene $refgenes \
        --features $features \
        --outdir . \
        --prefix ${dataset}_${exp_sim}_${full}_${abs_rel}
    """

    stub:
    def methods_zero   = method_replace_zero.join(' ')
    def methods_transf = method_transform_data.join(' ')
    def refgenes       = refgene.join(' ')
    def transformed2   = transformed.join(' ')
    """
    echo plot-original-vs-transformed.R \
        --original $original \
        --transformed $transformed2 \
        --method_zero $methods_zero \
        --method_transf $methods_transf \
        --refgene $refgenes \
        --features $features \
        --outdir . \
        --prefix ${dataset}_${exp_sim}_${full}_${abs_rel}
    touch ${dataset}_${exp_sim}_${full}_${abs_rel}.png
    touch ${dataset}_${exp_sim}_${full}_${abs_rel}-logx.png
    touch ${dataset}_${exp_sim}_${full}_${abs_rel}-logy.png
    touch ${dataset}_${exp_sim}_${full}_${abs_rel}-logxy.png
    """
}

/* plot log vs transformed data */

process PLOT_LOG_VS_TRANSFORMED {

    container 'suzannejin/scpropr:plot'
    tag "${dataset}-${exp_sim}-${full}-${abs_rel}"
    publishDir "${params.outdir}/${dataset}/plot/log-vs-transformed/${dataset}_${exp_sim}_${full}_${abs_rel}", mode: params.publish_dir_mode

    input:
    tuple val(dataset), 
          val(exp_sim), 
          val(full), 
          val(abs_rel), 
          val(method_replace_zero), 
          val(method_transform_data), 
          val(refgene),
          file(features), 
          file(barcodes),
          file(log2),
          file(transformed)

    output:
    path "${dataset}_${exp_sim}_${full}_${abs_rel}.png"
    path "${dataset}_${exp_sim}_${full}_${abs_rel}-logx.png"
    path "${dataset}_${exp_sim}_${full}_${abs_rel}-logy.png"
    path "${dataset}_${exp_sim}_${full}_${abs_rel}-logxy.png"
    path ".command.trace"
    path ".command.sh"
    path ".command.log"

    when:
    params.do_plot_log_vs_transf

    script:
    def methods_zero   = method_replace_zero.join(' ')
    def methods_transf = method_transform_data.join(' ')
    def refgenes       = refgene.join(' ')
    def transformed2   = transformed.join(' ')
    """
    plot-original-vs-transformed.R \
        --original $log2 \
        --transformed $transformed2 \
        --method_zero $methods_zero \
        --method_transf $methods_transf \
        --refgene $refgenes \
        --features $features \
        --outdir . \
        --prefix ${dataset}_${exp_sim}_${full}_${abs_rel} \
        --xlab log2
    """

    stub:
    def methods_zero   = method_replace_zero.join(' ')
    def methods_transf = method_transform_data.join(' ')
    def refgenes       = refgene.join(' ')
    def transformed2   = transformed.join(' ')
    """
    echo plot-original-vs-transformed.R \
        --original $log2 \
        --transformed $transformed2 \
        --method_zero $methods_zero \
        --method_transf $methods_transf \
        --refgene $refgenes \
        --features $features \
        --outdir . \
        --prefix ${dataset}_${exp_sim}_${full}_${abs_rel} \
        --xlab log2
    touch ${dataset}_${exp_sim}_${full}_${abs_rel}.png
    touch ${dataset}_${exp_sim}_${full}_${abs_rel}-logx.png
    touch ${dataset}_${exp_sim}_${full}_${abs_rel}-logy.png
    touch ${dataset}_${exp_sim}_${full}_${abs_rel}-logxy.png
    """
}