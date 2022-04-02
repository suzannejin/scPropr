/* plot original vs transformed data*/

process PLOT_ORIGINAL_VS_TRANSFORMED {

    container 'suzannejin/scpropr:plot'
    tag "${dataset}-${exp_sim}-${full}-${abs_rel}_${method_replace_zero}"
    publishDir "${params.outdir}/${dataset}/plot/original-vs-transformed/${dataset}_${exp_sim}_${full}_${abs_rel}_${method_replace_zero}", mode: params.publish_dir_mode

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
    file "${dataset}_${exp_sim}_${full}_${abs_rel}_${method_replace_zero}.png"
    file "${dataset}_${exp_sim}_${full}_${abs_rel}_${method_replace_zero}-logx.png"
    file "${dataset}_${exp_sim}_${full}_${abs_rel}_${method_replace_zero}-logy.png"
    file "${dataset}_${exp_sim}_${full}_${abs_rel}_${method_replace_zero}-logxy.png"
    file ".command.trace"
    file ".command.sh"

    script:
    def methods = method_transform_data.join(' ')
    def refgenes = refgene.join(' ')
    def transformed2 = transformed.join(' ')
    """
    plot-original-vs-transformed.R \
        --original $original \
        --transformed $transformed2 \
        --method $methods \
        --refgene $refgenes \
        --features $features \
        --outdir . \
        --prefix ${dataset}_${exp_sim}_${full}_${abs_rel}_${method_replace_zero}
    """

    stub:
    def methods = method_transform_data.join(' ')
    def refgenes = refgene.join(' ')
    def transformed2 = transformed.join(' ')
    """
    echo plot-original-vs-transformed.R \
        --original $original \
        --transformed $transformed2 \
        --method $methods \
        --refgene $refgenes \
        --features $features \
        --outdir . \
        --prefix ${dataset}_${exp_sim}_${full}_${abs_rel}_${method_replace_zero}
    touch ${dataset}_${exp_sim}_${full}_${abs_rel}_${method_replace_zero}.png
    touch ${dataset}_${exp_sim}_${full}_${abs_rel}_${method_replace_zero}-logx.png
    touch ${dataset}_${exp_sim}_${full}_${abs_rel}_${method_replace_zero}-logy.png
    touch ${dataset}_${exp_sim}_${full}_${abs_rel}_${method_replace_zero}-logxy.png
    """
}

/* plot log vs transformed data */

process PLOT_LOG_VS_TRANSFORMED {

    container 'suzannejin/scpropr:plot'
    tag "${dataset}-${exp_sim}-${full}-${abs_rel}_${method_replace_zero}"
    publishDir "${params.outdir}/${dataset}/plot/log-vs-transformed/${dataset}_${exp_sim}_${full}_${abs_rel}_${method_replace_zero}", mode: params.publish_dir_mode

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
    path "${dataset}_${exp_sim}_${full}_${abs_rel}_${method_replace_zero}.png"
    path "${dataset}_${exp_sim}_${full}_${abs_rel}_${method_replace_zero}-logx.png"
    path "${dataset}_${exp_sim}_${full}_${abs_rel}_${method_replace_zero}-logy.png"
    path "${dataset}_${exp_sim}_${full}_${abs_rel}_${method_replace_zero}-logxy.png"
    path ".command.trace"
    path ".command.sh"

    script:
    def methods = method_transform_data.join(' ')
    def refgenes = refgene.join(' ')
    def transformed2 = transformed.join(' ')
    """
    plot-original-vs-transformed.R \
        --original $log2 \
        --transformed $transformed2 \
        --method $methods \
        --refgene $refgenes \
        --features $features \
        --outdir . \
        --prefix ${dataset}_${exp_sim}_${full}_${abs_rel}_${method_replace_zero} \
        --xlab log2
    """

    stub:
    def methods = method_transform_data.join(' ')
    def refgenes = refgene.join(' ')
    def transformed2 = transformed.join(' ')
    """
    echo plot-original-vs-transformed.R \
        --original $log2 \
        --transformed $transformed2 \
        --method $methods \
        --refgene $refgenes \
        --features $features \
        --outdir . \
        --prefix ${dataset}_${exp_sim}_${full}_${abs_rel}_${method_replace_zero} \
        --xlab log2
    touch ${dataset}_${exp_sim}_${full}_${abs_rel}_${method_replace_zero}.png
    touch ${dataset}_${exp_sim}_${full}_${abs_rel}_${method_replace_zero}-logx.png
    touch ${dataset}_${exp_sim}_${full}_${abs_rel}_${method_replace_zero}-logy.png
    touch ${dataset}_${exp_sim}_${full}_${abs_rel}_${method_replace_zero}-logxy.png
    """
}