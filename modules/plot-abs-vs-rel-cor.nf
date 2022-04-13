
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

/* plot coefficients computed on log2 absolute data vs coefficients computed on transformed relative data - and colored by gene dropout */

process PLOT_LOG2ABS_VS_REL_COR_COLORED {
    label 'process_shigh'
    container 'suzannejin/scpropr:plot'
    tag "${dataset}-${exp_sim}-${full}-${method_replace_zero}-${color}"
    publishDir "${params.outdir}/${dataset}/plot/log2abs-vs-rel-cor-colored/${dataset}_${exp_sim}_${full}_${method_replace_zero}_${color}", mode: params.publish_dir_mode 

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
    each color

    output:
    file "${dataset}_${exp_sim}_${full}_${method_replace_zero}_${color}.png"
    file ".command.trace"
    file ".command.sh"

    script:
    def abs2 = abs.join(' ')
    def rel2 = rel.join(' ')
    """
    plot-log2abs-vs-rel-cor-colored.R \
        --ori $ori \
        --abs $abs2 \
        --rel $rel2 \
        --features $features \
        --color $color \
        --output ${dataset}_${exp_sim}_${full}_${method_replace_zero}_${color}.png
    """

    stub:
    def abs2 = abs.join(' ')
    def rel2 = rel.join(' ')
    """
    echo plot-log2abs-vs-rel-cor-colored.R \
        --ori $ori \
        --abs $abs2 \
        --rel $rel2 \
        --features $features \
        --color $color \
        --output ${dataset}_${exp_sim}_${full}_${method_replace_zero}_${color}.png
    touch ${dataset}_${exp_sim}_${full}_${method_replace_zero}_${color}.png
    """
}

process PLOT_LOG2ABS_VS_REL_COR_COLORED_INDIV {
    container 'suzannejin/scpropr:plot'
    tag "${dataset}-${exp_sim}-${full}-${method_replace_zero}-${method_transform_data}-${refgene}-${method_cor}-${method_eval}"
    storeDir "${params.outdir}/${dataset}/plot/log2abs-vs-rel-cor-colored-indiv/${dataset}_${exp_sim}_${full}_${method_replace_zero}_${method_transform_data}_${refgene}_${method_cor}"

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

    output:
    file "${dataset}_${exp_sim}_${full}_${method_replace_zero}_${method_transform_data}_${refgene}_${method_cor}.png"
    file ".command.trace"
    file ".command.sh"

    script:
    def refgene1 = refgene == 'NA' ? '' : "--refgene $refgene"
    """
    plot-log2abs-vs-rel-cor-colored-indiv.R \
        --ori $ori \
        --abs $abs \
        --rel $rel \
        --features $features \
        $refgene1 \
        --out ${dataset}_${exp_sim}_${full}_${method_replace_zero}_${method_transform_data}_${refgene}_${method_cor}.png 
    sleep 30
    """

    stub:
    """
    plot-log2abs-vs-rel-cor-colored-indiv.R \
        --ori $ori \
        --abs $abs \
        --rel $rel \
        --features $features \
        $refgene1 \
        --out ${dataset}_${exp_sim}_${full}_${method_replace_zero}_${method_transform_data}_${refgene}_${method_cor}.png 
    touch ${dataset}_${exp_sim}_${full}_${method_replace_zero}_${method_transform_data}_${refgene}_${method_cor}.png
    """
}