
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
    
    when:
    params.do_plot_abs_vs_rel_cor

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

    when:
    params.do_plot_log2abs_vs_rel_cor

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
    tag "${dataset}-${exp_sim}-${full}-${method_replace_zero}-${color}-seed${seed}"
    publishDir "${params.outdir}/${dataset}/plot/log2abs-vs-rel-cor-colored/${dataset}_${exp_sim}_${full}_${method_replace_zero}_${color}_seed${seed}", mode: params.publish_dir_mode 

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
    file "*.png"
    file ".command.trace"
    file ".command.sh"

    when:
    params.do_plot_log2abs_vs_rel_cor_colored

    script:
    def filter = 0.2
    def abs2 = abs.join(' ')
    def rel2 = rel.join(' ')
    def seed = params.random_seed
    """
    plot-log2abs-vs-rel-cor-colored.R \
        --ori $ori \
        --abs $abs2 \
        --rel $rel2 \
        --features $features \
        --color $color \
        --filter $filter \
        --seed $seed \
        --outdir .
    """

    stub:
    def filter = 0.2
    def abs2 = abs.join(' ')
    def rel2 = rel.join(' ')
    def seed = params.random_seed
    """
    echo plot-log2abs-vs-rel-cor-colored.R \
        --ori $ori \
        --abs $abs2 \
        --rel $rel2 \
        --features $features \
        --color $color \
        --filter $filter \
        --seed $seed \
        --outdir .
    touch log2abs-vs-rel-cor-colored.png log2abs-vs-rel-cor-black.png log2abs-vs-rel-cor-colored-filtered.png log2abs-vs-rel-cor-black-filtered.png
    """
}

process PLOT_LOG2ABS_VS_REL_COR_COLORED_INDIV {
    label 'process_high'
    container 'suzannejin/scpropr:plot'
    tag "${dataset}-${exp_sim}-${full}-${method_replace_zero}-${method_transform_data}-${refgene}-${method_cor}-seed${seed}"
    storeDir "${params.outdir}/${dataset}/plot/log2abs-vs-rel-cor-colored-indiv/${dataset}_${exp_sim}_${full}_${method_replace_zero}_${method_transform_data}_${refgene}_${method_cor}_seed${seed}"

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
    file "log2abs-vs-rel-cor-colored*.png"
    file ".command.trace"
    file ".command.sh"

    when:
    params.do_plot_log2abs_vs_rel_cor_colored_indiv

    script:
    def refgene1 = refgene == 'NA' ? '' : "--refgene $refgene"
    def filter = '--filter 0.2'
    def seed = params.random_seed
    """
    plot-log2abs-vs-rel-cor-colored-indiv.R \
        --ori $ori \
        --abs $abs \
        --rel $rel \
        --features $features \
        $refgene1 \
        $filter \
        --seed $seed \
        --outdir .
    sleep 30
    """

    stub:
    def refgene1 = refgene == 'NA' ? '' : "--refgene $refgene"
    def filter = '--filter 0.2'
    def seed = params.random_seed
    """
    plot-log2abs-vs-rel-cor-colored-indiv.R \
        --ori $ori \
        --abs $abs \
        --rel $rel \
        --features $features \
        $refgene1 \
        $filter \
        --seed $seed \
        --outdir .
    touch ${dataset}_${exp_sim}_${full}_${method_replace_zero}_${method_transform_data}_${refgene}_${method_cor}.png
    """
}