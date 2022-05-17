
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
    params.do_plot_abs_vs_rel_cor && params.do_transform_abs

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
    sleep 30
    """
}

/* plot coefficients computed on log2 absolute data vs coefficients computed on transformed relative data - and colored by gene dropout */

process PLOT_LOG2ABS_VS_REL_COR_COLORED {
    label 'process_shigh'
    container 'suzannejin/scpropr:plot'
    tag "${dataset}-${exp_sim}-${full}-${method_replace_zero}-${color}-seed${params.random_seed}"
    publishDir "${params.outdir}/${dataset}/plot/log2abs-vs-rel-cor-colored/${dataset}_${exp_sim}_${full}_${method_replace_zero}_${color}_seed${params.random_seed}", mode: params.publish_dir_mode 

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
    def filter = full == 'nozero' ? 0 : 0.2
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
    def filter = full == 'nozero' ? 0 : 0.2
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
    sleep 30
    """
}

process PLOT_LOG2ABS_VS_REL_COR_COLORED_INDIV {
    label 'process_high'
    container 'suzannejin/scpropr:plot'
    tag "${dataset}-${exp_sim}-${full}-${method_replace_zero}-${method_transform_data}-${refgene}-${method_cor}-seed${params.random_seed}"
    // publishDir "${params.outdir}/${dataset}/plot/log2abs-vs-rel-cor-colored-indiv/${dataset}_${exp_sim}_${full}_${method_replace_zero}_${method_transform_data}_${refgene}_${method_cor}_seed${params.random_seed}", publish_dir_mode: params.publish_dir_mode
    storeDir "${params.outdir}/${dataset}/plot/log2abs-vs-rel-cor-colored-indiv/${dataset}_${exp_sim}_${full}_${method_replace_zero}_${method_transform_data}_${refgene}_${method_cor}_seed${params.random_seed}"

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
    def filter = full == 'nozero' ? 0 : 0.2
    def seed = params.random_seed
    """
    plot-log2abs-vs-rel-cor-colored-indiv.R \
        --ori $ori \
        --abs $abs \
        --rel $rel \
        --features $features \
        $refgene1 \
        --filter $filter \
        --seed $seed \
        --outdir .
    sleep 30
    """

    stub:
    def refgene1 = refgene == 'NA' ? '' : "--refgene $refgene"
    def filter = full == 'nozero' ? 0 : 0.2
    def seed = params.random_seed
    """
    echo plot-log2abs-vs-rel-cor-colored-indiv.R \
        --ori $ori \
        --abs $abs \
        --rel $rel \
        --features $features \
        $refgene1 \
        --filter $filter \
        --seed $seed \
        --outdir .
    touch log2abs-vs-rel-cor-colored.png log2abs-vs-rel-cor-colored-filtered.png log2abs-vs-rel-cor-colored-filtered.png 
    sleep 30
    """
}

process PLOT_LOG2ABS_VS_REL_COR_COLORED_INDIV_GM {
    label 'process_high'
    container 'suzannejin/scpropr:plot'
    tag "${dataset}-${exp_sim}-${full}-${method_replace_zero}-${method_transform_data}-${refgene}-${method_cor}-seed${params.random_seed}"
    // publishDir "${params.outdir}/${dataset}/plot/log2abs-vs-rel-cor-colored-indiv-gm/${dataset}_${exp_sim}_${full}_${method_replace_zero}_${method_transform_data}_${refgene}_${method_cor}_seed${params.random_seed}", publish_dir_mode: params.publish_dir_mode
    storeDir "${params.outdir}/${dataset}/plot/log2abs-vs-rel-cor-colored-indiv-gm/${dataset}_${exp_sim}_${full}_${method_replace_zero}_${method_transform_data}_${refgene}_${method_cor}_seed${params.random_seed}"

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
    file "log2abs-vs-rel-cor-colored*-gm.png"
    file ".command.trace"
    file ".command.sh"

    when:
    params.do_plot_log2abs_vs_rel_cor_colored_indiv_gm

    script:
    def refgene1 = refgene == 'NA' ? '' : "--refgene $refgene"
    def filter = full == 'nozero' ? 0 : 0.2
    def seed = params.random_seed
    """
    plot-log2abs-vs-rel-cor-colored-indiv.R \
        --ori $ori \
        --abs $abs \
        --rel $rel \
        --features $features \
        $refgene1 \
        --filter $filter \
        --gm \
        --seed $seed \
        --outdir .
    sleep 30
    """

    stub:
    def refgene1 = refgene == 'NA' ? '' : "--refgene $refgene"
    def filter = full == 'nozero' ? 0 : 0.2
    def seed = params.random_seed
    """
    echo plot-log2abs-vs-rel-cor-colored-indiv.R \
        --ori $ori \
        --abs $abs \
        --rel $rel \
        --features $features \
        $refgene1 \
        --filter $filter \
        --gm \
        --seed $seed \
        --outdir .
    touch log2abs-vs-rel-cor-colored-gm.png log2abs-vs-rel-cor-colored-filtered-gm.png log2abs-vs-rel-cor-colored-filtered-gm.png 
    sleep 30
    """
}