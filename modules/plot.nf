process PLOT_ABS_VS_REL_CORR {
    cache false
    container 'suzannejin/scpropr:propr'
    tag "${cell_type}_${count_type2}"

    publishDir "${params.outdir}/${cell_type}/plot", mode: params.publish_dir_mode

    input:
    tuple val(cell_type), \
          val(count_type2), \
          file(corr)

    output:
    tuple val(cell_type), \
          val(count_type2), \
          file("*.png"), \
          emit: ch_plot_abs_vs_rel_corr

    script:
    def zero_method = params.show_methods_replace_zero
    def alr_ref = params.show_alr
    def transf_list1 = params.methods_transf_data
    def transf_list2 = params.show_methods_transf_data
    def corr_list1 = params.methods_corr
    def corr_list2 = params.show_methods_corr
    def random = params.nozero_mode ? 0 : 10000
    """
    Rscript ${baseDir}/bin/plot-abs-vs-rel-corr.R \
        ${cell_type} \
        ${count_type2} \
        ${zero_method} \
        ${alr_ref} \
        ${transf_list1} \
        ${corr_list1} \
        ${cell_type}_${count_type2}_${zero_method}_${alr_ref}_absolute_vs_relative_corr.png \
        ${random}

    Rscript ${baseDir}/bin/plot-abs-vs-rel-corr.R \
        ${cell_type} \
        ${count_type2} \
        ${zero_method} \
        ${alr_ref} \
        ${transf_list2} \
        ${corr_list2} \
        ${cell_type}_${count_type2}_${zero_method}_${alr_ref}_absolute_vs_relative_corr_show.png \
        ${random}
    """
}

process PLOT_ABS2_VS_REL_CORR {
    cache false
    container 'suzannejin/scpropr:propr'
    tag "${cell_type}_${count_type2}"

    publishDir "${params.outdir}/${cell_type}/plot", mode: params.publish_dir_mode

    input:
    tuple val(cell_type), \
          val(count_type2), \
          file(corr)

    output:
    tuple val(cell_type), \
          val(count_type2), \
          file("*.png"), \
          emit: ch_plot_abs2_vs_rel_corr

    script:
    def zero_method = params.show_methods_replace_zero
    def alr_ref = params.show_alr
    def transf_list1 = params.methods_transf_data
    def transf_list2 = params.show_methods_transf_data
    def corr_list1 = params.methods_corr
    def corr_list2 = params.show_methods_corr
    def random = params.nozero_mode ? 0 : 10000
    """
    Rscript ${baseDir}/bin/plot-abs2-vs-rel-corr.R \
        ${cell_type} \
        ${count_type2} \
        ${zero_method} \
        ${alr_ref} \
        ${transf_list1} \
        ${corr_list1} \
        ${cell_type}_${count_type2}_${zero_method}_${alr_ref}_absolute2_vs_relative_corr.png \
        ${random}

    Rscript ${baseDir}/bin/plot-abs2-vs-rel-corr.R \
        ${cell_type} \
        ${count_type2} \
        ${zero_method} \
        ${alr_ref} \
        ${transf_list2} \
        ${corr_list2} \
        ${cell_type}_${count_type2}_${zero_method}_${alr_ref}_absolute2_vs_relative_corr_show.png \
        ${random}
    """
}

process PLOT_ABS_VS_REL_TRANSF {
    cache false
    container 'suzannejin/scpropr:propr'
    tag "${cell_type}_${count_type2}"

    publishDir "${params.outdir}/${cell_type}/plot", mode: params.publish_dir_mode

    input:
    tuple val(cell_type), \
          val(count_type2), \
          file(transf)

    output:
    tuple val(cell_type), \
          val(count_type2), \
          file("*.png"), \
          emit: ch_plot_abs_vs_rel_transf

    script:
    def zero_method = params.show_methods_replace_zero
    def alr_ref = params.show_alr
    def transf_list1 = params.methods_transf_data
    def transf_list2 = params.show_methods_transf_data
    def random = 100000
    """
    Rscript ${baseDir}/bin/plot-abs-vs-rel-transf.R \
        ${cell_type} \
        ${count_type2} \
        ${zero_method} \
        ${alr_ref} \
        ${transf_list1} \
        ${cell_type}_${count_type2}_${zero_method}_${alr_ref}_absolute_vs_relative_transf.png \
        ${random}
    """
}

process PLOT_ABS2_VS_REL_TRANSF {
    cache false
    container 'suzannejin/scpropr:propr'
    tag "${cell_type}_${count_type2}"

    publishDir "${params.outdir}/${cell_type}/plot", mode: params.publish_dir_mode

    input:
    tuple val(cell_type), \
          val(count_type2), \
          file(transf)

    output:
    tuple val(cell_type), \
          val(count_type2), \
          file("*.png"), \
          emit: ch_plot_abs2_vs_rel_transf

    script:
    def zero_method = params.show_methods_replace_zero
    def alr_ref = params.show_alr
    def transf_list1 = params.methods_transf_data
    def transf_list2 = params.show_methods_transf_data
    def random = 100000
    """
    Rscript ${baseDir}/bin/plot-abs2-vs-rel-transf.R \
        ${cell_type} \
        ${count_type2} \
        ${zero_method} \
        ${alr_ref} \
        ${transf_list1} \
        ${cell_type}_${count_type2}_${zero_method}_${alr_ref}_absolute2_vs_relative_transf.png \
        ${random}
    """
}

// process PLOT_MODEL_VS_SIM {

// }

// process PLOT_EXP_VS_SIM {

// }

