process TRANSF {

    tag "${cell_type}-${method_zero}-${method_transf}-${count_type}"
    publishDir "${params.outdir}/${cell_type}/data", mode: params.publish_dir_mode

    input:
    tuple val(cell_type), \
          val(count_type), \
          file(count_file), \
          file(nozero_file), \
          val(method_zero), \
          val(method_transf)

    output:
    tuple val(cell_type), \
          val(count_type), \
          file("${count_type}_${method_zero}_${method_transf}.csv.gz"), \
          val(method_zero), \
          val(method_transf), \
          emit: ch_transf

    script:
    def nozero = params.nozero_mode ? params.nozero_genes : "NA"

    """
    Rscript ${baseDir}/bin/transf-data.R \
        ${count_file} \
        ${count_type}_${method_zero}_${method_transf}.csv.gz \
        ${method_zero} \
        ${method_transf} \
        ${nozero}
    """
}

process CORR {

    tag "${cell_type}-${method_zero}-${method_transf}-${method_corr}-${count_type}"
    publishDir "${params.outdir}/${cell_type}/corr", mode: params.publish_dir_mode

    input:
    tuple val(cell_type), \
          val(count_type), \
          file(count_file), \
          val(method_zero), \
          val(method_transf), \
          val(method_corr)

    output:
    tuple val(cell_type), \
          file("${count_type}_${method_zero}_${method_transf}_${method_corr}.csv.gz"), \
          emit: ch_corr

    script:
    """
    Rscript ${baseDir}/bin/get-corr.R \
        ${count_file} \
        ${method_corr} \
        ${count_type}_${method_zero}_${method_transf}_${method_corr}.csv.gz
    """
}