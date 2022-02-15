process RELATIVE {

    container 'suzannejin/scpropr:propr'
    tag "${cell_type}-${count_type2}"
    
    input:
    tuple val(cell_type), \
        val(absolute), \
        val(count_type2), \
        file(count_file)

    output:
    tuple val(cell_type), \
          val("relative"), \
          val(count_type2), \
          file("relative_${count_type2}.csv.gz"), \
          emit: ch_relative

    script:
    """
    Rscript ${baseDir}/bin/get-relative.R \
        ${count_file} \
        relative_${count_type2}.csv.gz
        
    """
}

process TRANSF {

    container 'suzannejin/scpropr:propr'

    tag "${cell_type}-${method_zero}-${method_transf}-${count_type1}-${count_type2}-${ref_gene}"
    publishDir "${params.outdir}/${cell_type}/data", mode: params.publish_dir_mode

    input:
    tuple val(cell_type), \
          val(count_type1), \
          val(count_type2), \
          file(count_file), \
          file(nozero_file), \
          val(method_zero), \
          val(method_transf), \
          val(ref_gene)

    output:
    tuple val(cell_type), \
          val(count_type1), \
          val(count_type2), \
          file("${count_type1}_${count_type2}_${method_zero}_${method_transf}_${ref_gene}.csv.gz"), \
          val(method_zero), \
          val(method_transf), \
          val(ref_gene), \
          emit: ch_transf

    script:
    def nozero = params.nozero_mode ? nozero_file : "NA"

    """
    Rscript ${baseDir}/bin/transf-data.R \
        ${count_file} \
        ${count_type1}_${count_type2}_${method_zero}_${method_transf}_${ref_gene}.csv.gz \
        ${method_zero} \
        ${method_transf} \
        ${nozero} \
        ${ref_gene}
    """
}

process CORR {

    container 'suzannejin/scpropr:propr'

    tag "${cell_type}-${method_zero}-${method_transf}-${method_corr}-${count_type1}-${count_type2}-${ref_gene}"
    publishDir "${params.outdir}/${cell_type}/corr", mode: params.publish_dir_mode

    input:
    tuple val(cell_type), \
          val(count_type1), \
          val(count_type2), \
          file(count_file), \
          val(method_zero), \
          val(method_transf), \
          val(ref_gene), \
          val(method_corr)

    output:
    tuple val(cell_type), \
          val(count_type1), \
          val(count_type2), \
          file("${count_type1}_${count_type2}_${method_zero}_${method_transf}_${ref_gene}_${method_corr}.csv.gz"), \
          val(method_zero), \
          val(method_transf), \
          val(ref_gene), \
          val(method_corr), \
          emit: ch_corr

    script:
    """
    Rscript ${baseDir}/bin/get-corr.R \
        ${count_file} \
        ${method_corr} \
        ${count_type1}_${count_type2}_${method_zero}_${method_transf}_${ref_gene}_${method_corr}.csv.gz \
        ${ref_gene}
    """
}