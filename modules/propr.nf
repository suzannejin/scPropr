process TRANSF {

    tag "${cell_type}-${method}-${count_type}"
    publishDir "${params.outdir}/${cell_type}/data", 
        mode: params.publish_dir_mode,
        saveAs: { filename -> if (filename.indexOf("_none.csv.gz") > 0) null else filename }

    input:
    tuple val(cell_type), \
          val(count_type), \
          file(count_file), \
          val(method)

    output:
    tuple val(cell_type), \
          val(count_type), \
          file("${count_type}_${method}.csv.gz"), \
          val(method), \
          emit: ch_clr

    script:
    if (method == "none")
    """
    mv ${count_file} ${count_type}_${method}.csv.gz
    """
    else
    """
    Rscript ${baseDir}/bin/transf-data.R \
        ${count_file} \
        ${method} \
        ${count_type}_${method}.csv.gz
    """
}

process CORR {

    tag "${cell_type}-${method1}-${method2}-${count_type}"
    publishDir "${params.outdir}/${cell_type}/corr", mode: params.publish_dir_mode

    input:
    tuple val(cell_type), \
          val(count_type), \
          file(count_file), \
          val(method1), \
          val(method2)

    output:
    tuple val(cell_type), \
          file("${count_type}_${method1}_${method2}.csv.gz"), \
          emit: ch_corr

    script:
    """
    Rscript ${baseDir}/bin/get-corr.R \
        ${count_file} \
        ${method2} \
        ${count_type}_${method1}_${method2}.csv.gz
    """
}