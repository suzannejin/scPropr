process CLR {

    tag "${cell_type}-${zerohandling}-${count_type}"
    publishDir "${params.outdir}/${cell_type}/clr", 
        mode: params.publish_dir_mode,
        saveAs: { filename -> if (filename.indexOf("_none.csv.gz") > 0) null else filename }

    input:
    tuple val(cell_type), \
          val(count_type), \
          file(count_file), \
          val(zerohandling)

    output:
    tuple val(cell_type), \
          val(count_type), \
          file("${count_type}_${zerohandling}.csv.gz"), \
          val(zerohandling), \
          emit: ch_clr
    file("${count_type}_zeroreplace_${zerohandling}.csv.gz")

    script:
    if (zerohandling == "none")
    """
    cp ${count_file} ${count_type}_zeroreplace_${zerohandling}.csv.gz
    mv ${count_file} ${count_type}_${zerohandling}.csv.gz
    """
    else
    """
    Rscript ${baseDir}/bin/clr-data.R \
        ${count_file} \
        ${zerohandling} \
        ${count_type}_${zerohandling}.csv.gz \
        ${count_type}_zeroreplace_${zerohandling}.csv.gz
    """
}

process CORR {

    tag "${cell_type}-${zerohandling}-${method}-${count_type}"
    publishDir "${params.outdir}/${cell_type}/corr", mode: params.publish_dir_mode

    input:
    tuple val(cell_type), \
          val(count_type), \
          file(count_file), \
          val(zerohandling), \
          val(method)

    output:
    tuple val(cell_type), \
          file("${count_type}_${zerohandling}_${method}.csv.gz"), \
          emit: ch_corr

    script:
    """
    Rscript ${baseDir}/bin/get-corr.R \
        ${count_file} \
        ${method} \
        ${count_type}_${zerohandling}_${method}.csv.gz
    """
}