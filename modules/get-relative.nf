
/* get relative data */

process GET_RELATIVE {

    label 'process_low_short'
    container 'suzannejin/scpropr:propr'
    tag "$dataset-$exp_sim-$full"

    input:
    tuple val(dataset),
          val(exp_sim),
          val(full),
          val(abs_rel),
          file(count),
          file(features),
          file(barcodes)

    output:
    tuple val(dataset),
          val(exp_sim),
          val(full),
          val('relative'),
          path("${dataset}_${exp_sim}_${full}_relative.csv.gz"),
          path(features),
          path(barcodes)

    script:
    """
    get-relative.R \
        $count \
        ${dataset}_${exp_sim}_${full}_relative.csv.gz
    """

    stub:
    """
    echo get-relative.R \
        $count \
        ${dataset}_${exp_sim}_${full}_relative.csv.gz
    touch ${dataset}_${exp_sim}_${full}_relative.csv.gz
    """
}

process GET_RELATIVE_FOR_PLOT {

    label 'process_low_short'
    container 'suzannejin/scpropr:propr'
    tag "$dataset-$exp_sim-$full"

    input:
    tuple val(dataset),
          val(exp_sim),
          val(full),
          val(abs_rel),
          file(count),
          file(features),
          file(barcodes)

    output:
    tuple val(dataset),
          val(exp_sim),
          val(full),
          val('relative'),
          path("${dataset}_${exp_sim}_${full}_relative.csv.gz"),
          path(features),
          path(barcodes)

    script:
    """
    get-relative.R \
        $count \
        ${dataset}_${exp_sim}_${full}_relative.csv.gz
    """

    stub:
    """
    echo get-relative.R \
        $count \
        ${dataset}_${exp_sim}_${full}_relative.csv.gz
    touch ${dataset}_${exp_sim}_${full}_relative.csv.gz
    """
}