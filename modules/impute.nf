process IMPUTE {
    label '${dataset}'

    input:
    tuple val(dataset),
          path(count),
          path(features),
          path(barcodes),
          val(imputation)
    
    output:
    tuple val(dataset),
          path(""),
          path(""),
          path(""),
          val(imputation)

    script:
    """
    impute_zeros.R \
        -i $count \
        -o

    """
}