
/* keep specified cells - for simulation */

process FILTER_GENES {
    
    label 'process_low_short'
    container 'suzannejin/scpropr:propr'
    tag "$dataset"

    input:
    tuple val(dataset),
          val(exp_sim),
          val(full),
          val(abs_rel),
          path(count),
          path(features),
          path(barcodes)

    output:
    tuple val(dataset),
          val(exp_sim),
          val("filtered${params.gene_dropout_threshold}"),
          val(abs_rel),
          path("${dataset}_${exp_sim}_filtered${params.gene_dropout_threshold}_${abs_rel}.csv.gz"),
          path("${dataset}_${exp_sim}_filtered${params.gene_dropout_threshold}_features.csv"),
          path("${dataset}_${exp_sim}_filtered${params.gene_dropout_threshold}_barcodes.csv")
    
    script:
    def threshold = params.gene_dropout_threshold
    """
    select-data.R \
        -i $count \
        -f $features \
        -b $barcodes \
        -o ${dataset}_${exp_sim}_filtered${threshold}_${abs_rel}.csv.gz \
        -o2 ${dataset}_${exp_sim}_filtered${threshold}_features.csv \
        -o3 ${dataset}_${exp_sim}_filtered${threshold}_barcodes.csv \
        --filter_gene $threhsold
    """

    stub:
    def threshold = params.gene_dropout_threshold
    """
    echo select-data.R \
        -i $count \
        -f $features \
        -b $barcodes \
        -o ${dataset}_${exp_sim}_filtered${threshold}_${abs_rel}.csv.gz \
        -o2 ${dataset}_${exp_sim}_filtered${threshold}_features.csv \
        -o3 ${dataset}_${exp_sim}_filtered${threshold}_barcodes.csv \
        --filter_gene $threhsold
    touch ${dataset}_${exp_sim}_filtered${threshold}_${abs_rel}.csv.gz
    touch ${dataset}_${exp_sim}_filtered${threshold}_features.csv
    touch ${dataset}_${exp_sim}_filtered${threshold}_barcodes.csv
    """

}


process SELECT_BARCODES_SIM {

    label 'process_low_short'
    container 'suzannejin/scpropr:propr'
    tag "$dataset"

    input:
    tuple val(dataset),
          val(exp_sim),
          val(full),
          val(abs_rel),
          path(count),
          path(features),
          path(barcodes),
          path(barcodes_sim)

    output:
    tuple val(dataset),
          val('simulation'),
          val(full),
          val(abs_rel),
          path("${dataset}_simulation_${full}_${abs_rel}.csv.gz"),
          path("${dataset}_simulation_${full}_features.csv"),
          path("${dataset}_simulation_${full}_barcodes.csv")
    
    script:
    """
    select-data.R \
        -i $count \
        -f $features \
        -b $barcodes \
        -o ${dataset}_simulation_${full}_${abs_rel}.csv.gz \
        -o2 ${dataset}_simulation_${full}_features.csv \
        -o3 ${dataset}_simulation_${full}_barcodes.csv \
        --select_cell $barcodes_sim
    """

    stub:
    """
    echo select-data.R \
        -i $count \
        -f $features \
        -b $barcodes \
        -o ${dataset}_simulation_${full}_${abs_rel}.csv.gz \
        -o2 ${dataset}_simulation_${full}_features.csv \
        -o3 ${dataset}_simulation_${full}_barcodes.csv \
        --select_cell $barcodes_sim
    touch ${dataset}_simulation_${full}_${abs_rel}.csv.gz 
    touch ${dataset}_simulation_${full}_features.csv 
    touch ${dataset}_simulation_${full}_barcodes.csv 
    """
}


process SELECT_NOZERO_GENES_LONG {

    label 'process_low_short'
    container 'suzannejin/scpropr:propr'
    tag "$dataset-$exp_sim"

    input:
    tuple val(dataset),
          val(exp_sim),
          val(full),
          val(abs_rel),
          path(count),
          path(features),
          path(barcodes)

    output:
    tuple val(dataset),
          val(exp_sim),
          val('nozerolong'),
          val(abs_rel),
          path("${dataset}_${exp_sim}_nozerolong_${abs_rel}.csv.gz"),
          path("${dataset}_${exp_sim}_nozerolong_features.csv"),
          path("${dataset}_${exp_sim}_nozerolong_barcodes.csv")
    
    script:
    """
    select-data.R \
        -i $count \
        -f $features \
        -b $barcodes \
        -o ${dataset}_${exp_sim}_nozerolong_${abs_rel}.csv.gz \
        -o2 ${dataset}_${exp_sim}_nozerolong_features.csv \
        -o3 ${dataset}_${exp_sim}_nozerolong_barcodes.csv \
        --filter_gene 0
    """

    stub:
    """
    echo select-data.R \
        -i $count \
        -f $features \
        -b $barcodes \
        -o ${dataset}_${exp_sim}_nozerolong_${abs_rel}.csv.gz \
        -o2 ${dataset}_${exp_sim}_nozerolong_features.csv \
        -o3 ${dataset}_${exp_sim}_nozerolong_barcodes.csv \
        --filter_gene 0
    touch ${dataset}_${exp_sim}_nozerolong_${abs_rel}.csv.gz
    touch ${dataset}_${exp_sim}_nozerolong_features.csv
    touch ${dataset}_${exp_sim}_nozerolong_barcodes.csv
    """
}

process SELECT_NOZERO_GENES_SQUARE {

    label 'process_low_short'
    container 'suzannejin/scpropr:propr'
    tag "$dataset-$exp_sim"

    input:
    tuple val(dataset),
          val(exp_sim),
          val(full),
          val(abs_rel),
          path(count),
          path(features),
          path(barcodes)

    output:
    tuple val(dataset),
          val(exp_sim),
          val('nozerosquare'),
          val(abs_rel),
          path("${dataset}_${exp_sim}_nozerosquare_${abs_rel}.csv.gz"),
          path("${dataset}_${exp_sim}_nozerosquare_features.csv"),
          path("${dataset}_${exp_sim}_nozerosquare_barcodes.csv")
    
    script:
    """
    get-reduced-dataset.R \
        -i $count \
        -f $features \
        -b $barcodes \
        -o ${dataset}_${exp_sim}_nozerosquare_${abs_rel}.csv.gz \
        -o2 ${dataset}_${exp_sim}_nozerosquare_features.csv \
        -o3 ${dataset}_${exp_sim}_nozerosquare_barcodes.csv 
    """

    stub:
    """
    echo get-reduced-dataset.R \
        -i $count \
        -f $features \
        -b $barcodes \
        -o ${dataset}_${exp_sim}_nozerosquare_${abs_rel}.csv.gz \
        -o2 ${dataset}_${exp_sim}_nozerosquare_features.csv \
        -o3 ${dataset}_${exp_sim}_nozerosquare_barcodes.csv 
    touch ${dataset}_${exp_sim}_nozerosquare_${abs_rel}.csv.gz
    touch ${dataset}_${exp_sim}_nozerosquare_features.csv
    touch ${dataset}_${exp_sim}_nozerosquare_barcodes.csv
    """
}


process SELECT_FIXED_GENES {

    label 'process_low_short'
    container 'suzannejin/scpropr:propr'
    tag "$dataset-$exp_sim"

    input:
    tuple val(dataset),
          val(exp_sim),
          val(full),
          val(abs_rel),
          path(count),
          path(features),
          path(barcodes),
          path(fixed_features)

    output:
    tuple val(dataset),
          val(exp_sim),
          val('fixed'),
          val(abs_rel),
          path("${dataset}_${exp_sim}_fixed_${abs_rel}.csv.gz"),
          path("${dataset}_${exp_sim}_fixed_features.csv"),
          path("${dataset}_${exp_sim}_fixed_barcodes.csv")
    
    script:
    """
    select-data.R \
        -i $count \
        -f $features \
        -b $barcodes \
        -o ${dataset}_${exp_sim}_fixed_${abs_rel}.csv.gz \
        -o2 ${dataset}_${exp_sim}_fixed_features.csv \
        -o3 ${dataset}_${exp_sim}_fixed_barcodes.csv \
        --select_gene $fixed_features
    """

    stub:
    """
    echo select-data.R \
        -i $count \
        -f $features \
        -b $barcodes \
        -o ${dataset}_${exp_sim}_fixed_${abs_rel}.csv.gz \
        -o2 ${dataset}_${exp_sim}_fixed_features.csv \
        -o3 ${dataset}_${exp_sim}_fixed_barcodes.csv \
        --select_gene $fixed_features
    touch ${dataset}_${exp_sim}_fixed_${abs_rel}.csv.gz
    touch ${dataset}_${exp_sim}_fixed_features.csv
    touch ${dataset}_${exp_sim}_fixed_barcodes.csv
    """
}
