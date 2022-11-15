
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
          val('filtered'),
          val(abs_rel),
          path("${dataset}_${exp_sim}_filtered_${abs_rel}.csv.gz"),
          path("${dataset}_${exp_sim}_filtered_features.csv"),
          path("${dataset}_${exp_sim}_filtered_barcodes.csv")
    
    script:
    """
    select-data.R \
        -i $count \
        -f $features \
        -b $barcodes \
        -o ${dataset}_${exp_sim}_filtered_${abs_rel}.csv.gz \
        -o2 ${dataset}_${exp_sim}_filtered_features.csv \
        -o3 ${dataset}_${exp_sim}_filtered_barcodes.csv \
        --filter_gene ${params.gene_dropout_threshold}
    """

    stub:
    """
    echo select-data.R \
        -i $count \
        -f $features \
        -b $barcodes \
        -o ${dataset}_${exp_sim}_filtered_${abs_rel}.csv.gz \
        -o2 ${dataset}_${exp_sim}_filtered_features.csv \
        -o3 ${dataset}_${exp_sim}_filtered_barcodes.csv \
        --filter_gene ${params.gene_dropout_threshold}
    touch ${dataset}_${exp_sim}_filtered_${abs_rel}.csv.gz
    touch ${dataset}_${exp_sim}_filtered_features.csv
    touch ${dataset}_${exp_sim}_filtered_barcodes.csv
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


/* keep only the genes with no zero components - for analysis */

process SELECT_NOZERO_GENES {

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
          val('nozero'),
          val(abs_rel),
          path("${dataset}_${exp_sim}_nozero_${abs_rel}.csv.gz"),
          path("${dataset}_${exp_sim}_nozero_features.csv"),
          path("${dataset}_${exp_sim}_nozero_barcodes.csv")
    
    script:
    """
    select-data.R \
        -i $count \
        -f $features \
        -b $barcodes \
        -o ${dataset}_${exp_sim}_nozero_${abs_rel}.csv.gz \
        -o2 ${dataset}_${exp_sim}_nozero_features.csv \
        -o3 ${dataset}_${exp_sim}_nozero_barcodes.csv \
        --filter_gene 0
    """

    stub:
    """
    echo select-data.R \
        -i $count \
        -f $features \
        -b $barcodes \
        -o ${dataset}_${exp_sim}_nozero_${abs_rel}.csv.gz \
        -o2 ${dataset}_${exp_sim}_nozero_features.csv \
        -o3 ${dataset}_${exp_sim}_nozero_barcodes.csv \
        --filter_gene 0
    touch ${dataset}_${exp_sim}_nozero_${abs_rel}.csv.gz
    touch ${dataset}_${exp_sim}_nozero_features.csv
    touch ${dataset}_${exp_sim}_nozero_barcodes.csv
    """
}

process GET_REDUCED_DATASET {

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
          val('nozero'),
          val(abs_rel),
          path("${dataset}_${exp_sim}_nozero_${abs_rel}.csv.gz"),
          path("${dataset}_${exp_sim}_nozero_features.csv"),
          path("${dataset}_${exp_sim}_nozero_barcodes.csv")
    
    script:
    """
    get-reduced-dataset.R \
        -i $count \
        -f $features \
        -b $barcodes \
        -o ${dataset}_${exp_sim}_nozero_${abs_rel}.csv.gz \
        -o2 ${dataset}_${exp_sim}_nozero_features.csv \
        -o3 ${dataset}_${exp_sim}_nozero_barcodes.csv \
        --ncell 1100 \
        --ngene 1100
    """

    stub:
    """
    echo get-reduced-dataset.R \
        -i $count \
        -f $features \
        -b $barcodes \
        -o ${dataset}_${exp_sim}_nozero_${abs_rel}.csv.gz \
        -o2 ${dataset}_${exp_sim}_nozero_features.csv \
        -o3 ${dataset}_${exp_sim}_nozero_barcodes.csv \
        --ncell 1100 \
        --ngene 1100
    touch ${dataset}_${exp_sim}_nozero_${abs_rel}.csv.gz
    touch ${dataset}_${exp_sim}_nozero_features.csv
    touch ${dataset}_${exp_sim}_nozero_barcodes.csv
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


process EXP_NOZERO_GENES_2_SELECT {

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
          path(features_experimental_nozero)

    output:
    tuple val(dataset),
          val(exp_sim),
          val('nozero'),
          val(abs_rel),
          path("${dataset}_${exp_sim}_nozero_${abs_rel}.csv.gz"),
          path("${dataset}_${exp_sim}_nozero_features.csv"),
          path("${dataset}_${exp_sim}_nozero_barcodes.csv")
    
    script:
    """
    select-data.R \
        -i $count \
        -f $features \
        -b $barcodes \
        -o ${dataset}_${exp_sim}_nozero_${abs_rel}.csv.gz \
        -o2 ${dataset}_${exp_sim}_nozero_features.csv \
        -o3 ${dataset}_${exp_sim}_nozero_barcodes.csv \
        --select_gene $features_experimental_nozero
    """

    stub:
    """
    echo select-data.R \
        -i $count \
        -f $features \
        -b $barcodes \
        -o ${dataset}_${exp_sim}_nozero_${abs_rel}.csv.gz \
        -o2 ${dataset}_${exp_sim}_nozero_features.csv \
        -o3 ${dataset}_${exp_sim}_nozero_barcodes.csv \
        --select_gene $features_experimental_nozero
    touch ${dataset}_${exp_sim}_nozero_${abs_rel}.csv.gz
    touch ${dataset}_${exp_sim}_nozero_features.csv
    touch ${dataset}_${exp_sim}_nozero_barcodes.csv
    """
}



