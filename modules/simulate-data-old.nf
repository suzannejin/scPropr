/* Simulate data using scDesign2 */

process SIMULATE_DATA_BYDEPTH {

    container 'suzannejin/scpropr:simulate'
    tag "${dataset}_simulate+${size_factor}"
    storeDir "${params.outdir}/${dataset}/simulated/${dataset}_simulate+${size_factor}"
    // publishDir "${params.outdir}/${dataset}/simulated/${dataset}_simulate+${size_factor}", mode: params.publish_dir_mode

    input:
    tuple val(dataset),
          val(exp_sim),
          val(full),
          val(abs_rel),
          file(model),
          file(features),
          file(barcodes)
    each size_factor 

    output:
    tuple file("${dataset}_simulate+${size_factor}_full_absolute.csv.gz"),
          file("${dataset}_simulate+${size_factor}_full_features.csv"),
          file("${dataset}_simulate+${size_factor}_full_barcodes.csv"),
          file(".command.trace"),
          file(".command.sh")

    script:
    """
    simulate-data-old.R \
        $model \
        ${dataset}_simulate+${size_factor}_full_absolute.csv.gz \
        --size_factor $size_factor
    mv $features ${dataset}_simulate+${size_factor}_full_features.csv
    mv $barcodes ${dataset}_simulate+${size_factor}_full_barcodes.csv
    """

    stub:
    """
    echo simulate-data-old.R \
        $model \
        ${dataset}_simulate+${size_factor}_full_absolute.csv.gz \
        --size_factor $size_factor
    touch ${dataset}_simulate+${size_factor}_full_absolute.csv.gz
    mv $features ${dataset}_simulate+${size_factor}_full_features.csv
    mv $barcodes ${dataset}_simulate+${size_factor}_full_barcodes.csv
    """
}

process SIMULATE_DATA_BYSLOPE {

    container 'suzannejin/scpropr:simulate'
    tag "${dataset}_simulate+s${slope}+n${ndata}+c${cell_factor}"
    storeDir "${params.outdir}/${dataset}/simulated/${dataset}_simulate+s${slope}+n${ndata}+c${cell_factor}+byslope"
    // publishDir "${params.outdir}/${dataset}/simulated/${dataset}_simulate+s${slope}+n${ndata}+c${cell_factor}+byslope", mode: params.publish_dir_mode

    input:
    tuple val(dataset),
          val(exp_sim),
          val(full),
          val(abs_rel),
          file(model),
          file(features),
          file(barcodes)
    each slope
    each ndata 
    each cell_factor

    output:
    tuple file("${dataset}_simulate+s${slope}+n${ndata}+c${cell_factor}+byslope_full_absolute.csv.gz"),
          file("${dataset}_simulate+s${slope}+n${ndata}+c${cell_factor}+byslope_full_features.csv"),
          file("${dataset}_simulate+s${slope}+n${ndata}+c${cell_factor}+byslope_full_barcodes.csv"),
          file(".command.trace"),
          file(".command.sh")

    script:
    """
    simulate-data-old.R \
        $model \
        ${dataset}_simulate+s${slope}+n${ndata}+c${cell_factor}+byslope_full_absolute.csv.gz \
        --slope $slope \
        --ndata $ndata \
        --cell_factor $cell_factor \
        --byslope
    mv $features ${dataset}_simulate+s${slope}+n${ndata}+c${cell_factor}+byslope_full_features.csv
    mv $barcodes ${dataset}_simulate+s${slope}+n${ndata}+c${cell_factor}+byslope_full_barcodes.csv
    """

    stub:
    """
    echo simulate-data-old.R \
        $model \
        ${dataset}_simulate+s${slope}+n${ndata}+c${cell_factor}+byslope_full_absolute.csv.gz \
        --slope $slope \
        --ndata $ndata \
        --cell_factor $cell_factor \
        --byslope
    touch ${dataset}_simulate+s${slope}+n${ndata}+c${cell_factor}+byslope_full_absolute.csv.gz
    mv $features ${dataset}_simulate+s${slope}+n${ndata}+c${cell_factor}+byslope_full_features.csv
    mv $barcodes ${dataset}_simulate+s${slope}+n${ndata}+c${cell_factor}+byslope_full_barcodes.csv
    """
}

process SIMULATE_DATA_BYSTEP {

    container 'suzannejin/scpropr:simulate'
    tag "${dataset}_simulate+s${slope}+n${ndata}+c${cell_factor}"
    storeDir "${params.outdir}/${dataset}/simulated/${dataset}_simulate+s${slope}+n${ndata}+c${cell_factor}+bystep"
    // publishDir "${params.outdir}/${dataset}/simulated/${dataset}_simulate+s${slope}+n${ndata}+c${cell_factor}+bystep", mode: params.publish_dir_mode

    input:
    tuple val(dataset),
          val(exp_sim),
          val(full),
          val(abs_rel),
          file(model),
          file(features),
          file(barcodes)
    each slope
    each ndata 
    each cell_factor

    output:
    tuple file("${dataset}_simulate+s${slope}+n${ndata}+c${cell_factor}+bystep_full_absolute.csv.gz"),
          file("${dataset}_simulate+s${slope}+n${ndata}+c${cell_factor}+bystep_full_features.csv"),
          file("${dataset}_simulate+s${slope}+n${ndata}+c${cell_factor}+bystep_full_barcodes.csv"),
          file(".command.trace"),
          file(".command.sh")

    script:
    """
    simulate-data-old.R \
        $model \
        ${dataset}_simulate+s${slope}+n${ndata}+c${cell_factor}+bystep_full_absolute.csv.gz \
        --slope $slope \
        --ndata $ndata \
        --cell_factor $cell_factor \
        --bystep
    mv $features ${dataset}_simulate+s${slope}+n${ndata}+c${cell_factor}+bystep_full_features.csv
    mv $barcodes ${dataset}_simulate+s${slope}+n${ndata}+c${cell_factor}+bystep_full_barcodes.csv
    """

    stub:
    """
    echo simulate-data-old.R \
        $model \
        ${dataset}_simulate+s${slope}+n${ndata}+c${cell_factor}+bystep_full_absolute.csv.gz \
        --slope $slope \
        --ndata $ndata \
        --cell_factor $cell_factor \
        --bystep
    touch ${dataset}_simulate+s${slope}+n${ndata}+c${cell_factor}+bystep_full_absolute.csv.gz
    mv $features ${dataset}_simulate+s${slope}+n${ndata}+c${cell_factor}+bystep_full_features.csv
    mv $barcodes ${dataset}_simulate+s${slope}+n${ndata}+c${cell_factor}+bystep_full_barcodes.csv
    """
}