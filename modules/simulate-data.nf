/* Simulate data using scDesign2 */

process SIMULATE_DATA_BYDEPTH {

    container 'suzannejin/scpropr:simulate'

    tag "${dataset}_simulate+${size_factor}"
    storeDir "${params.outdir}/${dataset}/data"

    input:
    val dataset
    path count
    path model
    val size_factor

    output:
    val dataset, emit: ch_dataset
    val "simulate+${size_factor}", emit: ch_type1
    path "${dataset}_simulate+${size_factor}_absolute.csv.gz", emit: ch_simulated

    script:
    """
    simulate-data.R \
        $model \
        ${dataset}_simulate+${size_factor}_absolute.csv.gz \
        --size_factor $size_factor
    """
}

process SIMULATE_DATA_BYSLOPE {

    container 'suzannejin/scpropr:simulate'

    tag "${dataset}_simulate+s${slope}n${ndata}"
    storeDir "${params.outdir}/${dataset}/data"

    input:
    val dataset
    path count
    path model
    val slope
    val ndata

    output:
    val dataset, emit: ch_dataset
    val "simulate+s${slope}n${ndata}byslope", emit: ch_type1
    path "${dataset}_simulate+s${slope}n${nadata}byslope_absolute.csv.gz", emit: ch_simulated

    script:
    """
    simulate-data.R \
        $model \
        ${dataset}_simulate+s${slope}n${ndata}byslope_absolute.csv.gz \
        --slope $slope \
        --ndata $ndata \
        --bymax
    """
}

process SIMULATE_DATA_BYSTEP {

    container 'suzannejin/scpropr:simulate'

    tag "${dataset}_simulate+s${slope}n${ndata}"
    storeDir "${params.outdir}/${dataset}/data"

    input:
    val dataset
    path count
    path model
    val slope
    val ndata

    output:
    val dataset, emit: ch_dataset
    val "simulate+s${slope}n${ndata}bystep", emit: ch_type1
    path "${dataset}_simulate+s${slope}n${nadata}bystep_absolute.csv.gz", emit: ch_simulated

    script:
    """
    simulate-data.R \
        $model \
        ${dataset}_simulate+s${slope}n${ndata}bystep_absolute.csv.gz \
        --slope $slope \
        --ndata $ndata \
        --bystep
    """
}