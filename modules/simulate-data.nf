/* Simulate data using scDesign2 */

process SIMULATE_DATA {

    container 'suzannejin/scpropr:simulate'
    tag "${dataset}_${exp_sim}+simulate+s${slope}+n${ndata}+c${cell_factor}"
    storeDir "${params.outdir}/${dataset}/simulated/${dataset}_${exp_sim}+simulate+s${slope}+n${ndata}+c${cell_factor}"

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
    tuple file("${dataset}_${exp_sim}+simulate+s${slope}+n${ndata}+c${cell_factor}_full_absolute.csv.gz"),
          file("${dataset}_${exp_sim}+simulate+s${slope}+n${ndata}+c${cell_factor}_full_features.csv"),
          file("${dataset}_${exp_sim}+simulate+s${slope}+n${ndata}+c${cell_factor}_full_barcodes.csv"),
          file(".command.trace"),
          file(".command.sh"),
          file(".command.log")

    script:
    """
    simulate-data.R \
        $model \
        ${dataset}_${exp_sim}+simulate+s${slope}+n${ndata}+c${cell_factor}_full_absolute.csv.gz \
        --slope $slope \
        --ndata $ndata \
        --cell_factor $cell_factor 
    mv $features ${dataset}_${exp_sim}+simulate+s${slope}+n${ndata}+c${cell_factor}_full_features.csv
    mv $barcodes ${dataset}_${exp_sim}+simulate+s${slope}+n${ndata}+c${cell_factor}_full_barcodes.csv
    """

    stub:
    """
    echo simulate-data.R \
        $model \
        ${dataset}_${exp_sim}+simulate+s${slope}+n${ndata}+c${cell_factor}_full_absolute.csv.gz \
        --slope $slope \
        --ndata $ndata \
        --cell_factor $cell_factor 
    touch ${dataset}_${exp_sim}+simulate+s${slope}+n${ndata}+c${cell_factor}_full_absolute.csv.gz
    mv $features ${dataset}_${exp_sim}+simulate+s${slope}+n${ndata}+c${cell_factor}_full_features.csv
    mv $barcodes ${dataset}_${exp_sim}+simulate+s${slope}+n${ndata}+c${cell_factor}_full_barcodes.csv
    """
}