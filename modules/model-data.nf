
/* model single cell data using scDesign2 */

process MODEL_DATA {

    container 'suzannejin/scpropr:simulate'

    tag "${cell_type}"
    storeDir "${params.outdir}/${dataset}/model"

    input:
    val dataset
    path count
    path features
    path barcodes_sim

    output:
    val dataset, emit: ch_dataset
    path "${dataset}_model.rds", emit: ch_model

    script:
    """
    model-data.R \
        $count \
        $features \
        $barcodes_sim \
        $dataset \
        ${dataset}_model.rds
    """
}