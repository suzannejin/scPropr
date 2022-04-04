
/* model single cell data using scDesign2 */

process MODEL_DATA {

    label 'process_high_long'
    container 'suzannejin/scpropr:simulate'
    tag "${dataset}"
    storeDir "${params.outdir}/${dataset}/model/${exp_sim}"
    // publishDir "${params.outdir}/${dataset}/model", mode: params.publish_dir_mode
    // TODO be able to model based on different homogeneous populations
    // for this I need to change the entire pipeline, to add another variable 'population'

    input:
    tuple val(dataset),
          val(exp_sim),
          val(full),
          val(abs_rel),
          file(count),
          file(features),
          file(barcodes)

    output:
    tuple file("${dataset}_${exp_sim}_model.rds"), 
          file(features), 
          file(barcodes),
          file(".command.trace"),
          file(".command.sh")

    script:
    """
    model-data.R \
        $count \
        $features \
        $dataset \
        ${dataset}_${exp_sim}_model.rds
    """

    stub:
    """
    echo model-data.R \
        $count \
        $features \
        $dataset \
        ${dataset}_model.rds
    touch ${dataset}_${exp_sim}_model.rds
    """
}