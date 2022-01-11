process MODEL {

    container 'suzannejin/scpropr:simulate'

    tag "${cell_type}"
    publishDir "${params.outdir}/${cell_type}/model", mode: params.publish_dir_mode

    input:
    tuple val(cell_type), \
          val(absolute), \
          val(experimental), \
          file(count_file), \
          file(features_file)

    output:
    tuple val(cell_type), \
          val(absolute), \
          val(experimental), \
          file("model.rds"), \
          emit: ch_model

    script:
    """
    Rscript ${baseDir}/bin/get-model.R \
        ${count_file} \
        ${features_file} \
        ${cell_type} \
        model.rds
    """
}

process SIMULATE {

    container 'suzannejin/scpropr:simulate'

    tag "${cell_type}_${size_factor}"
    publishDir "${params.outdir}/${cell_type}/simulate", mode: params.publish_dir_mode

    input:
    tuple val(cell_type), \
          val(absolute), \
          val(experimental), \
          file(model), \
          val(size_factor)

    output:
    tuple val(cell_type), \
          val(absolute), \
          val("simulate+${size_factor}"), \
          file("simulate+${size_factor}.csv.gz"), \
          emit: ch_simulated

    script:
    """
    Rscript ${baseDir}/bin/simulate-data.R \
        ${model} \
        ${size_factor} \
        simulate+${size_factor}.csv.gz
    """
}

process CHECK_DROPOUT {

    container 'suzannejin/scpropr:simulate'

    tag "${cell_type}_${simulated_size_factor}"
    publishDir "${params.outdir}/${cell_type}/simulate_dropout", mode: params.publish_dir_mode

    input:
    tuple val(cell_type), \
          val(absolute), \
          val(simulated_size_factor), \
          file(simulated_data), \
          file(features_pos)

    output:
    tuple val(cell_type), \
          file("${simulated_size_factor}_dropout.csv"), \
          file("${simulated_size_factor}_nonzerogenes_dropout.csv"), \
          emit: ch_dropout

    script:
    """
    Rscript ${baseDir}/bin/check-dropout.R \
        ${simulated_data} \
        NA \
        ${simulated_size_factor}_dropout.csv

    Rscript ${baseDir}/bin/check-dropout.R \
        ${simulated_data} \
        ${features_pos} \
        ${simulated_size_factor}_nonzerogenes_dropout.csv
    """
}