process MODEL {

    container 'suzannejin/scpropr:simulate'

    tag "${cell_type}"
    publishDir params.use_phase ? "${params.outdir}/${cell_type}/model_phase" : "${params.outdir}/${cell_type}/model", 
               mode: params.publish_dir_mode

    input:
    tuple val(cell_type), \
          val(absolute), \
          val(experimental), \
          file(count_file), \
          file(features_file), \
          file(phase_file)

    output:
    tuple val(cell_type), \
          val(absolute), \
          val(experimental), \
          file("model.rds"), \
          emit: ch_model

    script:
    def phase = params.use_phase ? phase_file : "NA" 

    """
    Rscript ${baseDir}/bin/get-model.R \
        ${count_file} \
        ${features_file} \
        ${phase} \
        ${cell_type} \
        model.rds
    """
}

process SIMULATE {

    container 'suzannejin/scpropr:simulate'

    tag "${cell_type}_simulate+${size_factor}"
    publishDir params.use_phase ? "${params.outdir}/${cell_type}/simulate_phase" : "${params.outdir}/${cell_type}/simulate", 
               mode: params.publish_dir_mode

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
        1 \
        simulate+${size_factor}.csv.gz
    """
}

process SIMULATE2 {

    container 'suzannejin/scpropr:simulate'

    tag "${cell_type}_merge+${cell_factor}"
    publishDir params.use_phase ? "${params.outdir}/${cell_type}/simulate_phase" : "${params.outdir}/${cell_type}/simulate", 
               mode: params.publish_dir_mode

    input:
    tuple val(cell_type), \
          val(absolute), \
          val(experimental), \
          file(model), \
          val(cell_factor)

    output:
    tuple val(cell_type), \
          val(absolute), \
          file("merged*.csv.gz"), \
          emit: ch_simulated2

    script:
    size_factors = params.size_factors
    """
    vals=\$( echo ${size_factors} | tr -d "[" | tr -d "]" | tr "," " " )
    out=merged+\$(echo \$vals | tr " " "+").csv.gz
    for i in \$vals; do
    Rscript ${baseDir}/bin/simulate-data.R \
        ${model} \
        \$i \
        ${cell_factor} \
        simulate+\$i.csv.gz
    cat simulate+\$i.csv.gz >> \$out
    done
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