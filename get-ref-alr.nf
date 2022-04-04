#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

input_count_file =  "/users/cn/sjin/projects/proportionality/run_scPropr/data/mouseStemCells/mouseStemCells_experimental_full_absolute.csv.gz"

// input count data
Channel
    .fromPath(input_count_file, checkIfExists:true)
    .map { it -> [ it.getSimpleName().tokenize('_')[0], it.getSimpleName().tokenize('_')[1], it.getSimpleName().tokenize('_')[2], it.getSimpleName().tokenize('_')[3], it ] }  
    .set { ch_count }  
// features
Channel
    .fromPath(params.input_features_file, checkIfExists:true)
    .map { it -> [ it.getSimpleName().tokenize('_')[0], it.getSimpleName().tokenize('_')[1], it ] }   
    .set { ch_features }
// barcodes
Channel
    .fromPath(params.input_barcodes_file, checkIfExists:true)
    .map { it -> [ it.getSimpleName().tokenize('_')[0], it.getSimpleName().tokenize('_')[1], it ] } 
    .set { ch_barcodes }
// main input channel
ch_count
    .combine( ch_features, by:[0,1] )
    .combine( ch_barcodes, by:[0,1] )
    .set{ ch_input }  // dataset, experimental, full, absolute, count, features, barcodes
procrustes_filter = 0  // % of zeros in a gene. Procrustes analysis will only be performed for genes below the threshold

// import module
include { COMPUTE_LRA            } from "${launchDir}/modules/compute-lra.nf"
include { PROCRUSTES_ANALYSIS    } from "${launchDir}/modules/procrustes-analysis.nf"
include { REPLACE_ZERO           } from "${launchDir}/modules/replace-zero.nf"
include { SELECT_GENE_PROCRUSTES } from "${launchDir}/modules/select-gene-procrustes.nf"


workflow {

    /* replace zero */
    REPLACE_ZERO(
        ch_input, 
        params.method_replace_zero.tokenize(',')
    )
    REPLACE_ZERO.out
        .map{ it -> [
            it[0].getSimpleName().tokenize('_')[0],
            it[0].getSimpleName().tokenize('_')[1],
            it[0].getSimpleName().tokenize('_')[2],
            it[0].getSimpleName().tokenize('_')[3],
            it[0].getSimpleName().tokenize('_')[4],
            it[0],
            it[1],
            it[2]
        ] } 
        .set{ ch_imputed }   // dataset, experimental, full, absolute, method_replace_zero, count, features, barcodes

    /* compute exact log ratio geometry */
    COMPUTE_LRA(ch_imputed)
    COMPUTE_LRA.out
        .map{ it -> [
            it[3].getSimpleName().tokenize('_')[0],
            it[3].getSimpleName().tokenize('_')[1],
            it[3].getSimpleName().tokenize('_')[2],
            it[3].getSimpleName().tokenize('_')[3],
            it[3].getSimpleName().tokenize('_')[4],
            it[0],
            it[1],
            it[2],
            it[3]
        ] }
        .set{ ch_lra }  // dataset, experimental, full, absolute, method_replace_zero, count, features, barcodes, lra

    /* select genes for procrustes analysis, given the limiter computational time and resources*/
    SELECT_GENE_PROCRUSTES(
        ch_input,
        procrustes_filter
    )
    ch_lra
        .combine(SELECT_GENE_PROCRUSTES.out, by:[0,1,2,3])
        .transpose()
        .set{ ch_input2procrustes }
    
    /* perform procrustes analysis for each gene */
    PROCRUSTES_ANALYSIS(ch_input2procrustes)
}