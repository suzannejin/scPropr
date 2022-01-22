#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


////////////////////////////////////////////////////
/* --          VALIDATE INPUTS                 -- */
////////////////////////////////////////////////////

// input count data
Channel
    .fromPath(params.count_file, checkIfExists:true)
    .map { it -> [ it.getParent().getParent().baseName, "absolute", "experimental", it ] }  
    .set { ch_count }
// cell type
ch_count
    .map { it -> [it[0]] }
    .set{ ch_cell }
// feature (gene) names
Channel
    .fromPath(params.features, checkIfExists:true)
    .map { it -> [ it.getParent().getParent().baseName, it ] }
    .set { ch_features }
// non-zero genes positions
Channel
    .fromPath(params.nozero_genes, checkIfExists:true)
    .map { it -> [ it.getParent().getParent().baseName, it ] }
    .set { ch_nozero }
// cells from a given cell cycle phase
Channel
    .fromPath(params.phase_pos, checkIfExists:true)
    .map { it -> [ it.getParent().getParent().baseName, it ] }
    .set { ch_phase }

// input for modeling step
ch_count
    .combine( ch_features, by:0 )
    .combine( ch_phase, by:0 )
    .set { ch_input2model }  // cell_type, count_type1 (absolute/relative), count_type2 (experimental/simulated), count_file, features_file, phase_file


////////////////////////////////////////////////////
/* --          IMPORT LOCAL MODULES            -- */
////////////////////////////////////////////////////

include { MODEL    } from "${baseDir}/modules/simulate.nf"
include { SIMULATE } from "${baseDir}/modules/simulate.nf"
include { SIMULATE2} from "${baseDir}/modules/simulate.nf"
include { RELATIVE } from "${baseDir}/modules/propr.nf"
include { TRANSF   } from "${baseDir}/modules/propr.nf"
include { CORR     } from "${baseDir}/modules/propr.nf"


////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////

workflow {

    if (params.model) {

    /* fit model */
    MODEL( ch_input2model )
    MODEL.out.ch_model
        .set{ ch_model2 }

    } else {

    Channel
        .fromPath(params.use_phase ? "${params.outdir}/*/model_phase/model.rds" : "${params.outdir}/*/model/model.rds")
        .map{ it -> [ it.getParent().getParent().baseName, "absolute", "experimental", it ] }
        .combine( ch_cell, by:0 )
        .set{ ch_model2 }
    }

    if (params.simulate) {

    ch_model2
        .combine( Channel.fromList(params.size_factors) )
        .set{ ch_model2simulate }
    ch_model2
        .combine( Channel.from(2..params.ncomb_max) )
        .set{ ch_model2simulate2 }

    /* simulate data */
    SIMULATE(ch_model2simulate)
    SIMULATE2(ch_model2simulate2)
    SIMULATE2.out.ch_simulated2
        .map{ it -> [ it[0], it[1], it[2].baseName.tokenize('.')[0], it[2] ] }
        .mix(SIMULATE.out.ch_simulated)
        .mix(ch_count)
        .set{ ch_absolute }

    }else{

    Channel
        .fromPath(params.use_phase ? "${params.outdir}/*/simulate_phase/{simulate,merged}*.csv.gz" : "${params.outdir}/*/simulate/{simulate,merged}*.csv.gz")
        .map{ it -> [ it.getParent().getParent().baseName, "absolute", it.baseName.tokenize(".")[0], it ] }
        .combine( ch_cell, by:0 )
        .mix(ch_count)
        .set{ ch_absolute }
    }

    /* reclose data - get relative data */
    RELATIVE(ch_absolute)
    RELATIVE.out.ch_relative
        .mix(ch_absolute)
        .combine(ch_nozero, by:0)
        .combine(Channel.fromList(params.methods_replace_zero.tokenize(",")))
        .combine(Channel.fromList(params.methods_transf_data.tokenize(",")))
        .combine(Channel.fromList(params.ref_genes))
        .map{ it -> it[6] == "alr" ? it : [it[0..6], "NA"].flatten() }
        .unique()
        .set { ch_count2transf }   // cell_type, count_type1, count_type2, count_file, nozero_file, zero_replacement_method, data_transf_method, ref_gene

    /* transform or normalize data */
    TRANSF(ch_count2transf)
    TRANSF.out.ch_transf
          .combine( Channel.fromList(params.methods_corr.tokenize(",")) )
        //   .filter{ (it[5] in ["log2", "clr", "alr"] ) || (it[7] == "cor") }
          .set { ch_transf2corr }

    /* compute association coefficients */
    CORR(ch_transf2corr)
}