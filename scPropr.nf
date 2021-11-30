#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


////////////////////////////////////////////////////
/* --          VALIDATE INPUTS                 -- */
////////////////////////////////////////////////////

Channel
    .fromPath(params.count_file, checkIfExists:true)
    .map { it -> [ it.getParent().getParent().baseName, it.baseName.tokenize('.')[0], it ] }  
    .set { ch_count }
Channel
    .fromPath(params.nozero_genes, checkIfExists:true)
    .map { it -> [ it.getParent().getParent().baseName, it ] }
    .set { ch_nozero }
ch_count
    .combine( ch_nozero, by:0 )
    .combine( Channel.fromList(params.methods_replace_zero.tokenize(",")) )
    .combine( Channel.fromList(params.methods_transf_data.tokenize(",")) )
    .filter{ (it[4] == "none" && it[5] == "tmm") || (it[4] != "none" && it[5] in ["log2","clr"]) }
    .set { ch_input }   // cell_type, count_type, count_file, nozero_file, zero_replacement_method, data_transf_method


////////////////////////////////////////////////////
/* --          IMPORT LOCAL MODULES            -- */
////////////////////////////////////////////////////

include { TRANSF } from "${baseDir}/modules/propr.nf"
include { CORR   } from "${baseDir}/modules/propr.nf"


////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////

workflow {
    /* 1st step: Transform or normalize data */
    TRANSF(ch_input)
    TRANSF.out.ch_transf
          .combine( Channel.fromList(params.methods_corr.tokenize(",")) )
          .filter{ (it[4] == "log2") || (it[4] == "clr") || (it[5] == "cor") }
          .set { ch_to_corr }

    /* 2nd step: Compute association coefficients */
    CORR(ch_to_corr)
}