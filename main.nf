#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


////////////////////////////////////////////////////
/* --          VALIDATE INPUTS                 -- */
////////////////////////////////////////////////////

Channel
    .fromPath(params.count_file)
    .map { it -> [ it.getParent().baseName, it.baseName.tokenize('.')[0], it ] }  
    .combine( Channel.fromList(params.zerohandling) )
    .set { ch_input }   // cell_type, count_type, count_file, zerohandling


////////////////////////////////////////////////////
/* --          IMPORT LOCAL MODULES            -- */
////////////////////////////////////////////////////

include { CLR  } from "${baseDir}/modules/propr.nf"
include { CORR } from "${baseDir}/modules/propr.nf"


////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////

workflow {
    /* 1st step: Compute CLR */
    CLR (ch_input)
    CLR.out.ch_clr
        .combine( Channel.fromList(params.methods) )
        .filter{ !(it[3] == "none" && it[4] == "rho") }
        .set { ch_to_corr }

    /* 2nd step: Compute association coefficients */
    CORR(ch_to_corr)
}