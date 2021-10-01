#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


////////////////////////////////////////////////////
/* --          VALIDATE INPUTS                 -- */
////////////////////////////////////////////////////

Channel
    .fromPath(params.count_file)
    .map { it -> [ it.getParent().baseName, it.baseName.tokenize('.')[0], it ] }  
    .combine( Channel.fromList(params.methods_transf) )
    .set { ch_input }   // cell_type, count_type, count_file, method_transf


////////////////////////////////////////////////////
/* --          IMPORT LOCAL MODULES            -- */
////////////////////////////////////////////////////

include { TRANSF } from "${baseDir}/modules/propr.nf"
include { CORR   } from "${baseDir}/modules/propr.nf"


////////////////////////////////////////////////////
/* --           RUN MAIN WORKFLOW              -- */
////////////////////////////////////////////////////

workflow {
    /* 1st step: Compute CLR */
    TRANSF(ch_input)
    TRANSF.out.ch_clr
          .combine( Channel.fromList(params.methods_corr) )
          .filter{ !(it[3] in params.no_rho && it[4] == "rho") }
          .set { ch_to_corr }

    /* 2nd step: Compute association coefficients */
    CORR(ch_to_corr)
}