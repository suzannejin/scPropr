
/* Process data
 *    - handle zero
 *    - normalize data
 *    - log transform data
 */

process PROCESS_DATA {
    
    label 'process_high_long'
    container 'suzannejin/scpropr:propr'
    tag "${dataset}-${exp_sim}-${full}-${abs_rel}-${method_replace_zero}-${method_transform_data}-${refgene}"
    storeDir "${params.outdir}/${dataset}/processed/${dataset}_${exp_sim}_${full}_${abs_rel}_${method_replace_zero}_${method_transform_data}_${refgene}"
    // publishDir "${params.outdir}/${dataset}/processed/${dataset}_${exp_sim}_${full}_${abs_rel}_${method_replace_zero}_${method_transform_data}_${refgene}", mode: params.publish_dir_mode
    
    input:
    tuple val(dataset), 
          val(exp_sim), 
          val(full), 
          val(abs_rel), 
          file(count), 
          file(features), 
          file(barcodes)
    each method_replace_zero
    each method_transform_data
    each refgene

    output:
    tuple file("${dataset}_${exp_sim}_${full}_${abs_rel}_${method_replace_zero}_${method_transform_data}_${refgene}.csv.gz"), 
          file("${dataset}_${exp_sim}_${full}_${abs_rel}_${method_replace_zero}_${method_transform_data}_${refgene}.list"),
          file("${dataset}_${exp_sim}_${full}_${abs_rel}_${method_replace_zero}_${method_transform_data}_${refgene}.lambda"),
          file(features),
          file(barcodes),
          file(".command.trace"),
          file(".command.sh"),
          file(".command.log")

    when:
    if (!params.do_transform_abs && abs_rel == 'absolute'){
        if (method_transform_data == 'log2' && refgene == 'NA'){true}
    }else{
        if (method_replace_zero == 'NA' && method_transform_data in ['log2','tmm','scran'] && refgene == 'NA'){true}
        else if (method_replace_zero != 'NA' && method_transform_data == 'clr' && refgene == 'NA' && full == 'full'){true}
        else if (method_replace_zero == 'NA' && method_transform_data == 'clr' && refgene == 'NA' && full == 'nozero'){true}
        else if (method_replace_zero != 'NA' && method_transform_data == 'alr' && refgene != 'NA' && full == 'full'){true}
        else if (method_replace_zero == 'NA' && method_transform_data == 'alr' && refgene in params.refgenes_nozero && full == 'nozero'){true}
    }

    script:
    if (full == 'full'){
        def ref_gene_cml = refgene == 'NA' ? " " : "--refgene $refgene" 
        def method_zero = method_replace_zero == 'NA' ? " " : "--method_zero $method_replace_zero"
        """
        process-data.R \
            -i $count \
            -o ${dataset}_${exp_sim}_${full}_${abs_rel}_${method_replace_zero}_${method_transform_data}_${refgene}.csv.gz \
            -o2 ${dataset}_${exp_sim}_${full}_${abs_rel}_${method_replace_zero}_${method_transform_data}_${refgene}.list \
            --lambda ${dataset}_${exp_sim}_${full}_${abs_rel}_${method_replace_zero}_${method_transform_data}_${refgene}.lambda \
            --features $features \
            $method_zero \
            --method_transf $method_transform_data \
            $ref_gene_cml 
        sleep 30
        """
    }else{
        def ref_gene_cml = refgene == 'NA' ? " " : "--refgene $refgene" 
        """
        process-data-nozero.R \
            -i $count \
            -o ${dataset}_${exp_sim}_${full}_${abs_rel}_${method_replace_zero}_${method_transform_data}_${refgene}.csv.gz \
            -o2 ${dataset}_${exp_sim}_${full}_${abs_rel}_${method_replace_zero}_${method_transform_data}_${refgene}.list \
            --features $features \
            --method_transf $method_transform_data \
            $ref_gene_cml 
        touch ${dataset}_${exp_sim}_${full}_${abs_rel}_${method_replace_zero}_${method_transform_data}_${refgene}.lambda
        sleep 30
        """
    }


    stub:
    if (full == 'full'){
        def ref_gene_cml = refgene == 'NA' ? " " : "--refgene $refgene" 
        def method_zero = method_replace_zero == 'NA' ? " " : "--method_zero $method_replace_zero"
        """
        echo process-data.R \
            -i $count \
            -o ${dataset}_${exp_sim}_${full}_${abs_rel}_${method_replace_zero}_${method_transform_data}_${refgene}.csv.gz \
            -o2 ${dataset}_${exp_sim}_${full}_${abs_rel}_${method_replace_zero}_${method_transform_data}_${refgene}.list \
            --lambda ${dataset}_${exp_sim}_${full}_${abs_rel}_${method_replace_zero}_${method_transform_data}_${refgene}.lambda \
            --features $features \
            $method_zero \
            --method_transf $method_transform_data \
            $ref_gene_cml 
        touch ${dataset}_${exp_sim}_${full}_${abs_rel}_${method_replace_zero}_${method_transform_data}_${refgene}.csv.gz
        touch ${dataset}_${exp_sim}_${full}_${abs_rel}_${method_replace_zero}_${method_transform_data}_${refgene}.list
        touch ${dataset}_${exp_sim}_${full}_${abs_rel}_${method_replace_zero}_${method_transform_data}_${refgene}.lambda
        sleep 30
        """
    }else{
        def ref_gene_cml = refgene == 'NA' ? " " : "--refgene $refgene" 
        """
        process-data-nozero.R \
            -i $count \
            -o ${dataset}_${exp_sim}_${full}_${abs_rel}_${method_replace_zero}_${method_transform_data}_${refgene}.csv.gz \
            -o2 ${dataset}_${exp_sim}_${full}_${abs_rel}_${method_replace_zero}_${method_transform_data}_${refgene}.list \
            --features $features \
            --method_transf $method_transform_data \
            $ref_gene_cml 
        touch ${dataset}_${exp_sim}_${full}_${abs_rel}_${method_replace_zero}_${method_transform_data}_${refgene}.csv.gz
        touch ${dataset}_${exp_sim}_${full}_${abs_rel}_${method_replace_zero}_${method_transform_data}_${refgene}.list
        touch ${dataset}_${exp_sim}_${full}_${abs_rel}_${method_replace_zero}_${method_transform_data}_${refgene}.lambda
        sleep 30
        """
    }
}