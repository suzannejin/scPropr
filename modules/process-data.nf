
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
          file(features),
          file(barcodes),
          file(".command.trace"),
          file(".command.sh")

    when:
    if (full == 'full'){
        (refgene == 'NA' && method_transform_data !in ['alr', 'alr2']) || (refgene != 'NA' && method_transform_data in ['alr', 'alr2'])
    } else if (full == 'nozero'){
        (refgene == 'NA' && method_transform_data !in ['alr', 'alr2']) || (refgene in params.refgenes_nozero && method_transform_data in ['alr', 'alr2'])
    }

    script:
    def ref_gene_cml = refgene == 'NA' ? " " : "--refgene $refgene" 
    """
    process-data.R \
        -i $count \
        -o ${dataset}_${exp_sim}_${full}_${abs_rel}_${method_replace_zero}_${method_transform_data}_${refgene}.csv.gz \
        -o2 ${dataset}_${exp_sim}_${full}_${abs_rel}_${method_replace_zero}_${method_transform_data}_${refgene}.list \
        --features $features \
        --method_zero $method_replace_zero \
        --method_transf $method_transform_data \
        $ref_gene_cml 
    sleep 30
    """

    stub:
    def ref_gene_cml = method_transform_data == 'alr' ? "--refgene $refgene" : " "
    """
    echo process-data.R \
        -i $count \
        -o ${dataset}_${exp_sim}_${full}_${abs_rel}_${method_replace_zero}_${method_transform_data}_${refgene}.csv.gz \
        -o2 ${dataset}_${exp_sim}_${full}_${abs_rel}_${method_replace_zero}_${method_transform_data}_${refgene}.list \
        --features $features \
        --method_zero $method_replace_zero \
        --method_transf $method_transform_data \
        $ref_gene_cml
    touch ${dataset}_${exp_sim}_${full}_${abs_rel}_${method_replace_zero}_${method_transform_data}_${refgene}.csv.gz
    touch ${dataset}_${exp_sim}_${full}_${abs_rel}_${method_replace_zero}_${method_transform_data}_${refgene}.list
    """
}