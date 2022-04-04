
process PROCRUSTES_ANALYSIS {
    label 'process_high'
    container 'suzannejin/scpropr:propr'
    tag "$dataset-$exp_sim-$full-$abs_rel-$method_replace_zero-$gene1"
    storeDir "${params.outdir}/${dataset}/procrustes/${dataset}_${exp_sim}_${full}_${abs_rel}_${method_replace_zero}_${gene1}"

    input:
    tuple val(dataset),
          val(exp_sim),
          val(full),
          val(abs_rel),
          val(method_replace_zero),
          file(count),
          file(features),
          file(barcodes),
          file(lra),
          val(gene)

    output:
    tuple file("${dataset}_${exp_sim}_${full}_${abs_rel}_${method_replace_zero}_${gene1}.stats"),
          file(".command.trace"),
          file(".command.sh")

    script:
    gene1 = gene.getSimpleName()
    """
    procrustes-analysis-parallel.R \
        --count $count \
        --lra $lra \
        --features $features \
        --gene $gene1 \
        --output ${dataset}_${exp_sim}_${full}_${abs_rel}_${method_replace_zero}_${gene1}.stats 
    """
}