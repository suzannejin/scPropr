
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
          file(features),
          file(barcodes),
          file(".command.trace"),
          file(".command.sh"),
          file(".command.log"),
          emit: count
    file("${dataset}_${exp_sim}_${full}_${abs_rel}_${method_replace_zero}_${method_transform_data}_${refgene}.list"),
        optional: true,
        emit: size_factor
    file("${dataset}_${exp_sim}_${full}_${abs_rel}_${method_replace_zero}_${method_transform_data}_${refgene}.errorbars.gz"),
        optional: true,
        emit: error_bars

    when:
    if (params.do_transform_abs || abs_rel == 'relative'){
        if (full.startsWith('nozero')){
            if (method_replace_zero == 'NA'){
                if (method_transform_data == 'alr' && refgene != 'NA'){true}
                else if (method_transform_data != 'alr' && refgene == 'NA'){true}
            }
        }else{
            if (method_replace_zero != 'NA' && method_transform_data != 'sanity'){
                if (method_transform_data == 'alr' && refgene != 'NA'){true}
                else if (method_transform_data != 'alr' && refgene == 'NA'){true}
            }else if (method_replace_zero == 'NA' && method_transform_data == 'sanity' && refgene == 'NA'){true}
        }
    }else{
        if (abs_rel == 'absolute'){
            if (method_transform_data == 'log2' && refgene == 'NA'){
                if (full.startsWith('nozero')){
                    if (method_replace_zero == 'NA'){true}
                } else{
                    if (method_replace_zero != 'NA'){true}
                }
            }
        }
    }

    script:
    if (method_transform_data != 'sanity'){
        def script = full.startsWith('nozero') ? "process-data-nozero.R" : "process-data.R"
        def ref_gene_cml = refgene == 'NA' ? " " : "--refgene $refgene" 
        def method_zero = method_replace_zero == 'NA' ? " " : "--method_zero $method_replace_zero"
        """
        $script \
            -i $count \
            -o ${dataset}_${exp_sim}_${full}_${abs_rel}_${method_replace_zero}_${method_transform_data}_${refgene}.csv.gz \
            -o2 ${dataset}_${exp_sim}_${full}_${abs_rel}_${method_replace_zero}_${method_transform_data}_${refgene}.list \
            --features $features \
            $method_zero \
            --method_transf $method_transform_data \
            $ref_gene_cml 
        """
    }else{
        """
        pack_count_data.R \
            -c $count \
            -f $features \
            -b $barcodes \
            -o tmp.csv
        Sanity -f tmp.csv
        transpose_and_compress_count.R \
            -i log_transcription_quotients.txt \
            -o ${dataset}_${exp_sim}_${full}_${abs_rel}_${method_replace_zero}_${method_transform_data}_${refgene}.csv.gz
        transpose_and_compress_count.R
            -i ltq_error_bars.txt \
            -o ${dataset}_${exp_sim}_${full}_${abs_rel}_${method_replace_zero}_${method_transform_data}_${refgene}.errorbars.gz
        """
    }
    

    stub:
    if (method_transform_data != 'sanity'){
        def script = full.startsWith('nozero') ? "process-data-nozero.R" : "process-data.R"
        def ref_gene_cml = refgene == 'NA' ? " " : "--refgene $refgene" 
        def method_zero = method_replace_zero == 'NA' ? " " : "--method_zero $method_replace_zero"
        """
        echo $script \
            -i $count \
            -o ${dataset}_${exp_sim}_${full}_${abs_rel}_${method_replace_zero}_${method_transform_data}_${refgene}.csv.gz \
            -o2 ${dataset}_${exp_sim}_${full}_${abs_rel}_${method_replace_zero}_${method_transform_data}_${refgene}.list \
            --features $features \
            $method_zero \
            --method_transf $method_transform_data \
            $ref_gene_cml 
        touch ${dataset}_${exp_sim}_${full}_${abs_rel}_${method_replace_zero}_${method_transform_data}_${refgene}.csv.gz
        touch ${dataset}_${exp_sim}_${full}_${abs_rel}_${method_replace_zero}_${method_transform_data}_${refgene}.list
        """
    }else{
        """
        echo pack_count_data.R \
            -c $count \
            -f $features \
            -b $barcodes \
            -o tmp.csv
        echo Sanity -f tmp.csv
        echo transpose_and_compress_count.R \
            -i log_transcription_quotients.txt \
            -o ${dataset}_${exp_sim}_${full}_${abs_rel}_${method_replace_zero}_${method_transform_data}_${refgene}.csv.gz
        echo transpose_and_compress_count.R
            -i ltq_error_bars.txt \
            -o ${dataset}_${exp_sim}_${full}_${abs_rel}_${method_replace_zero}_${method_transform_data}_${refgene}.errorbars.gz
        touch ${dataset}_${exp_sim}_${full}_${abs_rel}_${method_replace_zero}_${method_transform_data}_${refgene}.csv.gz
        touch ${dataset}_${exp_sim}_${full}_${abs_rel}_${method_replace_zero}_${method_transform_data}_${refgene}.errorbars.gz
        """
    }
}