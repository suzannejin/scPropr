## usage, eg: bash bin/evaluate-log2abs-vs-rel-cor-table-full.sh > ../run_scPropr/results-zcompositions/mouseStemCells/evaluate/rho/evaluate.stats

stat="rho"
full="fixed2"
filt="evaluate"
basedir="/users/cn/sjin/projects/proportionality/run_scPropr/results-corrected-zcompositions-15112022/mouseStemCells/${filt}/${stat}"
type="phaseS+simulate+s1+n5+c1+d1 phaseS+simulate+s2+n5+c1+d1 phaseS+simulate+s4+n5+c1+d1 phaseS+simulate+s6+n5+c1+d1" # phaseS+simulate+s8+n5+c1+d1"
method_transf="NA_log2_NA NA_tmm_NA NA_scran_NA NA_clr_NA NA_alr_MT-RNR2 NA_alr_RPL8"
# method_transf="NA_log2_NA NA_tmm_NA NA_scran_NA bgm_log2_NA bgm_tmm_NA bgm_scran_NA bgm_clr_NA bgm_alr_MT-RNR2 bgm_alr_RPL8 czm_log2_NA czm_tmm_NA czm_scran_NA czm_clr_NA czm_alr_MT-RNR2 czm_alr_RPL8 pseudocount_log2_NA pseudocount_tmm_NA pseudocount_scran_NA pseudocount_clr_NA pseudocount_alr_MT-RNR2 pseudocount_alr_RPL8" #  zcompositions_alr_EIF4A1  zcompositions_alr_HIST1H1B zcompositions_alr_UBE2C"
method_cor="cor rho pcor pcor.shrink vlr"

echo "type,transf,cor,val"
for t in $type; do
    for i in $method_cor; do
        for j in $method_transf; do
            filename=${basedir}/mouseStemCells_${t}_${full}_${j}_${i}/mouseStemCells_${t}_${full}_${j}_${i}_${stat}.csv
            # j=${j##NA_}
            # j=${j##zcompositions_}
            echo -n $t,$j,$i
            if [ -f $filename ]; then
                awk -F',' 'NR==2{printf ",%f\n", $2}' $filename
            else
                echo ","
                echo $filename
                exit 1
            fi
        done
    done
done

