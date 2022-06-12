## usage, eg: bash bin/evaluate-log2abs-vs-rel-cor-table-full.sh > ../run_scPropr/results-zcompositions/mouseStemCells/evaluate/rho/evaluate.stats

stat="rho"
full="nozero"
filt="evaluate"
basedir="/users/cn/sjin/projects/proportionality/run_scPropr/results-default-pseudocount-withnozero-10062022/mouseStemCells/${filt}/${stat}"
type="experimental phaseS+simulate+s2+n5+c3+d3 phaseS+simulate+s4+n5+c3+d3 phaseS+simulate+s6+n5+c3+d3 phaseS+simulate+s8+n5+c3+d3"
# method_transf="NA_log2_NA NA_tmm_NA NA_scran_NA pseudocount_clr_NA pseudocount_alr_EIF4A1 pseudocount_alr_RPL8 pseudocount_alr_MT-RNR2 pseudocount_alr_HIST1H1B pseudocount_alr_UBE2C"
method_transf="NA_log2_NA NA_tmm_NA NA_scran_NA NA_clr_NA NA_alr_RPL8 NA_alr_MT-RNR2"
method_cor="cor rho pcor pcor.shrink vlr"

echo "type,transf,cor,val"
for t in $type; do
    for i in $method_cor; do
        for j in $method_transf; do
            filename=${basedir}/mouseStemCells_${t}_${full}_${j}_${i}/mouseStemCells_${t}_${full}_${j}_${i}_${stat}.csv
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

