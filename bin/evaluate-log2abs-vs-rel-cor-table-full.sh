stat="rho"
method_zero="pseudocount"
basedir="/users/cn/sjin/projects/proportionality/run_scPropr/results-${method_zero}/mouseStemCells/evaluate/${stat}"
type="phaseS+rel phaseS+rel+simulate+s2+n5+c1 phaseS+rel+simulate+s4+n5+c1 phaseS+rel+simulate+s6+n5+c1 phaseS+rel+simulate+s8+n5+c1"
method_transf="log2_NA tmm_NA tmm2_NA scran_NA clr_NA alr_EIF4A1 alr_RPL8 alr_MT-RNR2 alr_TOP2A alr_UBE2C"
method_cor="cor rho pcor pcor.shrink vlr"

echo "type,transf,cor,val"
for t in $type; do
    for i in $method_cor; do
        for j in $method_transf; do
            filename=${basedir}/mouseStemCells_${t}_full_${method_zero}_${j}_${i}/mouseStemCells_${t}_full_${method_zero}_${j}_${i}_${stat}.csv
            echo -n $t,$j,$i
            if [ -f $filename ]; then
                awk -F',' 'NR==2{printf ",%f\n", $2}' $filename
            else
                echo ","
            fi
        done
    done
done

