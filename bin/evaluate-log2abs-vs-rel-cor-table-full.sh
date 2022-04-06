basedir="/users/cn/sjin/projects/proportionality/run_scPropr/results/mouseStemCells/evaluate/pearson"
type="phaseS simulate+s2+n5+c1+byslope simulate+s4+n5+c1+byslope simulate+s6+n5+c1+byslope simulate+s8+n5+c1+byslope"
method_transf="log2_NA tmm_NA scran_NA clr_NA alr_EIF4A1 alr_LACTB2"
method_cor="cor rho pcor pcor.shrink vlr"

echo "type,transf,cor,val"
for t in $type; do
    for i in $method_cor; do
        for j in $method_transf; do
            filename=${basedir}/mouseStemCells_${t}_full_zcompositions_${j}_${i}/mouseStemCells_${t}_full_zcompositions_${j}_${i}_pearson.csv
            echo -n $t,$j,$i
            if [ -f $filename ]; then
                awk -F',' 'NR==2{printf ",%f\n", $2}' $filename
            else
                echo ","
            fi
        done
    done
done

