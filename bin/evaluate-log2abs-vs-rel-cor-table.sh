basedir="/users/cn/sjin/projects/proportionality/run_scPropr/results/mouseStemCells/evaluate/pearson"
type="simulate+s8+n5+c1"
method_transf="log2_NA tmm_NA scran_NA clr_NA alr_EIF4A1 alr_LACTB2 alr_RPL8"
method_cor="cor rho pcor pcor.shrink vlr"

echo -n "method"
for i in $method_transf; do echo -n ","$i; done
echo

for i in $method_cor; do
    echo -n $i
    for j in $method_transf; do
        filename=${basedir}/mouseStemCells_${type}_full_zcompositions_${j}_${i}/mouseStemCells_${type}_full_zcompositions_${j}_${i}_pearson.csv
        if [ -f $filename ]; then
            awk -F',' 'NR==2{printf ",%f", $2}' $filename
        else
            echo -n ","
        fi
    done
    echo
done

