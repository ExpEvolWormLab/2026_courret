mkdir -p output

perm=1

for i in $(seq 1 10); do
    echo "Boucle $i"

    for j in $(seq 1 100); do
        echo "Permutation $perm"

        # Randomiser phénotype
        shuf $pheno_file > ${data/pheno_perm.txt}_${perm}

        # Lancer GEMMA en arrière-plan
        gemma \
            -g data/genotypes_for_gemma.txt \
            -k data/kinship.txt \
            -lmm 4 \
            -p data/gemma_pheno.txt_${perm} \
            -o perm_$perm \
            > output/perm_$perm.txt &

        perm=$((perm + 1))
    done
    wait
    echo "Boucle $i terminé"

done

echo "Toutes les permutations sont terminées"