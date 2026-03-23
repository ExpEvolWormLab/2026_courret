#!/bin/bash

mkdir -p output

perm=1

for i in $(seq 1 50); do
    echo "Boucle $i"

    for j in $(seq 1 20); do
        echo "Permutation $perm"

        # Randomiser phénotype
        shuf data/gemma_pheno_log_CV.txt > data/pheno_log_CV_perm_${perm}.txt

        # Lancer GEMMA en arrière-plan
        gemma \
            -g data/genotypes_for_gemma.txt \
            -k data/kinship.txt \
            -lmm 4 \
            -p data/pheno_log_CV_perm_${perm}.txt \
            -o perm_log_CV_$perm \
            > output/perm_log_CV_$perm.txt &

        perm=$((perm + 1))
    done
    wait
    echo "Boucle $i terminé"

done

echo "Toutes les permutations sont terminées"