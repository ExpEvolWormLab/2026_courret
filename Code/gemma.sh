nohup gemma \
  -g data/genotypes_for_gemma__log_var.txt \
  -k data/kinship.txt \
  -lmm 4 \
  -p data/gemma_pheno_var.txt \
  -o gwas_log_var_results \
  > gemma_output_txt/output_gemma_log_var.txt &