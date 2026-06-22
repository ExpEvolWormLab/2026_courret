library(data.table)
library(qqman)

rm(list = ls())

chr_I   <- fread("GridLMM_results_mixed_perCHR/gwas_results_area.F_CHR_I_NGM.txt")
chr_II  <- fread("GridLMM_results_mixed_perCHR/gwas_results_area.F_CHR_II_NGM.txt")
chr_III <- fread("GridLMM_results_mixed_perCHR/gwas_results_area.F_CHR_III_NGM.txt")
chr_IV  <- fread("GridLMM_results_mixed_perCHR/gwas_results_area.F_CHR_IV_NGM.txt")
chr_V   <- fread("GridLMM_results_mixed_perCHR/gwas_results_area.F_CHR_V_NGM.txt")
chr_X   <- fread("GridLMM_results_mixed_perCHR/gwas_results_area.F_CHR_X_NGM.txt")
chr_MtDNA <- fread("GridLMM_results_mixed_perCHR/gwas_results_area.F_CHR_MtDNA_NGM.txt")

chr_I_var   <- fread("GridLMM_results_mixed_perCHR/gwas_results_log_area.F.var_CHR_I_NGM.txt")
chr_II_var  <- fread("GridLMM_results_mixed_perCHR/gwas_results_log_area.F.var_CHR_II_NGM.txt")
chr_III_var <- fread("GridLMM_results_mixed_perCHR/gwas_results_log_area.F.var_CHR_III_NGM.txt")
chr_IV_var  <- fread("GridLMM_results_mixed_perCHR/gwas_results_log_area.F.var_CHR_IV_NGM.txt")
chr_V_var   <- fread("GridLMM_results_mixed_perCHR/gwas_results_log_area.F.var_CHR_V_NGM.txt")
chr_X_var   <- fread("GridLMM_results_mixed_perCHR/gwas_results_log_area.F.var_CHR_X_NGM.txt")
chr_MtDNA_var <- fread("GridLMM_results_mixed_perCHR/gwas_results_log_area.F.var_CHR_MtDNA_NGM.txt")

chr_I_CV   <- fread("GridLMM_results_mixed_perCHR/gwas_results_log_CV_areaF_CHR_I_NGM.txt")
chr_II_CV  <- fread("GridLMM_results_mixed_perCHR/gwas_results_log_CV_areaF_CHR_II_NGM.txt")
chr_III_CV <- fread("GridLMM_results_mixed_perCHR/gwas_results_log_CV_areaF_CHR_III_NGM.txt")
chr_IV_CV  <- fread("GridLMM_results_mixed_perCHR/gwas_results_log_CV_areaF_CHR_IV_NGM.txt")
chr_V_CV   <- fread("GridLMM_results_mixed_perCHR/gwas_results_log_CV_areaF_CHR_V_NGM.txt")
chr_X_CV   <- fread("GridLMM_results_mixed_perCHR/gwas_results_log_CV_areaF_CHR_X_NGM.txt")
chr_MtDNA_CV <- fread("GridLMM_results_mixed_perCHR/gwas_results_log_CV_areaF_CHR_MtDNA_NGM.txt")

GWAS.mean <- rbind(chr_I, chr_II, chr_III, chr_IV, chr_V, chr_X)
GWAS.var <- rbind(chr_I_var, chr_II_var, chr_III_var, chr_IV_var, chr_V_var, chr_X_var)
GWAS.cv <- rbind(chr_I_CV, chr_II_CV, chr_III_CV, chr_IV_CV, chr_V_CV, chr_X_CV)


str(GWAS.mean)

# =========================================================
# LOAD GWAS RESULTS NaCl
# =========================================================

chr_I.NaCl   <- fread("GridLMM_results_mixed_perCHR_NaCl/gwas_results_area.F_CHR_I_NaCl.txt")
chr_II.NaCl  <- fread("GridLMM_results_mixed_perCHR_NaCl/gwas_results_area.F_CHR_II_NaCl.txt")
chr_III.NaCl <- fread("GridLMM_results_mixed_perCHR_NaCl/gwas_results_area.F_CHR_III_NaCl.txt")
chr_IV.NaCl  <- fread("GridLMM_results_mixed_perCHR_NaCl/gwas_results_area.F_CHR_IV_NaCl.txt")
chr_V.NaCl   <- fread("GridLMM_results_mixed_perCHR_NaCl/gwas_results_area.F_CHR_V_NaCl.txt")
chr_X.NaCl   <- fread("GridLMM_results_mixed_perCHR_NaCl/gwas_results_area.F_CHR_X_NaCl.txt")
chr_MtDNA.NaCl <- fread("GridLMM_results_mixed_perCHR_NaCl/gwas_results_area.F_CHR_MtDNA_NaCl.txt")

chr_I_var.NaCl   <- fread("GridLMM_results_mixed_perCHR_NaCl/gwas_results_log_area.F.var_CHR_I_NaCl.txt")
chr_II_var.NaCl  <- fread("GridLMM_results_mixed_perCHR_NaCl/gwas_results_log_area.F.var_CHR_II_NaCl.txt")
chr_III_var.NaCl <- fread("GridLMM_results_mixed_perCHR_NaCl/gwas_results_log_area.F.var_CHR_III_NaCl.txt")
chr_IV_var.NaCl  <- fread("GridLMM_results_mixed_perCHR_NaCl/gwas_results_log_area.F.var_CHR_IV_NaCl.txt")
chr_V_var.NaCl   <- fread("GridLMM_results_mixed_perCHR_NaCl/gwas_results_log_area.F.var_CHR_V_NaCl.txt")
chr_X_var.NaCl   <- fread("GridLMM_results_mixed_perCHR_NaCl/gwas_results_log_area.F.var_CHR_X_NaCl.txt")
chr_MtDNA_var.NaCl <- fread("GridLMM_results_mixed_perCHR_NaCl/gwas_results_log_area.F.var_CHR_MtDNA_NaCl.txt")

chr_I_CV.NaCl   <- fread("GridLMM_results_mixed_perCHR_NaCl/gwas_results_log_CV_areaF_CHR_I_NaCl.txt")
chr_II_CV.NaCl  <- fread("GridLMM_results_mixed_perCHR_NaCl/gwas_results_log_CV_areaF_CHR_II_NaCl.txt")
chr_III_CV.NaCl <- fread("GridLMM_results_mixed_perCHR_NaCl/gwas_results_log_CV_areaF_CHR_III_NaCl.txt")
chr_IV_CV.NaCl  <- fread("GridLMM_results_mixed_perCHR_NaCl/gwas_results_log_CV_areaF_CHR_IV_NaCl.txt")
chr_V_CV.NaCl   <- fread("GridLMM_results_mixed_perCHR_NaCl/gwas_results_log_CV_areaF_CHR_V_NaCl.txt")
chr_X_CV.NaCl   <- fread("GridLMM_results_mixed_perCHR_NaCl/gwas_results_log_CV_areaF_CHR_X_NaCl.txt")
chr_MtDNA_CV.NaCl <- fread("GridLMM_results_mixed_perCHR_NaCl/gwas_results_log_CV_areaF_CHR_MtDNA_NaCl.txt")

GWAS.mean.NaCl <- rbind(chr_I.NaCl, chr_II.NaCl, chr_III.NaCl, chr_IV.NaCl, chr_V.NaCl, chr_X.NaCl)
GWAS.var.NaCl <- rbind(chr_I_var.NaCl, chr_II_var.NaCl, chr_III_var.NaCl, chr_IV_var.NaCl, chr_V_var.NaCl, chr_X_var.NaCl)
GWAS.cv.NaCl <- rbind(chr_I_CV.NaCl, chr_II_CV.NaCl, chr_III_CV.NaCl, chr_IV_CV.NaCl, chr_V_CV.NaCl, chr_X_CV.NaCl)

# =========================================================
# Pruning
# =========================================================


system("mkdir -p plink_pruning/pruning_08")
system("/opt/homebrew/bin/plink --vcf data/Cemee2025.imputed.SNP.filtered.final_renamed.vcf.gz --indep-pairwise 200 10 0.8 --out plink_pruning/pruning_08/CeMEE.pruned.0.8 --allow-extra-chr --double-id")

# SNPs conservés après pruning
snps_pruned <- fread(
  "plink_pruning/pruning_08/CeMEE.pruned.0.8.prune.in",
  header = FALSE
)

setnames(snps_pruned, "V1", "X_ID")
head(snps_pruned)

GWAS.mean.pruned <- GWAS.mean[X_ID %in% snps_pruned$X_ID]
GWAS.var.pruned  <- GWAS.var[X_ID %in% snps_pruned$X_ID]
GWAS.cv.pruned   <- GWAS.cv[X_ID %in% snps_pruned$X_ID]

GWAS.mean.NaCl.pruned <- GWAS.mean.NaCl[X_ID %in% snps_pruned$X_ID]
GWAS.var.NaCl.pruned  <- GWAS.var.NaCl[X_ID %in% snps_pruned$X_ID]
GWAS.cv.NaCl.pruned   <- GWAS.cv.NaCl[X_ID %in% snps_pruned$X_ID]


write.table(GWAS.mean.pruned, file = "data/GWAS.mean.pruned.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(GWAS.var.pruned, file = "data/GWAS.var.pruned.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(GWAS.cv.pruned, file = "data/GWAS.cv.pruned.txt", sep = "\t", quote = FALSE, row.names = FALSE)

write.table(GWAS.mean.NaCl.pruned, file = "data/GWAS.mean.NaCl.pruned.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(GWAS.var.NaCl.pruned, file = "data/GWAS.var.NaCl.pruned.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(GWAS.cv.NaCl.pruned, file = "data/GWAS.cv.NaCl.pruned.txt", sep = "\t", quote = FALSE, row.names = FALSE)
