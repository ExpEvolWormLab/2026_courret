

## My results
a=read.table("~/Documents/GitHub/2025_Breeding_values/GridLMM_results_mixed_perCHR/gwas_results_area.F_CHR_I_NGM.txt",h=TRUE)
head(a)
#b=read.table("data_FM/GWAS.mean.txt",h=TRUE)
b=read.table("GridLMM_results_mixed_perCHR_FM/gwas_results_area.F_CHR_I_NGM.txt",h=TRUE)
#paste0("GridLMM_results_mixed_perCHR_FM/gwas_results_",traitID,"_CHR_",chrID,"_NGM.txt")
head(b)

a2=a[,c("X_ID","beta.5","p_value_REML")]
names(a2)[2:3]=c("beta_FM","pval_FM")

b2=b[,c("X_ID","beta.5","p_value_REML")]
names(b2)[2:3]=c("beta_EC","pval_EC")

c=merge(a2,b2)
dim(c)
plot(c$beta_FM,c$beta_EC)
plot(log10(c$pval_FM),log10(c$pval_EC))

FM.env=new.env()
load('~/Documents/GitHub/2025_Breeding_values/RData/Preparatory_GLMM.RData',envir=FM.env)
ls(FM.env)



dim(FM.env$all_BLUEs_genotyped_NGM)
head(all_BLUEs_genotyped_NGM)








