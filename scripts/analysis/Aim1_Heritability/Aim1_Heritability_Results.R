# rm(list = ls())
# 
# library(data.table)
# 
# sumstats <- fread("/data/williamsjacr/Aim1_Overall/GWAS_SumStats/AA_Overall_HM3_01_2024.txt")
# sumstats <- as.data.frame(sumstats)
# 
# Z <- as.numeric(sumstats$BETA/sumstats$SE)
# P <- as.numeric(sumstats$P)
# A1 <- sumstats$non_effect_allele
# A2 <- sumstats$effect_allele
# N_eff <- as.numeric(sumstats$N_eff)
# 
# SNP_ID <- sumstats$SNPID
# sumstats <- data.frame(SNP = SNP_ID,A1 = A1,A2 = A2, P = P, Z = Z,N_eff = N_eff)
# fwrite(sumstats,file = "/data/williamsjacr/Aim1_Overall/Intermediate/sumstats_Overall",row.names = F,col.names = T,quote = F,sep = " ")
# 
# system(paste0("munge_sumstats.py",
#               " --sumstats /data/williamsjacr/Aim1_Overall/Intermediate/sumstats_Overall",
#               " --out /data/williamsjacr/Aim1_Overall/Intermediate/sum_align",
#               " --merge-alleles /data/williamsjacr/1KG_ldsc_HM3/AFR/AFR_hm3_ldsc/w_hm3_ldsc.snplist",
#               " --chunksize 500000 --N-col N_eff --a2 A2 --a1 A1 --p P --signed-sumstats Z,0 --snp SNP",
#               " --maf-min 0.01"))
# 
# system(paste0("ldsc.py",
#               " --h2 /data/williamsjacr/Aim1_Overall/Intermediate/sum_align.sumstats.gz",
#               " --ref-ld-chr /data/williamsjacr/1KG_ldsc_HM3/AFR/AFR_hm3_ldsc/",
#               " --w-ld-chr /data/williamsjacr/1KG_ldsc_HM3/AFR/AFR_hm3_ldsc/",
#               " --out /data/williamsjacr/Aim1_Overall/Results/AFR_HM3_h2_observed_01"))
# 
# rm(list = ls())
# 
# library(data.table)
# 
# sumstats <- fread("/data/williamsjacr/Aim1_Overall/GWAS_SumStats/AA_Overall_HM3_MEGA_01_2024.txt")
# sumstats <- as.data.frame(sumstats)
# 
# Z <- as.numeric(sumstats$BETA/sumstats$SE)
# P <- as.numeric(sumstats$P)
# A1 <- sumstats$non_effect_allele
# A2 <- sumstats$effect_allele
# N_eff <- as.numeric(sumstats$N_eff)
# 
# SNP_ID <- sumstats$SNPID
# sumstats <- data.frame(SNP = SNP_ID,A1 = A1,A2 = A2, P = P, Z = Z,N_eff = N_eff)
# fwrite(sumstats,file = "/data/williamsjacr/Aim1_Overall/Intermediate/sumstats_Overall",row.names = F,col.names = T,quote = F,sep = " ")
# 
# system(paste0("munge_sumstats.py",
#               " --sumstats /data/williamsjacr/Aim1_Overall/Intermediate/sumstats_Overall",
#               " --out /data/williamsjacr/Aim1_Overall/Intermediate/sum_align",
#               " --merge-alleles /data/williamsjacr/1KG_ldsc_HM3_MEGA/AFR/AFR_hm3_mega_ldsc/hm3_mega_ldsc.snplist",
#               " --chunksize 500000 --N-col N_eff --a2 A2 --a1 A1 --p P --signed-sumstats Z,0 --snp SNP",
#               " --maf-min 0.01"))
# 
# system(paste0("ldsc.py",
#               " --h2 /data/williamsjacr/Aim1_Overall/Intermediate/sum_align.sumstats.gz",
#               " --ref-ld-chr /data/williamsjacr/1KG_ldsc_HM3_MEGA/AFR/AFR_hm3_mega_ldsc/",
#               " --w-ld-chr /data/williamsjacr/1KG_ldsc_HM3_MEGA/AFR/AFR_hm3_mega_ldsc/",
#               " --out /data/williamsjacr/Aim1_Overall/Results/AFR_HM3_MEGA_h2_observed_01"))
# 
# rm(list = ls())
# 
# library(data.table)
# 
# sumstats <- fread("/data/williamsjacr/Aim1_Overall/GWAS_SumStats/AA_Overall_HM3_05_2024.txt")
# sumstats <- as.data.frame(sumstats)
# 
# Z <- as.numeric(sumstats$BETA/sumstats$SE)
# P <- as.numeric(sumstats$P)
# A1 <- sumstats$non_effect_allele
# A2 <- sumstats$effect_allele
# N_eff <- as.numeric(sumstats$N_eff)
# 
# SNP_ID <- sumstats$SNPID
# sumstats <- data.frame(SNP = SNP_ID,A1 = A1,A2 = A2, P = P, Z = Z,N_eff = N_eff)
# fwrite(sumstats,file = "/data/williamsjacr/Aim1_Overall/Intermediate/sumstats_Overall",row.names = F,col.names = T,quote = F,sep = " ")
# 
# system(paste0("munge_sumstats.py",
#               " --sumstats /data/williamsjacr/Aim1_Overall/Intermediate/sumstats_Overall",
#               " --out /data/williamsjacr/Aim1_Overall/Intermediate/sum_align",
#               " --merge-alleles /data/williamsjacr/1KG_ldsc_HM3/AFR/AFR_hm3_ldsc/w_hm3_ldsc.snplist",
#               " --chunksize 500000 --N-col N_eff --a2 A2 --a1 A1 --p P --signed-sumstats Z,0 --snp SNP",
#               " --maf-min 0.01"))
# 
# system(paste0("ldsc.py",
#               " --h2 /data/williamsjacr/Aim1_Overall/Intermediate/sum_align.sumstats.gz",
#               " --ref-ld-chr /data/williamsjacr/1KG_ldsc_HM3/AFR/AFR_hm3_ldsc/",
#               " --w-ld-chr /data/williamsjacr/1KG_ldsc_HM3/AFR/AFR_hm3_ldsc/",
#               " --out /data/williamsjacr/Aim1_Overall/Results/AFR_HM3_h2_observed_05"))
# 
# rm(list = ls())
# 
# library(data.table)
# 
# sumstats <- fread("/data/williamsjacr/Aim1_Overall/GWAS_SumStats/AA_Overall_HM3_MEGA_05_2024.txt")
# sumstats <- as.data.frame(sumstats)
# 
# Z <- as.numeric(sumstats$BETA/sumstats$SE)
# P <- as.numeric(sumstats$P)
# A1 <- sumstats$non_effect_allele
# A2 <- sumstats$effect_allele
# N_eff <- as.numeric(sumstats$N_eff)
# 
# SNP_ID <- sumstats$SNPID
# sumstats <- data.frame(SNP = SNP_ID,A1 = A1,A2 = A2, P = P, Z = Z,N_eff = N_eff)
# fwrite(sumstats,file = "/data/williamsjacr/Aim1_Overall/Intermediate/sumstats_Overall",row.names = F,col.names = T,quote = F,sep = " ")
# 
# system(paste0("munge_sumstats.py",
#               " --sumstats /data/williamsjacr/Aim1_Overall/Intermediate/sumstats_Overall",
#               " --out /data/williamsjacr/Aim1_Overall/Intermediate/sum_align",
#               " --merge-alleles /data/williamsjacr/1KG_ldsc_HM3_MEGA/AFR/AFR_hm3_mega_ldsc/hm3_mega_ldsc.snplist",
#               " --chunksize 500000 --N-col N_eff --a2 A2 --a1 A1 --p P --signed-sumstats Z,0 --snp SNP",
#               " --maf-min 0.01"))
# 
# system(paste0("ldsc.py",
#               " --h2 /data/williamsjacr/Aim1_Overall/Intermediate/sum_align.sumstats.gz",
#               " --ref-ld-chr /data/williamsjacr/1KG_ldsc_HM3_MEGA/AFR/AFR_hm3_mega_ldsc/",
#               " --w-ld-chr /data/williamsjacr/1KG_ldsc_HM3_MEGA/AFR/AFR_hm3_mega_ldsc/",
#               " --out /data/williamsjacr/Aim1_Overall/Results/AFR_HM3_MEGA_h2_observed_05"))





















rm(list = ls())

library(data.table)

sumstats <- fread("/data/williamsjacr/Aim1_Overall/GWAS_SumStats/AA_Overall_HM3_01_2024.txt")
sumstats <- as.data.frame(sumstats)

Z <- as.numeric(sumstats$BETA/sumstats$SE)
P <- as.numeric(sumstats$P)
A1 <- sumstats$non_effect_allele
A2 <- sumstats$effect_allele
N_eff <- as.numeric(sumstats$N_eff)

SNP_ID <- sumstats$SNPID
sumstats <- data.frame(SNP = SNP_ID,A1 = A1,A2 = A2, P = P, Z = Z,N_eff = N_eff)
fwrite(sumstats,file = "/data/williamsjacr/Aim1_Overall/Intermediate/sumstats_Overall",row.names = F,col.names = T,quote = F,sep = " ")

system(paste0("munge_sumstats.py",
              " --sumstats /data/williamsjacr/Aim1_Overall/Intermediate/sumstats_Overall",
              " --out /data/williamsjacr/Aim1_Overall/Intermediate/sum_align",
              " --merge-alleles /data/williamsjacr/1KG_ldsc_HM3/AFR/AFR_hm3_ldsc/w_hm3_ldsc.snplist",
              " --chunksize 500000 --N-col N_eff --a2 A2 --a1 A1 --p P --signed-sumstats Z,0 --snp SNP",
              " --maf-min 0.01"))

system(paste0("ldsc.py",
              " --h2 /data/williamsjacr/Aim1_Overall/Intermediate/sum_align.sumstats.gz",
              " --ref-ld-chr /data/williamsjacr/AFR_bd38_MEGA_HM3_Jia_rsid_July2024/",
              " --w-ld-chr /data/williamsjacr/AFR_bd38_MEGA_HM3_Jia_rsid_July2024/",
              " --out /data/williamsjacr/Aim1_Overall/Results/AFR_HM3_AABCG_h2_observed_01"))

rm(list = ls())

library(data.table)

sumstats <- fread("/data/williamsjacr/Aim1_Overall/GWAS_SumStats/AA_Overall_HM3_MEGA_01_2024.txt")
sumstats <- as.data.frame(sumstats)

Z <- as.numeric(sumstats$BETA/sumstats$SE)
P <- as.numeric(sumstats$P)
A1 <- sumstats$non_effect_allele
A2 <- sumstats$effect_allele
N_eff <- as.numeric(sumstats$N_eff)

SNP_ID <- sumstats$SNPID
sumstats <- data.frame(SNP = SNP_ID,A1 = A1,A2 = A2, P = P, Z = Z,N_eff = N_eff)
fwrite(sumstats,file = "/data/williamsjacr/Aim1_Overall/Intermediate/sumstats_Overall",row.names = F,col.names = T,quote = F,sep = " ")

system(paste0("munge_sumstats.py",
              " --sumstats /data/williamsjacr/Aim1_Overall/Intermediate/sumstats_Overall",
              " --out /data/williamsjacr/Aim1_Overall/Intermediate/sum_align",
              " --merge-alleles /data/williamsjacr/1KG_ldsc_HM3_MEGA/AFR/AFR_hm3_mega_ldsc/hm3_mega_ldsc.snplist",
              " --chunksize 500000 --N-col N_eff --a2 A2 --a1 A1 --p P --signed-sumstats Z,0 --snp SNP",
              " --maf-min 0.01"))

system(paste0("ldsc.py",
              " --h2 /data/williamsjacr/Aim1_Overall/Intermediate/sum_align.sumstats.gz",
              " --ref-ld-chr /data/williamsjacr/AFR_bd38_MEGA_HM3_Jia_rsid_July2024/",
              " --w-ld-chr /data/williamsjacr/AFR_bd38_MEGA_HM3_Jia_rsid_July2024/",
              " --out /data/williamsjacr/Aim1_Overall/Results/AFR_HM3_MEGA_AABCG_h2_observed_01"))

rm(list = ls())

library(data.table)

sumstats <- fread("/data/williamsjacr/Aim1_Overall/GWAS_SumStats/AA_Overall_HM3_05_2024.txt")
sumstats <- as.data.frame(sumstats)

Z <- as.numeric(sumstats$BETA/sumstats$SE)
P <- as.numeric(sumstats$P)
A1 <- sumstats$non_effect_allele
A2 <- sumstats$effect_allele
N_eff <- as.numeric(sumstats$N_eff)

SNP_ID <- sumstats$SNPID
sumstats <- data.frame(SNP = SNP_ID,A1 = A1,A2 = A2, P = P, Z = Z,N_eff = N_eff)
fwrite(sumstats,file = "/data/williamsjacr/Aim1_Overall/Intermediate/sumstats_Overall",row.names = F,col.names = T,quote = F,sep = " ")

system(paste0("munge_sumstats.py",
              " --sumstats /data/williamsjacr/Aim1_Overall/Intermediate/sumstats_Overall",
              " --out /data/williamsjacr/Aim1_Overall/Intermediate/sum_align",
              " --merge-alleles /data/williamsjacr/1KG_ldsc_HM3/AFR/AFR_hm3_ldsc/w_hm3_ldsc.snplist",
              " --chunksize 500000 --N-col N_eff --a2 A2 --a1 A1 --p P --signed-sumstats Z,0 --snp SNP",
              " --maf-min 0.01"))

system(paste0("ldsc.py",
              " --h2 /data/williamsjacr/Aim1_Overall/Intermediate/sum_align.sumstats.gz",
              " --ref-ld-chr /data/williamsjacr/AFR_bd38_MEGA_HM3_Jia_rsid_July2024/",
              " --w-ld-chr /data/williamsjacr/AFR_bd38_MEGA_HM3_Jia_rsid_July2024/",
              " --out /data/williamsjacr/Aim1_Overall/Results/AFR_HM3_AABCG_h2_observed_05"))

rm(list = ls())

library(data.table)

sumstats <- fread("/data/williamsjacr/Aim1_Overall/GWAS_SumStats/AA_Overall_HM3_MEGA_05_2024.txt")
sumstats <- as.data.frame(sumstats)

Z <- as.numeric(sumstats$BETA/sumstats$SE)
P <- as.numeric(sumstats$P)
A1 <- sumstats$non_effect_allele
A2 <- sumstats$effect_allele
N_eff <- as.numeric(sumstats$N_eff)

SNP_ID <- sumstats$SNPID
sumstats <- data.frame(SNP = SNP_ID,A1 = A1,A2 = A2, P = P, Z = Z,N_eff = N_eff)
fwrite(sumstats,file = "/data/williamsjacr/Aim1_Overall/Intermediate/sumstats_Overall",row.names = F,col.names = T,quote = F,sep = " ")

system(paste0("munge_sumstats.py",
              " --sumstats /data/williamsjacr/Aim1_Overall/Intermediate/sumstats_Overall",
              " --out /data/williamsjacr/Aim1_Overall/Intermediate/sum_align",
              " --merge-alleles /data/williamsjacr/1KG_ldsc_HM3_MEGA/AFR/AFR_hm3_mega_ldsc/hm3_mega_ldsc.snplist",
              " --chunksize 500000 --N-col N_eff --a2 A2 --a1 A1 --p P --signed-sumstats Z,0 --snp SNP",
              " --maf-min 0.01"))

system(paste0("ldsc.py",
              " --h2 /data/williamsjacr/Aim1_Overall/Intermediate/sum_align.sumstats.gz",
              " --ref-ld-chr /data/williamsjacr/AFR_bd38_MEGA_HM3_Jia_rsid_July2024/",
              " --w-ld-chr /data/williamsjacr/AFR_bd38_MEGA_HM3_Jia_rsid_July2024/",
              " --out /data/williamsjacr/Aim1_Overall/Results/AFR_HM3_MEGA_AABCG_h2_observed_05"))





















rm(list = ls())

library(data.table)

sumstats <- fread("/data/williamsjacr/Aim1_Overall/GWAS_SumStats/EAS_Overall_HM3_01.txt")
sumstats <- as.data.frame(sumstats)

Z <- as.numeric(sumstats$BETA_meta/sumstats$SE_meta)
P <- as.numeric(sumstats$P_meta)
A1 <- sumstats$non_effect_allele_meta
A2 <- sumstats$effect_allele_meta
N_eff <- as.numeric(sumstats$N_eff_meta)

SNP_ID <- sumstats$SNPID
sumstats <- data.frame(SNP = SNP_ID,A1 = A1,A2 = A2, P = P, Z = Z,N_eff = N_eff)
fwrite(sumstats,file = "/data/williamsjacr/Aim1_Overall/Intermediate/sumstats_Overall",row.names = F,col.names = T,quote = F,sep = " ")

system(paste0("munge_sumstats.py",
              " --sumstats /data/williamsjacr/Aim1_Overall/Intermediate/sumstats_Overall",
              " --out /data/williamsjacr/Aim1_Overall/Intermediate/sum_align",
              " --merge-alleles /data/williamsjacr/1KG_ldsc_HM3/EAS/EAS_hm3_ldsc/w_hm3_ldsc.snplist",
              " --chunksize 500000 --N-col N_eff --a2 A2 --a1 A1 --p P --signed-sumstats Z,0 --snp SNP",
              " --maf-min 0.01"))

system(paste0("ldsc.py",
              " --h2 /data/williamsjacr/Aim1_Overall/Intermediate/sum_align.sumstats.gz",
              " --ref-ld-chr /data/williamsjacr/1KG_ldsc_HM3/EAS/EAS_hm3_ldsc/",
              " --w-ld-chr /data/williamsjacr/1KG_ldsc_HM3/EAS/EAS_hm3_ldsc/",
              " --out /data/williamsjacr/Aim1_Overall/Results/EAS_HM3_h2_observed_01"))

rm(list = ls())

library(data.table)

sumstats <- fread("/data/williamsjacr/Aim1_Overall/GWAS_SumStats/EAS_Overall_HM3_MEGA_01.txt")
sumstats <- as.data.frame(sumstats)

Z <- as.numeric(sumstats$BETA_meta/sumstats$SE_meta)
P <- as.numeric(sumstats$P_meta)
A1 <- sumstats$non_effect_allele_meta
A2 <- sumstats$effect_allele_meta
N_eff <- as.numeric(sumstats$N_eff_meta)

SNP_ID <- sumstats$SNPID
sumstats <- data.frame(SNP = SNP_ID,A1 = A1,A2 = A2, P = P, Z = Z,N_eff = N_eff)
fwrite(sumstats,file = "/data/williamsjacr/Aim1_Overall/Intermediate/sumstats_Overall",row.names = F,col.names = T,quote = F,sep = " ")

system(paste0("munge_sumstats.py",
              " --sumstats /data/williamsjacr/Aim1_Overall/Intermediate/sumstats_Overall",
              " --out /data/williamsjacr/Aim1_Overall/Intermediate/sum_align",
              " --merge-alleles /data/williamsjacr/1KG_ldsc_HM3_MEGA/EAS/EAS_hm3_mega_ldsc/hm3_mega_ldsc.snplist",
              " --chunksize 500000 --N-col N_eff --a2 A2 --a1 A1 --p P --signed-sumstats Z,0 --snp SNP",
              " --maf-min 0.01"))

system(paste0("ldsc.py",
              " --h2 /data/williamsjacr/Aim1_Overall/Intermediate/sum_align.sumstats.gz",
              " --ref-ld-chr /data/williamsjacr/1KG_ldsc_HM3_MEGA/EAS/EAS_hm3_mega_ldsc/",
              " --w-ld-chr /data/williamsjacr/1KG_ldsc_HM3_MEGA/EAS/EAS_hm3_mega_ldsc/",
              " --out /data/williamsjacr/Aim1_Overall/Results/EAS_HM3_MEGA_h2_observed_01"))


rm(list = ls())

library(data.table)

sumstats <- fread("/data/williamsjacr/Aim1_Overall/GWAS_SumStats/EAS_Overall_HM3_05.txt")
sumstats <- as.data.frame(sumstats)

Z <- as.numeric(sumstats$BETA_meta/sumstats$SE_meta)
P <- as.numeric(sumstats$P_meta)
A1 <- sumstats$non_effect_allele_meta
A2 <- sumstats$effect_allele_meta
N_eff <- as.numeric(sumstats$N_eff_meta)

SNP_ID <- sumstats$SNPID
sumstats <- data.frame(SNP = SNP_ID,A1 = A1,A2 = A2, P = P, Z = Z,N_eff = N_eff)
fwrite(sumstats,file = "/data/williamsjacr/Aim1_Overall/Intermediate/sumstats_Overall",row.names = F,col.names = T,quote = F,sep = " ")

system(paste0("munge_sumstats.py",
              " --sumstats /data/williamsjacr/Aim1_Overall/Intermediate/sumstats_Overall",
              " --out /data/williamsjacr/Aim1_Overall/Intermediate/sum_align",
              " --merge-alleles /data/williamsjacr/1KG_ldsc_HM3/EAS/EAS_hm3_ldsc/w_hm3_ldsc.snplist",
              " --chunksize 500000 --N-col N_eff --a2 A2 --a1 A1 --p P --signed-sumstats Z,0 --snp SNP",
              " --maf-min 0.01"))

system(paste0("ldsc.py",
              " --h2 /data/williamsjacr/Aim1_Overall/Intermediate/sum_align.sumstats.gz",
              " --ref-ld-chr /data/williamsjacr/1KG_ldsc_HM3/EAS/EAS_hm3_ldsc/",
              " --w-ld-chr /data/williamsjacr/1KG_ldsc_HM3/EAS/EAS_hm3_ldsc/",
              " --out /data/williamsjacr/Aim1_Overall/Results/EAS_HM3_h2_observed_05"))

rm(list = ls())

library(data.table)

sumstats <- fread("/data/williamsjacr/Aim1_Overall/GWAS_SumStats/EAS_Overall_HM3_MEGA_05.txt")
sumstats <- as.data.frame(sumstats)

Z <- as.numeric(sumstats$BETA_meta/sumstats$SE_meta)
P <- as.numeric(sumstats$P_meta)
A1 <- sumstats$non_effect_allele_meta
A2 <- sumstats$effect_allele_meta
N_eff <- as.numeric(sumstats$N_eff_meta)

SNP_ID <- sumstats$SNPID
sumstats <- data.frame(SNP = SNP_ID,A1 = A1,A2 = A2, P = P, Z = Z,N_eff = N_eff)
fwrite(sumstats,file = "/data/williamsjacr/Aim1_Overall/Intermediate/sumstats_Overall",row.names = F,col.names = T,quote = F,sep = " ")

system(paste0("munge_sumstats.py",
              " --sumstats /data/williamsjacr/Aim1_Overall/Intermediate/sumstats_Overall",
              " --out /data/williamsjacr/Aim1_Overall/Intermediate/sum_align",
              " --merge-alleles /data/williamsjacr/1KG_ldsc_HM3_MEGA/EAS/EAS_hm3_mega_ldsc/hm3_mega_ldsc.snplist",
              " --chunksize 500000 --N-col N_eff --a2 A2 --a1 A1 --p P --signed-sumstats Z,0 --snp SNP",
              " --maf-min 0.01"))

system(paste0("ldsc.py",
              " --h2 /data/williamsjacr/Aim1_Overall/Intermediate/sum_align.sumstats.gz",
              " --ref-ld-chr /data/williamsjacr/1KG_ldsc_HM3_MEGA/EAS/EAS_hm3_mega_ldsc/",
              " --w-ld-chr /data/williamsjacr/1KG_ldsc_HM3_MEGA/EAS/EAS_hm3_mega_ldsc/",
              " --out /data/williamsjacr/Aim1_Overall/Results/EAS_HM3_MEGA_h2_observed_05"))



























rm(list = ls())

library(data.table)

sumstats <- fread("/data/williamsjacr/Aim1_Overall/GWAS_SumStats/Latina_Overall_HM3_01.txt")
sumstats <- as.data.frame(sumstats)

Z <- as.numeric(sumstats$BETA/sumstats$SE)
P <- as.numeric(sumstats$P)
A1 <- sumstats$non_effect_allele
A2 <- sumstats$effect_allele
N_eff <- as.numeric(sumstats$N_eff)

SNP_ID <- sumstats$SNP
sumstats <- data.frame(SNP = SNP_ID,A1 = A1,A2 = A2, P = P, Z = Z,N_eff = N_eff)
fwrite(sumstats,file = "/data/williamsjacr/Aim1_Overall/Intermediate/sumstats_Overall",row.names = F,col.names = T,quote = F,sep = " ")

system(paste0("munge_sumstats.py",
              " --sumstats /data/williamsjacr/Aim1_Overall/Intermediate/sumstats_Overall",
              " --out /data/williamsjacr/Aim1_Overall/Intermediate/sum_align",
              " --merge-alleles /data/williamsjacr/1KG_ldsc_HM3/AMR/AMR_hm3_ldsc/w_hm3_ldsc.snplist",
              " --chunksize 500000 --N-col N_eff --a2 A2 --a1 A1 --p P --signed-sumstats Z,0 --snp SNP",
              " --maf-min 0.01"))

system(paste0("ldsc.py",
              " --h2 /data/williamsjacr/Aim1_Overall/Intermediate/sum_align.sumstats.gz",
              " --ref-ld-chr /data/williamsjacr/1KG_ldsc_HM3/AMR/AMR_hm3_ldsc/",
              " --w-ld-chr /data/williamsjacr/1KG_ldsc_HM3/AMR/AMR_hm3_ldsc/",
              " --out /data/williamsjacr/Aim1_Overall/Results/AMR_HM3_h2_observed_01"))

rm(list = ls())

library(data.table)

sumstats <- fread("/data/williamsjacr/Aim1_Overall/GWAS_SumStats/Latina_Overall_HM3_MEGA_01.txt")
sumstats <- as.data.frame(sumstats)

Z <- as.numeric(sumstats$BETA/sumstats$SE)
P <- as.numeric(sumstats$P)
A1 <- sumstats$non_effect_allele
A2 <- sumstats$effect_allele
N_eff <- as.numeric(sumstats$N_eff)

SNP_ID <- sumstats$SNP
sumstats <- data.frame(SNP = SNP_ID,A1 = A1,A2 = A2, P = P, Z = Z,N_eff = N_eff)
fwrite(sumstats,file = "/data/williamsjacr/Aim1_Overall/Intermediate/sumstats_Overall",row.names = F,col.names = T,quote = F,sep = " ")

system(paste0("munge_sumstats.py",
              " --sumstats /data/williamsjacr/Aim1_Overall/Intermediate/sumstats_Overall",
              " --out /data/williamsjacr/Aim1_Overall/Intermediate/sum_align",
              " --merge-alleles /data/williamsjacr/1KG_ldsc_HM3_MEGA/AMR/AMR_hm3_mega_ldsc/hm3_mega_ldsc.snplist",
              " --chunksize 500000 --N-col N_eff --a2 A2 --a1 A1 --p P --signed-sumstats Z,0 --snp SNP",
              " --maf-min 0.01"))

system(paste0("ldsc.py",
              " --h2 /data/williamsjacr/Aim1_Overall/Intermediate/sum_align.sumstats.gz",
              " --ref-ld-chr /data/williamsjacr/1KG_ldsc_HM3_MEGA/AMR/AMR_hm3_mega_ldsc/",
              " --w-ld-chr /data/williamsjacr/1KG_ldsc_HM3_MEGA/AMR/AMR_hm3_mega_ldsc/",
              " --out /data/williamsjacr/Aim1_Overall/Results/AMR_HM3_MEGA_h2_observed_01"))

rm(list = ls())

library(data.table)

sumstats <- fread("/data/williamsjacr/Aim1_Overall/GWAS_SumStats/Latina_Overall_HM3_05.txt")
sumstats <- as.data.frame(sumstats)

Z <- as.numeric(sumstats$BETA/sumstats$SE)
P <- as.numeric(sumstats$P)
A1 <- sumstats$non_effect_allele
A2 <- sumstats$effect_allele
N_eff <- as.numeric(sumstats$N_eff)

SNP_ID <- sumstats$SNP
sumstats <- data.frame(SNP = SNP_ID,A1 = A1,A2 = A2, P = P, Z = Z,N_eff = N_eff)
fwrite(sumstats,file = "/data/williamsjacr/Aim1_Overall/Intermediate/sumstats_Overall",row.names = F,col.names = T,quote = F,sep = " ")

system(paste0("munge_sumstats.py",
              " --sumstats /data/williamsjacr/Aim1_Overall/Intermediate/sumstats_Overall",
              " --out /data/williamsjacr/Aim1_Overall/Intermediate/sum_align",
              " --merge-alleles /data/williamsjacr/1KG_ldsc_HM3/AMR/AMR_hm3_ldsc/w_hm3_ldsc.snplist",
              " --chunksize 500000 --N-col N_eff --a2 A2 --a1 A1 --p P --signed-sumstats Z,0 --snp SNP",
              " --maf-min 0.01"))

system(paste0("ldsc.py",
              " --h2 /data/williamsjacr/Aim1_Overall/Intermediate/sum_align.sumstats.gz",
              " --ref-ld-chr /data/williamsjacr/1KG_ldsc_HM3/AMR/AMR_hm3_ldsc/",
              " --w-ld-chr /data/williamsjacr/1KG_ldsc_HM3/AMR/AMR_hm3_ldsc/",
              " --out /data/williamsjacr/Aim1_Overall/Results/AMR_HM3_h2_observed_05"))

rm(list = ls())

library(data.table)

sumstats <- fread("/data/williamsjacr/Aim1_Overall/GWAS_SumStats/Latina_Overall_HM3_MEGA_05.txt")
sumstats <- as.data.frame(sumstats)

Z <- as.numeric(sumstats$BETA/sumstats$SE)
P <- as.numeric(sumstats$P)
A1 <- sumstats$non_effect_allele
A2 <- sumstats$effect_allele
N_eff <- as.numeric(sumstats$N_eff)

SNP_ID <- sumstats$SNP
sumstats <- data.frame(SNP = SNP_ID,A1 = A1,A2 = A2, P = P, Z = Z,N_eff = N_eff)
fwrite(sumstats,file = "/data/williamsjacr/Aim1_Overall/Intermediate/sumstats_Overall",row.names = F,col.names = T,quote = F,sep = " ")

system(paste0("munge_sumstats.py",
              " --sumstats /data/williamsjacr/Aim1_Overall/Intermediate/sumstats_Overall",
              " --out /data/williamsjacr/Aim1_Overall/Intermediate/sum_align",
              " --merge-alleles /data/williamsjacr/1KG_ldsc_HM3_MEGA/AMR/AMR_hm3_mega_ldsc/hm3_mega_ldsc.snplist",
              " --chunksize 500000 --N-col N_eff --a2 A2 --a1 A1 --p P --signed-sumstats Z,0 --snp SNP",
              " --maf-min 0.01"))

system(paste0("ldsc.py",
              " --h2 /data/williamsjacr/Aim1_Overall/Intermediate/sum_align.sumstats.gz",
              " --ref-ld-chr /data/williamsjacr/1KG_ldsc_HM3_MEGA/AMR/AMR_hm3_mega_ldsc/",
              " --w-ld-chr /data/williamsjacr/1KG_ldsc_HM3_MEGA/AMR/AMR_hm3_mega_ldsc/",
              " --out /data/williamsjacr/Aim1_Overall/Results/AMR_HM3_MEGA_h2_observed_05"))























rm(list = ls())

library(data.table)

sumstats <- fread("/data/williamsjacr/Aim1_Overall/GWAS_SumStats/EUR_Overall_HM3_01.txt")
sumstats <- as.data.frame(sumstats)

Z <- as.numeric(sumstats$BETA_meta/sumstats$SE_meta)
P <- as.numeric(sumstats$P_meta)
A1 <- sumstats$non_effect_allele_meta
A2 <- sumstats$effect_allele_meta
N_eff <- as.numeric(sumstats$N_eff_meta)

SNP_ID <- sumstats$SNP_ID
sumstats <- data.frame(SNP = SNP_ID,A1 = A1,A2 = A2, P = P, Z = Z,N_eff = N_eff)
fwrite(sumstats,file = "/data/williamsjacr/Aim1_Overall/Intermediate/sumstats_Overall",row.names = F,col.names = T,quote = F,sep = " ")

system(paste0("munge_sumstats.py",
              " --sumstats /data/williamsjacr/Aim1_Overall/Intermediate/sumstats_Overall",
              " --out /data/williamsjacr/Aim1_Overall/Intermediate/sum_align",
              " --merge-alleles /data/williamsjacr/1KG_ldsc_HM3/EUR/EUR_hm3_ldsc/w_hm3_ldsc.snplist",
              " --chunksize 500000 --N-col N_eff --a2 A2 --a1 A1 --p P --signed-sumstats Z,0 --snp SNP",
              " --maf-min 0.01"))

system(paste0("ldsc.py",
              " --h2 /data/williamsjacr/Aim1_Overall/Intermediate/sum_align.sumstats.gz",
              " --ref-ld-chr /data/williamsjacr/1KG_ldsc_HM3/EUR/EUR_hm3_ldsc/",
              " --w-ld-chr /data/williamsjacr/1KG_ldsc_HM3/EUR/EUR_hm3_ldsc/",
              " --out /data/williamsjacr/Aim1_Overall/Results/EUR_HM3_h2_observed_01"))

rm(list = ls())

library(data.table)

sumstats <- fread("/data/williamsjacr/Aim1_Overall/GWAS_SumStats/EUR_Overall_HM3_MEGA_01.txt")
sumstats <- as.data.frame(sumstats)

Z <- as.numeric(sumstats$BETA_meta/sumstats$SE_meta)
P <- as.numeric(sumstats$P_meta)
A1 <- sumstats$non_effect_allele_meta
A2 <- sumstats$effect_allele_meta
N_eff <- as.numeric(sumstats$N_eff_meta)

SNP_ID <- sumstats$SNP_ID
sumstats <- data.frame(SNP = SNP_ID,A1 = A1,A2 = A2, P = P, Z = Z,N_eff = N_eff)
fwrite(sumstats,file = "/data/williamsjacr/Aim1_Overall/Intermediate/sumstats_Overall",row.names = F,col.names = T,quote = F,sep = " ")

system(paste0("munge_sumstats.py",
              " --sumstats /data/williamsjacr/Aim1_Overall/Intermediate/sumstats_Overall",
              " --out /data/williamsjacr/Aim1_Overall/Intermediate/sum_align",
              " --merge-alleles /data/williamsjacr/1KG_ldsc_HM3_MEGA/EUR/EUR_hm3_mega_ldsc/hm3_mega_ldsc.snplist",
              " --chunksize 500000 --N-col N_eff --a2 A2 --a1 A1 --p P --signed-sumstats Z,0 --snp SNP",
              " --maf-min 0.01"))

system(paste0("ldsc.py",
              " --h2 /data/williamsjacr/Aim1_Overall/Intermediate/sum_align.sumstats.gz",
              " --ref-ld-chr /data/williamsjacr/1KG_ldsc_HM3_MEGA/EUR/EUR_hm3_mega_ldsc/",
              " --w-ld-chr /data/williamsjacr/1KG_ldsc_HM3_MEGA/EUR/EUR_hm3_mega_ldsc/",
              " --out /data/williamsjacr/Aim1_Overall/Results/EUR_HM3_MEGA_h2_observed_01"))


rm(list = ls())

library(data.table)

sumstats <- fread("/data/williamsjacr/Aim1_Overall/GWAS_SumStats/EUR_Overall_HM3_05.txt")
sumstats <- as.data.frame(sumstats)

Z <- as.numeric(sumstats$BETA_meta/sumstats$SE_meta)
P <- as.numeric(sumstats$P_meta)
A1 <- sumstats$non_effect_allele_meta
A2 <- sumstats$effect_allele_meta
N_eff <- as.numeric(sumstats$N_eff_meta)

SNP_ID <- sumstats$SNP_ID
sumstats <- data.frame(SNP = SNP_ID,A1 = A1,A2 = A2, P = P, Z = Z,N_eff = N_eff)
fwrite(sumstats,file = "/data/williamsjacr/Aim1_Overall/Intermediate/sumstats_Overall",row.names = F,col.names = T,quote = F,sep = " ")

system(paste0("munge_sumstats.py",
              " --sumstats /data/williamsjacr/Aim1_Overall/Intermediate/sumstats_Overall",
              " --out /data/williamsjacr/Aim1_Overall/Intermediate/sum_align",
              " --merge-alleles /data/williamsjacr/1KG_ldsc_HM3/EUR/EUR_hm3_ldsc/w_hm3_ldsc.snplist",
              " --chunksize 500000 --N-col N_eff --a2 A2 --a1 A1 --p P --signed-sumstats Z,0 --snp SNP",
              " --maf-min 0.01"))

system(paste0("ldsc.py",
              " --h2 /data/williamsjacr/Aim1_Overall/Intermediate/sum_align.sumstats.gz",
              " --ref-ld-chr /data/williamsjacr/1KG_ldsc_HM3/EUR/EUR_hm3_ldsc/",
              " --w-ld-chr /data/williamsjacr/1KG_ldsc_HM3/EUR/EUR_hm3_ldsc/",
              " --out /data/williamsjacr/Aim1_Overall/Results/EUR_HM3_h2_observed_05"))

rm(list = ls())

library(data.table)

sumstats <- fread("/data/williamsjacr/Aim1_Overall/GWAS_SumStats/EUR_Overall_HM3_MEGA_05.txt")
sumstats <- as.data.frame(sumstats)

Z <- as.numeric(sumstats$BETA_meta/sumstats$SE_meta)
P <- as.numeric(sumstats$P_meta)
A1 <- sumstats$non_effect_allele_meta
A2 <- sumstats$effect_allele_meta
N_eff <- as.numeric(sumstats$N_eff_meta)

SNP_ID <- sumstats$SNP_ID
sumstats <- data.frame(SNP = SNP_ID,A1 = A1,A2 = A2, P = P, Z = Z,N_eff = N_eff)
fwrite(sumstats,file = "/data/williamsjacr/Aim1_Overall/Intermediate/sumstats_Overall",row.names = F,col.names = T,quote = F,sep = " ")

system(paste0("munge_sumstats.py",
              " --sumstats /data/williamsjacr/Aim1_Overall/Intermediate/sumstats_Overall",
              " --out /data/williamsjacr/Aim1_Overall/Intermediate/sum_align",
              " --merge-alleles /data/williamsjacr/1KG_ldsc_HM3_MEGA/EUR/EUR_hm3_mega_ldsc/hm3_mega_ldsc.snplist",
              " --chunksize 500000 --N-col N_eff --a2 A2 --a1 A1 --p P --signed-sumstats Z,0 --snp SNP",
              " --maf-min 0.01"))

system(paste0("ldsc.py",
              " --h2 /data/williamsjacr/Aim1_Overall/Intermediate/sum_align.sumstats.gz",
              " --ref-ld-chr /data/williamsjacr/1KG_ldsc_HM3_MEGA/EUR/EUR_hm3_mega_ldsc/",
              " --w-ld-chr /data/williamsjacr/1KG_ldsc_HM3_MEGA/EUR/EUR_hm3_mega_ldsc/",
              " --out /data/williamsjacr/Aim1_Overall/Results/EUR_HM3_MEGA_h2_observed_05"))
