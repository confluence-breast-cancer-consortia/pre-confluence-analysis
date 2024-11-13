rm(list = ls())

# all files retrieved from https://zenodo.org/records/10515792 version 4

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
fwrite(sumstats,file = "/data/williamsjacr/Aim4_PartitionedHeritability/Intermediate_Files/sumstats_Overall",row.names = F,col.names = T,quote = F,sep = " ")

#---------Munge Sumstats file

system(paste0("munge_sumstats.py", 
              " --sumstats /data/williamsjacr/Aim4_PartitionedHeritability/Intermediate_Files/sumstats_Overall",
              " --out /data/williamsjacr/Aim4_PartitionedHeritability/Intermediate_Files/sum_align",
              " --merge-alleles /data/williamsjacr/1KG_ldsc_HM3/AFR/AFR_hm3_ldsc/w_hm3_ldsc.snplist",
              " --chunksize 500000 --N-col N_eff --a2 A2 --a1 A1 --p P --signed-sumstats Z,0 --snp SNP",
              " --maf-min 0.01"))

#---------ldsc h2 script

system(paste0("ldsc.py",
              " --h2 /data/williamsjacr/Aim4_PartitionedHeritability/Intermediate_Files/sum_align.sumstats.gz",
              " --w-ld-chr /data/williamsjacr/1KG_ldsc_HM3/AFR/AFR_hm3_ldsc/",
              " --ref-ld-chr /data/williamsjacr/1KG_Partitioned_LDSC_HM3/AFR/AFR_hm3_ldsc/baselineLD.", 
              " --frqfile-chr /data/williamsjacr/1KG_Partitioned_LDSC_HM3/AFR/chr_clean",  
              " --overlap-annot",
              " --out /data/williamsjacr/Aim4_PartitionedHeritability/Results/AFR_05_HM3_My_LDSC",  
              " --print-cov --print-coefficients --print-delete-vals"))

system(paste0("ldsc.py",
              " --h2 /data/williamsjacr/Aim4_PartitionedHeritability/Intermediate_Files/sum_align.sumstats.gz",
              " --w-ld-chr /data/williamsjacr/AFR_bd38_MEGA_HM3_Jia_rsid_July2024/",
              " --ref-ld-chr /data/williamsjacr/AABCG_baseline_LD/baselineLD.", 
              " --frqfile-chr /data/williamsjacr/AABCG_MAF/chr_clean",  
              " --overlap-annot",
              " --out /data/williamsjacr/Aim4_PartitionedHeritability/Results/AFR_05_HM3_AABCG_LDSC",  
              " --print-cov --print-coefficients --print-delete-vals"))