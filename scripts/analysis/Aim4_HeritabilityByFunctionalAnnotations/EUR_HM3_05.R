rm(list = ls())

# all files retrieved from https://zenodo.org/records/10515792 version 4

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
fwrite(sumstats,file = "/data/williamsjacr/Aim4_PartitionedHeritability/Intermediate_Files/sumstats_Overall",row.names = F,col.names = T,quote = F,sep = " ")

#---------Munge Sumstats file

system(paste0("munge_sumstats.py", 
              " --sumstats /data/williamsjacr/Aim4_PartitionedHeritability/Intermediate_Files/sumstats_Overall",
              " --out /data/williamsjacr/Aim4_PartitionedHeritability/Intermediate_Files/sum_align",
              " --merge-alleles /data/williamsjacr/1KG_ldsc_HM3/EUR/EUR_hm3_ldsc/w_hm3_ldsc.snplist",
              " --chunksize 500000 --N-col N_eff --a2 A2 --a1 A1 --p P --signed-sumstats Z,0 --snp SNP",
              " --maf-min 0.01"))

#---------ldsc h2 script

system(paste0("ldsc.py",
              " --h2 /data/williamsjacr/Aim4_PartitionedHeritability/Intermediate_Files/sum_align.sumstats.gz",
              " --w-ld-chr /data/williamsjacr/1KG_ldsc_HM3/EUR/EUR_hm3_ldsc/",
              " --ref-ld-chr /data/williamsjacr/Aim4_PartitionedHeritability/LDSC_Files/EUR_Baseline/baselineLD.", 
              " --frqfile-chr /data/williamsjacr/Aim4_PartitionedHeritability/LDSC_Files/1000G_Phase3_frq/1000G.EUR.QC.",  
              " --overlap-annot",
              " --out /data/williamsjacr/Aim4_PartitionedHeritability/Results/EUR_05_HM3_Their_LDSC",  
              " --print-cov --print-coefficients --print-delete-vals"))

system(paste0("ldsc.py",
              " --h2 /data/williamsjacr/Aim4_PartitionedHeritability/Intermediate_Files/sum_align.sumstats.gz",
              " --w-ld-chr /data/williamsjacr/1KG_ldsc_HM3/EUR/EUR_hm3_ldsc/",
              " --ref-ld-chr /data/williamsjacr/1KG_Partitioned_LDSC_HM3/EUR/EUR_hm3_ldsc/baselineLD.", 
              " --frqfile-chr /data/williamsjacr/1KG_Partitioned_LDSC_HM3/EUR/chr_clean",  
              " --overlap-annot",
              " --out /data/williamsjacr/Aim4_PartitionedHeritability/Results/EUR_05_HM3_My_LDSC",  
              " --print-cov --print-coefficients --print-delete-vals"))

system(paste0("ldsc.py",
              " --h2 /data/williamsjacr/Aim4_PartitionedHeritability/Intermediate_Files/sum_align.sumstats.gz",
              " --w-ld-chr /data/williamsjacr/1KG_ldsc_w_ld_chr_HM3/EUR/",
              " --ref-ld-chr /data/williamsjacr/1KG_Partitioned_LDSC_HM3/EUR/EUR_hm3_ldsc/baselineLD.", 
              " --frqfile-chr /data/williamsjacr/1KG_Partitioned_LDSC_HM3/EUR/chr_clean",  
              " --overlap-annot",
              " --out /data/williamsjacr/Aim4_PartitionedHeritability/Results/EUR_05_HM3_w_ld_chr",  
              " --print-cov --print-coefficients --print-delete-vals"))
