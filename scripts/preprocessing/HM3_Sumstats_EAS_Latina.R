# dx run app-swiss-army-knife -iin=PreConfluenceAnalysis:Aim1_Heritability/JacobWilliams/Clean_Data_Scripts/r_with_plink.tar.gz -iin=PreConfluenceAnalysis:Aim1_Heritability/JacobWilliams/Clean_Data_Scripts/HM3_Sumstats_EAS_Latina.R -iin=PreConfluenceAnalysis:Aim1_Heritability/JacobWilliams/Clean_Data_Scripts/HM3_Sumstats_EAS_Latina.sh -iin=PreConfluenceData:BCAC_Asian/bcac_asian_icogs_onco_erneg_meta1.txt.gz -iin=PreConfluenceData:Latina_Summary_Results/Result_merged_ucsf_kp_mec_cama_allchr.txt -iin=PreConfluenceData:BCAC_Asian/bcac_asian_icogs_onco_erpos_meta1.txt.gz -icmd="bash HM3_Sumstats_EAS_Latina.sh" -y --destination PreConfluenceAnalysis:Aim2_Polygenicity/HZ/Results/Clean_summary_data_update/ --instance-type mem1_ssd1_v2_x36

rm(list = ls())
## Load HM3 SNPs
w_hm3 <- read.csv("w_hm3_updated.snplist")

unique_id <- paste0(w_hm3$chr,"_",w_hm3$pos,"_",
                    w_hm3$allele1,"_",w_hm3$allele2)

merge_data <- data.frame(rsid = w_hm3$rsid, unique_id = unique_id)

# SNP IDs
library(data.table)
library(stringr)
library(dplyr)

rm(list=setdiff(ls(), c("w_hm3","merge_data")))

# bcac_asian erneg
bcac_erneg <- as.data.frame(fread("bcac_asian_icogs_onco_erneg_meta1.txt"))

bcac_erneg_update <- bcac_erneg %>% 
  mutate(unique_allele1 = toupper(ifelse(Allele1 < Allele2, Allele1,Allele2)),
         unique_allele2 = toupper(ifelse(Allele1 < Allele2, Allele2,Allele1)),
         chr_pos = sub("^([^_]+_[^_]+)_.*", "\\1", MarkerName),
         unique_SNP_id = paste0(chr_pos,"_",unique_allele1,"_",unique_allele2)) %>% 
  mutate(effect_allele = toupper(Allele1),
         non_effect_allele = toupper(Allele2),
         Freq_effect = Freq1,
         BETA = Effect,
         SE = StdErr,
         P = `P-value`,
         N_eff = 1/(2*Freq_effect*(1-Freq_effect)*SE^2)) %>% 
  select(unique_SNP_id, effect_allele, non_effect_allele, Freq_effect, BETA, SE, P, N_eff)

names_data <- colnames(bcac_erneg_update)

bcac_erneg_update$CHR_POS <- unlist(lapply(strsplit(bcac_erneg_update$unique_SNP_id,"_"),function(x){ifelse(sum(is.na(x)) == 1,"40_1",paste0(x[1],"_",x[2]))}))

bcac_erneg_update$unique_id1 <- paste0(bcac_erneg_update$CHR_POS,"_",
                                       bcac_erneg_update$effect_allele,"_",bcac_erneg_update$non_effect_allele)
bcac_erneg_update$unique_id2 <- paste0(bcac_erneg_update$CHR_POS,"_",
                                       bcac_erneg_update$non_effect_allele,"_",bcac_erneg_update$effect_allele)

bcac_erneg_update <- left_join(bcac_erneg_update,merge_data,by = c("unique_id1" = "unique_id"))
bcac_erneg_update$unique_SNP_id[!is.na(bcac_erneg_update$rsid)] <- bcac_erneg_update$rsid[!is.na(bcac_erneg_update$rsid)]

bcac_erneg_update <- subset(bcac_erneg_update,select = -c(rsid))

bcac_erneg_update <- left_join(bcac_erneg_update,merge_data,by = c("unique_id2" = "unique_id"))
bcac_erneg_update$unique_SNP_id[!is.na(bcac_erneg_update$rsid)] <- bcac_erneg_update$rsid[!is.na(bcac_erneg_update$rsid)]

bcac_erneg_update_hm3 <-  bcac_erneg_update %>% 
  filter(Freq_effect >= 0.01 & Freq_effect <= 0.99) %>%
  filter(unique_SNP_id %in% w_hm3$rsid) 

N_eff <- bcac_erneg_update_hm3$N_eff
N_cut_low <- quantile(N_eff,probs = 0.1,na.rm = TRUE)

bcac_erneg_sub_quan <-  bcac_erneg_update %>% 
  filter(Freq_effect >= 0.01 & Freq_effect <= 0.99) %>%
  filter(unique_SNP_id %in% w_hm3$rsid) %>% 
  filter(N_eff > N_cut_low)

bcac_erneg_sub_quan <- bcac_erneg_sub_quan[,names_data]

write.table(bcac_erneg_sub_quan, file = "EAS_ERNeg_HM3_01.txt", row.names = F, col.names = T, quote = F)

bcac_erneg_update_hm3 <-  bcac_erneg_update %>% 
  filter(Freq_effect >= 0.05 & Freq_effect <= 0.95) %>%
  filter(unique_SNP_id %in% w_hm3$rsid) 

N_eff <- bcac_erneg_update_hm3$N_eff
N_cut_low <- quantile(N_eff,probs = 0.1,na.rm = TRUE)

bcac_erneg_sub_quan <-  bcac_erneg_update %>% 
  filter(Freq_effect >= 0.05 & Freq_effect <= 0.95) %>%
  filter(unique_SNP_id %in% w_hm3$rsid) %>% 
  filter(N_eff > N_cut_low)

bcac_erneg_sub_quan <- bcac_erneg_sub_quan[,names_data]

write.table(bcac_erneg_sub_quan, file = "EAS_ERNeg_HM3_05.txt", row.names = F, col.names = T, quote = F)

rm(list=setdiff(ls(), c("w_hm3","merge_data")))

# bcac_asian erpos
bcac_erpos <- as.data.frame(fread("bcac_asian_icogs_onco_erpos_meta1.txt"))

bcac_erpos_update <- bcac_erpos %>% 
  mutate(unique_allele1 = toupper(ifelse(Allele1 < Allele2, Allele1,Allele2)),
         unique_allele2 = toupper(ifelse(Allele1 < Allele2, Allele2,Allele1)),
         chr_pos = sub("^([^_]+_[^_]+)_.*", "\\1", MarkerName),
         unique_SNP_id = paste0(chr_pos,"_",unique_allele1,"_",unique_allele2)) %>% 
  mutate(effect_allele = toupper(Allele1),
         non_effect_allele = toupper(Allele2),
         Freq_effect = Freq1,
         BETA = Effect,
         SE = StdErr,
         P = `P-value`,
         N_eff = 1/(2*Freq_effect*(1-Freq_effect)*SE^2)) %>% 
  select(unique_SNP_id, effect_allele, non_effect_allele, Freq_effect, BETA, SE, P, N_eff)

names_data <- colnames(bcac_erpos_update)

bcac_erpos_update$CHR_POS <- unlist(lapply(strsplit(bcac_erpos_update$unique_SNP_id,"_"),function(x){ifelse(sum(is.na(x)) == 1,"40_1",paste0(x[1],"_",x[2]))}))

bcac_erpos_update$unique_id1 <- paste0(bcac_erpos_update$CHR_POS,"_",
                                       bcac_erpos_update$effect_allele,"_",bcac_erpos_update$non_effect_allele)
bcac_erpos_update$unique_id2 <- paste0(bcac_erpos_update$CHR_POS,"_",
                                       bcac_erpos_update$non_effect_allele,"_",bcac_erpos_update$effect_allele)

bcac_erpos_update <- left_join(bcac_erpos_update,merge_data,by = c("unique_id1" = "unique_id"))
bcac_erpos_update$unique_SNP_id[!is.na(bcac_erpos_update$rsid)] <- bcac_erpos_update$rsid[!is.na(bcac_erpos_update$rsid)]

bcac_erpos_update <- subset(bcac_erpos_update,select = -c(rsid))

bcac_erpos_update <- left_join(bcac_erpos_update,merge_data,by = c("unique_id2" = "unique_id"))
bcac_erpos_update$unique_SNP_id[!is.na(bcac_erpos_update$rsid)] <- bcac_erpos_update$rsid[!is.na(bcac_erpos_update$rsid)]

bcac_erpos_update_hm3 <-  bcac_erpos_update %>% 
  filter(Freq_effect >= 0.01 & Freq_effect <= 0.99) %>%
  filter(unique_SNP_id %in% w_hm3$rsid) 

N_eff <- bcac_erpos_update_hm3$N_eff
N_cut_low <- quantile(N_eff,probs = 0.1,na.rm = TRUE)

bcac_erpos_sub_quan <-  bcac_erpos_update %>% 
  filter(Freq_effect >= 0.01 & Freq_effect <= 0.99) %>%
  filter(unique_SNP_id %in% w_hm3$rsid) %>% 
  filter(N_eff > N_cut_low)

bcac_erpos_sub_quan <- bcac_erpos_sub_quan[,names_data]

write.table(bcac_erpos_sub_quan, file = "EAS_ERPos_HM3_01.txt", row.names = F, col.names = T, quote = F)

bcac_erpos_update_hm3 <-  bcac_erpos_update %>% 
  filter(Freq_effect >= 0.05 & Freq_effect <= 0.95) %>%
  filter(unique_SNP_id %in% w_hm3$rsid) 

N_eff <- bcac_erpos_update_hm3$N_eff
N_cut_low <- quantile(N_eff,probs = 0.1,na.rm = TRUE)

bcac_erpos_sub_quan <-  bcac_erpos_update %>% 
  filter(Freq_effect >= 0.05 & Freq_effect <= 0.95) %>%
  filter(unique_SNP_id %in% w_hm3$rsid) %>% 
  filter(N_eff > N_cut_low)

bcac_erpos_sub_quan <- bcac_erpos_sub_quan[,names_data]

write.table(bcac_erpos_sub_quan, file = "EAS_ERPos_HM3_05.txt", row.names = F, col.names = T, quote = F)

rm(list=setdiff(ls(), c("w_hm3","merge_data")))

# latina
latina_overall <- as.data.frame(fread("Result_merged_ucsf_kp_mec_cama_allchr.txt"))
latina_overall_update <- latina_overall %>% 
  mutate(unique_allele1 = ifelse(A2 < A1, A2,A1),
         unique_allele2 = ifelse(A2 < A1, A1,A2),
         chr_pos = paste0(CHR,"_",BP),
         unique_SNP_id = paste0(chr_pos,"_",unique_allele1,"_",unique_allele2)) %>% 
  mutate(effect_allele = toupper(A1),
         non_effect_allele = toupper(A2),
         #for effect sample size calculation, using the Allele frequency in controls
         Freq_effect = Freq,
         P = P,
         BETA = log(OR),
         N_eff = as.integer(1/(2*Freq_effect*(1-Freq_effect)*SE^2))) %>% 
  select(unique_SNP_id, effect_allele, non_effect_allele, Freq_effect, BETA, SE, P, N_eff, SNP)

names_data <- colnames(latina_overall_update)

latina_overall_update$CHR_POS <- unlist(lapply(strsplit(latina_overall_update$unique_SNP_id,"_"),function(x){ifelse(sum(is.na(x)) == 1,"40_1",paste0(x[1],"_",x[2]))}))

latina_overall_update$unique_id1 <- paste0(latina_overall_update$CHR_POS,"_",
                                           latina_overall_update$effect_allele,"_",latina_overall_update$non_effect_allele)
latina_overall_update$unique_id2 <- paste0(latina_overall_update$CHR_POS,"_",
                                           latina_overall_update$non_effect_allele,"_",latina_overall_update$effect_allele)

latina_overall_update <- left_join(latina_overall_update,merge_data,by = c("unique_id1" = "unique_id"))
latina_overall_update$SNP[!is.na(latina_overall_update$rsid)] <- latina_overall_update$rsid[!is.na(latina_overall_update$rsid)]

latina_overall_update <- subset(latina_overall_update,select = -c(rsid))

latina_overall_update <- left_join(latina_overall_update,merge_data,by = c("unique_id2" = "unique_id"))
latina_overall_update$SNP[!is.na(latina_overall_update$rsid)] <- latina_overall_update$rsid[!is.na(latina_overall_update$rsid)]

latina_overall_update_hm3 <-  latina_overall_update %>%
  filter(Freq_effect >= 0.01 & Freq_effect <= 0.99) %>%
  filter(SNP %in% w_hm3$rsid) 

N_eff <- latina_overall_update_hm3$N_eff
N_cut_low <- quantile(N_eff,probs = 0.1,na.rm = TRUE)

latina_overall_sub_quan <-  latina_overall_update %>%
  filter(Freq_effect >= 0.01 & Freq_effect <= 0.99) %>%
  filter(SNP %in% w_hm3$rsid) %>% 
  filter(N_eff > N_cut_low)

latina_overall_sub_quan <- latina_overall_sub_quan[,names_data]

print(nrow(latina_overall_sub_quan))

write.table(latina_overall_sub_quan, file = "Latina_Overall_HM3_01.txt", row.names = F, col.names = T, quote = F)

latina_overall_update_hm3 <-  latina_overall_update %>%
  filter(Freq_effect >= 0.05 & Freq_effect <= 0.95) %>%
  filter(SNP %in% w_hm3$rsid) 

N_eff <- latina_overall_update_hm3$N_eff
N_cut_low <- quantile(N_eff,probs = 0.1,na.rm = TRUE)

latina_overall_sub_quan <-  latina_overall_update %>%
  filter(Freq_effect >= 0.05 & Freq_effect <= 0.95) %>%
  filter(SNP %in% w_hm3$rsid) %>% 
  filter(N_eff > N_cut_low)

latina_overall_sub_quan <- latina_overall_sub_quan[,names_data]

print(nrow(latina_overall_sub_quan))

write.table(latina_overall_sub_quan, file = "Latina_Overall_HM3_05.txt", row.names = F, col.names = T, quote = F)

rm(list=setdiff(ls(), c("w_hm3","merge_data")))

system("rm bcac_asian_icogs_onco_erneg_meta1.txt")
system("rm bcac_asian_icogs_onco_erpos_meta1.txt")
system("rm Result_merged_ucsf_kp_mec_cama_allchr.txt")
system("rm w_hm3_updated.snplist")