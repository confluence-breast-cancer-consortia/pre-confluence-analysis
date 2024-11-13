rm(list = ls())

# dx run app-swiss-army-knife -iin=PreConfluenceAnalysis:Aim1_Heritability/JacobWilliams/Clean_Data_Scripts/r_with_plink.tar.gz -iin=PreConfluenceAnalysis:Aim1_Heritability/JacobWilliams/Clean_Data_Scripts/HM3_Sumstats_EUR_EAS.R -iin=PreConfluenceAnalysis:Aim1_Heritability/JacobWilliams/Clean_Data_Scripts/HM3_Sumstats_EUR_EAS.sh -icmd="bash HM3_Sumstats_EUR_EAS.sh" -y --destination PreConfluenceAnalysis:Aim2_Polygenicity/HZ/Results/Clean_summary_data_update/ --instance-type mem1_ssd1_v2_x36

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

# EAS
EAS_data <- fread("Clean_summary_data/EAS_BCAC_BBJ_meta_sumdata.txt")
EAS_data <- as.data.frame(EAS_data)
names_data <- colnames(EAS_data)

EAS_data$CHR_POS <- unlist(lapply(strsplit(EAS_data$unique_SNP_id,"_"),function(x){ifelse(sum(is.na(x)) == 1,"40_1",paste0(x[1],"_",x[2]))}))

EAS_data$unique_id1 <- paste0(EAS_data$CHR_POS,"_",
                              EAS_data$effect_allele_meta,"_",EAS_data$non_effect_allele_meta)
EAS_data$unique_id2 <- paste0(EAS_data$CHR_POS,"_",
                              EAS_data$non_effect_allele_meta,"_",EAS_data$effect_allele_meta)

EAS_data <- left_join(EAS_data,merge_data,by = c("unique_id1" = "unique_id"))
EAS_data$SNPID[!is.na(EAS_data$rsid)] <- EAS_data$rsid[!is.na(EAS_data$rsid)]

EAS_data <- subset(EAS_data,select = -c(rsid))

EAS_data <- left_join(EAS_data,merge_data,by = c("unique_id2" = "unique_id"))
EAS_data$SNPID[!is.na(EAS_data$rsid)] <- EAS_data$rsid[!is.na(EAS_data$rsid)]

EAS_data_sub <- EAS_data %>% 
  filter(!is.na(effect_allele_BCAC) & !is.na(effect_allele_BBJ)) %>% 
  filter(Freq_effect_BCAC >= 0.01 & Freq_effect_BCAC <= 0.99) %>%
  filter(SNPID %in% w_hm3$rsid)

N_eff <- EAS_data_sub$N_eff_meta
N_cut_low <- quantile(N_eff,probs = 0.1,na.rm = TRUE)

EAS_data_sub_quan <-  EAS_data %>% 
  filter(!is.na(effect_allele_BCAC) & !is.na(effect_allele_BBJ)) %>% 
  filter(Freq_effect_BCAC >= 0.01 & Freq_effect_BCAC <= 0.99) %>%
  filter(SNPID %in% w_hm3$rsid) %>% 
  filter(N_eff_meta > N_cut_low)

EAS_data_sub_quan <- EAS_data_sub_quan[,names_data]

write.table(EAS_data_sub_quan, file = "EAS_Overall_HM3_01.txt", row.names = F, col.names = T, quote = F)

EAS_data_sub <- EAS_data %>% 
  filter(!is.na(effect_allele_BCAC) & !is.na(effect_allele_BBJ)) %>% 
  filter(Freq_effect_BCAC >= 0.05 & Freq_effect_BCAC <= 0.95) %>%
  filter(SNPID %in% w_hm3$rsid)

N_eff <- EAS_data_sub$N_eff_meta
N_cut_low <- quantile(N_eff,probs = 0.1,na.rm = TRUE)

EAS_data_sub_quan <-  EAS_data %>% 
  filter(!is.na(effect_allele_BCAC) & !is.na(effect_allele_BBJ)) %>% 
  filter(Freq_effect_BCAC >= 0.05 & Freq_effect_BCAC <= 0.95) %>%
  filter(SNPID %in% w_hm3$rsid) %>% 
  filter(N_eff_meta > N_cut_low)

EAS_data_sub_quan <- EAS_data_sub_quan[,names_data]

write.table(EAS_data_sub_quan, file = "EAS_Overall_HM3_05.txt", row.names = F, col.names = T, quote = F)

rm(list=setdiff(ls(), c("w_hm3","merge_data")))

# # AA Overall
# AA_overall <- fread("Clean_summary_data/AA_overall_sumdata.txt")
# AA_overall <- as.data.frame(AA_overall)
# names_data <- colnames(AA_overall)
# 
# AA_overall$unique_id1 <- paste0(AA_overall$CHR,"_",AA_overall$POS,"_",
#                               AA_overall$effect_allele,"_",AA_overall$non_effect_allele)
# AA_overall$unique_id2 <- paste0(AA_overall$CHR,"_",AA_overall$POS,"_",
#                               AA_overall$non_effect_allele,"_",AA_overall$effect_allele)
# 
# AA_overall <- left_join(AA_overall,merge_data,by = c("unique_id1" = "unique_id"))
# AA_overall$ID[!is.na(AA_overall$rsid)] <- AA_overall$rsid[!is.na(AA_overall$rsid)]
# 
# AA_overall <- subset(AA_overall,select = -c(rsid))
# 
# AA_overall <- left_join(AA_overall,merge_data,by = c("unique_id2" = "unique_id"))
# AA_overall$ID[!is.na(AA_overall$rsid)] <- AA_overall$rsid[!is.na(AA_overall$rsid)]
# 
# AA_overall_update <-  AA_overall %>% 
#   filter(Freq_effect >= 0.01 & Freq_effect <= 0.99) %>%
#   filter(ID %in% w_hm3$rsid) 
# 
# N_eff <- AA_overall_update$N_eff
# N_cut_low <- quantile(N_eff,probs = 0.1,na.rm = TRUE)
# 
# AA_overall_sub_quan <-  AA_overall %>% 
#   filter(Freq_effect >= 0.01 & Freq_effect <= 0.99) %>%
#   filter(ID %in% w_hm3$rsid) %>% 
#   filter(N_eff > N_cut_low)
# 
# AA_overall_sub_quan <- AA_overall_sub_quan[,names_data]
# 
# write.table(AA_overall_sub_quan, file = "AA_overall_sumdata_01.txt", row.names = F, col.names = T, quote = F)
# 
# AA_overall_update <-  AA_overall %>% 
#   filter(Freq_effect >= 0.05 & Freq_effect <= 0.95) %>%
#   filter(ID %in% w_hm3$rsid) 
# 
# N_eff <- AA_overall_update$N_eff
# N_cut_low <- quantile(N_eff,probs = 0.1,na.rm = TRUE)
# 
# AA_overall_sub_quan <-  AA_overall %>% 
#   filter(Freq_effect >= 0.05 & Freq_effect <= 0.95) %>%
#   filter(ID %in% w_hm3$rsid) %>% 
#   filter(N_eff > N_cut_low)
# 
# AA_overall_sub_quan <- AA_overall_sub_quan[,names_data]
# 
# write.table(AA_overall_sub_quan, file = "AA_overall_sumdata_05.txt", row.names = F, col.names = T, quote = F)
# 
# rm(list=setdiff(ls(), c("w_hm3","merge_data")))
# 
# # AA_ERpos
# AA_ERpos <- fread("Clean_summary_data/AA_ERpos_sumdata.txt")
# AA_ERpos <- as.data.frame(AA_ERpos)
# names_data <- colnames(AA_ERpos)
# 
# AA_ERpos$unique_id1 <- paste0(AA_ERpos$CHR,"_",AA_ERpos$POS,"_",
#                               AA_ERpos$effect_allele,"_",AA_ERpos$non_effect_allele)
# AA_ERpos$unique_id2 <- paste0(AA_ERpos$CHR,"_",AA_ERpos$POS,"_",
#                               AA_ERpos$non_effect_allele,"_",AA_ERpos$effect_allele)
# 
# AA_ERpos <- left_join(AA_ERpos,merge_data,by = c("unique_id1" = "unique_id"))
# AA_ERpos$ID[!is.na(AA_ERpos$rsid)] <- AA_ERpos$rsid[!is.na(AA_ERpos$rsid)]
# 
# AA_ERpos <- subset(AA_ERpos,select = -c(rsid))
# 
# AA_ERpos <- left_join(AA_ERpos,merge_data,by = c("unique_id2" = "unique_id"))
# AA_ERpos$ID[!is.na(AA_ERpos$rsid)] <- AA_ERpos$rsid[!is.na(AA_ERpos$rsid)]
# 
# AA_ERpos_update <-  AA_ERpos %>%
#   filter(Freq_effect >= 0.01 & Freq_effect <= 0.99) %>%
#   filter(ID %in% w_hm3$rsid)
# 
# N_eff <- AA_ERpos_update$N_eff
# N_cut_low <- quantile(N_eff,probs = 0.1,na.rm = TRUE)
# 
# AA_ERpos_sub_quan <-  AA_ERpos %>% 
#   filter(Freq_effect >= 0.01 & Freq_effect <= 0.99) %>%
#   filter(ID %in% w_hm3$rsid) %>% 
#   filter(N_eff > N_cut_low)
# 
# AA_ERpos_sub_quan <- AA_ERpos_sub_quan[,names_data]
# 
# write.table(AA_ERpos_sub_quan, file = "AA_ERpos_sumdata_01.txt", row.names = F, col.names = T, quote = F)
# 
# AA_ERpos_update <-  AA_ERpos %>%
#   filter(Freq_effect >= 0.05 & Freq_effect <= 0.95) %>%
#   filter(ID %in% w_hm3$rsid)
# 
# N_eff <- AA_ERpos_update$N_eff
# N_cut_low <- quantile(N_eff,probs = 0.1,na.rm = TRUE)
# 
# AA_ERpos_sub_quan <-  AA_ERpos %>% 
#   filter(Freq_effect >= 0.05 & Freq_effect <= 0.95) %>%
#   filter(ID %in% w_hm3$rsid) %>% 
#   filter(N_eff > N_cut_low)
# 
# AA_ERpos_sub_quan <- AA_ERpos_sub_quan[,names_data]
# 
# write.table(AA_ERpos_sub_quan, file = "AA_ERpos_sumdata_05.txt", row.names = F, col.names = T, quote = F)
# 
# rm(list=setdiff(ls(), c("w_hm3","merge_data")))
# 
# # AA_ERneg
# AA_ERneg <- fread("Clean_summary_data/AA_ERneg_sumdata.txt")
# AA_ERneg <- as.data.frame(AA_ERneg)
# names_data <- colnames(AA_ERneg)
# 
# AA_ERneg$unique_id1 <- paste0(AA_ERneg$CHR,"_",AA_ERneg$POS,"_",
#                               AA_ERneg$effect_allele,"_",AA_ERneg$non_effect_allele)
# AA_ERneg$unique_id2 <- paste0(AA_ERneg$CHR,"_",AA_ERneg$POS,"_",
#                               AA_ERneg$non_effect_allele,"_",AA_ERneg$effect_allele)
# 
# AA_ERneg <- left_join(AA_ERneg,merge_data,by = c("unique_id1" = "unique_id"))
# AA_ERneg$ID[!is.na(AA_ERneg$rsid)] <- AA_ERneg$rsid[!is.na(AA_ERneg$rsid)]
# 
# AA_ERneg <- subset(AA_ERneg,select = -c(rsid))
# 
# AA_ERneg <- left_join(AA_ERneg,merge_data,by = c("unique_id2" = "unique_id"))
# AA_ERneg$ID[!is.na(AA_ERneg$rsid)] <- AA_ERneg$rsid[!is.na(AA_ERneg$rsid)]
# 
# AA_ERneg_update <-  AA_ERneg %>% 
#   filter(Freq_effect >= 0.01 & Freq_effect <= 0.99) %>%
#   filter(ID %in% w_hm3$rsid) 
# 
# N_eff <- AA_ERneg_update$N_eff
# N_cut_low <- quantile(N_eff,probs = 0.1,na.rm = TRUE)
# 
# AA_ERneg_sub_quan <-  AA_ERneg %>% 
#   filter(Freq_effect >= 0.01 & Freq_effect <= 0.99) %>%
#   filter(ID %in% w_hm3$rsid) %>% 
#   filter(N_eff > N_cut_low)
# 
# AA_ERneg_sub_quan <- AA_ERneg_sub_quan[,names_data]
# 
# write.table(AA_ERneg_sub_quan, file = "AA_ERneg_sumdata_01.txt", row.names = F, col.names = T, quote = F)
# 
# AA_ERneg_update <-  AA_ERneg %>% 
#   filter(Freq_effect >= 0.05 & Freq_effect <= 0.95) %>%
#   filter(ID %in% w_hm3$rsid) 
# 
# N_eff <- AA_ERneg_update$N_eff
# N_cut_low <- quantile(N_eff,probs = 0.1,na.rm = TRUE)
# 
# AA_ERneg_sub_quan <-  AA_ERneg %>% 
#   filter(Freq_effect >= 0.05 & Freq_effect <= 0.95) %>%
#   filter(ID %in% w_hm3$rsid) %>% 
#   filter(N_eff > N_cut_low)
# 
# AA_ERneg_sub_quan <- AA_ERneg_sub_quan[,names_data]
# 
# write.table(AA_ERneg_sub_quan, file = "AA_ERneg_sumdata_05.txt", row.names = F, col.names = T, quote = F)
# 
# rm(list=setdiff(ls(), c("w_hm3","merge_data")))

# Eur_Overall
Eur_Overall <- fread("Clean_summary_data/European_BCAC_icogs_onco_sumdata.txt")
Eur_Overall <- as.data.frame(Eur_Overall)
names_data <- colnames(Eur_Overall)

Eur_Overall$unique_id1 <- paste0(Eur_Overall$CHR,"_",Eur_Overall$POS,"_",
                              Eur_Overall$effect_allele_meta,"_",Eur_Overall$non_effect_allele_meta)
Eur_Overall$unique_id2 <- paste0(Eur_Overall$CHR,"_",Eur_Overall$POS,"_",
                              Eur_Overall$non_effect_allele_meta,"_",Eur_Overall$effect_allele_meta)

Eur_Overall <- left_join(Eur_Overall,merge_data,by = c("unique_id1" = "unique_id"))
Eur_Overall$SNP_ID[!is.na(Eur_Overall$rsid)] <- Eur_Overall$rsid[!is.na(Eur_Overall$rsid)]

Eur_Overall <- subset(Eur_Overall,select = -c(rsid))

Eur_Overall <- left_join(Eur_Overall,merge_data,by = c("unique_id2" = "unique_id"))
Eur_Overall$SNP_ID[!is.na(Eur_Overall$rsid)] <- Eur_Overall$rsid[!is.na(Eur_Overall$rsid)]

Eur_Overall_sub <- Eur_Overall %>% 
  filter(!is.na(effect_allele_iCOGs) & !is.na(effect_allele_Onco)) %>%
  filter(Freq_effect_Onco >= 0.01 & Freq_effect_Onco <= 0.99) %>%
  filter(SNP_ID %in% w_hm3$rsid)

N_eff <- Eur_Overall_sub$N_eff_meta
N_cut_low <- quantile(N_eff,probs = 0.1,na.rm = TRUE)

Eur_Overall_sub_quan <-  Eur_Overall %>% 
  filter(!is.na(effect_allele_iCOGs) & !is.na(effect_allele_Onco)) %>%
  filter(Freq_effect_Onco >= 0.01 & Freq_effect_Onco <= 0.99) %>%
  filter(SNP_ID %in% w_hm3$rsid) %>% 
  filter(N_eff_meta > N_cut_low)

Eur_Overall_sub_quan <- Eur_Overall_sub_quan[,names_data]

write.table(Eur_Overall_sub_quan, file = "EUR_Overall_HM3_01.txt", row.names = F, col.names = T, quote = F)

Eur_Overall_sub <- Eur_Overall %>% 
  filter(!is.na(effect_allele_iCOGs) & !is.na(effect_allele_Onco)) %>%
  filter(Freq_effect_Onco >= 0.05 & Freq_effect_Onco <= 0.95) %>%
  filter(SNP_ID %in% w_hm3$rsid)

N_eff <- Eur_Overall_sub$N_eff_meta
N_cut_low <- quantile(N_eff,probs = 0.1,na.rm = TRUE)

Eur_Overall_sub_quan <-  Eur_Overall %>% 
  filter(!is.na(effect_allele_iCOGs) & !is.na(effect_allele_Onco)) %>%
  filter(Freq_effect_Onco >= 0.05 & Freq_effect_Onco <= 0.95) %>%
  filter(SNP_ID %in% w_hm3$rsid) %>% 
  filter(N_eff_meta > N_cut_low)

Eur_Overall_sub_quan <- Eur_Overall_sub_quan[,names_data]

write.table(Eur_Overall_sub_quan, file = "EUR_Overall_HM3_05.txt", row.names = F, col.names = T, quote = F)

rm(list=setdiff(ls(), c("w_hm3","merge_data")))

# Eur_ERpos
Eur_ERpos_Overall <- fread("Clean_summary_data/European_ERpos_BCAC_icogs_onco_sumdata.txt")
Eur_ERpos_Overall <- as.data.frame(Eur_ERpos_Overall)
names_data <- colnames(Eur_ERpos_Overall)

Eur_ERpos_Overall$phase3_1kg_id[is.na(Eur_ERpos_Overall$phase3_1kg_id)] <- paste0(Eur_ERpos_Overall$chr[is.na(Eur_ERpos_Overall$phase3_1kg_id)],":",Eur_ERpos_Overall$position_b37[is.na(Eur_ERpos_Overall$phase3_1kg_id)])

Eur_ERpos_Overall$phase3_1kg_id[str_detect(Eur_ERpos_Overall$phase3_1kg_id,"rs")] <- unlist(lapply( strsplit(Eur_ERpos_Overall$phase3_1kg_id[str_detect(Eur_ERpos_Overall$phase3_1kg_id,"rs")],":")  ,function(x){x[1]}))
Eur_ERpos_Overall$phase3_1kg_id[!str_detect(Eur_ERpos_Overall$phase3_1kg_id,"rs")] <- paste0(Eur_ERpos_Overall$chr[!str_detect(Eur_ERpos_Overall$phase3_1kg_id,"rs")],":",Eur_ERpos_Overall$position_b37[!str_detect(Eur_ERpos_Overall$phase3_1kg_id,"rs")])

Eur_ERpos_Overall$unique_id1 <- paste0(Eur_ERpos_Overall$chr,"_",Eur_ERpos_Overall$position_b37,"_",
                              Eur_ERpos_Overall$effect_allele_meta,"_",Eur_ERpos_Overall$non_effect_allele_meta)
Eur_ERpos_Overall$unique_id2 <- paste0(Eur_ERpos_Overall$chr,"_",Eur_ERpos_Overall$position_b37,"_",
                              Eur_ERpos_Overall$non_effect_allele_meta,"_",Eur_ERpos_Overall$effect_allele_meta)

Eur_ERpos_Overall <- left_join(Eur_ERpos_Overall,merge_data,by = c("unique_id1" = "unique_id"))
Eur_ERpos_Overall$phase3_1kg_id[!is.na(Eur_ERpos_Overall$rsid)] <- Eur_ERpos_Overall$rsid[!is.na(Eur_ERpos_Overall$rsid)]

Eur_ERpos_Overall <- subset(Eur_ERpos_Overall,select = -c(rsid))

Eur_ERpos_Overall <- left_join(Eur_ERpos_Overall,merge_data,by = c("unique_id2" = "unique_id"))
Eur_ERpos_Overall$phase3_1kg_id[!is.na(Eur_ERpos_Overall$rsid)] <- Eur_ERpos_Overall$rsid[!is.na(Eur_ERpos_Overall$rsid)]

Eur_ERpos_Overall_sub <- Eur_ERpos_Overall %>% 
  filter(!is.na(BETA_iCOGs) & !is.na(BETA_Onco)) %>% 
  filter(Freq_effect_Onco >= 0.01 & Freq_effect_Onco <= 0.99) %>%
  filter(phase3_1kg_id %in% w_hm3$rsid)

N_eff <- Eur_ERpos_Overall_sub$N_eff_meta
N_cut_low <- quantile(N_eff,probs = 0.1,na.rm = TRUE)

Eur_ERpos_Overall_sub_quan <-  Eur_ERpos_Overall %>% 
  filter(!is.na(BETA_iCOGs) & !is.na(BETA_Onco)) %>% 
  filter(Freq_effect_Onco >= 0.01 & Freq_effect_Onco <= 0.99) %>%
  filter(phase3_1kg_id %in% w_hm3$rsid) %>% 
  filter(N_eff_meta > N_cut_low)

Eur_ERpos_Overall_sub_quan <- Eur_ERpos_Overall_sub_quan[,names_data]

write.table(Eur_ERpos_Overall_sub_quan, file = "EUR_ERPos_HM3_01.txt", row.names = F, col.names = T, quote = F)

Eur_ERpos_Overall_sub <- Eur_ERpos_Overall %>% 
  filter(!is.na(BETA_iCOGs) & !is.na(BETA_Onco)) %>% 
  filter(Freq_effect_Onco >= 0.05 & Freq_effect_Onco <= 0.95) %>%
  filter(phase3_1kg_id %in% w_hm3$rsid)

N_eff <- Eur_ERpos_Overall_sub$N_eff_meta
N_cut_low <- quantile(N_eff,probs = 0.1,na.rm = TRUE)

Eur_ERpos_Overall_sub_quan <-  Eur_ERpos_Overall %>% 
  filter(!is.na(BETA_iCOGs) & !is.na(BETA_Onco)) %>% 
  filter(Freq_effect_Onco >= 0.05 & Freq_effect_Onco <= 0.95) %>%
  filter(phase3_1kg_id %in% w_hm3$rsid) %>% 
  filter(N_eff_meta > N_cut_low)

Eur_ERpos_Overall_sub_quan <- Eur_ERpos_Overall_sub_quan[,names_data]

write.table(Eur_ERpos_Overall_sub_quan, file = "EUR_ERPos_HM3_05.txt", row.names = F, col.names = T, quote = F)

rm(list=setdiff(ls(), c("w_hm3","merge_data")))

# Eur_ERneg
Eur_ERneg_Overall <- fread("Clean_summary_data/European_ERneg_BCAC_icogs_onco_sumdata.txt")
Eur_ERneg_Overall <- as.data.frame(Eur_ERneg_Overall)
names_data <- colnames(Eur_ERneg_Overall)

sum(is.na(Eur_ERneg_Overall$phase3_1kg_id)) 
Eur_ERneg_Overall$phase3_1kg_id[is.na(Eur_ERneg_Overall$phase3_1kg_id)] <- paste0(Eur_ERneg_Overall$chr[is.na(Eur_ERneg_Overall$phase3_1kg_id)],":",Eur_ERneg_Overall$position_b37[is.na(Eur_ERneg_Overall$phase3_1kg_id)])

Eur_ERneg_Overall$phase3_1kg_id[str_detect(Eur_ERneg_Overall$phase3_1kg_id,"rs")] <- unlist(lapply( strsplit(Eur_ERneg_Overall$phase3_1kg_id[str_detect(Eur_ERneg_Overall$phase3_1kg_id,"rs")],":")  ,function(x){x[1]}))
Eur_ERneg_Overall$phase3_1kg_id[!str_detect(Eur_ERneg_Overall$phase3_1kg_id,"rs")] <- paste0(Eur_ERneg_Overall$chr[!str_detect(Eur_ERneg_Overall$phase3_1kg_id,"rs")],":",Eur_ERneg_Overall$position_b37[!str_detect(Eur_ERneg_Overall$phase3_1kg_id,"rs")])

Eur_ERneg_Overall$unique_id1 <- paste0(Eur_ERneg_Overall$chr,"_",Eur_ERneg_Overall$position_b37,"_",
                              Eur_ERneg_Overall$effect_allele_meta,"_",Eur_ERneg_Overall$non_effect_allele_meta)
Eur_ERneg_Overall$unique_id2 <- paste0(Eur_ERneg_Overall$chr,"_",Eur_ERneg_Overall$position_b37,"_",
                              Eur_ERneg_Overall$non_effect_allele_meta,"_",Eur_ERneg_Overall$effect_allele_meta)

Eur_ERneg_Overall <- left_join(Eur_ERneg_Overall,merge_data,by = c("unique_id1" = "unique_id"))
Eur_ERneg_Overall$phase3_1kg_id[!is.na(Eur_ERneg_Overall$rsid)] <- Eur_ERneg_Overall$rsid[!is.na(Eur_ERneg_Overall$rsid)]

Eur_ERneg_Overall <- subset(Eur_ERneg_Overall,select = -c(rsid))

Eur_ERneg_Overall <- left_join(Eur_ERneg_Overall,merge_data,by = c("unique_id2" = "unique_id"))
Eur_ERneg_Overall$phase3_1kg_id[!is.na(Eur_ERneg_Overall$rsid)] <- Eur_ERneg_Overall$rsid[!is.na(Eur_ERneg_Overall$rsid)]

Eur_ERneg_Overall_sub <- Eur_ERneg_Overall %>% 
  filter(!is.na(BETA_iCOGs) & !is.na(BETA_Onco)) %>% 
  filter(Freq_effect_Onco >= 0.01 & Freq_effect_Onco <= 0.99) %>%
  filter(phase3_1kg_id %in% w_hm3$rsid)

N_eff <- Eur_ERneg_Overall_sub$N_eff_meta
N_cut_low <- quantile(N_eff,probs = 0.1,na.rm = TRUE)

Eur_ERneg_Overall_sub_quan <-  Eur_ERneg_Overall %>% 
  filter(!is.na(BETA_iCOGs) & !is.na(BETA_Onco)) %>%
  filter(Freq_effect_Onco >= 0.01 & Freq_effect_Onco <= 0.99) %>%
  filter(phase3_1kg_id %in% w_hm3$rsid) %>% 
  filter(N_eff_meta > N_cut_low)

Eur_ERneg_Overall_sub_quan <- Eur_ERneg_Overall_sub_quan[,names_data]

write.table(Eur_ERneg_Overall_sub_quan, file = "EUR_ERNeg_HM3_01.txt", row.names = F, col.names = T, quote = F)

Eur_ERneg_Overall_sub <- Eur_ERneg_Overall %>% 
  filter(!is.na(BETA_iCOGs) & !is.na(BETA_Onco)) %>% 
  filter(Freq_effect_Onco >= 0.05 & Freq_effect_Onco <= 0.95) %>%
  filter(phase3_1kg_id %in% w_hm3$rsid)

N_eff <- Eur_ERneg_Overall_sub$N_eff_meta
N_cut_low <- quantile(N_eff,probs = 0.1,na.rm = TRUE)

Eur_ERneg_Overall_sub_quan <-  Eur_ERneg_Overall %>% 
  filter(!is.na(BETA_iCOGs) & !is.na(BETA_Onco)) %>%
  filter(Freq_effect_Onco >= 0.05 & Freq_effect_Onco <= 0.95) %>%
  filter(phase3_1kg_id %in% w_hm3$rsid) %>% 
  filter(N_eff_meta > N_cut_low)

Eur_ERneg_Overall_sub_quan <- Eur_ERneg_Overall_sub_quan[,names_data]

write.table(Eur_ERneg_Overall_sub_quan, file = "EUR_ERNeg_HM3_05.txt", row.names = F, col.names = T, quote = F)

rm(list=setdiff(ls(), c("w_hm3","merge_data")))

system("rm -r Clean_summary_data/")
system("rm w_hm3_updated.snplist")