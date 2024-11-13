rm(list = ls())
library(data.table)
library(stringr)
library(dplyr)

# dx run app-swiss-army-knife -iin=PreConfluenceAnalysis:Aim1_Heritability/JacobWilliams/Clean_Data_Scripts/r_with_plink.tar.gz -iin=PreConfluenceAnalysis:Aim1_Heritability/JacobWilliams/Clean_Data_Scripts/AA_Overall_2024.R -iin=PreConfluenceAnalysis:Aim1_Heritability/JacobWilliams/Clean_Data_Scripts/AA_Overall_2024.sh -iin=PreConfluenceAnalysis:Aim1_Heritability/JacobWilliams/AA_Overall_Zheng2024/AABCG_sum_GRCh37.csv -iin=PreConfluenceAnalysis:Aim1_Heritability/JacobWilliams/w_hm3_updated.snplist -iin=PreConfluenceAnalysis:Aim1_Heritability/JacobWilliams/snpinfo_hm3_mega.csv  -icmd="bash AA_Overall_2024.sh" -y --destination PreConfluenceAnalysis:Aim2_Polygenicity/HZ/Results/Clean_summary_data_update/ --instance-type mem1_ssd1_v2_x36

w_hm3 <- read.csv("w_hm3_updated.snplist")

unique_id <- paste0(w_hm3$chr,"_",w_hm3$pos,"_",
                    w_hm3$allele1,"_",w_hm3$allele2)

merge_data <- data.frame(rsid = w_hm3$rsid, unique_id = unique_id)


AA_Updated <- read.csv("AABCG_sum_GRCh37.csv")

AA_Updated$unique_id1 <- paste0(AA_Updated$CHR,"_",AA_Updated$POS,"_",
                                AA_Updated$effect_allele,"_",AA_Updated$non_effect_allele)
AA_Updated$unique_id2 <- paste0(AA_Updated$CHR,"_",AA_Updated$POS,"_",
                                AA_Updated$non_effect_allele,"_",AA_Updated$effect_allele)

AA_Updated <- left_join(AA_Updated,merge_data,by = c("unique_id1" = "unique_id"))
AA_Updated$SNPID[!is.na(AA_Updated$rsid)] <- AA_Updated$rsid[!is.na(AA_Updated$rsid)]

AA_Updated <- subset(AA_Updated,select = -c(rsid))

AA_Updated <- left_join(AA_Updated,merge_data,by = c("unique_id2" = "unique_id"))
AA_Updated$SNPID[!is.na(AA_Updated$rsid)] <- AA_Updated$rsid[!is.na(AA_Updated$rsid)]

AA_Updated_sub <- AA_Updated %>% 
  filter(!is.na(effect_allele) & !is.na(effect_allele)) %>% 
  filter(Freq_effect >= 0.05 & Freq_effect <= 0.95) %>%
  filter(SNPID %in% w_hm3$rsid)

N_eff <- AA_Updated_sub$N_eff
N_cut_low <- quantile(N_eff,probs = 0.1,na.rm = TRUE)

AA_Updated_sub_quan <-  AA_Updated %>% 
  filter(!is.na(effect_allele) & !is.na(effect_allele)) %>% 
  filter(Freq_effect >= 0.05 & Freq_effect <= 0.95) %>%
  filter(SNPID %in% w_hm3$rsid) %>% 
  filter(N_eff > N_cut_low)

AA_Updated_sub_quan <- AA_Updated_sub_quan[,c("SNPID","CHR","POS","effect_allele","non_effect_allele","Freq_effect","BETA","SE","P","N_eff")]

write.table(AA_Updated_sub_quan, file = "AA_Overall_HM3_05_2024.txt", row.names = F, col.names = T, quote = F)

rm(list = ls())

w_hm3 <- read.csv("w_hm3_updated.snplist")

unique_id <- paste0(w_hm3$chr,"_",w_hm3$pos,"_",
                    w_hm3$allele1,"_",w_hm3$allele2)

merge_data <- data.frame(rsid = w_hm3$rsid, unique_id = unique_id)


AA_Updated <- read.csv("AABCG_sum_GRCh37.csv")

AA_Updated$unique_id1 <- paste0(AA_Updated$CHR,"_",AA_Updated$POS,"_",
                                AA_Updated$effect_allele,"_",AA_Updated$non_effect_allele)
AA_Updated$unique_id2 <- paste0(AA_Updated$CHR,"_",AA_Updated$POS,"_",
                                AA_Updated$non_effect_allele,"_",AA_Updated$effect_allele)

AA_Updated <- left_join(AA_Updated,merge_data,by = c("unique_id1" = "unique_id"))
AA_Updated$SNPID[!is.na(AA_Updated$rsid)] <- AA_Updated$rsid[!is.na(AA_Updated$rsid)]

AA_Updated <- subset(AA_Updated,select = -c(rsid))

AA_Updated <- left_join(AA_Updated,merge_data,by = c("unique_id2" = "unique_id"))
AA_Updated$SNPID[!is.na(AA_Updated$rsid)] <- AA_Updated$rsid[!is.na(AA_Updated$rsid)]

AA_Updated_sub <- AA_Updated %>% 
  filter(!is.na(effect_allele) & !is.na(effect_allele)) %>% 
  filter(Freq_effect >= 0.01 & Freq_effect <= 0.99) %>%
  filter(SNPID %in% w_hm3$rsid)

N_eff <- AA_Updated_sub$N_eff
N_cut_low <- quantile(N_eff,probs = 0.1,na.rm = TRUE)

AA_Updated_sub_quan <-  AA_Updated %>% 
  filter(!is.na(effect_allele) & !is.na(effect_allele)) %>% 
  filter(Freq_effect >= 0.01 & Freq_effect <= 0.99) %>%
  filter(SNPID %in% w_hm3$rsid) %>% 
  filter(N_eff > N_cut_low)

AA_Updated_sub_quan <- AA_Updated_sub_quan[,c("SNPID","CHR","POS","effect_allele","non_effect_allele","Freq_effect","BETA","SE","P","N_eff")]

write.table(AA_Updated_sub_quan, file = "AA_Overall_HM3_01_2024.txt", row.names = F, col.names = T, quote = F)


rm(list = ls())

hm3_mega <- read.csv("snpinfo_hm3_mega.csv")

unique_id <- paste0(hm3_mega$chr,"_",hm3_mega$pos,"_",
                    hm3_mega$allele1,"_",hm3_mega$allele2)

merge_data <- data.frame(rsid = hm3_mega$rsid, unique_id = unique_id)


AA_Updated <- read.csv("AABCG_sum_GRCh37.csv")

AA_Updated$unique_id1 <- paste0(AA_Updated$CHR,"_",AA_Updated$POS,"_",
                                AA_Updated$effect_allele,"_",AA_Updated$non_effect_allele)
AA_Updated$unique_id2 <- paste0(AA_Updated$CHR,"_",AA_Updated$POS,"_",
                                AA_Updated$non_effect_allele,"_",AA_Updated$effect_allele)

AA_Updated <- left_join(AA_Updated,merge_data,by = c("unique_id1" = "unique_id"))
AA_Updated$SNPID[!is.na(AA_Updated$rsid)] <- AA_Updated$rsid[!is.na(AA_Updated$rsid)]

AA_Updated <- subset(AA_Updated,select = -c(rsid))

AA_Updated <- left_join(AA_Updated,merge_data,by = c("unique_id2" = "unique_id"))
AA_Updated$SNPID[!is.na(AA_Updated$rsid)] <- AA_Updated$rsid[!is.na(AA_Updated$rsid)]

AA_Updated_sub <- AA_Updated %>% 
  filter(!is.na(effect_allele) & !is.na(effect_allele)) %>% 
  filter(Freq_effect >= 0.05 & Freq_effect <= 0.95) %>%
  filter(SNPID %in% hm3_mega$rsid)

N_eff <- AA_Updated_sub$N_eff
N_cut_low <- quantile(N_eff,probs = 0.1,na.rm = TRUE)

AA_Updated_sub_quan <-  AA_Updated %>% 
  filter(!is.na(effect_allele) & !is.na(effect_allele)) %>% 
  filter(Freq_effect >= 0.05 & Freq_effect <= 0.95) %>%
  filter(SNPID %in% hm3_mega$rsid) %>% 
  filter(N_eff > N_cut_low)

AA_Updated_sub_quan <- AA_Updated_sub_quan[,c("SNPID","CHR","POS","effect_allele","non_effect_allele","Freq_effect","BETA","SE","P","N_eff")]

write.table(AA_Updated_sub_quan, file = "AA_Overall_HM3_MEGA_05_2024.txt", row.names = F, col.names = T, quote = F)

rm(list = ls())

hm3_mega <- read.csv("snpinfo_hm3_mega.csv")

unique_id <- paste0(hm3_mega$chr,"_",hm3_mega$pos,"_",
                    hm3_mega$allele1,"_",hm3_mega$allele2)

merge_data <- data.frame(rsid = hm3_mega$rsid, unique_id = unique_id)


AA_Updated <- read.csv("AABCG_sum_GRCh37.csv")

AA_Updated$unique_id1 <- paste0(AA_Updated$CHR,"_",AA_Updated$POS,"_",
                                AA_Updated$effect_allele,"_",AA_Updated$non_effect_allele)
AA_Updated$unique_id2 <- paste0(AA_Updated$CHR,"_",AA_Updated$POS,"_",
                                AA_Updated$non_effect_allele,"_",AA_Updated$effect_allele)

AA_Updated <- left_join(AA_Updated,merge_data,by = c("unique_id1" = "unique_id"))
AA_Updated$SNPID[!is.na(AA_Updated$rsid)] <- AA_Updated$rsid[!is.na(AA_Updated$rsid)]

AA_Updated <- subset(AA_Updated,select = -c(rsid))

AA_Updated <- left_join(AA_Updated,merge_data,by = c("unique_id2" = "unique_id"))
AA_Updated$SNPID[!is.na(AA_Updated$rsid)] <- AA_Updated$rsid[!is.na(AA_Updated$rsid)]

AA_Updated_sub <- AA_Updated %>% 
  filter(!is.na(effect_allele) & !is.na(effect_allele)) %>% 
  filter(Freq_effect >= 0.01 & Freq_effect <= 0.99) %>%
  filter(SNPID %in% hm3_mega$rsid)

N_eff <- AA_Updated_sub$N_eff
N_cut_low <- quantile(N_eff,probs = 0.1,na.rm = TRUE)

AA_Updated_sub_quan <-  AA_Updated %>% 
  filter(!is.na(effect_allele) & !is.na(effect_allele)) %>% 
  filter(Freq_effect >= 0.01 & Freq_effect <= 0.99) %>%
  filter(SNPID %in% hm3_mega$rsid) %>% 
  filter(N_eff > N_cut_low)

AA_Updated_sub_quan <- AA_Updated_sub_quan[,c("SNPID","CHR","POS","effect_allele","non_effect_allele","Freq_effect","BETA","SE","P","N_eff")]

write.table(AA_Updated_sub_quan, file = "AA_Overall_HM3_MEGA_01_2024.txt", row.names = F, col.names = T, quote = F)




