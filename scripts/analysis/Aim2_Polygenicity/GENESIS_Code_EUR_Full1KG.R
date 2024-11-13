rm(list = ls())

# devtools::install("/data/williamsjacr/software/GENESIS/")
library(GENESIS)
library(data.table)
library(tidyverse)

anc <- "EUR"
MAF <- "05"

load("/data/williamsjacr/software/GENESIS/data/LDwindow1MB_cutoff0.1.RData")
load("/data/williamsjacr/software/GENESIS/data/w_hm3.noMHC.snplist.RData")
w_hm3_ldsc <- w_hm3.noMHC.snplist

sumstats <- fread("/data/BB_Bioinformatics/ProjectData/Preconfluence/GWAS_SumStats/EUR_Overall_HM3_05.txt")
sumstats <- sumstats[sumstats$P_meta > 0,]
dat <- inner_join(sumstats,w_hm3_ldsc, by=c("SNP_ID"="SNP"))

dat_s <- dat %>% mutate(t1 = effect_allele_meta == A1 & non_effect_allele_meta == A2) %>%
  mutate(t2 = non_effect_allele_meta == A1 & effect_allele_meta == A2)

dat_s <- dat_s %>% mutate(
  Z = case_when(
    t1 ~ BETA_meta / SE_meta,
    t2 ~ - BETA_meta / SE_meta,
    !t1 & !t2 ~ NA
  )
)

dat_s <- dat_s %>% select(SNP_ID, Z, N_eff_meta) %>% filter(!is.na(Z)) %>%
  rename(SNP=SNP_ID, N=N_eff_meta) # 960217

dat_s1 <- preprocessing(dat_s,filter = T,dataLD = dataLD)

fit2 <- genesis(dat_s1[,1:3], filter=F, modelcomponents=2, cores=20, LDcutoff=0.1, LDwindow=1, c0=10, startingpic=0.005,qqplot.name="/data/williamsjacr/Aim2_Genesis/EUR_HM3_05_PreExtract/Results/qqplot_fit2_Full1KG",dataLD = dataLD)

est2 <- fit2$estimates$`Parameter (pic, sigmasq, a) estimates` # the model parameter estimates
v2 <- fit2$estimates$`Covariance matrix of parameter estimates` # the covariance matrix of model parameter estimtaes

starting <- rep(0,5)
starting[1] <- est2[1]
starting[2] <- 1/9
starting[3] <- est2[2]*5
starting[4] <- starting[3]/10
starting[5] <- est2[3]

fit3 <- genesis(dat_s1[,1:3], filter=F, modelcomponents=3, cores=20, LDcutoff=0.1, LDwindow=1, c0=10,starting=starting,qqplot.name="/data/williamsjacr/Aim2_Genesis/EUR_HM3_05_PreExtract/Results/qqplot_fit3_Full1KG",dataLD = dataLD)

est <- fit3$estimates$`Parameter (pic, p1, sigmasq1, sigmasq2, a) estimates` # the model parameter estimates
v <- fit3$estimates$`Covariance matrix of parameter estimates` # the covariance matrix of model parameter estimtaes

save(fit2,fit3,file="/data/williamsjacr/Aim2_Genesis/EUR_HM3_05_PreExtract/Results/fit2_fit3_Full1KG.RData")