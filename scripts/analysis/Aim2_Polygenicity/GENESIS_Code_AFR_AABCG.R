rm(list = ls())

# devtools::install("/data/williamsjacr/software/GENESIS/")
library(GENESIS)
library(data.table)
library(tidyverse)

anc <- "AFR"
MAF <- "05"

load("/data/williamsjacr/Aim2_Genesis/AFR_HM3_05_AABCG/Final_CompiledResult_AFR_AABCG/AFR_MAF05_HM3_LDwindow_1MB_cutoff0.1.RData")
w_hm3_ldsc <- read.delim(paste0("/data/williamsjacr/1KG_ldsc_HM3/",anc,"/",anc,"_hm3_ldsc/w_hm3_ldsc.snplist"))

sumstats <- fread("/data/BB_Bioinformatics/ProjectData/Preconfluence/GWAS_SumStats/AA_Overall_HM3_05_2024.txt")
sumstats <- sumstats[sumstats$P > 0,]
dat <- inner_join(sumstats,w_hm3_ldsc,by=c("SNPID"="SNP"))

dat_s <- dat %>% mutate(t1 = effect_allele == A1 & non_effect_allele == A2) %>%
  mutate(t2 = non_effect_allele == A1 & effect_allele == A2)

dat_s <- dat_s %>% mutate(
  Z = case_when(
    t1 ~ BETA / SE,
    t2 ~ - BETA / SE,
    !t1 & !t2 ~ NA
  )
)

tolerance = c(1e-6,1e-8,1e-9,1e-6,10e3,2)

dat_s <- dat_s %>% select(SNPID, Z, N_eff) %>% filter(!is.na(Z)) %>%
  rename(SNP=SNPID, N=N_eff) # 960217

dat_s1 <- preprocessing(dat_s,filter = T,dataLD = dataLD)

fit2 <- genesis(dat_s1[,1:3], filter=F, modelcomponents=2, cores=30, LDcutoff=0.1, LDwindow=1, c0=10, startingpic=0.005,qqplot.name="/data/williamsjacr/Aim2_Genesis/AFR_HM3_05_AABCG/Results/qqplot_fit2_unrelated",dataLD = dataLD,tolerance = tolerance)

est2 <- fit2$estimates$`Parameter (pic, sigmasq, a) estimates` # the model parameter estimates
v2 <- fit2$estimates$`Covariance matrix of parameter estimates` # the covariance matrix of model parameter estimtaes

starting <- rep(0,5)
starting[1] <- est2[1]
starting[2] <- 1/9
starting[3] <- est2[2]*5
starting[4] <- starting[3]/10
starting[5] <- est2[3]

tolerance = c(1e-6,1e-5,1e-8,1e-8,1e-9,1e-6,10e3,2)

fit3 <- genesis(dat_s1[,1:3], filter=F, modelcomponents=3, cores=30, LDcutoff=0.1, LDwindow=1, c0=10,starting=starting,qqplot.name="/data/williamsjacr/Aim2_Genesis/AFR_HM3_05_AABCG/Results/qqplot_fit3_unrelated",dataLD = dataLD,tolerance = tolerance)

est <- fit3$estimates$`Parameter (pic, p1, sigmasq1, sigmasq2, a) estimates` # the model parameter estimates
v <- fit3$estimates$`Covariance matrix of parameter estimates` # the covariance matrix of model parameter estimtaes

save(fit2,fit3,file="/data/williamsjacr/Aim2_Genesis/AFR_HM3_05_AABCG/Results/fit2_fit3_unrelated.RData")