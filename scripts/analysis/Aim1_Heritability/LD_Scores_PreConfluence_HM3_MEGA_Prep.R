rm(list = ls())

library(dplyr)
library(stringr)
library(readxl)

eth <- c("EUR","AFR","EAS","AMR","SAS")

hm3_mega <- read.csv("/data/williamsjacr/ldsc/snpinfo_hm3_mega.csv")
hm3_mega_ldsc <- hm3_mega[,c("rsid","allele1","allele2")]
colnames(hm3_mega_ldsc) <- c("SNP","A1","A2")
write.table(hm3_mega_ldsc,file = "/data/williamsjacr/ldsc/hm3_mega_ldsc.snplist",row.names = FALSE,quote=F,sep = '\t')
write.table(hm3_mega_ldsc[,1,drop = FALSE],file = "/data/williamsjacr/1KG_ldsc_HM3_MEGA/hm3_mega_ldsc_printsnps.snplist",row.names = FALSE,quote=F,sep = '\t',col.names = FALSE)
write.table(hm3_mega_ldsc[,c(1,2)],file = "/data/williamsjacr/ldsc/ref_allele_hm3_mega.txt",row.names = FALSE,quote=F,col.names = FALSE)

Gazal_Inbreeding_amp_Related <- read_excel("/data/williamsjacr/Gazal_Inbreeding &amp; Related.xls")
colnames(Gazal_Inbreeding_amp_Related) <- c("IID","Super_Pop","Pop","Gender","Q_Score","f","p_value","Mating_Type","TGP2457","TGP2261")
IDS_BAD <- Gazal_Inbreeding_amp_Related$IID[Gazal_Inbreeding_amp_Related$TGP2457 != "YES" | Gazal_Inbreeding_amp_Related$TGP2261 != "YES"]

for(i in 1:length(eth)){
  
  kg.dir <- paste0("/data/BB_Bioinformatics/ProjectData/1000G_full_data/GRCh37/",eth[i],"/")
  local_dir <- paste0("/data/williamsjacr/1KG_ldsc_HM3_MEGA/",eth[i],"/")
  system(paste0("mkdir ",local_dir))
  system(paste0("mkdir ",local_dir,eth[i],"_hm3_mega_ldsc/"))
  system(paste0("cp /data/williamsjacr/ldsc/hm3_mega_ldsc.snplist ",local_dir,eth[i],"_hm3_mega_ldsc/"))
  
  for(j in 1:22){
    fam <- read.table(paste0(kg.dir,"chr",j,".fam"))
    fam <- fam[!(fam$V2 %in% IDS_BAD),]
    write.table(cbind(0,fam$V2),file = paste0(local_dir,"keeplist"),row.names = F,col.names = F,quote=F)
    
    system(paste0("/data/williamsjacr/software/plink2 --bfile ",kg.dir,"chr",j," --keep ",local_dir,"keeplist --ref-allele /data/williamsjacr/ldsc/ref_allele_hm3_mega.txt --make-bed --out ",local_dir,"chr_clean",j))
    system(paste0("rm ",paste0(local_dir,"keeplist")))
  }
  
  if(j == 6){
    bim <- read.table(paste0(local_dir,"chr_clean",j,".bim"))
    bim <- bim[bim$V4 >= 25002566 & bim$V4 <= 34994414,]
    write.table(bim$V2,file = paste0(local_dir,"remove_region_chr6"),row.names = F,col.names = F,quote=F)
    system(paste0("/data/williamsjacr/software/plink2 --bfile ",local_dir,"chr_clean",j," --exclude ",local_dir,"remove_region_chr6 --make-bed --out ",local_dir,"chr_clean",j))
    
    system(paste0("rm ",paste0(local_dir,"chr_clean",j,".bim~")))
    system(paste0("rm ",paste0(local_dir,"chr_clean",j,".fam~")))
    system(paste0("rm ",paste0(local_dir,"chr_clean",j,".bed~")))
    system(paste0("rm ",paste0(local_dir,"remove_region_chr6")))
  }
}