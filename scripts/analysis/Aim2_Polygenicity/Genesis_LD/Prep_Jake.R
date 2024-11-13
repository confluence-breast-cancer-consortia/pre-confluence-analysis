rm(list =ls())

for(extract in c("PreExtract")){
  for(MAF in c("05")){
    for(anc in c("EUR","AMR","AFR","SAS","EAS")){
      for(hm3 in c("HM3")){
        dir <- paste0("/data/williamsjacr/Aim2_Genesis/",anc,"_",hm3,"_",MAF,"_",extract,"/")
        print(dir)
        system(paste0("mkdir ",dir))
        system(paste0("mkdir ",dir,"BEDFiles/"))
        system(paste0("mkdir ",dir,"FreqFiles/"))
        system(paste0("mkdir ",dir,"PairwiseLDFiles/"))
        system(paste0("mkdir ",dir,"Final_CompiledResult/"))
        system(paste0("mkdir ",dir,"Results/"))
        if(hm3 == "HM3"){
          w_hm3_ldsc <- read.delim(paste0("/data/williamsjacr/1KG_ldsc_HM3/",anc,"/",anc,"_hm3_ldsc/w_hm3_ldsc.snplist"))
          write.table(w_hm3_ldsc[,1],paste0(dir,"BEDFiles/extract.snplist"),row.names = F,col.names = F,quote=F)
        }else{
          hm3_mega_ldsc <- read.delim(paste0("/data/williamsjacr/1KG_ldsc_HM3_MEGA/",anc,"/",anc,"_hm3_mega_ldsc/hm3_mega_ldsc.snplist"))
          write.table(hm3_mega_ldsc[,1],paste0(dir,"BEDFiles/extract.snplist"),row.names = F,col.names = F,quote=F)
        } 
      }
    }
  }
}