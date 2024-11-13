rm(list = ls())

# for anc in EUR AMR EAS SAS AFR;
# do
# sbatch /spin1/home/linux/williamsjacr/Aim2_Genesis/Genesis_LD/CompileFinal_Jake.sh ${anc}
# done

temp <- commandArgs(TRUE)
idx <- as.numeric(temp[1]) # 1-9

anc <- temp[2]
extract <- temp[3]
MAF <- temp[4]
hm3 <- temp[5]

if(idx < 4){
  ldwindowkb <- 500
}else if (idx < 7){
  ldwindowkb <- 1000
}else{
  ldwindowkb <- 2000
}

cutoff <- c(0.05, 0.1, 0.2)[(idx %% 3) + 1]

dir <- paste0("/data/williamsjacr/Aim2_Genesis/",anc,"_",hm3,"_",MAF,"_",extract,"/")
filename <- paste0(anc,"_MAF",MAF,"_",hm3)

ldscorefileprefix <- paste0(dir,"PairwiseLDFiles/",filename,"_ldwindow_",ldwindowkb/1000,"MB.") 
ldscorefilepostfix <- paste0(".cutoff",cutoff,".ldscore")

# %--------------------------------------------------------------------
# merge 22 chromosomes LD score
# %--------------------------------------------------------------------
# Combine the 22 chromosomes LDscore into one file.
linuxcommand <- paste(paste0("awk 'FNR >1' ",ldscorefileprefix,"*",ldscorefilepostfix),
                      paste0("> ",  ldscorefileprefix, ldscorefilepostfix,".all"))
system(linuxcommand)
linuxcommand <- paste0("rm ",ldscorefileprefix,"*",ldscorefilepostfix)
system(linuxcommand) 

library(data.table)
# read LD file
dataLD <- fread(paste0(ldscorefileprefix, ldscorefilepostfix,".all"))
colnames(dataLD) = c("SNPname" , "LD.score", "Nstar", "TaggingSNPs",
                     "pairwiseLD" ,"LD.score.correct", "CHR" ,"BP","MAF" ,"LLD")

dataLD <- dataLD[,-c(2,5,9,10)]

save(dataLD, file=paste0(dir,"Final_CompiledResult/",filename,"_LDwindow_",ldwindowkb/1000,"MB_cutoff",cutoff,".RData"),
     compress="xz",compression_level=9)
