rm(list = ls())

# for anc in EUR AMR EAS SAS AFR;
# do
# sbatch /spin1/home/linux/williamsjacr/Aim2_Genesis/Genesis_LD/Per_Chr_Jake.sh ${anc}
# done

temp <- commandArgs(TRUE)
chr <- as.numeric(temp[1]) # 1-22

plinkpath <- "/data/williamsjacr/software/plink"

anc <- temp[2]
extract <- temp[3]
MAF <- temp[4]
hm3 <- temp[5]

bfile <- paste0("/data/BB_Bioinformatics/ProjectData/1000G_unrelated/1000G_full_data/GRCh37/",anc,"/chr")

dir <- paste0("/data/williamsjacr/Aim2_Genesis/",anc,"_",hm3,"_",MAF,"_",extract,"/")
filename <- paste0(anc,"_MAF",MAF,"_",hm3)

##--------------#--------------#--------------
# Filter bed files 
##--------------#--------------#--------------
plinkcode = paste(plinkpath,
                  "--bfile", paste0(bfile, chr),
                  "--extract ",paste0(dir, "BEDFiles/extract.snplist"),
                  "--maf", 0.05,
                  "--make-bed",
                  "--out",paste0(dir,"BEDFiles/",filename,chr))
system(plinkcode)
##--------------#--------------#--------------
# 01. Calculate Allele frequency in 1000 Genome.
##--------------#--------------#--------------
bfile <- paste0("",dir,"BEDFiles/",filename,chr)
plinkcode = paste(plinkpath,
                  "--bfile", bfile,
                  "--freq",
                  "--out",paste0(dir,"FreqFiles/",filename,chr))
system(plinkcode) 


##--------------#--------------#--------------
# Calculate pairwise LD for all SNPs in w_hm3.snplist from 1000G_EUR_Phase3_plink (MAF>5%)
##--------------#--------------#--------------


for(ldwindowkb in c(500,1000,2000)){
  plinkcode = paste(plinkpath,
                    "--noweb",
                    "--bfile", bfile,
                    "--r2",
                    "--ld-window-kb", ldwindowkb,
                    "--ld-window", 99999,
                    "--ld-window-r2", 0, #minimum bound
                    "--out",paste0("",dir,"PairwiseLDFiles/",filename,"_ldwindow_",ldwindowkb/1000,"MB.",chr))
  system(plinkcode)
}

frqfile <- paste0(dir,"FreqFiles/",filename,chr,".frq")
for(ldwindowkb in c(500,1000,2000)){
  pairwiseLDfile <- paste0(dir,"PairwiseLDFiles/",filename,"_ldwindow_",ldwindowkb/1000,"MB.",chr,".ld")
  for(cutoff in c(0.05, 0.1, 0.2)){
    ##--------------#--------------#--------------
    # Calculate LD score
    ##--------------#--------------#--------------
    # ------------------------------------------------
    # no --ld-snp-list flag in plink,
    # .ld file doesn't calculate LD with itself (1.0)
    # .ld file is upper triangular.
    # ------------------------------------------------
    library(data.table)
    # read the pairwise LD data
    fam <- fread(paste0(dir,"BEDFiles/",filename,chr,".fam"))
    N <- nrow(fam)
    df <- fread(pairwiseLDfile)
    colnames(df) <- c("CHR_A",	"BP_A",	"SNP_A",	"CHR_B",	"BP_B",	"SNP_B",	"R2")
    
    
    snp <- NULL; maf <- NULL;
    tem <- fread(frqfile)
    snp <- c(snp, tem$SNP)
    maf <- c(maf, tem$MAF)
    
    inx <- which((maf>0.05|maf==0.05) & (maf<0.5|maf==0.5))
    commonSNP <- snp[inx]
    commonMAF <- maf[inx]
    K <- length(commonSNP)
    df_commonSNP <- cbind(commonSNP, commonMAF)
    df_commonSNP <- data.frame(df_commonSNP)
    
    df$inx = rep(1,nrow(df))
    inx_extract <- which(df$SNP_A %in% commonSNP)
    df <- df[inx_extract,]
    
    inx_extract0 <- which(df$SNP_B %in% commonSNP)
    df <- df[inx_extract0,]
    
    
    #---------------------------------
    dftemA <- unique(df[,c(1,2,3)])
    dftemB <- unique(df[,c(4,5,6)]);
    colnames(dftemB) <- colnames(dftemA)
    dftem <- unique(rbind(dftemA,dftemB))
    snptem <- unique(c(df$SNP_A, df$SNP_B))
    if(length(snptem) < length(snp)){
      bim <- fread(paste0(bfile,".bim"))
      s0 <- setdiff(snp, snptem)
      df0 <- matrix(0,length(s0),3)
      for(i in 1:length(s0)){
        snpid <- s0[i]
        inx <- which(bim$V2==snpid)
        df0[i,] <- c(chr, bim$V4[inx], snpid)
      }
      df0 <- data.frame(df0)
      colnames(df0) <- colnames(dftem)
      dftem <- unique(rbind(dftem, df0))
    }
    dftem$CHR_B <- dftem$CHR_A
    dftem$BP_B <- dftem$BP_A
    dftem$SNP_B <- dftem$SNP_A
    dftem$R2 <- 1 
    dftem$inx <- 1
    #---------------------------------
    
    
    #---------------------------------
    inx_extract1 <- which(df$R2>=cutoff)
    df <- df[inx_extract1,]
    #---------------------------------
    df.rev <- df
    df.rev$SNP_A <- df$SNP_B
    df.rev$SNP_B <- df$SNP_A
    df.rev$BP_A <- df$BP_B
    df.rev$BP_B <- df$BP_A
    
    
    dfnew <- rbind(df, df.rev)
    dfnew <- rbind(dfnew, dftem)
    
    ldscore <- aggregate(dfnew$R2, by=list(dfnew$SNP_A),FUN=sum)
    Nstar <- aggregate(dfnew$inx,by=list(dfnew$SNP_A), FUN=sum)
    result <- cbind(ldscore,Nstar[,2])
    colnames(result) <- c("SNPname","LD score","Nstar")
    result <- data.frame(result)
    
    # # ------------------------------------------------
    # get the index list for each SNP that the pairwise r^2 >threshold.
    TaggingSNPs <- aggregate(dfnew$SNP_B,by=list(dfnew$SNP_A),FUN= function(t) list(as.character(t)))
    result$TaggingSNPs <- TaggingSNPs[,2]
    result$TaggingSNPs <- vapply(result$TaggingSNPs, paste, collapse = ",", character(1L))
    # # ------------------------------------------------
    
    # # ------------------------------------------------
    # get the index list for each SNP that the pairwise r^2 >threshold.
    pairwiseLD <- aggregate(dfnew$R2,by=list(dfnew$SNP_A),FUN= function(t) list(as.character(t)))
    result$pairwiseLD <- pairwiseLD[,2]
    result$pairwiseLD <- vapply(result$pairwiseLD, paste, collapse = ",", character(1L))
    # # ------------------------------------------------
    
    result$LD.score.correct <- (result$LD.score * (N-1) -  result$Nstar)/(N-2)
    
    bp_snp <- unique(dfnew[,1:3])
    result <- merge(result, bp_snp,by.x="SNPname",by.y="SNP_A")
    result <- merge(result,df_commonSNP, by.x="SNPname",by.y="commonSNP",sort=F)
    
    binbound <- c(4.9,7.05, 9.82, 13.09, 17.08, 21.47, 26.48, 31.9, 37.73, 43.87, 51)/100
    # quantile normalized to a normal distribution N(0,1)
    quantNorm <- function(x){qnorm(rank(x,ties.method = "average")/(length(x)+1))}
    
    LLD <- rep(0, nrow(result));   maf <- as.numeric(as.character(result$commonMAF))
    for(i in 1:(length(binbound)-1)){
      inx <- which(maf>binbound[i] & maf<binbound[i+1])
      tem <- as.numeric(as.character(result$LD.score.correct))[inx]
      LLD[inx] <- quantNorm(tem)
    }
    result$LLD <- LLD
    
    write.table(result,file=paste0(dir,"PairwiseLDFiles/",filename,"_ldwindow_",ldwindowkb/1000,"MB.",chr,".cutoff",cutoff,".ldscore"),
                quote=F,sep = "\t",row.names=F,col.names = T)
  }
}


