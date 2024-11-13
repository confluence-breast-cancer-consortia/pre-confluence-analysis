rm(list = ls())

i <- 1

binary_continuous <- read.delim("/data/williamsjacr/Aim4_PartitionedHeritability/LDSC_Files/binary_continuous.annot")

annotation_names <- colnames(binary_continuous)[6:101]
binary_annotation_names <- names(which(apply(binary_continuous[,6:101],2,function(x){sum(x %in% c(0,1))}) == nrow(binary_continuous)))
binary_annotation_names <- gsub(".bed","",binary_annotation_names)

for(i in 1:22){
  system(paste0("cp /data/williamsjacr/Aim4_PartitionedHeritability/LDSC_Files/EAS_Baseline/baselineLD.",i,".annot.gz /data/williamsjacr/Aim4_PartitionedHeritability/LDSC_Files/AFR_AMR_annot_gz/"))
  system(paste0("gunzip /data/williamsjacr/Aim4_PartitionedHeritability/LDSC_Files/AFR_AMR_annot_gz/baselineLD.",i,".annot.gz"))
  EAS_Baseline <- read.delim(paste0("/data/williamsjacr/Aim4_PartitionedHeritability/LDSC_Files/AFR_AMR_annot_gz/baselineLD.",i,".annot"))
  colnames(EAS_Baseline) <- gsub(".bed","",colnames(EAS_Baseline))
  system(paste0("rm /data/williamsjacr/Aim4_PartitionedHeritability/LDSC_Files/AFR_AMR_annot_gz/baselineLD.",i,".annot"))
  
  system(paste0("cp /data/williamsjacr/Aim4_PartitionedHeritability/LDSC_Files/EUR_Baseline/baselineLD.",i,".annot.gz /data/williamsjacr/Aim4_PartitionedHeritability/LDSC_Files/AFR_AMR_annot_gz/"))
  system(paste0("gunzip /data/williamsjacr/Aim4_PartitionedHeritability/LDSC_Files/AFR_AMR_annot_gz/baselineLD.",i,".annot.gz"))
  EUR_Baseline <- read.delim(paste0("/data/williamsjacr/Aim4_PartitionedHeritability/LDSC_Files/AFR_AMR_annot_gz/baselineLD.",i,".annot"))
  colnames(EUR_Baseline) <- gsub(".bed","",colnames(EUR_Baseline))
  system(paste0("rm /data/williamsjacr/Aim4_PartitionedHeritability/LDSC_Files/AFR_AMR_annot_gz/baselineLD.",i,".annot"))
  
  All_Baseline <- unique(rbind(EAS_Baseline,EUR_Baseline))
  
  All_Baseline <- All_Baseline[,colnames(All_Baseline) %in% c("CHR","BP","SNP","CM","base",binary_annotation_names)]
  All_Baseline <- All_Baseline[,!(colnames(All_Baseline) %in% c("MAFbin1","MAFbin2","MAFbin3","MAFbin4","MAFbin5","MAFbin6","MAFbin7","MAFbin8","MAFbin9","MAFbin10"))]
  All_Baseline <- unique(All_Baseline)
  
  write.table(x = All_Baseline,file = paste0("/data/williamsjacr/Aim4_PartitionedHeritability/LDSC_Files/AFR_AMR_annot_gz/baselineLD.",i,".annot"), row.names = F, col.names = T, quote = F, sep = '\t')
  system(paste0("gzip /data/williamsjacr/Aim4_PartitionedHeritability/LDSC_Files/AFR_AMR_annot_gz/baselineLD.",i,".annot"))
}