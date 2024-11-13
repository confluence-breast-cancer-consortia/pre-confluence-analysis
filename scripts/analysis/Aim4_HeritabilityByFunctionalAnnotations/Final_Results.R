rm(list = ls())

library(ggplot2)

AFR_Results <- read.delim("/data/williamsjacr/Aim4_PartitionedHeritability/Results/AFR_05_HM3_My_LDSC.results")
AFR_Results$Prop._h2[AFR_Results$Prop._h2 < 0] <- 0
AFR_Results$Enrichment <- AFR_Results$Prop._h2/AFR_Results$Prop._SNPs
AFR_Results$Category <- gsub("L2_0","",AFR_Results$Category)
AFR_Results$Category[AFR_Results$Category == "H3K27ac_Hnisz"] <- "H3K27ac"
AFR_Results$Category[AFR_Results$Category == "SuperEnhancer_Hnisz"] <- "SuperEnhancer"
AFR_Results$Category[AFR_Results$Category == "Ancient_Sequence_Age_Human_Promoter"] <- "Ancient Promoter"
AFR_Results$Category[AFR_Results$Category == "TFBS_ENCODE"] <- "TFBS"
AFR_Results$Category[AFR_Results$Category == "H3K4me3_Trynka"] <- "H3K4me3"
AFR_Results$Category[AFR_Results$Category == "H3K4me1_Trynka"] <- "H3K4me1"
AFR_Results <- AFR_Results[AFR_Results$Category!= "base",]
AFR_Table <- AFR_Results[AFR_Results$Enrichment_p < 0.05/nrow(AFR_Results),]

AFR_AABCG_Results <- read.delim("/data/williamsjacr/Aim4_PartitionedHeritability/Results/AFR_05_HM3_AABCG_LDSC.results")
AFR_AABCG_Results$Prop._h2[AFR_AABCG_Results$Prop._h2 < 0] <- 0
AFR_AABCG_Results$Enrichment <- AFR_AABCG_Results$Prop._h2/AFR_AABCG_Results$Prop._SNPs
AFR_AABCG_Results$Category <- gsub("L2_0","",AFR_AABCG_Results$Category)
AFR_AABCG_Results$Category[AFR_AABCG_Results$Category == "H3K27ac_Hnisz"] <- "H3K27ac"
AFR_AABCG_Results$Category[AFR_AABCG_Results$Category == "SuperEnhancer_Hnisz"] <- "SuperEnhancer"
AFR_AABCG_Results$Category[AFR_AABCG_Results$Category == "Ancient_Sequence_Age_Human_Promoter"] <- "Ancient Promoter"
AFR_AABCG_Results$Category[AFR_AABCG_Results$Category == "TFBS_ENCODE"] <- "TFBS"
AFR_AABCG_Results$Category[AFR_AABCG_Results$Category == "H3K4me3_Trynka"] <- "H3K4me3"
AFR_AABCG_Results$Category[AFR_AABCG_Results$Category == "H3K4me1_Trynka"] <- "H3K4me1"
AFR_AABCG_Results <- AFR_AABCG_Results[AFR_AABCG_Results$Category!= "base",]
AFR_AABCG_Table <- AFR_AABCG_Results[AFR_AABCG_Results$Enrichment_p < 0.05/nrow(AFR_AABCG_Results),]

AMR_Results <- read.delim("/data/williamsjacr/Aim4_PartitionedHeritability/Results/AMR_05_HM3_My_LDSC.results")
AMR_Results$Prop._h2[AMR_Results$Prop._h2 < 0] <- 0
AMR_Results$Enrichment <- AMR_Results$Prop._h2/AMR_Results$Prop._SNPs
AMR_Results$Category <- gsub("L2_0","",AMR_Results$Category)
AMR_Results$Category[AMR_Results$Category == "H3K27ac_Hnisz"] <- "H3K27ac"
AMR_Results$Category[AMR_Results$Category == "SuperEnhancer_Hnisz"] <- "SuperEnhancer"
AMR_Results$Category[AMR_Results$Category == "Ancient_Sequence_Age_Human_Promoter"] <- "Ancient Promoter"
AMR_Results$Category[AMR_Results$Category == "TFBS_ENCODE"] <- "TFBS"
AMR_Results$Category[AMR_Results$Category == "H3K4me3_Trynka"] <- "H3K4me3"
AMR_Results$Category[AMR_Results$Category == "H3K4me1_Trynka"] <- "H3K4me1"
AMR_Results <- AMR_Results[AMR_Results$Category!= "base",]
AMR_Table <- AMR_Results[AMR_Results$Enrichment_p < 0.05/nrow(AMR_Results),]

EUR_Mine_Results <- read.delim("/data/williamsjacr/Aim4_PartitionedHeritability/Results/EUR_05_HM3_My_LDSC.results")
EUR_Mine_Results$Prop._h2[EUR_Mine_Results$Prop._h2 < 0] <- 0
EUR_Mine_Results$Enrichment <- EUR_Mine_Results$Prop._h2/EUR_Mine_Results$Prop._SNPs
EUR_Mine_Results$Category <- gsub("L2_0","",EUR_Mine_Results$Category)
EUR_Mine_Results$Category[EUR_Mine_Results$Category == "H3K27ac_Hnisz"] <- "H3K27ac"
EUR_Mine_Results$Category[EUR_Mine_Results$Category == "SuperEnhancer_Hnisz"] <- "SuperEnhancer"
EUR_Mine_Results$Category[EUR_Mine_Results$Category == "Ancient_Sequence_Age_Human_Promoter"] <- "Ancient Promoter"
EUR_Mine_Results$Category[EUR_Mine_Results$Category == "TFBS_ENCODE"] <- "TFBS"
EUR_Mine_Results$Category[EUR_Mine_Results$Category == "H3K4me3_Trynka"] <- "H3K4me3"
EUR_Mine_Results$Category[EUR_Mine_Results$Category == "H3K4me1_Trynka"] <- "H3K4me1"
EUR_Mine_Results <- EUR_Mine_Results[EUR_Mine_Results$Category!= "base",]
EUR_Mine_Table <- EUR_Mine_Results[EUR_Mine_Results$Enrichment_p < 0.05/nrow(EUR_Mine_Results),]

EUR_Theirs_Results <- read.delim("/data/williamsjacr/Aim4_PartitionedHeritability/Results/EUR_05_HM3_Their_LDSC.results")
EUR_Theirs_Results$Prop._h2[EUR_Theirs_Results$Prop._h2 < 0] <- 0
EUR_Theirs_Results$Enrichment <- EUR_Theirs_Results$Prop._h2/EUR_Theirs_Results$Prop._SNPs
EUR_Theirs_Results$Category <- gsub("L2_0","",EUR_Theirs_Results$Category)
EUR_Theirs_Results$Category[EUR_Theirs_Results$Category == "H3K27ac_Hnisz"] <- "H3K27ac"
EUR_Theirs_Results$Category[EUR_Theirs_Results$Category == "SuperEnhancer_Hnisz"] <- "SuperEnhancer"
EUR_Theirs_Results$Category[EUR_Theirs_Results$Category == "Ancient_Sequence_Age_Human_Promoter"] <- "Ancient Promoter"
EUR_Theirs_Results$Category[EUR_Theirs_Results$Category == "TFBS_ENCODE"] <- "TFBS"
EUR_Theirs_Results$Category[EUR_Theirs_Results$Category == "H3K4me3_Trynka"] <- "H3K4me3"
EUR_Theirs_Results$Category[EUR_Theirs_Results$Category == "H3K4me1_Trynka"] <- "H3K4me1"
EUR_Theirs_Results <- EUR_Theirs_Results[EUR_Theirs_Results$Category!= "base",]
EUR_Theirs_Table <- EUR_Theirs_Results[EUR_Theirs_Results$Enrichment_p < 0.05/73,]

EAS_Mine_Results <- read.delim("/data/williamsjacr/Aim4_PartitionedHeritability/Results/EAS_05_HM3_My_LDSC.results")
EAS_Mine_Results$Prop._h2[EAS_Mine_Results$Prop._h2 < 0] <- 0
EAS_Mine_Results$Enrichment <- EAS_Mine_Results$Prop._h2/EAS_Mine_Results$Prop._SNPs
EAS_Mine_Results$Category <- gsub("L2_0","",EAS_Mine_Results$Category)
EAS_Mine_Results$Category[EAS_Mine_Results$Category == "H3K27ac_Hnisz"] <- "H3K27ac"
EAS_Mine_Results$Category[EAS_Mine_Results$Category == "SuperEnhancer_Hnisz"] <- "SuperEnhancer"
EAS_Mine_Results$Category[EAS_Mine_Results$Category == "Ancient_Sequence_Age_Human_Promoter"] <- "Ancient Promoter"
EAS_Mine_Results$Category[EAS_Mine_Results$Category == "TFBS_ENCODE"] <- "TFBS"
EAS_Mine_Results$Category[EAS_Mine_Results$Category == "H3K4me3_Trynka"] <- "H3K4me3"
EAS_Mine_Results$Category[EAS_Mine_Results$Category == "H3K4me1_Trynka"] <- "H3K4me1"
EAS_Mine_Results <- EAS_Mine_Results[EAS_Mine_Results$Category!= "base",]
EAS_Mine_Table <- EAS_Mine_Results[EAS_Mine_Results$Enrichment_p < 0.05/nrow(EAS_Mine_Results),]

EAS_Theirs_Results <- read.delim("/data/williamsjacr/Aim4_PartitionedHeritability/Results/EAS_05_HM3_Their_LDSC.results")
EAS_Theirs_Results$Prop._h2[EAS_Theirs_Results$Prop._h2 < 0] <- 0
EAS_Theirs_Results$Enrichment <- EAS_Theirs_Results$Prop._h2/EAS_Theirs_Results$Prop._SNPs
EAS_Theirs_Results$Category <- gsub("L2_0","",EAS_Theirs_Results$Category)
EAS_Theirs_Results$Category[EAS_Theirs_Results$Category == "H3K27ac_Hnisz"] <- "H3K27ac"
EAS_Theirs_Results$Category[EAS_Theirs_Results$Category == "SuperEnhancer_Hnisz"] <- "SuperEnhancer"
EAS_Theirs_Results$Category[EAS_Theirs_Results$Category == "Ancient_Sequence_Age_Human_Promoter"] <- "Ancient Promoter"
EAS_Theirs_Results$Category[EAS_Theirs_Results$Category == "TFBS_ENCODE"] <- "TFBS"
EAS_Theirs_Results$Category[EAS_Theirs_Results$Category == "H3K4me3_Trynka"] <- "H3K4me3"
EAS_Theirs_Results$Category[EAS_Theirs_Results$Category == "H3K4me1_Trynka"] <- "H3K4me1"
EAS_Theirs_Results <- EAS_Theirs_Results[EAS_Theirs_Results$Category!= "base",]
EAS_Theirs_Table <- EAS_Theirs_Results[EAS_Theirs_Results$Enrichment_p < 0.05/73,]



All_EUR_Table <- data.frame(Ancestry = "EUR",Name = EUR_Mine_Results$Category,Percent_SNPs = round(EUR_Mine_Results$Prop._SNPs,digits = 3),Percent_H2 = paste0(round(EUR_Mine_Results$Prop._h2,digits = 3)," (",round(EUR_Mine_Results$Prop._h2_std_error,digits = 3),")"),
                            Enrichment = paste0(round(EUR_Mine_Results$Enrichment,digits = 3)," (",round(EUR_Mine_Results$Enrichment_std_error,digits = 3),")"),Enrichment_PValue = EUR_Mine_Results$Enrichment_p)
All_EAS_Table <- data.frame(Ancestry = "EAS",Name = EAS_Mine_Results$Category,Percent_SNPs = round(EAS_Mine_Results$Prop._SNPs,digits = 3),Percent_H2 = paste0(round(EAS_Mine_Results$Prop._h2,digits = 3)," (",round(EAS_Mine_Results$Prop._h2_std_error,digits = 3),")"),
                            Enrichment = paste0(round(EAS_Mine_Results$Enrichment,digits = 3)," (",round(EAS_Mine_Results$Enrichment_std_error,digits = 3),")"),Enrichment_PValue = EAS_Mine_Results$Enrichment_p)
All_AMR_Table <- data.frame(Ancestry = "AMR",Name = AMR_Results$Category,Percent_SNPs = round(AMR_Results$Prop._SNPs,digits = 3),Percent_H2 = paste0(round(AMR_Results$Prop._h2,digits = 3)," (",round(AMR_Results$Prop._h2_std_error,digits = 3),")"),
                            Enrichment = paste0(round(AMR_Results$Enrichment,digits = 3)," (",round(AMR_Results$Enrichment_std_error,digits = 3),")"),Enrichment_PValue = AMR_Results$Enrichment_p)
All_AA_Table <- data.frame(Ancestry = "AA",Name = AFR_AABCG_Results$Category,Percent_SNPs = round(AFR_AABCG_Results$Prop._SNPs,digits = 3),Percent_H2 = paste0(round(AFR_AABCG_Results$Prop._h2,digits = 3)," (",round(AFR_AABCG_Results$Prop._h2_std_error,digits = 3),")"),
                           Enrichment = paste0(round(AFR_AABCG_Results$Enrichment,digits = 3)," (",round(AFR_AABCG_Results$Enrichment_std_error,digits = 3),")"),Enrichment_PValue = AFR_AABCG_Results$Enrichment_p)

All <- rbind(All_EUR_Table,All_EAS_Table,All_AMR_Table,All_AA_Table)
write.csv(All,file = "All.csv",row.names = FALSE)

EUR_Results <- data.frame(Ancestry = "EUR",Name = EUR_Mine_Table$Category,Percent_SNPs = round(EUR_Mine_Table$Prop._SNPs,digits = 3),Percent_H2 = paste0(round(EUR_Mine_Table$Prop._h2,digits = 3)," (",round(EUR_Mine_Table$Prop._h2_std_error,digits = 3),")"),
                          Enrichment = paste0(round(EUR_Mine_Table$Enrichment,digits = 3)," (",round(EUR_Mine_Table$Enrichment_std_error,digits = 3),")"),Enrichment_PValue = EUR_Mine_Table$Enrichment_p)
EAS_Results <- data.frame(Ancestry = "EAS",Name = EAS_Mine_Table$Category,Percent_SNPs = round(EAS_Mine_Table$Prop._SNPs,digits = 3),Percent_H2 = paste0(round(EAS_Mine_Table$Prop._h2,digits = 3)," (",round(EAS_Mine_Table$Prop._h2_std_error,digits = 3),")"),
                          Enrichment = paste0(round(EAS_Mine_Table$Enrichment,digits = 3)," (",round(EAS_Mine_Table$Enrichment_std_error,digits = 3),")"),Enrichment_PValue = EAS_Mine_Table$Enrichment_p)
AABCG_Results <- data.frame(Ancestry = "AA",Name = AFR_AABCG_Table$Category,Percent_SNPs = round(AFR_AABCG_Table$Prop._SNPs,digits = 3),Percent_H2 = paste0(round(AFR_AABCG_Table$Prop._h2,digits = 3)," (",round(AFR_AABCG_Table$Prop._h2_std_error,digits = 3),")"),
                            Enrichment = paste0(round(AFR_AABCG_Table$Enrichment,digits = 3)," (",round(AFR_AABCG_Table$Enrichment_std_error,digits = 3),")"),Enrichment_PValue = AFR_AABCG_Table$Enrichment_p)

write.csv(rbind(EUR_Results,EAS_Results,AABCG_Results),file = "tmp.csv",row.names = FALSE)

Figure_Data <- rbind(data.frame(Ancestry = "AFR",AFR_AABCG_Results),data.frame(Ancestry = "AMR",AMR_Results),data.frame(Ancestry = "EUR",EUR_Mine_Results),data.frame(Ancestry = "EAS",EAS_Mine_Results))
Figure_Data <- Figure_Data[Figure_Data$Category %in% EUR_Results$Name,]
Figure_Data <- Figure_Data[Figure_Data$Category != "Ancient Promoter",]
Figure_Data$Ancestry[Figure_Data$Ancestry == "AFR"] <- "AA"

# Figure_Data <- Figure_Data[Figure_Data$Ancestry != "AMR",]

C <- matrix(0,nrow = 3,ncol = 4)
C[,1] <- 1; C[1,2] <- -1;C[2,3] <- -1;C[3,4] <- -1

library(expm)

heterogenity_dat <- NULL

for(cat in unique(Figure_Data$Category)){
  tmp <- Figure_Data[Figure_Data$Category == cat,]
  Beta <- matrix(c(tmp$Enrichment[tmp$Ancestry == "EUR"],tmp$Enrichment[tmp$Ancestry == "AA"],tmp$Enrichment[tmp$Ancestry == "AMR"],tmp$Enrichment[tmp$Ancestry == "EAS"]),ncol = 1)
  Sigma <- diag(c(tmp$Enrichment_std_error[tmp$Ancestry == "EUR"]^2,tmp$Enrichment_std_error[tmp$Ancestry == "AA"]^2,tmp$Enrichment_std_error[tmp$Ancestry == "AMR"]^2,tmp$Enrichment_std_error[tmp$Ancestry == "EAS"]^2))
  Q <- C%*%Beta
  Q_Stand <- solve(sqrtm(C%*%Sigma%*%t(C)))%*%Q
  test_statistic <- as.numeric(t(Q_Stand)%*%Q_Stand)
  # test_statistic <- t(C%*%Beta)%*%solve(C%*%Sigma%*%t(C))%*%C%*%Beta
  p_value <- pchisq(test_statistic,df = 4,lower.tail = FALSE) 
  heterogenity_dat <- rbind(heterogenity_dat,data.frame(Category = cat,Chi_Squared_Stat = test_statistic,Chi_Squared_Pvalue = p_value))
}

write.csv(heterogenity_dat,file = "/data/williamsjacr/Aim4_PartitionedHeritability/Results/Heterogenity_Test.csv",row.names = FALSE)

library(ggplot2)

theme_Publication <- function(base_size=12) {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, )
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.1), hjust = 0.5),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            title = element_text(size = 24),
            axis.title.y = element_text(angle=90,vjust =2,size = 24,face = "bold"),
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            axis.text.y = element_text(size = 20),
            axis.line = element_line(colour="black",size=2),
            axis.ticks = element_line(),
            strip.text.x = element_text(size = 20),
            legend.text=element_text(size=20),
            # panel.grid.major = element_line(colour="#f0f0f0"),
            # panel.grid.minor = element_line(colour="#f0f0f0"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            #legend.position = "bottom",
            #legend.direction = "horizontal",
            #legend.key.size= unit(0.2, "cm"),
            #legend.margin = unit(0, "cm"),
            legend.title = element_text(face="bold.italic", size =24),
            #legend.text = element_text(face ="bold"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

scale_fill_Publication <- function(...){   
  library(scales)   
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#EF7E3D","#ffd558","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)    
} 

scale_colour_Publication <- function(...){   
  library(scales)   
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#EF7E3D","#ffd558","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
}

Figure_Data$Category <- factor(Figure_Data$Category,levels = EUR_Results$Name[order(EUR_Results$Enrichment,decreasing = TRUE)])
p<-ggplot(data=Figure_Data, aes(x=Ancestry, y=Enrichment,fill = Ancestry)) +
  geom_bar(stat="identity",color="black") + geom_errorbar(aes(ymin=Enrichment, ymax=Enrichment+Enrichment_std_error), width=.2,position=position_dodge(.9)) + 
  theme_Publication() + ggtitle("Significantly Enriched Annotations") + facet_wrap(Category ~ ., scales = "free_y", ncol = 2) + 
  scale_fill_Publication()

pdf(paste0("/data/williamsjacr/Aim4_PartitionedHeritability/Results/Enrichment.pdf"), width=20, height=20)
p + geom_hline(yintercept=1, linetype="dashed", color = "black")
dev.off()
