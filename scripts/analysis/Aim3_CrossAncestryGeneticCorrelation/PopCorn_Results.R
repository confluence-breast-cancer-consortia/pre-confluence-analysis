rm(list = ls())

library(corrplot)

Overall_HM3_results <- read.delim("/data/williamsjacr/Aim3_GeneticCorrelation/PopCorn_Results/Overall_HM3_results.txt")

hm3_cor_mat <- diag(4)
colnames(hm3_cor_mat) <- c("AFR","EAS","EUR","H/L")
rownames(hm3_cor_mat) <- c("AFR","EAS","EUR","H/L")
hm3_cor_mat[1,2] <- Overall_HM3_results$Val[Overall_HM3_results$X == "AFR_EAS"]
hm3_cor_mat[2,1] <- Overall_HM3_results$Val[Overall_HM3_results$X == "AFR_EAS"]
hm3_cor_mat[1,3] <- Overall_HM3_results$Val[Overall_HM3_results$X == "AFR_EUR"]
hm3_cor_mat[3,1] <- Overall_HM3_results$Val[Overall_HM3_results$X == "AFR_EUR"]
hm3_cor_mat[1,4] <- Overall_HM3_results$Val[Overall_HM3_results$X == "AFR_AMR"]
hm3_cor_mat[4,1] <- Overall_HM3_results$Val[Overall_HM3_results$X == "AFR_AMR"]
hm3_cor_mat[2,3] <- Overall_HM3_results$Val[Overall_HM3_results$X == "EAS_EUR"]
hm3_cor_mat[3,2] <- Overall_HM3_results$Val[Overall_HM3_results$X == "EAS_EUR"]
hm3_cor_mat[2,4] <- Overall_HM3_results$Val[Overall_HM3_results$X == "AMR_EAS"]
hm3_cor_mat[4,2] <- Overall_HM3_results$Val[Overall_HM3_results$X == "AMR_EAS"]
hm3_cor_mat[3,4] <- Overall_HM3_results$Val[Overall_HM3_results$X == "AMR_EUR"]
hm3_cor_mat[4,3] <- Overall_HM3_results$Val[Overall_HM3_results$X == "AMR_EUR"]

pdf(paste0("/data/williamsjacr/Aim3_GeneticCorrelation/PopCorn_Results/HM3_CorPlot.pdf"), width=4, height=4)
# corrplot.mixed(hm3_cor_mat, order = 'AOE',number.cex=2.5,cl.cex = 1.5,tl.cex = 2.5)
corrplot.mixed(hm3_cor_mat, order = 'AOE')
dev.off()


Overall_HM3_MEGA_results <- read.delim("/data/williamsjacr/Aim3_GeneticCorrelation/PopCorn_Results/Overall_HM3_MEGA_results.txt")

hm3_cor_mat <- diag(4)
colnames(hm3_cor_mat) <- c("AFR","EAS","EUR","H/L")
rownames(hm3_cor_mat) <- c("AFR","EAS","EUR","H/L")
hm3_cor_mat[1,2] <- Overall_HM3_MEGA_results$Val[Overall_HM3_MEGA_results$X == "AFR_EAS"]
hm3_cor_mat[2,1] <- Overall_HM3_MEGA_results$Val[Overall_HM3_MEGA_results$X == "AFR_EAS"]
hm3_cor_mat[1,3] <- Overall_HM3_MEGA_results$Val[Overall_HM3_MEGA_results$X == "AFR_EUR"]
hm3_cor_mat[3,1] <- Overall_HM3_MEGA_results$Val[Overall_HM3_MEGA_results$X == "AFR_EUR"]
hm3_cor_mat[1,4] <- Overall_HM3_MEGA_results$Val[Overall_HM3_MEGA_results$X == "AFR_AMR"]
hm3_cor_mat[4,1] <- Overall_HM3_MEGA_results$Val[Overall_HM3_MEGA_results$X == "AFR_AMR"]
hm3_cor_mat[2,3] <- Overall_HM3_MEGA_results$Val[Overall_HM3_MEGA_results$X == "EAS_EUR"]
hm3_cor_mat[3,2] <- Overall_HM3_MEGA_results$Val[Overall_HM3_MEGA_results$X == "EAS_EUR"]
hm3_cor_mat[2,4] <- Overall_HM3_MEGA_results$Val[Overall_HM3_MEGA_results$X == "AMR_EAS"]
hm3_cor_mat[4,2] <- Overall_HM3_MEGA_results$Val[Overall_HM3_MEGA_results$X == "AMR_EAS"]
hm3_cor_mat[3,4] <- Overall_HM3_MEGA_results$Val[Overall_HM3_MEGA_results$X == "AMR_EUR"]
hm3_cor_mat[4,3] <- Overall_HM3_MEGA_results$Val[Overall_HM3_MEGA_results$X == "AMR_EUR"]

pdf(paste0("/data/williamsjacr/Aim3_GeneticCorrelation/PopCorn_Results/HM3_Mega_CorPlot.pdf"), width=4, height=4)
# corrplot.mixed(hm3_cor_mat, order = 'AOE',number.cex=2.5,cl.cex = 1.5,tl.cex = 2.5)
corrplot.mixed(hm3_cor_mat, order = 'AOE')
dev.off()
