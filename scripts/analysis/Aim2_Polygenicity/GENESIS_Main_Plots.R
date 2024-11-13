rm(list = ls())

library(ggplot2)

ff <- function(x,y){
  x*y/(x+y)
}

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
            axis.title.y = element_text(angle=90,vjust =2,size = 16,face = "bold"),
            axis.title.x = element_text(size = 16,face = "bold"),
            axis.text.x = element_text(size = 20),
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
            legend.position=c(0.7,0.3),
            #legend.direction = "horizontal",
            #legend.key.size= unit(0.2, "cm"),
            #legend.margin = unit(0, "cm"),
            legend.title=element_blank(),
            #legend.text = element_text(face ="bold"),
            plot.margin=unit(c(0.5,0.5,0.5,0.5), "cm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

load("/data/williamsjacr/Aim2_Genesis/EUR_HM3_05_PreExtract/Results/EUR_Full1KG_Phenotypic_Variance.RData")

plot_data <- results[,c("n_seq","pheno.bestGWAS","AUC_GWAS","Genetic_Variance_GWAS")]
plot_data <- data.frame(Ancestry = "EUR",plot_data)

load("/data/williamsjacr/Aim2_Genesis/AFR_HM3_05_AABCG/Results/AABCG_Phenotypic_Variance.RData")

plot_data <- rbind(plot_data,data.frame(Ancestry = "AA",results[,c("n_seq","pheno.bestGWAS","AUC_GWAS","Genetic_Variance_GWAS")]))

load("/data/williamsjacr/Aim2_Genesis/EAS_HM3_05_PreExtract/Results/EAS_Full1KG_Phenotypic_Variance.RData")

plot_data <- rbind(plot_data,data.frame(Ancestry = "EAS",results[,c("n_seq","pheno.bestGWAS","AUC_GWAS","Genetic_Variance_GWAS")]))

p1 <- ggplot(plot_data,aes(x = n_seq/1000, y = pheno.bestGWAS,color = Ancestry)) + geom_line(size = 1) + theme_Publication() + ggtitle("Phenotypic Variance") + 
  ylab("Percentage of phenotypic variance explained") + xlab("Total sample size assuming 1:1 case:control ratio (in thousands)") + 
  scale_color_manual(values=c("#386cb0","#EF7E3D","#ffd558","#7fc97f","#ef3b2c","#973999")) + geom_vline(xintercept = ff(100000,100000)*4/1000,color = "#ef3b2c",linetype = "dashed",size = 1)

ggsave(p1,filename="/data/williamsjacr/Aim2_Genesis/Phenotypic_Variance.pdf",width = 10,height = 6.1804697157)

p2 <- ggplot(plot_data,aes(x = n_seq/1000, y = AUC_GWAS,color = Ancestry)) + geom_line(size = 1) + theme_Publication() + ggtitle("AUC") +
  ylab("AUC associated with the PRS") + xlab("Total sample size assuming 1:1 case:control ratio (in thousands)") + 
  scale_color_manual(values=c("#386cb0","#EF7E3D","#ffd558","#7fc97f","#ef3b2c","#973999")) + geom_vline(xintercept = ff(100000,100000)*4/1000,color = "#ef3b2c",linetype = "dashed",size = 1)

ggsave(p2,filename="/data/williamsjacr/Aim2_Genesis/AUC.pdf",width = 10,height = 6.1804697157)

p3 <- ggplot(plot_data,aes(x = n_seq/1000, y = Genetic_Variance_GWAS,color = Ancestry)) + geom_line(size = 1) + theme_Publication() + ggtitle("Genetic Variance") +
  ylab("Percentage of genetic variance explained") + xlab("Total sample size assuming 1:1 case:control ratio (in thousands)") + 
  scale_color_manual(values=c("#386cb0","#EF7E3D","#ffd558","#7fc97f","#ef3b2c","#973999")) + geom_vline(xintercept = ff(100000,100000)*4/1000,color = "#ef3b2c",linetype = "dashed",size = 1)

ggsave(p3,filename="/data/williamsjacr/Aim2_Genesis/Genetic_Variance.pdf",width = 10,height = 6.1804697157)
