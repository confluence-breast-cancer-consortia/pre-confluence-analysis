rm(list = ls())
library(devtools)
# devtools::install("/data/williamsjacr/software/GENESIS/")
# library(GENESIS)
library(gridExtra)
library(grid)
library(ggplot2)

Table1 <- NULL
plot_data <- NULL

dmixssnp <- function(x,est){
  
  if(length(est)==5){
    pic = est[1]
    p0 = est[2]
    s1 = sqrt(est[3])
    s2 = sqrt(est[4])
    den <- function(x){return((p0 * dnorm(x/s1)/s1 + (1-p0)*dnorm(x/s2) /s2))}
  }
  
  if(length(est)==3){
    pic = est[1]
    s1 = sqrt(est[2])
    den <- function(x){return(dnorm(x/s1)/s1)}
  }
  return(den(x))
}

### EUR Full 1KG

load("/data/williamsjacr/Aim2_Genesis/EUR_HM3_05_PreExtract/Results/fit2_fit3_Full1KG.RData")

Table1 <- rbind(Table1,
                data.frame(Data = "EUR_HM3_Full1KG",Number_of_SNPs = fit3$estimates$`Number of sSNPs (sd)`, 
                           Total_Heritability_logOR = fit3$estimates$`Total heritability in log-odds-ratio scale (sd)`,
                           Number_of_SNPs_LargestVariance = fit3$estimates$`Number of sSNPs in the cluster with larger variance component (sd)`,
                           Heritability_LargestVariance = fit3$estimates$`Heritability explained by the cluster with larger variance component (sd)`,
                           pic = fit3$estimates$`Parameter (pic, p1, sigmasq1, sigmasq2, a) estimates`[1]))

est <- fit3$estimates$`Parameter (pic, p1, sigmasq1, sigmasq2, a) estimates` # the model parameter estimates
v <- fit3$estimates$`Covariance matrix of parameter estimates` # the covariance matrix of model parameter estimates

x_seq <- seq(-0.02,0.02,length.out = 1000); 
y_seq <- apply(matrix(x_seq,ncol=1),1,function(t) dmixssnp(t,est))

plot_data <- rbind(plot_data,data.frame(x = x_seq,y = y_seq, group = "EUR Full 1KG"))

### EUR Unrelated

load("/data/williamsjacr/Aim2_Genesis/EUR_HM3_05_PreExtract/Results/fit2_fit3_unrelated.RData")

Table1 <- rbind(Table1,
                data.frame(Data = "EUR_HM3_05_Unrelated",Number_of_SNPs = fit3$estimates$`Number of sSNPs (sd)`, 
                           Total_Heritability_logOR = fit3$estimates$`Total heritability in log-odds-ratio scale (sd)`,
                           Number_of_SNPs_LargestVariance = fit3$estimates$`Number of sSNPs in the cluster with larger variance component (sd)`,
                           Heritability_LargestVariance = fit3$estimates$`Heritability explained by the cluster with larger variance component (sd)`,
                           pic = fit3$estimates$`Parameter (pic, p1, sigmasq1, sigmasq2, a) estimates`[1]))

est <- fit3$estimates$`Parameter (pic, p1, sigmasq1, sigmasq2, a) estimates` # the model parameter estimates
v <- fit3$estimates$`Covariance matrix of parameter estimates` # the covariance matrix of model parameter estimates

x_seq <- seq(-0.02,0.02,length.out = 1000); 
y_seq <- apply(matrix(x_seq,ncol=1),1,function(t) dmixssnp(t,est))

plot_data <- rbind(plot_data,data.frame(x = x_seq,y = y_seq, group = "EUR Unrelated"))

### AFR Unrelated

load("/data/williamsjacr/Aim2_Genesis/AFR_HM3_05_PreExtract/Results/fit2_fit3_unrelated.RData")

Table1 <- rbind(Table1,
                data.frame(Data = "AFR_HM3_05_Unrelated",Number_of_SNPs = fit3$estimates$`Number of sSNPs (sd)`, 
                           Total_Heritability_logOR = fit3$estimates$`Total heritability in log-odds-ratio scale (sd)`,
                           Number_of_SNPs_LargestVariance = fit3$estimates$`Number of sSNPs in the cluster with larger variance component (sd)`,
                           Heritability_LargestVariance = fit3$estimates$`Heritability explained by the cluster with larger variance component (sd)`,
                           pic = fit3$estimates$`Parameter (pic, p1, sigmasq1, sigmasq2, a) estimates`[1]))

est <- fit3$estimates$`Parameter (pic, p1, sigmasq1, sigmasq2, a) estimates` # the model parameter estimates
v <- fit3$estimates$`Covariance matrix of parameter estimates` # the covariance matrix of model parameter estimates

x_seq <- seq(-0.02,0.02,length.out = 1000); 
y_seq <- apply(matrix(x_seq,ncol=1),1,function(t) dmixssnp(t,est))

plot_data <- rbind(plot_data,data.frame(x = x_seq,y = y_seq, group = "AFR Unrelated"))

### AFR AABCG

load("/data/williamsjacr/Aim2_Genesis/AFR_HM3_05_AABCG/Results/fit2_fit3_unrelated.RData")

Table1 <- rbind(Table1,
                data.frame(Data = "AFR_HM3_05_AABCG",Number_of_SNPs = fit3$estimates$`Number of sSNPs (sd)`, 
                           Total_Heritability_logOR = fit3$estimates$`Total heritability in log-odds-ratio scale (sd)`,
                           Number_of_SNPs_LargestVariance = fit3$estimates$`Number of sSNPs in the cluster with larger variance component (sd)`,
                           Heritability_LargestVariance = fit3$estimates$`Heritability explained by the cluster with larger variance component (sd)`,
                           pic = fit3$estimates$`Parameter (pic, p1, sigmasq1, sigmasq2, a) estimates`[1]))

est <- fit3$estimates$`Parameter (pic, p1, sigmasq1, sigmasq2, a) estimates` # the model parameter estimates
v <- fit3$estimates$`Covariance matrix of parameter estimates` # the covariance matrix of model parameter estimates

x_seq <- seq(-0.02,0.02,length.out = 1000); 
y_seq <- apply(matrix(x_seq,ncol=1),1,function(t) dmixssnp(t,est))

plot_data <- rbind(plot_data,data.frame(x = x_seq,y = y_seq, group = "AFR AABCG"))

### EAS Unrelated

load("/data/williamsjacr/Aim2_Genesis/EAS_HM3_05_PreExtract/Results/fit2_fit3_Full1KG.RData")

Table1 <- rbind(Table1,
                data.frame(Data = "EAS_HM3_05_Unrelated",Number_of_SNPs = fit3$estimates$`Number of sSNPs (sd)`, 
                           Total_Heritability_logOR = fit3$estimates$`Total heritability in log-odds-ratio scale (sd)`,
                           Number_of_SNPs_LargestVariance = fit3$estimates$`Number of sSNPs in the cluster with larger variance component (sd)`,
                           Heritability_LargestVariance = fit3$estimates$`Heritability explained by the cluster with larger variance component (sd)`,
                           pic = fit3$estimates$`Parameter (pic, p1, sigmasq1, sigmasq2, a) estimates`[1]))

est <- fit3$estimates$`Parameter (pic, p1, sigmasq1, sigmasq2, a) estimates` # the model parameter estimates
v <- fit3$estimates$`Covariance matrix of parameter estimates` # the covariance matrix of model parameter estimates

x_seq <- seq(-0.02,0.02,length.out = 1000); 
y_seq <- apply(matrix(x_seq,ncol=1),1,function(t) dmixssnp(t,est))

plot_data <- rbind(plot_data,data.frame(x = x_seq,y = y_seq, group = "EAS Full 1KG"))

### EAS Full 1KG

load("/data/williamsjacr/Aim2_Genesis/EAS_HM3_05_PreExtract/Results/fit2_fit3_unrelated.RData")

Table1 <- rbind(Table1,
                data.frame(Data = "EAS_HM3_05_Full1KG",Number_of_SNPs = fit3$estimates$`Number of sSNPs (sd)`, 
                           Total_Heritability_logOR = fit3$estimates$`Total heritability in log-odds-ratio scale (sd)`,
                           Number_of_SNPs_LargestVariance = fit3$estimates$`Number of sSNPs in the cluster with larger variance component (sd)`,
                           Heritability_LargestVariance = fit3$estimates$`Heritability explained by the cluster with larger variance component (sd)`,
                           pic = fit3$estimates$`Parameter (pic, p1, sigmasq1, sigmasq2, a) estimates`[1]))

est <- fit3$estimates$`Parameter (pic, p1, sigmasq1, sigmasq2, a) estimates` # the model parameter estimates
v <- fit3$estimates$`Covariance matrix of parameter estimates` # the covariance matrix of model parameter estimates

x_seq <- seq(-0.02,0.02,length.out = 1000); 
y_seq <- apply(matrix(x_seq,ncol=1),1,function(t) dmixssnp(t,est))

plot_data <- rbind(plot_data,data.frame(x = x_seq,y = y_seq, group = "EAS Unrelated"))

pdf(paste0("/data/williamsjacr/Aim2_Genesis/Density.pdf"), width=15, height=15)
# do.call(grid.arrange,c(lgv,ncol=4,as.table=T))

plot_data <- plot_data[plot_data$group %in% c("EAS Full 1KG","AFR AABCG","EUR Full 1KG"),]
plot_data$group[plot_data$group == "EAS Full 1KG"] <- "EAS"
plot_data$group[plot_data$group == "EUR Full 1KG"] <- "EUR"
plot_data$group[plot_data$group == "AFR AABCG"] <- "AFR"

ggplot(plot_data,aes(x = x,y = y,color = group)) +
  geom_line(size=0.8) + theme_bw() + labs(x = "Odds Ratio", y = "Probability density") +
  scale_color_manual(values=c("#386cb0","#EF7E3D","#ffd558","#7fc97f","#ef3b2c","#973999")) +
  theme(text = element_text(size=24), 
        # panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        plot.title=element_text(size=24,face="italic",hjust=0.5),
        axis.text = element_text(size=24), 
        legend.key.width = unit(6, "line"),
        legend.position=c(.85,.9), #"none",
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm")) + theme(legend.title=element_blank())

dev.off()

### AMR Unrelated

load("/data/williamsjacr/Aim2_Genesis/AMR_HM3_05_PreExtract/Results/fit2_fit3_unrelated.RData")

Table1 <- rbind(Table1,
                data.frame(Data = "AMR_Fit3",Number_of_SNPs = fit3$estimates$`Number of sSNPs (sd)`, 
                           Total_Heritability_logOR = fit3$estimates$`Total heritability in log-odds-ratio scale (sd)`,
                           Number_of_SNPs_LargestVariance = fit3$estimates$`Number of sSNPs in the cluster with larger variance component (sd)`,
                           Heritability_LargestVariance = fit3$estimates$`Heritability explained by the cluster with larger variance component (sd)`,
                           pic = fit3$estimates$`Parameter (pic, p1, sigmasq1, sigmasq2, a) estimates`[1]))

### AMR Unrelated

load("/data/williamsjacr/Aim2_Genesis/AMR_HM3_05_PreExtract/Results/fit2_fit3_unrelated.RData")

Table1 <- rbind(Table1,
                data.frame(Data = "AMR_Fit2",Number_of_SNPs = fit2$estimates$`Number of sSNPs (sd)`, 
                           Total_Heritability_logOR = fit2$estimates$`Total heritability in log-odds-ratio scale (sd)`,
                           Number_of_SNPs_LargestVariance = NA,
                           Heritability_LargestVariance = NA,
                           pic = fit2$estimates$`Parameter (pic, sigmasq, a) estimates`[1]))


######### Table 1
herit <- unlist(lapply(strsplit(Table1$Total_Heritability_logOR, " "),function(x){as.numeric(x[1])}))
sd_herit <- unlist(lapply(strsplit(Table1$Total_Heritability_logOR, " "),function(x){as.numeric(gsub(".*\\((.*)\\).*", "\\1",x[2], " "))})) 

Table1$AUC <- paste0(round(pnorm(sqrt(herit/2)),digits = 3)," (",round(sd_herit * dnorm(sqrt(herit/2))/(sqrt(2)*2*sqrt(herit)),digits =5),")")

Table1

######### Figure 2

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
            axis.title.y = element_text(angle=90,vjust =2,size = 24,face = "bold"),
            axis.title.x = element_blank(),
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

p1 <- ggplot() +
  geom_line(data=results, aes(x=n_seq/1000, y=pheno.best,linetype="dotted"), size=2.5,color="grey38") +
  geom_line(data=results, aes(x=n_seq/1000, y=pheno.bestGWAS,linetype="solid"), size=2.5,color="grey38") +
  geom_point(data=data.frame(n = ff(118474,96201)/1000,pheno.best = current_sample_optimum$pheno.best),aes(x=n, y=pheno.best, color="lightgoldenrod3"), size=6)+    
  geom_point(data=data.frame(n = ff(118474,96201)/1000,pheno.best = current_sample_gwas$pheno.best),aes(x=n, y=pheno.best, color="lightgoldenrod3"), size=6)+    
  geom_point(data=data.frame(n = 2*ff(118474,96201)/1000,pheno.best = double_sample_optimum$pheno.best),aes(x=n, y=pheno.best,color="palevioletred"), size=6)+    
  geom_point(data=data.frame(n = 2*ff(118474,96201)/1000,pheno.best = double_sample_gwas$pheno.best),aes(x=n, y=pheno.best,color="palevioletred"), size=6)+    
  geom_point(data=data.frame(n = 4*ff(118474,96201)/1000,pheno.best = quadruple_sample_optimum$pheno.best),aes(x=n, y=pheno.best,color="royalblue1"), size=6)+    
  geom_point(data=data.frame(n = 4*ff(118474,96201)/1000,pheno.best = quadruple_sample_gwas$pheno.best),aes(x=n, y=pheno.best,color="royalblue1"), size=6)+    
  theme_bw() + scale_color_manual(labels = c("Current", "Double","Quadruple"),values=c("royalblue1","palevioletred","lightgoldenrod3")) +
  scale_linetype_manual(values = c('solid', 'dotted'),labels = c("Optimized p","5e-8")) + 
  labs(x = paste0("\n"),y = "\n") +
  ggtitle("EUR") +
  theme_Publication()

load("/data/williamsjacr/Aim2_Genesis/EUR_HM3_05_PreExtract/Results/EUR_Unrelated_Phenotypic_Variance.RData")

p2 <- ggplot() +
  geom_line(data=results, aes(x=n_seq/1000, y=pheno.best,linetype="dotted"), size=2.5,color="grey38") +
  geom_line(data=results, aes(x=n_seq/1000, y=pheno.bestGWAS,linetype="solid"), size=2.5,color="grey38") +
  geom_point(data=data.frame(n = ff(118474,96201)/1000,pheno.best = current_sample_optimum$pheno.best),aes(x=n, y=pheno.best, color="lightgoldenrod3"), size=6)+    
  geom_point(data=data.frame(n = ff(118474,96201)/1000,pheno.best = current_sample_gwas$pheno.best),aes(x=n, y=pheno.best, color="lightgoldenrod3"), size=6)+    
  geom_point(data=data.frame(n = 2*ff(118474,96201)/1000,pheno.best = double_sample_optimum$pheno.best),aes(x=n, y=pheno.best,color="palevioletred"), size=6)+    
  geom_point(data=data.frame(n = 2*ff(118474,96201)/1000,pheno.best = double_sample_gwas$pheno.best),aes(x=n, y=pheno.best,color="palevioletred"), size=6)+    
  geom_point(data=data.frame(n = 4*ff(118474,96201)/1000,pheno.best = quadruple_sample_optimum$pheno.best),aes(x=n, y=pheno.best,color="royalblue1"), size=6)+    
  geom_point(data=data.frame(n = 4*ff(118474,96201)/1000,pheno.best = quadruple_sample_gwas$pheno.best),aes(x=n, y=pheno.best,color="royalblue1"), size=6)+    
  theme_bw() + scale_color_manual(labels = c("Current", "Double","Quadruple"),values=c("royalblue1","palevioletred","lightgoldenrod3")) +
  scale_linetype_manual(values = c('solid', 'dotted'),labels = c("Optimized p","5e-8")) + 
  labs(x = paste0("\n"),y = "\n") +
  ggtitle("EUR Unrelated") +
  theme_Publication() + 
  theme(legend.position="none")

load("/data/williamsjacr/Aim2_Genesis/AFR_HM3_05_PreExtract/Results/AFR_Unrelated_Phenotypic_Variance.RData")

p3 <- ggplot() +
  geom_line(data=results, aes(x=n_seq/1000, y=pheno.best,linetype="dotted"), size=2.5,color="grey38") +
  geom_line(data=results, aes(x=n_seq/1000, y=pheno.bestGWAS,linetype="solid"), size=2.5,color="grey38") +
  geom_point(data=data.frame(n = ff(9235,10184)/1000,pheno.best = current_sample_optimum$pheno.best),aes(x=n, y=pheno.best, color="lightgoldenrod3"), size=6)+    
  geom_point(data=data.frame(n = ff(9235,10184)/1000,pheno.best = current_sample_gwas$pheno.best),aes(x=n, y=pheno.best, color="lightgoldenrod3"), size=6)+    
  geom_point(data=data.frame(n = 2*ff(9235,10184)/1000,pheno.best = double_sample_optimum$pheno.best),aes(x=n, y=pheno.best,color="palevioletred"), size=6)+    
  geom_point(data=data.frame(n = 2*ff(9235,10184)/1000,pheno.best = double_sample_gwas$pheno.best),aes(x=n, y=pheno.best,color="palevioletred"), size=6)+    
  geom_point(data=data.frame(n = 4*ff(9235,10184)/1000,pheno.best = quadruple_sample_optimum$pheno.best),aes(x=n, y=pheno.best,color="royalblue1"), size=6)+    
  geom_point(data=data.frame(n = 4*ff(9235,10184)/1000,pheno.best = quadruple_sample_gwas$pheno.best),aes(x=n, y=pheno.best,color="royalblue1"), size=6)+    
  theme_bw() + scale_color_manual(labels = c("Current", "Double","Quadruple"),values=c("royalblue1","palevioletred","lightgoldenrod3")) +
  scale_linetype_manual(values = c('solid', 'dotted'),labels = c("Optimized p","5e-8")) + 
  labs(x = paste0("\n"),y = "\n") +
  ggtitle("AA Unrelated") +
  theme_Publication() + 
  theme(legend.position="none")

load("/data/williamsjacr/Aim2_Genesis/AFR_HM3_05_AABCG/Results/AABCG_Phenotypic_Variance.RData")

p4 <- ggplot() +
  geom_line(data=results, aes(x=n_seq/1000, y=pheno.best,linetype="dotted"), size=2.5,color="grey38") +
  geom_line(data=results, aes(x=n_seq/1000, y=pheno.bestGWAS,linetype="solid"), size=2.5,color="grey38") +
  geom_point(data=data.frame(n = ff(9235,10184)/1000,pheno.best = current_sample_optimum$pheno.best),aes(x=n, y=pheno.best, color="lightgoldenrod3"), size=6)+    
  geom_point(data=data.frame(n = ff(9235,10184)/1000,pheno.best = current_sample_gwas$pheno.best),aes(x=n, y=pheno.best, color="lightgoldenrod3"), size=6)+    
  geom_point(data=data.frame(n = 2*ff(9235,10184)/1000,pheno.best = double_sample_optimum$pheno.best),aes(x=n, y=pheno.best,color="palevioletred"), size=6)+    
  geom_point(data=data.frame(n = 2*ff(9235,10184)/1000,pheno.best = double_sample_gwas$pheno.best),aes(x=n, y=pheno.best,color="palevioletred"), size=6)+    
  geom_point(data=data.frame(n = 4*ff(9235,10184)/1000,pheno.best = quadruple_sample_optimum$pheno.best),aes(x=n, y=pheno.best,color="royalblue1"), size=6)+    
  geom_point(data=data.frame(n = 4*ff(9235,10184)/1000,pheno.best = quadruple_sample_gwas$pheno.best),aes(x=n, y=pheno.best,color="royalblue1"), size=6)+    
  theme_bw() + scale_color_manual(labels = c("Current", "Double","Quadruple"),values=c("royalblue1","palevioletred","lightgoldenrod3")) +
  scale_linetype_manual(values = c('solid', 'dotted'),labels = c("Optimized p","5e-8")) + 
  labs(x = paste0("\n"),y = "\n") +
  ggtitle("AA") +
  theme_Publication() + 
  theme(legend.position="none")

load("/data/williamsjacr/Aim2_Genesis/EAS_HM3_05_PreExtract/Results/EAS_Full1KG_Phenotypic_Variance.RData")

p5 <- ggplot() +
  geom_line(data=results, aes(x=n_seq/1000, y=pheno.best,linetype="dotted"), size=2.5,color="grey38") +
  geom_line(data=results, aes(x=n_seq/1000, y=pheno.bestGWAS,linetype="solid"), size=2.5,color="grey38") +
  geom_point(data=data.frame(n = ff(20393,86329)/1000,pheno.best = current_sample_optimum$pheno.best),aes(x=n, y=pheno.best, color="lightgoldenrod3"), size=6)+    
  geom_point(data=data.frame(n = ff(20393,86329)/1000,pheno.best = current_sample_gwas$pheno.best),aes(x=n, y=pheno.best, color="lightgoldenrod3"), size=6)+    
  geom_point(data=data.frame(n = 2*ff(20393,86329)/1000,pheno.best = double_sample_optimum$pheno.best),aes(x=n, y=pheno.best,color="palevioletred"), size=6)+    
  geom_point(data=data.frame(n = 2*ff(20393,86329)/1000,pheno.best = double_sample_gwas$pheno.best),aes(x=n, y=pheno.best,color="palevioletred"), size=6)+    
  geom_point(data=data.frame(n = 4*ff(20393,86329)/1000,pheno.best = quadruple_sample_optimum$pheno.best),aes(x=n, y=pheno.best,color="royalblue1"), size=6)+    
  geom_point(data=data.frame(n = 4*ff(20393,86329)/1000,pheno.best = quadruple_sample_gwas$pheno.best),aes(x=n, y=pheno.best,color="royalblue1"), size=6)+    
  theme_bw() + scale_color_manual(labels = c("Current", "Double","Quadruple"),values=c("royalblue1","palevioletred","lightgoldenrod3")) +
  scale_linetype_manual(values = c('solid', 'dotted'),labels = c("Optimized p","5e-8")) + 
  labs(x = paste0("\n"),y = "\n") +
  ggtitle("EAS") +
  theme_Publication() + 
  theme(legend.position="none")

load("/data/williamsjacr/Aim2_Genesis/EAS_HM3_05_PreExtract/Results/EAS_Unrelated_Phenotypic_Variance.RData")

p6 <- ggplot() +
  geom_line(data=results, aes(x=n_seq/1000, y=pheno.best,linetype="dotted"), size=2.5,color="grey38") +
  geom_line(data=results, aes(x=n_seq/1000, y=pheno.bestGWAS,linetype="solid"), size=2.5,color="grey38") +
  geom_point(data=data.frame(n = ff(20393,86329)/1000,pheno.best = current_sample_optimum$pheno.best),aes(x=n, y=pheno.best, color="lightgoldenrod3"), size=6)+    
  geom_point(data=data.frame(n = ff(20393,86329)/1000,pheno.best = current_sample_gwas$pheno.best),aes(x=n, y=pheno.best, color="lightgoldenrod3"), size=6)+    
  geom_point(data=data.frame(n = 2*ff(20393,86329)/1000,pheno.best = double_sample_optimum$pheno.best),aes(x=n, y=pheno.best,color="palevioletred"), size=6)+    
  geom_point(data=data.frame(n = 2*ff(20393,86329)/1000,pheno.best = double_sample_gwas$pheno.best),aes(x=n, y=pheno.best,color="palevioletred"), size=6)+    
  geom_point(data=data.frame(n = 4*ff(20393,86329)/1000,pheno.best = quadruple_sample_optimum$pheno.best),aes(x=n, y=pheno.best,color="royalblue1"), size=6)+    
  geom_point(data=data.frame(n = 4*ff(20393,86329)/1000,pheno.best = quadruple_sample_gwas$pheno.best),aes(x=n, y=pheno.best,color="royalblue1"), size=6)+    
  theme_bw() + scale_color_manual(labels = c("Current", "Double","Quadruple"),values=c("royalblue1","palevioletred","lightgoldenrod3")) +
  scale_linetype_manual(values = c('solid', 'dotted'),labels = c("Optimized p","5e-8")) + 
  labs(x = paste0("\n"),y = "\n") +
  ggtitle("EAS Unrelated") +
  theme_Publication() + 
  theme(legend.position="none")

load("/data/williamsjacr/Aim2_Genesis/AMR_HM3_05_PreExtract/Results/AMR_Unrelated_Phenotypic_Variance.RData")

p7 <- ggplot() +
  geom_line(data=results, aes(x=n_seq/1000, y=pheno.best,linetype="dotted"), size=2.5,color="grey38") +
  geom_line(data=results, aes(x=n_seq/1000, y=pheno.bestGWAS,linetype="solid"), size=2.5,color="grey38") +
  geom_point(data=data.frame(n = ff(2396,7468)/1000,pheno.best = current_sample_optimum$pheno.best),aes(x=n, y=pheno.best, color="lightgoldenrod3"), size=6)+    
  geom_point(data=data.frame(n = ff(2396,7468)/1000,pheno.best = current_sample_gwas$pheno.best),aes(x=n, y=pheno.best, color="lightgoldenrod3"), size=6)+    
  geom_point(data=data.frame(n = 2*ff(2396,7468)/1000,pheno.best = double_sample_optimum$pheno.best),aes(x=n, y=pheno.best,color="palevioletred"), size=6)+    
  geom_point(data=data.frame(n = 2*ff(2396,7468)/1000,pheno.best = double_sample_gwas$pheno.best),aes(x=n, y=pheno.best,color="palevioletred"), size=6)+    
  geom_point(data=data.frame(n = 4*ff(2396,7468)/1000,pheno.best = quadruple_sample_optimum$pheno.best),aes(x=n, y=pheno.best,color="royalblue1"), size=6)+    
  geom_point(data=data.frame(n = 4*ff(2396,7468)/1000,pheno.best = quadruple_sample_gwas$pheno.best),aes(x=n, y=pheno.best,color="royalblue1"), size=6)+    
  theme_bw() + scale_color_manual(labels = c("Current", "Double","Quadruple"),values=c("royalblue1","palevioletred","lightgoldenrod3")) +
  scale_linetype_manual(values = c('solid', 'dotted'),labels = c("Optimized p","5e-8")) + 
  labs(x = paste0("\n"),y = "\n") +
  ggtitle("AMR") +
  theme_Publication() + 
  theme(legend.position="none")

load("/data/williamsjacr/Aim2_Genesis/AMR_HM3_05_PreExtract/Results/AMR_Unrelated_Phenotypic_VarianceFit2.RData")

p8 <- ggplot() +
  geom_line(data=results, aes(x=n_seq/1000, y=pheno.best,linetype="dotted"), size=2.5,color="grey38") +
  geom_line(data=results, aes(x=n_seq/1000, y=pheno.bestGWAS,linetype="solid"), size=2.5,color="grey38") +
  geom_point(data=data.frame(n = ff(2396,7468)/1000,pheno.best = current_sample_optimum$pheno.best),aes(x=n, y=pheno.best, color="lightgoldenrod3"), size=6)+    
  geom_point(data=data.frame(n = ff(2396,7468)/1000,pheno.best = current_sample_gwas$pheno.best),aes(x=n, y=pheno.best, color="lightgoldenrod3"), size=6)+    
  geom_point(data=data.frame(n = 2*ff(2396,7468)/1000,pheno.best = double_sample_optimum$pheno.best),aes(x=n, y=pheno.best,color="palevioletred"), size=6)+    
  geom_point(data=data.frame(n = 2*ff(2396,7468)/1000,pheno.best = double_sample_gwas$pheno.best),aes(x=n, y=pheno.best,color="palevioletred"), size=6)+    
  geom_point(data=data.frame(n = 4*ff(2396,7468)/1000,pheno.best = quadruple_sample_optimum$pheno.best),aes(x=n, y=pheno.best,color="royalblue1"), size=6)+    
  geom_point(data=data.frame(n = 4*ff(2396,7468)/1000,pheno.best = quadruple_sample_gwas$pheno.best),aes(x=n, y=pheno.best,color="royalblue1"), size=6)+    
  theme_bw() + scale_color_manual(labels = c("Current", "Double","Quadruple"),values=c("royalblue1","palevioletred","lightgoldenrod3")) +
  scale_linetype_manual(values = c('solid', 'dotted'),labels = c("Optimized p","5e-8")) + 
  labs(x = paste0("\n"),y = "\n") +
  ggtitle("AMR") +
  theme_Publication() + 
  theme(legend.position="none")


pdf(paste0("/data/williamsjacr/Aim2_Genesis/Phenotypic_Variance.pdf"), width=2*10, height=2*6.1804697157)
# do.call(grid.arrange,c(lgv,ncol=4,as.table=T))
grid.arrange(p1,p4,p5,p7, nrow = 2,bottom=textGrob("x: Total sample size assuming 1:1 case:control ratio (in thousands)", vjust=0,gp=gpar(fontsize=25)),
             left=textGrob("y: Percentage of phenotypic variance explained", just=c(0.5,1.5), gp=gpar(fontsize=25), rot=90))

dev.off()



load("/data/williamsjacr/Aim2_Genesis/EUR_HM3_05_PreExtract/Results/EUR_Full1KG_Phenotypic_Variance.RData")

p1 <- ggplot() +
  geom_line(data=results, aes(x=n_seq/1000, y=Genetic_Variance_Optimum,linetype="dotted"), size=2.5,color="grey38") +
  geom_line(data=results, aes(x=n_seq/1000, y=Genetic_Variance_GWAS,linetype="solid"), size=2.5,color="grey38") +
  geom_point(data=data.frame(n = ff(118474,96201)/1000,gv.best = current_sample_optimum$gv.best),aes(x=n, y=gv.best,color="lightgoldenrod3"), size=6)+    
  geom_point(data=data.frame(n = ff(118474,96201)/1000,gv.best = current_sample_gwas$gv.best),aes(x=n, y=gv.best,color="lightgoldenrod3"), size=6)+    
  geom_point(data=data.frame(n = 2*ff(118474,96201)/1000,gv.best = double_sample_optimum$gv.best),aes(x=n, y=gv.best,color="palevioletred"), size=6)+    
  geom_point(data=data.frame(n = 2*ff(118474,96201)/1000,gv.best = double_sample_gwas$gv.best),aes(x=n, y=gv.best,color="palevioletred"), size=6)+    
  geom_point(data=data.frame(n = 4*ff(118474,96201)/1000,gv.best = quadruple_sample_optimum$gv.best),aes(x=n, y=gv.best,color="royalblue1"), size=6)+    
  geom_point(data=data.frame(n = 4*ff(118474,96201)/1000,gv.best = quadruple_sample_gwas$gv.best),aes(x=n, y=gv.best,color="royalblue1"), size=6)+    
  theme_bw() + scale_color_manual(labels = c("Current", "Double","Quadruple"),values=c("royalblue1","palevioletred","lightgoldenrod3")) +
  scale_linetype_manual(values = c('solid', 'dotted'),labels = c("Optimized p","5e-8")) + 
  labs(x = paste0("\n"),y = "\n") +
  ggtitle("EUR") +
  theme_Publication()

load("/data/williamsjacr/Aim2_Genesis/EUR_HM3_05_PreExtract/Results/EUR_Unrelated_Phenotypic_Variance.RData")

p2 <- ggplot() +
  geom_line(data=results, aes(x=n_seq/1000, y=Genetic_Variance_Optimum,linetype="dotted"), size=2.5,color="grey38") +
  geom_line(data=results, aes(x=n_seq/1000, y=Genetic_Variance_GWAS,linetype="solid"), size=2.5,color="grey38") +
  geom_point(data=data.frame(n = ff(118474,96201)/1000,gv.best = current_sample_optimum$gv.best),aes(x=n, y=gv.best,color="lightgoldenrod3"), size=6)+    
  geom_point(data=data.frame(n = ff(118474,96201)/1000,gv.best = current_sample_gwas$gv.best),aes(x=n, y=gv.best,color="lightgoldenrod3"), size=6)+    
  geom_point(data=data.frame(n = 2*ff(118474,96201)/1000,gv.best = double_sample_optimum$gv.best),aes(x=n, y=gv.best,color="palevioletred"), size=6)+    
  geom_point(data=data.frame(n = 2*ff(118474,96201)/1000,gv.best = double_sample_gwas$gv.best),aes(x=n, y=gv.best,color="palevioletred"), size=6)+    
  geom_point(data=data.frame(n = 4*ff(118474,96201)/1000,gv.best = quadruple_sample_optimum$gv.best),aes(x=n, y=gv.best,color="royalblue1"), size=6)+    
  geom_point(data=data.frame(n = 4*ff(118474,96201)/1000,gv.best = quadruple_sample_gwas$gv.best),aes(x=n, y=gv.best,color="royalblue1"), size=6)+    
  theme_bw() + scale_color_manual(labels = c("Current", "Double","Quadruple"),values=c("royalblue1","palevioletred","lightgoldenrod3")) +
  scale_linetype_manual(values = c('solid', 'dotted'),labels = c("Optimized p","5e-8")) + 
  labs(x = paste0("\n"),y = "\n") +
  ggtitle("EUR Unrelated") +
  theme_Publication() + 
  theme(legend.position="none")

load("/data/williamsjacr/Aim2_Genesis/AFR_HM3_05_PreExtract/Results/AFR_Unrelated_Phenotypic_Variance.RData")

p3 <- ggplot() +
  geom_line(data=results, aes(x=n_seq/1000, y=Genetic_Variance_Optimum,linetype="dotted"), size=2.5,color="grey38") +
  geom_line(data=results, aes(x=n_seq/1000, y=Genetic_Variance_GWAS,linetype="solid"), size=2.5,color="grey38") +
  geom_point(data=data.frame(n = ff(9235,10184)/1000,gv.best = current_sample_optimum$gv.best),aes(x=n, y=gv.best,color="lightgoldenrod3"), size=6)+    
  geom_point(data=data.frame(n = ff(9235,10184)/1000,gv.best = current_sample_gwas$gv.best),aes(x=n, y=gv.best,color="lightgoldenrod3"), size=6)+    
  geom_point(data=data.frame(n = 2*ff(9235,10184)/1000,gv.best = double_sample_optimum$gv.best),aes(x=n, y=gv.best,color="palevioletred"), size=6)+    
  geom_point(data=data.frame(n = 2*ff(9235,10184)/1000,gv.best = double_sample_gwas$gv.best),aes(x=n, y=gv.best,color="palevioletred"), size=6)+    
  geom_point(data=data.frame(n = 4*ff(9235,10184)/1000,gv.best = quadruple_sample_optimum$gv.best),aes(x=n, y=gv.best,color="royalblue1"), size=6)+    
  geom_point(data=data.frame(n = 4*ff(9235,10184)/1000,gv.best = quadruple_sample_gwas$gv.best),aes(x=n, y=gv.best,color="royalblue1"), size=6)+    
  theme_bw() + scale_color_manual(labels = c("Current", "Double","Quadruple"),values=c("royalblue1","palevioletred","lightgoldenrod3")) +
  scale_linetype_manual(values = c('solid', 'dotted'),labels = c("Optimized p","5e-8")) + 
  labs(x = paste0("\n"),y = "\n") +
  ggtitle("AA Unrelated") +
  theme_Publication() + 
  theme(legend.position="none")

load("/data/williamsjacr/Aim2_Genesis/AFR_HM3_05_AABCG/Results/AABCG_Phenotypic_Variance.RData")

p4 <- ggplot() +
  geom_line(data=results, aes(x=n_seq/1000, y=Genetic_Variance_Optimum,linetype="dotted"), size=2.5,color="grey38") +
  geom_line(data=results, aes(x=n_seq/1000, y=Genetic_Variance_GWAS,linetype="solid"), size=2.5,color="grey38") +
  geom_point(data=data.frame(n = ff(9235,10184)/1000,gv.best = current_sample_optimum$gv.best),aes(x=n, y=gv.best,color="lightgoldenrod3"), size=6)+    
  geom_point(data=data.frame(n = ff(9235,10184)/1000,gv.best = current_sample_gwas$gv.best),aes(x=n, y=gv.best,color="lightgoldenrod3"), size=6)+    
  geom_point(data=data.frame(n = 2*ff(9235,10184)/1000,gv.best = double_sample_optimum$gv.best),aes(x=n, y=gv.best,color="palevioletred"), size=6)+    
  geom_point(data=data.frame(n = 2*ff(9235,10184)/1000,gv.best = double_sample_gwas$gv.best),aes(x=n, y=gv.best,color="palevioletred"), size=6)+    
  geom_point(data=data.frame(n = 4*ff(9235,10184)/1000,gv.best = quadruple_sample_optimum$gv.best),aes(x=n, y=gv.best,color="royalblue1"), size=6)+    
  geom_point(data=data.frame(n = 4*ff(9235,10184)/1000,gv.best = quadruple_sample_gwas$gv.best),aes(x=n, y=gv.best,color="royalblue1"), size=6)+    
  theme_bw() + scale_color_manual(labels = c("Current", "Double","Quadruple"),values=c("royalblue1","palevioletred","lightgoldenrod3")) +
  scale_linetype_manual(values = c('solid', 'dotted'),labels = c("Optimized p","5e-8")) + 
  labs(x = paste0("\n"),y = "\n") +
  ggtitle("AA") +
  theme_Publication() + 
  theme(legend.position="none")

load("/data/williamsjacr/Aim2_Genesis/EAS_HM3_05_PreExtract/Results/EAS_Full1KG_Phenotypic_Variance.RData")

p5 <- ggplot() +
  geom_line(data=results, aes(x=n_seq/1000, y=Genetic_Variance_Optimum,linetype="dotted"), size=2.5,color="grey38") +
  geom_line(data=results, aes(x=n_seq/1000, y=Genetic_Variance_GWAS,linetype="solid"), size=2.5,color="grey38") +
  geom_point(data=data.frame(n = ff(20393,86329)/1000,gv.best = current_sample_optimum$gv.best),aes(x=n, y=gv.best,color="lightgoldenrod3"), size=6)+    
  geom_point(data=data.frame(n = ff(20393,86329)/1000,gv.best = current_sample_gwas$gv.best),aes(x=n, y=gv.best,color="lightgoldenrod3"), size=6)+    
  geom_point(data=data.frame(n = 2*ff(20393,86329)/1000,gv.best = double_sample_optimum$gv.best),aes(x=n, y=gv.best,color="palevioletred"), size=6)+    
  geom_point(data=data.frame(n = 2*ff(20393,86329)/1000,gv.best = double_sample_gwas$gv.best),aes(x=n, y=gv.best,color="palevioletred"), size=6)+    
  geom_point(data=data.frame(n = 4*ff(20393,86329)/1000,gv.best = quadruple_sample_optimum$gv.best),aes(x=n, y=gv.best,color="royalblue1"), size=6)+    
  geom_point(data=data.frame(n = 4*ff(20393,86329)/1000,gv.best = quadruple_sample_gwas$gv.best),aes(x=n, y=gv.best,color="royalblue1"), size=6)+    
  theme_bw() + scale_color_manual(labels = c("Current", "Double","Quadruple"),values=c("royalblue1","palevioletred","lightgoldenrod3")) +
  scale_linetype_manual(values = c('solid', 'dotted'),labels = c("Optimized p","5e-8")) + 
  labs(x = paste0("\n"),y = "\n") +
  ggtitle("EAS") +
  theme_Publication() + 
  theme(legend.position="none")

load("/data/williamsjacr/Aim2_Genesis/EAS_HM3_05_PreExtract/Results/EAS_Unrelated_Phenotypic_Variance.RData")

p6 <- ggplot() +
  geom_line(data=results, aes(x=n_seq/1000, y=Genetic_Variance_Optimum,linetype="dotted"), size=2.5,color="grey38") +
  geom_line(data=results, aes(x=n_seq/1000, y=Genetic_Variance_GWAS,linetype="solid"), size=2.5,color="grey38") +
  geom_point(data=data.frame(n = ff(20393,86329)/1000,gv.best = current_sample_optimum$gv.best),aes(x=n, y=gv.best,color="lightgoldenrod3"), size=6)+    
  geom_point(data=data.frame(n = ff(20393,86329)/1000,gv.best = current_sample_gwas$gv.best),aes(x=n, y=gv.best,color="lightgoldenrod3"), size=6)+    
  geom_point(data=data.frame(n = 2*ff(20393,86329)/1000,gv.best = double_sample_optimum$gv.best),aes(x=n, y=gv.best,color="palevioletred"), size=6)+    
  geom_point(data=data.frame(n = 2*ff(20393,86329)/1000,gv.best = double_sample_gwas$gv.best),aes(x=n, y=gv.best,color="palevioletred"), size=6)+    
  geom_point(data=data.frame(n = 4*ff(20393,86329)/1000,gv.best = quadruple_sample_optimum$gv.best),aes(x=n, y=gv.best,color="royalblue1"), size=6)+    
  geom_point(data=data.frame(n = 4*ff(20393,86329)/1000,gv.best = quadruple_sample_gwas$gv.best),aes(x=n, y=gv.best,color="royalblue1"), size=6)+    
  theme_bw() + scale_color_manual(labels = c("Current", "Double","Quadruple"),values=c("royalblue1","palevioletred","lightgoldenrod3")) +
  scale_linetype_manual(values = c('solid', 'dotted'),labels = c("Optimized p","5e-8")) + 
  labs(x = paste0("\n"),y = "\n") +
  ggtitle("EAS Unrelated") +
  theme_Publication() + 
  theme(legend.position="none")

load("/data/williamsjacr/Aim2_Genesis/AMR_HM3_05_PreExtract/Results/AMR_Unrelated_Phenotypic_Variance.RData")

p7 <- ggplot() +
  geom_line(data=results, aes(x=n_seq/1000, y=Genetic_Variance_Optimum,linetype="dotted"), size=2.5,color="grey38") +
  geom_line(data=results, aes(x=n_seq/1000, y=Genetic_Variance_GWAS,linetype="solid"), size=2.5,color="grey38") +
  geom_point(data=data.frame(n = ff(2396,7468)/1000,gv.best = current_sample_optimum$gv.best),aes(x=n, y=gv.best, color="lightgoldenrod3"), size=6)+    
  geom_point(data=data.frame(n = ff(2396,7468)/1000,gv.best = current_sample_gwas$gv.best),aes(x=n, y=gv.best, color="lightgoldenrod3"), size=6)+    
  geom_point(data=data.frame(n = 2*ff(2396,7468)/1000,gv.best = double_sample_optimum$gv.best),aes(x=n, y=gv.best,color="palevioletred"), size=6)+    
  geom_point(data=data.frame(n = 2*ff(2396,7468)/1000,gv.best = double_sample_gwas$gv.best),aes(x=n, y=gv.best,color="palevioletred"), size=6)+    
  geom_point(data=data.frame(n = 4*ff(2396,7468)/1000,gv.best = quadruple_sample_optimum$gv.best),aes(x=n, y=gv.best,color="royalblue1"), size=6)+    
  geom_point(data=data.frame(n = 4*ff(2396,7468)/1000,gv.best = quadruple_sample_gwas$gv.best),aes(x=n, y=gv.best,color="royalblue1"), size=6)+    
  theme_bw() + scale_color_manual(labels = c("Current", "Double","Quadruple"),values=c("royalblue1","palevioletred","lightgoldenrod3")) +
  scale_linetype_manual(values = c('solid', 'dotted'),labels = c("Optimized p","5e-8")) + 
  labs(x = paste0("\n"),y = "\n") +
  ggtitle("AMR") +
  theme_Publication() + 
  theme(legend.position="none")

load("/data/williamsjacr/Aim2_Genesis/AMR_HM3_05_PreExtract/Results/AMR_Unrelated_Phenotypic_VarianceFit2.RData")

p8 <- ggplot() +
  geom_line(data=results, aes(x=n_seq/1000, y=Genetic_Variance_Optimum,linetype="dotted"), size=2.5,color="grey38") +
  geom_line(data=results, aes(x=n_seq/1000, y=Genetic_Variance_GWAS,linetype="solid"), size=2.5,color="grey38") +
  geom_point(data=data.frame(n = ff(2396,7468)/1000,gv.best = current_sample_optimum$gv.best),aes(x=n, y=gv.best, color="lightgoldenrod3"), size=6)+    
  geom_point(data=data.frame(n = ff(2396,7468)/1000,gv.best = current_sample_gwas$gv.best),aes(x=n, y=gv.best, color="lightgoldenrod3"), size=6)+    
  geom_point(data=data.frame(n = 2*ff(2396,7468)/1000,gv.best = double_sample_optimum$gv.best),aes(x=n, y=gv.best,color="palevioletred"), size=6)+    
  geom_point(data=data.frame(n = 2*ff(2396,7468)/1000,gv.best = double_sample_gwas$gv.best),aes(x=n, y=gv.best,color="palevioletred"), size=6)+    
  geom_point(data=data.frame(n = 4*ff(2396,7468)/1000,gv.best = quadruple_sample_optimum$gv.best),aes(x=n, y=gv.best,color="royalblue1"), size=6)+    
  geom_point(data=data.frame(n = 4*ff(2396,7468)/1000,gv.best = quadruple_sample_gwas$gv.best),aes(x=n, y=gv.best,color="royalblue1"), size=6)+    
  theme_bw() + scale_color_manual(labels = c("Current", "Double","Quadruple"),values=c("royalblue1","palevioletred","lightgoldenrod3")) +
  scale_linetype_manual(values = c('solid', 'dotted'),labels = c("Optimized p","5e-8")) + 
  labs(x = paste0("\n"),y = "\n") +
  ggtitle("AMR") +
  theme_Publication() + 
  theme(legend.position="none")

pdf(paste0("/data/williamsjacr/Aim2_Genesis/Genetic_Variance.pdf"), width=2*10, height=2*6.1804697157)
# do.call(grid.arrange,c(lgv,ncol=4,as.table=T))
grid.arrange(p1,p4,p5,p7, nrow = 2,bottom=textGrob("x: Total sample size assuming 1:1 case:control ratio (in thousands)", vjust=0,gp=gpar(fontsize=25)),
             left=textGrob("y: Percentage of genetic variance explained", just=c(0.5,1.5), gp=gpar(fontsize=25), rot=90))

dev.off()


######### Figure 3

load("/data/williamsjacr/Aim2_Genesis/EUR_HM3_05_PreExtract/Results/EUR_Full1KG_Phenotypic_Variance.RData")

p1 <- ggplot() +
  geom_line(data=results, aes(x=n_seq/1000, y=AUC_Optimum,linetype="dotted"), size=2.5,color="grey38") +
  geom_line(data=results, aes(x=n_seq/1000, y=AUC_GWAS,linetype="solid"), size=2.5,color="grey38") +
  geom_point(data=data.frame(n = ff(118474,96201)/1000,auc.best = current_sample_optimum$auc.best),aes(x=n, y=auc.best,color="lightgoldenrod3"), size=6)+    
  geom_point(data=data.frame(n = ff(118474,96201)/1000,auc.best = current_sample_gwas$auc.best),aes(x=n, y=auc.best,color="lightgoldenrod3"), size=6)+    
  geom_point(data=data.frame(n = 2*ff(118474,96201)/1000,auc.best = double_sample_optimum$auc.best),aes(x=n, y=auc.best,color="palevioletred"), size=6)+    
  geom_point(data=data.frame(n = 2*ff(118474,96201)/1000,auc.best = double_sample_gwas$auc.best),aes(x=n, y=auc.best,color="palevioletred"), size=6)+    
  geom_point(data=data.frame(n = 4*ff(118474,96201)/1000,auc.best = quadruple_sample_optimum$auc.best),aes(x=n, y=auc.best,color="royalblue1"), size=6)+    
  geom_point(data=data.frame(n = 4*ff(118474,96201)/1000,auc.best = quadruple_sample_gwas$auc.best),aes(x=n, y=auc.best,color="royalblue1"), size=6)+    
  theme_bw() + scale_color_manual(labels = c("Current", "Double","Quadruple"),values=c("royalblue1","palevioletred","lightgoldenrod3")) +
  scale_linetype_manual(values = c('solid', 'dotted'),labels = c("Optimized p","5e-8")) + 
  labs(x = paste0("\n"),y = "\n") +
  ggtitle("EUR") +
  theme_Publication()

load("/data/williamsjacr/Aim2_Genesis/EUR_HM3_05_PreExtract/Results/EUR_Unrelated_Phenotypic_Variance.RData")

p2 <- ggplot() +
  geom_line(data=results, aes(x=n_seq/1000, y=AUC_Optimum,linetype="dotted"), size=2.5,color="grey38") +
  geom_line(data=results, aes(x=n_seq/1000, y=AUC_GWAS,linetype="solid"), size=2.5,color="grey38") +
  geom_point(data=data.frame(n = ff(118474,96201)/1000,auc.best = current_sample_optimum$auc.best),aes(x=n, y=auc.best,color="lightgoldenrod3"), size=6)+    
  geom_point(data=data.frame(n = ff(118474,96201)/1000,auc.best = current_sample_gwas$auc.best),aes(x=n, y=auc.best,color="lightgoldenrod3"), size=6)+    
  geom_point(data=data.frame(n = 2*ff(118474,96201)/1000,auc.best = double_sample_optimum$auc.best),aes(x=n, y=auc.best,color="palevioletred"), size=6)+    
  geom_point(data=data.frame(n = 2*ff(118474,96201)/1000,auc.best = double_sample_gwas$auc.best),aes(x=n, y=auc.best,color="palevioletred"), size=6)+    
  geom_point(data=data.frame(n = 4*ff(118474,96201)/1000,auc.best = quadruple_sample_optimum$auc.best),aes(x=n, y=auc.best,color="royalblue1"), size=6)+    
  geom_point(data=data.frame(n = 4*ff(118474,96201)/1000,auc.best = quadruple_sample_gwas$auc.best),aes(x=n, y=auc.best,color="royalblue1"), size=6)+    
  theme_bw() + scale_color_manual(labels = c("Current", "Double","Quadruple"),values=c("royalblue1","palevioletred","lightgoldenrod3")) +
  scale_linetype_manual(values = c('solid', 'dotted'),labels = c("Optimized p","5e-8")) + 
  labs(x = paste0("\n"),y = "\n") +
  ggtitle("EUR Unrelated") +
  theme_Publication() + 
  theme(legend.position="none")

load("/data/williamsjacr/Aim2_Genesis/AFR_HM3_05_PreExtract/Results/AFR_Unrelated_Phenotypic_Variance.RData")

p3 <- ggplot() +
  geom_line(data=results, aes(x=n_seq/1000, y=AUC_Optimum,linetype="dotted"), size=2.5,color="grey38") +
  geom_line(data=results, aes(x=n_seq/1000, y=AUC_GWAS,linetype="solid"), size=2.5,color="grey38") +
  geom_point(data=data.frame(n = ff(9235,10184)/1000,auc.best = current_sample_optimum$auc.best),aes(x=n, y=auc.best,color="lightgoldenrod3"), size=6)+    
  geom_point(data=data.frame(n = ff(9235,10184)/1000,auc.best = current_sample_gwas$auc.best),aes(x=n, y=auc.best,color="lightgoldenrod3"), size=6)+    
  geom_point(data=data.frame(n = 2*ff(9235,10184)/1000,auc.best = double_sample_optimum$auc.best),aes(x=n, y=auc.best,color="palevioletred"), size=6)+    
  geom_point(data=data.frame(n = 2*ff(9235,10184)/1000,auc.best = double_sample_gwas$auc.best),aes(x=n, y=auc.best,color="palevioletred"), size=6)+    
  geom_point(data=data.frame(n = 4*ff(9235,10184)/1000,auc.best = quadruple_sample_optimum$auc.best),aes(x=n, y=auc.best,color="royalblue1"), size=6)+    
  geom_point(data=data.frame(n = 4*ff(9235,10184)/1000,auc.best = quadruple_sample_gwas$auc.best),aes(x=n, y=auc.best,color="royalblue1"), size=6)+    
  theme_bw() + scale_color_manual(labels = c("Current", "Double","Quadruple"),values=c("royalblue1","palevioletred","lightgoldenrod3")) +
  scale_linetype_manual(values = c('solid', 'dotted'),labels = c("Optimized p","5e-8")) + 
  labs(x = paste0("\n"),y = "\n") +
  ggtitle("AA Unrelated") +
  theme_Publication() + 
  theme(legend.position="none")

load("/data/williamsjacr/Aim2_Genesis/AFR_HM3_05_AABCG/Results/AABCG_Phenotypic_Variance.RData")

p4 <- ggplot() +
  geom_line(data=results, aes(x=n_seq/1000, y=AUC_Optimum,linetype="dotted"), size=2.5,color="grey38") +
  geom_line(data=results, aes(x=n_seq/1000, y=AUC_GWAS,linetype="solid"), size=2.5,color="grey38") +
  geom_point(data=data.frame(n = ff(9235,10184)/1000,auc.best = current_sample_optimum$auc.best),aes(x=n, y=auc.best,color="lightgoldenrod3"), size=6)+    
  geom_point(data=data.frame(n = ff(9235,10184)/1000,auc.best = current_sample_gwas$auc.best),aes(x=n, y=auc.best,color="lightgoldenrod3"), size=6)+    
  geom_point(data=data.frame(n = 2*ff(9235,10184)/1000,auc.best = double_sample_optimum$auc.best),aes(x=n, y=auc.best,color="palevioletred"), size=6)+    
  geom_point(data=data.frame(n = 2*ff(9235,10184)/1000,auc.best = double_sample_gwas$auc.best),aes(x=n, y=auc.best,color="palevioletred"), size=6)+    
  geom_point(data=data.frame(n = 4*ff(9235,10184)/1000,auc.best = quadruple_sample_optimum$auc.best),aes(x=n, y=auc.best,color="royalblue1"), size=6)+    
  geom_point(data=data.frame(n = 4*ff(9235,10184)/1000,auc.best = quadruple_sample_gwas$auc.best),aes(x=n, y=auc.best,color="royalblue1"), size=6)+    
  theme_bw() + scale_color_manual(labels = c("Current", "Double","Quadruple"),values=c("royalblue1","palevioletred","lightgoldenrod3")) +
  scale_linetype_manual(values = c('solid', 'dotted'),labels = c("Optimized p","5e-8")) + 
  labs(x = paste0("\n"),y = "\n") +
  ggtitle("AA") +
  theme_Publication() + 
  theme(legend.position="none")

load("/data/williamsjacr/Aim2_Genesis/EAS_HM3_05_PreExtract/Results/EAS_Full1KG_Phenotypic_Variance.RData")

p5 <- ggplot() +
  geom_line(data=results, aes(x=n_seq/1000, y=AUC_Optimum,linetype="dotted"), size=2.5,color="grey38") +
  geom_line(data=results, aes(x=n_seq/1000, y=AUC_GWAS,linetype="solid"), size=2.5,color="grey38") +
  geom_point(data=data.frame(n = ff(20393,86329)/1000,auc.best = current_sample_optimum$auc.best),aes(x=n, y=auc.best,color="lightgoldenrod3"), size=6)+    
  geom_point(data=data.frame(n = ff(20393,86329)/1000,auc.best = current_sample_gwas$auc.best),aes(x=n, y=auc.best,color="lightgoldenrod3"), size=6)+    
  geom_point(data=data.frame(n = 2*ff(20393,86329)/1000,auc.best = double_sample_optimum$auc.best),aes(x=n, y=auc.best,color="palevioletred"), size=6)+    
  geom_point(data=data.frame(n = 2*ff(20393,86329)/1000,auc.best = double_sample_gwas$auc.best),aes(x=n, y=auc.best,color="palevioletred"), size=6)+    
  geom_point(data=data.frame(n = 4*ff(20393,86329)/1000,auc.best = quadruple_sample_optimum$auc.best),aes(x=n, y=auc.best,color="royalblue1"), size=6)+    
  geom_point(data=data.frame(n = 4*ff(20393,86329)/1000,auc.best = quadruple_sample_gwas$auc.best),aes(x=n, y=auc.best,color="royalblue1"), size=6)+    
  theme_bw() + scale_color_manual(labels = c("Current", "Double","Quadruple"),values=c("royalblue1","palevioletred","lightgoldenrod3")) +
  scale_linetype_manual(values = c('solid', 'dotted'),labels = c("Optimized p","5e-8")) + 
  labs(x = paste0("\n"),y = "\n") +
  ggtitle("EAS") +
  theme_Publication() + 
  theme(legend.position="none")

load("/data/williamsjacr/Aim2_Genesis/EAS_HM3_05_PreExtract/Results/EAS_Unrelated_Phenotypic_Variance.RData")

p6 <- ggplot() +
  geom_line(data=results, aes(x=n_seq/1000, y=AUC_Optimum,linetype="dotted"), size=2.5,color="grey38") +
  geom_line(data=results, aes(x=n_seq/1000, y=AUC_GWAS,linetype="solid"), size=2.5,color="grey38") +
  geom_point(data=data.frame(n = ff(20393,86329)/1000,auc.best = current_sample_optimum$auc.best),aes(x=n, y=auc.best,color="lightgoldenrod3"), size=6)+    
  geom_point(data=data.frame(n = ff(20393,86329)/1000,auc.best = current_sample_gwas$auc.best),aes(x=n, y=auc.best,color="lightgoldenrod3"), size=6)+    
  geom_point(data=data.frame(n = 2*ff(20393,86329)/1000,auc.best = double_sample_optimum$auc.best),aes(x=n, y=auc.best,color="palevioletred"), size=6)+    
  geom_point(data=data.frame(n = 2*ff(20393,86329)/1000,auc.best = double_sample_gwas$auc.best),aes(x=n, y=auc.best,color="palevioletred"), size=6)+    
  geom_point(data=data.frame(n = 4*ff(20393,86329)/1000,auc.best = quadruple_sample_optimum$auc.best),aes(x=n, y=auc.best,color="royalblue1"), size=6)+    
  geom_point(data=data.frame(n = 4*ff(20393,86329)/1000,auc.best = quadruple_sample_gwas$auc.best),aes(x=n, y=auc.best,color="royalblue1"), size=6)+    
  theme_bw() + scale_color_manual(labels = c("Current", "Double","Quadruple"),values=c("royalblue1","palevioletred","lightgoldenrod3")) +
  scale_linetype_manual(values = c('solid', 'dotted'),labels = c("Optimized p","5e-8")) + 
  labs(x = paste0("\n"),y = "\n") +
  ggtitle("EAS Unrelated") +
  theme_Publication() + 
  theme(legend.position="none")

load("/data/williamsjacr/Aim2_Genesis/AMR_HM3_05_PreExtract/Results/AMR_Unrelated_Phenotypic_Variance.RData")

p7 <- ggplot() +
  geom_line(data=results, aes(x=n_seq/1000, y=AUC_Optimum,linetype="dotted"), size=2.5,color="grey38") +
  geom_line(data=results, aes(x=n_seq/1000, y=AUC_GWAS,linetype="solid"), size=2.5,color="grey38") +
  geom_point(data=data.frame(n = ff(2396,7468)/1000,auc.best = current_sample_optimum$auc.best),aes(x=n, y=auc.best, color="lightgoldenrod3"), size=6)+    
  geom_point(data=data.frame(n = ff(2396,7468)/1000,auc.best = current_sample_gwas$auc.best),aes(x=n, y=auc.best, color="lightgoldenrod3"), size=6)+    
  geom_point(data=data.frame(n = 2*ff(2396,7468)/1000,auc.best = double_sample_optimum$auc.best),aes(x=n, y=auc.best,color="palevioletred"), size=6)+    
  geom_point(data=data.frame(n = 2*ff(2396,7468)/1000,auc.best = double_sample_gwas$auc.best),aes(x=n, y=auc.best,color="palevioletred"), size=6)+    
  geom_point(data=data.frame(n = 4*ff(2396,7468)/1000,auc.best = quadruple_sample_optimum$auc.best),aes(x=n, y=auc.best,color="royalblue1"), size=6)+    
  geom_point(data=data.frame(n = 4*ff(2396,7468)/1000,auc.best = quadruple_sample_gwas$auc.best),aes(x=n, y=auc.best,color="royalblue1"), size=6)+    
  theme_bw() + scale_color_manual(labels = c("Current", "Double","Quadruple"),values=c("royalblue1","palevioletred","lightgoldenrod3")) +
  scale_linetype_manual(values = c('solid', 'dotted'),labels = c("Optimized p","5e-8")) + 
  labs(x = paste0("\n"),y = "\n") +
  ggtitle("AMR") +
  theme_Publication() + 
  theme(legend.position="none")

pdf(paste0("/data/williamsjacr/Aim2_Genesis/AUC.pdf"), width=2*10, height=2*6.1804697157)
# do.call(grid.arrange,c(lgv,ncol=4,as.table=T))
grid.arrange(p1,p4,p5,p7, nrow = 2,bottom=textGrob("x: Total sample size assuming 1:1 case:control ratio (in thousands)", vjust=0,gp=gpar(fontsize=25)),
             left=textGrob("y: AUC associated with the PRS", just=c(0.5,1.5), gp=gpar(fontsize=25), rot=90))

dev.off()
