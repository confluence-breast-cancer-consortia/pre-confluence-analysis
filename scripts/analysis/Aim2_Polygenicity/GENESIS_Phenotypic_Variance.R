rm(list = ls())
tmp <- as.numeric(commandArgs(TRUE)[1])

library(GENESIS)

if(tmp == 1){
  load("/data/williamsjacr/Aim2_Genesis/EUR_HM3_05_PreExtract/Results/fit2_fit3_Full1KG.RData")
  save_file <- "/data/williamsjacr/Aim2_Genesis/EUR_HM3_05_PreExtract/Results/EUR_Full1KG_Phenotypic_Variance.RData"
  num.cases <- 118474
  num.controls <- 96201
}else if(tmp == 2){
  load("/data/williamsjacr/Aim2_Genesis/EUR_HM3_05_PreExtract/Results/fit2_fit3_unrelated.RData")
  save_file <- "/data/williamsjacr/Aim2_Genesis/EUR_HM3_05_PreExtract/Results/EUR_Unrelated_Phenotypic_Variance.RData"
  num.cases <- 118474
  num.controls <- 96201
}else if(tmp == 3){
  load("/data/williamsjacr/Aim2_Genesis/AFR_HM3_05_PreExtract/Results/fit2_fit3_unrelated.RData")
  save_file <- "/data/williamsjacr/Aim2_Genesis/AFR_HM3_05_PreExtract/Results/AFR_Unrelated_Phenotypic_Variance.RData"
  num.cases <- 9235
  num.controls <- 10184
}else if(tmp == 4){
  load("/data/williamsjacr/Aim2_Genesis/AMR_HM3_05_PreExtract/Results/fit2_fit3_unrelated.RData")
  save_file <- "/data/williamsjacr/Aim2_Genesis/AMR_HM3_05_PreExtract/Results/AMR_Unrelated_Phenotypic_VarianceFit3.RData"
  num.cases <- 2396
  num.controls <- 7468
}else if(tmp == 5){
  load("/data/williamsjacr/Aim2_Genesis/EAS_HM3_05_PreExtract/Results/fit2_fit3_Full1KG.RData")
  save_file <- "/data/williamsjacr/Aim2_Genesis/EAS_HM3_05_PreExtract/Results/EAS_Full1KG_Phenotypic_Variance.RData"
  num.cases <- 20393
  num.controls <- 86329
}else if(tmp == 6){
  load("/data/williamsjacr/Aim2_Genesis/EAS_HM3_05_PreExtract/Results/fit2_fit3_unrelated.RData")
  save_file <- "/data/williamsjacr/Aim2_Genesis/EAS_HM3_05_PreExtract/Results/EAS_Unrelated_Phenotypic_Variance.RData"
  num.cases <- 20393
  num.controls <- 86329
}else if(tmp == 7){
  load("/data/williamsjacr/Aim2_Genesis/AFR_HM3_05_AABCG/Results/fit2_fit3_unrelated.RData")
  save_file <- "/data/williamsjacr/Aim2_Genesis/AFR_HM3_05_AABCG/Results/AABCG_Phenotypic_Variance.RData"
  num.cases <- 9235
  num.controls <- 10184
}else{
  load("/data/williamsjacr/Aim2_Genesis/AMR_HM3_05_PreExtract/Results/fit2_fit3_unrelated.RData")
  save_file <- "/data/williamsjacr/Aim2_Genesis/AMR_HM3_05_PreExtract/Results/AMR_Unrelated_Phenotypic_VarianceFit2.RData"
  num.cases <- 2396
  num.controls <- 7468
  
  ff <- function(x,y){
    x*y/(x+y)
  }
  
  est <- fit2$estimates$`Parameter (pic, sigmasq, a) estimates` # the model parameter estimates
  v <- fit2$estimates$`Covariance matrix of parameter estimates` # the covariance matrix of model parameter estimtaes
  herit0 <- as.numeric(strsplit(fit2$estimates$`Total heritability in log-odds-ratio scale (sd)`, " ")[[1]][1])
  M <- fit2$estimates$`Total number of SNPs in the GWAS study after quality control`
  
  pheno_add <- function(est,n,summarydata_OutlierIndep=NULL,gwas.significance=5e-8){
    
    if(length(est)==3)components=2
    if(length(est)==5)components=3
    
    if(components==2){
      pic = est[1]
      sig = sqrt(est[2])
      den <- function(x){return(dnorm(x/sig)/sig )}
      herit0 <- pic*M*sig^2
    }
    
    if(components==3){
      pic = est[1]
      p1 = est[2]
      s1 = sqrt(est[3])
      s2 = sqrt(est[4])
      den <- function(x){return(p1 * dnorm(x/s1)/s1 + (1-p1)*dnorm(x/s2) /s2)}
      herit0 <- pic*M*(p1*est[3] + (1-p1)*est[4])
    }
    
    tem0 <- function(x){return(x^2*den(x))}
    c_gwsignificance = abs(qnorm(gwas.significance/2))
    pow <- function(x){return(1 - pnorm(c_gwsignificance - sqrt(n)*x) + pnorm(-c_gwsignificance - sqrt(n)*x) )}
    
    
    if(is.null(summarydata_OutlierIndep)){
      pheno.var1 = 0
    }
    
    if(!is.null(summarydata_OutlierIndep)){
      beta = as.numeric(as.character(summarydata_OutlierIndep$z))/sqrt(as.numeric(as.character(summarydata_OutlierIndep$n)))
      n_OutlierIndep = length(beta)
      tau = sqrt(1/as.numeric(as.character(summarydata_OutlierIndep$n)))
      pheno.var1 = sum( (beta^2 - tau^2) * pow(beta) )
    }
    
    return(pheno.var1)
  }
  
  
  #---------------------------------------------------------------------------------------#
  # 1. calculate at current sample size---------------#
  n.effect <- ff(num.cases, num.controls)
  n.plot <- 4*n.effect
  
  obj <- polyriskpredict(N=n.effect, Ps=c(.5,.5), Sig2s=c(est[2],est[2]), M=M, M1=M*est[1], type="optimum")
  p.tem <-  as.numeric(obj$alpha)
  pheno.tem1 <- pheno_add(est,n.effect,summarydata_OutlierIndep = NULL, gwas.significance=p.tem)
  p.best <- p.tem
  pheno.best <- 2*(qnorm(as.numeric(obj$AUC)))^2 + pheno.tem1
  gv.best <- pheno.best/herit0*100
  auc.best <- pnorm(sqrt(pheno.best/2))
  
  #-----at genomewide significance levle---------------#
  objGWAS <- polyriskpredict(N=n.effect, Ps=c(.5,.5), Sig2s=c(est[2],est[2]), M=M, M1=M*est[1], type="GWAS")
  p.temGWAS <- as.numeric(objGWAS$alpha)
  pheno.tem1GWAS <- pheno_add(est,n.effect,summarydata_OutlierIndep = NULL, gwas.significance=p.temGWAS)
  p.bestGWAS <- p.temGWAS
  pheno.bestGWAS <- 2*(qnorm(as.numeric(objGWAS$AUC)))^2 + pheno.tem1GWAS
  gv.bestGWAS <- pheno.bestGWAS/herit0*100
  
  auc.bestGWAS <- pnorm(sqrt(pheno.bestGWAS/2))
  
  current_sample_optimum <- as.data.frame(cbind(n.plot/1e3, pheno.best, gv.best, auc.best))
  current_sample_gwas <- as.data.frame(cbind(n.plot/1e3,pheno.bestGWAS, gv.bestGWAS, auc.bestGWAS))
  
  #---------------------------------------------------------------------------------------#
  #---------------------------------------------------------------------------------------#
  # 2 calculate at 2*current sample size---------------#
  n.effect <- 2*ff(num.cases, num.controls)
  n.plot <- 4*n.effect
  
  obj <- polyriskpredict(N=n.effect, Ps=c(.5,.5), Sig2s=c(est[2],est[2]), M=M, M1=M*est[1], type="optimum")
  p.tem <- as.numeric(obj$alpha)
  pheno.tem1 <- pheno_add(est,n.effect,summarydata_OutlierIndep = NULL, gwas.significance=p.tem)
  p.best <- p.tem
  pheno.best <- 2*(qnorm(as.numeric(obj$AUC)))^2 + pheno.tem1
  gv.best <- pheno.best/herit0*100
  
  auc.best <- pnorm(sqrt(pheno.best/2))
  
  #-----at genomewide significance levle---------------#
  objGWAS <- polyriskpredict(N=n.effect, Ps=c(.5,.5), Sig2s=c(est[2],est[2]), M=M, M1=M*est[1], type="GWAS")
  p.temGWAS <- as.numeric(objGWAS$alpha)
  pheno.tem1GWAS <- pheno_add(est,n.effect,summarydata_OutlierIndep = NULL, gwas.significance=p.temGWAS)
  p.bestGWAS <- p.temGWAS
  pheno.bestGWAS <- 2*(qnorm(as.numeric(objGWAS$AUC)))^2 + pheno.tem1GWAS
  gv.bestGWAS <- pheno.bestGWAS/herit0*100
  
  auc.bestGWAS <- pnorm(sqrt(pheno.bestGWAS/2))
  
  double_sample_optimum <- as.data.frame(cbind(n.plot/1e3, pheno.best, gv.best, auc.best))
  double_sample_gwas <- as.data.frame(cbind(n.plot/1e3, pheno.bestGWAS, gv.bestGWAS, auc.bestGWAS))
  
  #---------------------------------------------------------------------------------------#
  #---------------------------------------------------------------------------------------#
  # 3. calculate at 4*current sample size---------------#
  n.effect <- 4*ff(num.cases, num.controls)
  n.plot <- 4*n.effect
  
  obj <- polyriskpredict(N=n.effect, Ps=c(.5,.5), Sig2s=c(est[2],est[2]), M=M, M1=M*est[1], type="optimum")
  p.tem <- as.numeric(obj$alpha)
  pheno.tem1 <- pheno_add(est,n.effect,summarydata_OutlierIndep = NULL, gwas.significance=p.tem)
  p.best <- p.tem
  pheno.best <- 2*(qnorm(as.numeric(obj$AUC)))^2 + pheno.tem1
  gv.best <- pheno.best/herit0*100
  auc.best <- pnorm(sqrt(pheno.best/2))
  
  #-----at genomewide significance levle---------------#
  objGWAS <- polyriskpredict(N=n.effect, Ps=c(.5,.5), Sig2s=c(est[2],est[2]), M=M, M1=M*est[1], type="GWAS")
  p.temGWAS <- as.numeric(objGWAS$alpha)
  pheno.tem1GWAS <- pheno_add(est,n.effect,summarydata_OutlierIndep = NULL, gwas.significance=p.temGWAS)
  p.bestGWAS <- p.temGWAS
  pheno.bestGWAS <- 2*(qnorm(as.numeric(objGWAS$AUC)))^2 + pheno.tem1GWAS
  gv.bestGWAS <- pheno.bestGWAS/herit0*100
  auc.bestGWAS <- pnorm(sqrt(pheno.bestGWAS/2))
  
  quadruple_sample_optimum <- as.data.frame(cbind(n.plot/1e3, pheno.best, gv.best, auc.best))
  quadruple_sample_gwas <- as.data.frame(cbind(n.plot/1e3, pheno.bestGWAS, gv.bestGWAS, auc.bestGWAS))
  
  phenotypic_variance <- vector()
  p.best <- vector()
  pheno.best <- vector()
  Genetic_Variance_Optimum <- vector()
  
  p.bestGWAS <- vector()
  pheno.bestGWAS <- vector()
  AUC_GWAS <- vector()
  Genetic_Variance_GWAS <- vector()
  
  n_seq <- round(seq(1000,1000000,length.out = 4000))
  count <- 1
  
  for(i in 1:length(n_seq)){
    # PRS is calculated with SNPs included at optimum p-value threshold
    predict <- polyriskpredict(N=n_seq[count], Ps=c(.5,.5), Sig2s=c(est[2],est[2]), M=M, M1=M*est[1], type="optimum")
    p.tem <-  as.numeric(predict$alpha)
    pheno.tem1 = pheno_add(est,n_seq[count],summarydata_OutlierIndep = NULL, gwas.significance=p.tem)
    p.best[count] <- p.tem
    pheno.best[count] = 2*(qnorm(as.numeric(predict$AUC)))^2 + pheno.tem1
    Genetic_Variance_Optimum[count] = pheno.best[count]/herit0*100
    
    #-----at genomewide significance levle
    predict <- polyriskpredict(N=n_seq[count], Ps=c(.5,.5), Sig2s=c(est[2],est[2]), M=M, M1=M*est[1], type="GWAS",alp.GWAS=5e-8)
    p.temGWAS <- as.numeric(predict$alpha)
    pheno.tem1GWAS <- pheno_add(est,n_seq[count],summarydata_OutlierIndep = NULL, gwas.significance=p.temGWAS)
    p.bestGWAS[count] <- p.temGWAS
    pheno.bestGWAS[count] <- 2*(qnorm(as.numeric(predict$AUC)))^2 + pheno.tem1GWAS
    Genetic_Variance_GWAS[count] <- pheno.bestGWAS[count]/herit0*100
    
    count <- count + 1
  }
  
  AUC_Optimum <- pnorm(sqrt(pheno.best/2))
  AUC_GWAS <- pnorm(sqrt(pheno.bestGWAS/2))
  
  results <- data.frame(n_seq = n_seq,pheno.best = pheno.best, pheno.bestGWAS = pheno.bestGWAS,
                        AUC_Optimum = AUC_Optimum, Genetic_Variance_Optimum = Genetic_Variance_Optimum,
                        AUC_GWAS = AUC_GWAS,Genetic_Variance_GWAS = Genetic_Variance_GWAS)
  
  save(results,current_sample_optimum,double_sample_optimum,quadruple_sample_optimum,current_sample_gwas,double_sample_gwas,quadruple_sample_gwas,file = save_file)
  
  
}

ff <- function(x,y){
  x*y/(x+y)
}

est <- fit3$estimates$`Parameter (pic, p1, sigmasq1, sigmasq2, a) estimates` # the model parameter estimates
v <- fit3$estimates$`Covariance matrix of parameter estimates` # the covariance matrix of model parameter estimtaes
herit0 <- as.numeric(strsplit(fit3$estimates$`Total heritability in log-odds-ratio scale (sd)`, " ")[[1]][1])
M <- fit3$estimates$`Total number of SNPs in the GWAS study after quality control`

pheno_add <- function(est,n,summarydata_OutlierIndep=NULL,gwas.significance=5e-8){
  
  if(length(est)==3)components=2
  if(length(est)==5)components=3
  
  if(components==2){
    pic = est[1]
    sig = sqrt(est[2])
    den <- function(x){return(dnorm(x/sig)/sig )}
    herit0 <- pic*M*sig^2
  }
  
  if(components==3){
    pic = est[1]
    p1 = est[2]
    s1 = sqrt(est[3])
    s2 = sqrt(est[4])
    den <- function(x){return(p1 * dnorm(x/s1)/s1 + (1-p1)*dnorm(x/s2) /s2)}
    herit0 <- pic*M*(p1*est[3] + (1-p1)*est[4])
  }
  
  tem0 <- function(x){return(x^2*den(x))}
  c_gwsignificance = abs(qnorm(gwas.significance/2))
  pow <- function(x){return(1 - pnorm(c_gwsignificance - sqrt(n)*x) + pnorm(-c_gwsignificance - sqrt(n)*x) )}
  
  
  if(is.null(summarydata_OutlierIndep)){
    pheno.var1 = 0
  }
  
  if(!is.null(summarydata_OutlierIndep)){
    beta = as.numeric(as.character(summarydata_OutlierIndep$z))/sqrt(as.numeric(as.character(summarydata_OutlierIndep$n)))
    n_OutlierIndep = length(beta)
    tau = sqrt(1/as.numeric(as.character(summarydata_OutlierIndep$n)))
    pheno.var1 = sum( (beta^2 - tau^2) * pow(beta) )
  }
  
  return(pheno.var1)
}


#---------------------------------------------------------------------------------------#
# 1. calculate at current sample size---------------#
n.effect <- ff(num.cases, num.controls)
n.plot <- 4*n.effect

obj <- polyriskpredict(N=n.effect, Ps=c(est[2],1-est[2]), Sig2s=c(est[3],est[4]), M=M, M1=M*est[1], type="optimum")
p.tem <-  as.numeric(obj$alpha)
pheno.tem1 <- pheno_add(est,n.effect,summarydata_OutlierIndep = NULL, gwas.significance=p.tem)
p.best <- p.tem
pheno.best <- 2*(qnorm(as.numeric(obj$AUC)))^2 + pheno.tem1
gv.best <- pheno.best/herit0*100
auc.best <- pnorm(sqrt(pheno.best/2))

#-----at genomewide significance levle---------------#
objGWAS <- polyriskpredict(N=n.effect, Ps=c(est[2],1-est[2]), Sig2s=c(est[3],est[4]), M=M, M1=M*est[1], type="GWAS")
p.temGWAS <- as.numeric(objGWAS$alpha)
pheno.tem1GWAS <- pheno_add(est,n.effect,summarydata_OutlierIndep = NULL, gwas.significance=p.temGWAS)
p.bestGWAS <- p.temGWAS
pheno.bestGWAS <- 2*(qnorm(as.numeric(objGWAS$AUC)))^2 + pheno.tem1GWAS
gv.bestGWAS <- pheno.bestGWAS/herit0*100

auc.bestGWAS <- pnorm(sqrt(pheno.bestGWAS/2))

current_sample_optimum <- as.data.frame(cbind(n.plot/1e3, pheno.best, gv.best, auc.best))
current_sample_gwas <- as.data.frame(cbind(n.plot/1e3,pheno.bestGWAS, gv.bestGWAS, auc.bestGWAS))

#---------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------#
# 2 calculate at 2*current sample size---------------#
n.effect <- 2*ff(num.cases, num.controls)
n.plot <- 4*n.effect

obj <- polyriskpredict(N=n.effect, Ps=c(est[2],1-est[2]), Sig2s=c(est[3],est[4]), M=M, M1=M*est[1], type="optimum")
p.tem <- as.numeric(obj$alpha)
pheno.tem1 <- pheno_add(est,n.effect,summarydata_OutlierIndep = NULL, gwas.significance=p.tem)
p.best <- p.tem
pheno.best <- 2*(qnorm(as.numeric(obj$AUC)))^2 + pheno.tem1
gv.best <- pheno.best/herit0*100

auc.best <- pnorm(sqrt(pheno.best/2))

#-----at genomewide significance levle---------------#
objGWAS <- polyriskpredict(N=n.effect, Ps=c(est[2],1-est[2]), Sig2s=c(est[3],est[4]), M=M, M1=M*est[1], type="GWAS")
p.temGWAS <- as.numeric(objGWAS$alpha)
pheno.tem1GWAS <- pheno_add(est,n.effect,summarydata_OutlierIndep = NULL, gwas.significance=p.temGWAS)
p.bestGWAS <- p.temGWAS
pheno.bestGWAS <- 2*(qnorm(as.numeric(objGWAS$AUC)))^2 + pheno.tem1GWAS
gv.bestGWAS <- pheno.bestGWAS/herit0*100

auc.bestGWAS <- pnorm(sqrt(pheno.bestGWAS/2))

double_sample_optimum <- as.data.frame(cbind(n.plot/1e3, pheno.best, gv.best, auc.best))
double_sample_gwas <- as.data.frame(cbind(n.plot/1e3, pheno.bestGWAS, gv.bestGWAS, auc.bestGWAS))

#---------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------#
# 3. calculate at 4*current sample size---------------#
n.effect <- 4*ff(num.cases, num.controls)
n.plot <- 4*n.effect

obj <- polyriskpredict(N=n.effect, Ps=c(est[2],1-est[2]), Sig2s=c(est[3],est[4]), M=M, M1=M*est[1], type="optimum")
p.tem <- as.numeric(obj$alpha)
pheno.tem1 <- pheno_add(est,n.effect,summarydata_OutlierIndep = NULL, gwas.significance=p.tem)
p.best <- p.tem
pheno.best <- 2*(qnorm(as.numeric(obj$AUC)))^2 + pheno.tem1
gv.best <- pheno.best/herit0*100
auc.best <- pnorm(sqrt(pheno.best/2))

#-----at genomewide significance levle---------------#
objGWAS <- polyriskpredict(N=n.effect, Ps=c(est[2],1-est[2]), Sig2s=c(est[3],est[4]), M=M, M1=M*est[1], type="GWAS")
p.temGWAS <- as.numeric(objGWAS$alpha)
pheno.tem1GWAS <- pheno_add(est,n.effect,summarydata_OutlierIndep = NULL, gwas.significance=p.temGWAS)
p.bestGWAS <- p.temGWAS
pheno.bestGWAS <- 2*(qnorm(as.numeric(objGWAS$AUC)))^2 + pheno.tem1GWAS
gv.bestGWAS <- pheno.bestGWAS/herit0*100
auc.bestGWAS <- pnorm(sqrt(pheno.bestGWAS/2))

quadruple_sample_optimum <- as.data.frame(cbind(n.plot/1e3, pheno.best, gv.best, auc.best))
quadruple_sample_gwas <- as.data.frame(cbind(n.plot/1e3, pheno.bestGWAS, gv.bestGWAS, auc.bestGWAS))







phenotypic_variance <- vector()
p.best <- vector()
pheno.best <- vector()
Genetic_Variance_Optimum <- vector()

p.bestGWAS <- vector()
pheno.bestGWAS <- vector()
AUC_GWAS <- vector()
Genetic_Variance_GWAS <- vector()

n_seq <- round(seq(1000,1000000,length.out = 4000))
n_effect_seq <- n_seq/4
count <- 1

for(i in 1:length(n_effect_seq)){
  # PRS is calculated with SNPs included at optimum p-value threshold
  predict <- polyriskpredict(N=n_effect_seq[count], Ps=c(est[2],1-est[2]), Sig2s=c(est[3],est[4]), M=M, M1=M*est[1], type="optimum")
  p.tem <-  as.numeric(predict$alpha)
  pheno.tem1 = pheno_add(est,n_effect_seq[count],summarydata_OutlierIndep = NULL, gwas.significance=p.tem)
  p.best[count] <- p.tem
  pheno.best[count] = 2*(qnorm(as.numeric(predict$AUC)))^2 + pheno.tem1
  Genetic_Variance_Optimum[count] = pheno.best[count]/herit0*100
  
  #-----at genomewide significance levle
  predict <- polyriskpredict(N=n_effect_seq[count], Ps=c(est[2],1-est[2]), Sig2s=c(est[3],est[4]), M=M, M1=M*est[1], type="GWAS",alp.GWAS=5e-8)
  p.temGWAS <- as.numeric(predict$alpha)
  pheno.tem1GWAS <- pheno_add(est,n_effect_seq[count],summarydata_OutlierIndep = NULL, gwas.significance=p.temGWAS)
  p.bestGWAS[count] <- p.temGWAS
  pheno.bestGWAS[count] <- 2*(qnorm(as.numeric(predict$AUC)))^2 + pheno.tem1GWAS
  Genetic_Variance_GWAS[count] <- pheno.bestGWAS[count]/herit0*100
  
  count <- count + 1
}

AUC_Optimum <- pnorm(sqrt(pheno.best/2))
AUC_GWAS <- pnorm(sqrt(pheno.bestGWAS/2))

results <- data.frame(n_seq = n_seq,pheno.best = pheno.best, pheno.bestGWAS = pheno.bestGWAS,
                      AUC_Optimum = AUC_Optimum, Genetic_Variance_Optimum = Genetic_Variance_Optimum,
                      AUC_GWAS = AUC_GWAS,Genetic_Variance_GWAS = Genetic_Variance_GWAS)

save(results,current_sample_optimum,double_sample_optimum,quadruple_sample_optimum,current_sample_gwas,double_sample_gwas,quadruple_sample_gwas,file = save_file)