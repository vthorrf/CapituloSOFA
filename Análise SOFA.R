###===--- Modelo Semi-Paramétrico da Abordagem SOFA ---===###
## Baseado em Franco et al. (2021)

# Clear old objects and analysis====
rm(list=ls())
dev.off()
cat("\014")

### Load all packages and functions====
require(jagsUI)
DCompare <- function(dics, rounding=4) {
  ddics <- dics - min(dics)  # Difference of the DICs with the best fitting
                             # -min(Dics)-model.
  
  LRs <- exp( (-.5) * ddics) # Return the Likelihood Ratio.
  LRs <- round(LRs,rounding) # Rounding for aesthetics
  
  wM <- LRs/sum(LRs)         # Return the Bayesian model weights. This is the
                             # posterior probability of the models given the
                             # data, assuming exaustiveness of the models.
  wM <- round(wM, rounding)  # Rounding for aesthetics
  
  Result <- matrix(cbind(dics,ddics,LRs,wM),ncol=4)
  colnames(Result) <- c("DIC","dDIC","LR","w")
  row.names(Result) <- sapply(1:length(dics), function(g) 
    if(g < 10) { paste("Model_",0,g,sep="") } 
    else { paste("Model_",g,sep="") })
  return(Result)
}

### Read data====
df <- read.csv("Análise SOFA.csv")
score <- c(rowMeans(df[,4:11]))
test <- data.frame(cond=df$cond,score)
y <- test$score[-2]
x <- test$cond[-2]
CCS <- df[-2,12:22]; scale <- CCS
cols <- which(sign(psych::irt.fa(CCS, 1)$fa$loadings) == 1)
scale[ ,cols] = 6 - CCS[ ,cols]
scale <- scale[,-11]

### Combinar dados para rodar a análise====
n     = length(x)
nk    = tryCatch(nk, error=function(e) (round(sqrt(n-2),0) + 1) )
knots = as.vector(quantile(x, seq(0,.95,len=nk)))
X <- splines::bs(x, knots=knots, degree=3)
scale = if(min(scale) == 0) {scale} else{scale - 1}
dataList = list( y=y, X=X, n=n, nk=nk, scale=scale, Rmat=diag(2), Nitem=ncol(scale) )
# Define some MCMC parameters for JAGS
nthin    = 1                           # How Much Thinning?
nchains  = parallel::detectCores() - 1 # How Many Chains?
nburnin  = 7500                        # How Many Burn-in Samples?
nsamples = 27500                       # How Many Recorded Samples?
nadapt   = 45000                       # How Many Adaptive Samples?

### Modelo de Traço====
params <- c("b0", "Beta", "err", "Um", "Mu", "S",
            "Disc", "Diff", "muT", "tau", "pred")
set.seed(1234)
fitT <- jagsUI::jags(dataList, inits=NULL, params,
                     model.file="./SOFATraco.txt",
                     n.chains=nchains, n.adapt=nadapt,
                     n.iter=nsamples,  n.burnin=nburnin,
                     n.thin=nthin, modules=c('glm'),
                     factories=NULL, parallel=TRUE,
                     n.cores=nchains, DIC=TRUE,
                     store.data=FALSE, codaOnly=FALSE,
                     bugs.format=FALSE, verbose=TRUE)

### Modelo de Traço Independente====
set.seed(1234)
fitTI <- jagsUI::jags(dataList, inits=NULL, params,
                      model.file="./SOFATracoInd.txt",
                      n.chains=nchains, n.adapt=nadapt,
                      n.iter=nsamples,  n.burnin=nburnin,
                      n.thin=nthin, modules=c('glm'),
                      factories=NULL, parallel=TRUE,
                      n.cores=nchains, DIC=TRUE,
                      store.data=FALSE, codaOnly=FALSE,
                      bugs.format=FALSE, verbose=TRUE)

### Modelo de Estado====
paramE <- c("b0", "Beta", "err", "Delta1", "Delta2", "Um",
            "Mu", "S", "Disc", "Diff", "muT", "tau", "pred")
set.seed(1234)
fitE <- jagsUI::jags(dataList, inits=NULL, paramE,
                     model.file="./SOFAEstado.txt",
                     n.chains=nchains, n.adapt=nadapt,
                     n.iter=nsamples,  n.burnin=nburnin,
                     n.thin=nthin, modules=c('glm'),
                     factories=NULL, parallel=TRUE,
                     n.cores=nchains, DIC=TRUE,
                     store.data=FALSE, codaOnly=FALSE,
                     bugs.format=FALSE, verbose=TRUE)

### Modelo de Estado Independente SOFA====
set.seed(1234)
fitEI <- jagsUI::jags(dataList, inits=NULL, paramE,
                      model.file="./SOFAEstadoInd.txt",
                      n.chains=nchains, n.adapt=nadapt,
                      n.iter=nsamples,  n.burnin=nburnin,
                      n.thin=nthin, modules=c('glm'),
                      factories=NULL, parallel=TRUE,
                      n.cores=nchains, DIC=TRUE,
                      store.data=FALSE, codaOnly=FALSE,
                      bugs.format=FALSE, verbose=TRUE)

### Modelo de Estado Independente Coop-Comp====
set.seed(1234)
fitCC <- jagsUI::jags(dataList, inits=NULL, paramE,
                      model.file="./SOFAEstadoCCp.txt",
                      n.chains=nchains, n.adapt=nadapt,
                      n.iter=nsamples,  n.burnin=nburnin,
                      n.thin=nthin, modules=c('glm'),
                      factories=NULL, parallel=TRUE,
                      n.cores=nchains, DIC=TRUE,
                      store.data=FALSE, codaOnly=FALSE,
                      bugs.format=FALSE, verbose=TRUE)

### Compare====
# HDI 95% Correlação
COR1 <- sapply(1:nrow(fitT$sims.list$Um), function(g) cor(fitT$sims.list$Um[g,,])[1,2])
COR2 <- sapply(1:nrow(fitTI$sims.list$Um), function(g) cor(fitTI$sims.list$Um[g,,])[1,2])
COR3 <- sapply(1:nrow(fitE$sims.list$Um), function(g) cor(fitE$sims.list$Um[g,,])[1,2])
COR4 <- sapply(1:nrow(fitEI$sims.list$Um), function(g) cor(fitEI$sims.list$Um[g,,])[1,2])
COR5 <- sapply(1:nrow(fitCC$sims.list$Um), function(g) cor(fitCC$sims.list$Um[g,,])[1,2])
CORS <- cbind(COR1, COR2, COR3, COR4, COR5)
round(apply(CORS, 2, quantile, c(.025, .975)), 2)

# DICs
round(DCompare(c(fitT$DIC, fitTI$DIC, fitE$DIC, fitEI$DIC, fitCC$DIC)), 3)

####====---- The End ----====####