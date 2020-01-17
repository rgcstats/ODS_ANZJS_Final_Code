##################################################
# ODS Simulation Study Continuous Response
# Final version for ANZJS Paper
##################################################

#################################
# Load packages
#################################

devtools::install_github("rgcstats/ODS") # note: this is different from the CRAN ODS package (different authors)
library(ODS)
library(sampling)

#################################
# scenarios
#################################

R <- 500
scenarios <- expand.grid( design=c("EPS","strat1","strat2","strat3",
                                   "ushaped0.10","ushaped0.30","ushaped0.60","ushaped1.00"),
                          N=5000, n=c(300,1000) , beta1=c(0,0.7,0.85) , beta0=0 , beta2=-0.5,
                          stringsAsFactors=FALSE)
beta0 <- 1 ; beta2 <- -0.5 ; # sd <- 1
estimator.names <- c("OLS","WLS","MSL","ML2x")
quantity.names <- c("b0","b1","b2","sd","prob.x2","mu0.x1","mu1.x1","sd.x1")
t1 <- proc.time()[3]
all.simresults <- NULL
for(k in c(1:nrow(scenarios))){
  cat("************************ \n")
  cat("STARTING SCENARIO ",k,"\n")
  cat(paste(colnames(scenarios),as.character(scenarios[k,]),sep="=",collapse=" ; "),"\n")
  cat("************************ \n")
  N <- scenarios[k,"N"] ; n <- scenarios[k,"n"] ; design <- scenarios[k,"design"]
  beta1 <- scenarios[k,"beta1"] ; beta2 <- scenarios[k,"beta2"] ; beta0 <- scenarios[k,"beta0"]
  prob.x2 <- 0.5 ; mu0.x1 <- 0 ; mu1.x1 <- 0 ; sd.x1 <- 1
  sd <- sqrt( 1 - beta1^2 - beta2^2*prob.x2*(1-prob.x2) )
  # so E[Y]=0 and var[Y]=1 
  sim.estimators <- array(data=NA,dim=c(length(quantity.names),R,length(estimator.names)),
                          dimnames=list(quantity=quantity.names,r=as.character(c(1:R)),est=estimator.names))
  set.seed(50155)
  for(r in c(1:R)){
    if(r%%25==0) cat("starting simulation ",r,"\n")
    # generate population
    popn <- data.frame(z=1*(runif(N)<=prob.x2))
    popn$x <- mu0.x1*(1-popn$z) + mu1.x1*popn$z + sd.x1*rnorm(N)
    popn$y <- beta0 + beta1*popn$x + beta2*(popn$z-prob.x2) + sd*rnorm(N)
    # select sample according to specified design, and create popn$pi, pi.fn, p.s and specs
    if(design=="EPS"){
      s <- sample(N,n,replace=FALSE)
      popn$pi <- n/N
      meta.pi.fn <- function(u){function(y){rep(u,length(y))}}
      pi.fn <- meta.pi.fn(n/N)
    }
    if(design=="strat1"){
      cutoffs <- c(-1,1) 
      popn$h <- cut(x=popn$y,breaks=c(-Inf,cutoffs,Inf),include.lowest=TRUE)
      Nh <- table(popn$h)
      nh <- round(2/3*n*Nh/N + c(n/6,0,n/6))
    }
    if(design=="strat2"){
      cutoffs <- 0
      popn$h <- cut(x=popn$y,breaks=c(-Inf,cutoffs,Inf),include.lowest=TRUE)
      Nh <- table(popn$h)
      nh <- round(c(0.2*n,0.8*n))
    }
    if(design=="strat3"){
      cutoffs <- qnorm(p=c(0.25,0.75))
      popn$h <- cut(x=popn$y,breaks=c(-Inf,cutoffs,Inf),include.lowest=TRUE)
      Nh <- table(popn$h)
      nh <- rep(round(n/3),3)
    }
    if(design %in% c("strat1","strat2","strat3")){
      popn$pi <- (nh/Nh)[popn$h]
      meta.pi.fn <- function(fh,cutoffs){
        function(y){
          h <- cut(x=y,breaks=c(-Inf,cutoffs,Inf),include.lowest=TRUE)
          return(as.vector(fh[h]))
        }
      }
      pi.fn <- meta.pi.fn(nh/Nh,cutoffs)
    }
    u.tuner <- NA
    if(design=="ushaped0.10") u.tuner <- 0.10
    if(design=="ushaped0.25") u.tuner <- 0.25
    if(design=="ushaped0.30") u.tuner <- 0.3
    if(design=="ushaped0.50") u.tuner <- 0.5
    if(design=="ushaped0.60") u.tuner <- 0.6
    if(design=="ushaped0.75") u.tuner <- 0.75
    if(design=="ushaped1.00") u.tuner <- 1
    if(design=="ushaped0.30") u.tuner <- 0.3
    if(design=="ushaped0.60") u.tuner <- 0.6
    if(design %in% paste0("ushaped",c("0.10","0.25","0.30","0.50","0.60","0.75","1.00"))){
      integrand <- function(z,tuner) dnorm(z)*sqrt(tuner+(1-tuner)*z^2)
      A <- integrate(integrand,-Inf,Inf,tuner=u.tuner)$value
      popn$pi <- sqrt(u.tuner+(1-u.tuner)*popn$y^2)/A*n
      meta.pi.fn <- function(A,n,tuner){
        pi.fn <- function(y){
          sqrt(tuner+(1-tuner)*y^2)/A*n
        }
        return(pi.fn)
      }
      pi.fn <- meta.pi.fn(A=A,n=n,tuner=u.tuner)
    }
    s <- c(1:N)[runif(N)<=popn$pi]
    samp <- popn[s,]
    # Calculate linear model coefficient estimators
    ols.fit <- lm(y~x+z,data=samp)
    sim.estimators[c("b0","b1","b2","sd"),r,"OLS"] <- c(ols.fit$coef,summary(ols.fit)$sigma)
    if(design!="EPS"){
      wls.fit <- lm(y~x+z,data=samp,weights=1/samp$pi)
      sigma.wls <- sum(wls.fit$resid^2/samp$pi)/sum(1/samp$pi)
      sim.estimators[c("b0","b1","b2","sd"),r,"WLS"] <- c(wls.fit$coef,sigma.wls)
      if(design %in% c("strat1","strat2","strat3")){
        sim.estimators[c("b0","b1","b2","sd"),r,"MSL"] <- lm.samplik.2x.poisson.strat(ys=samp$y,xs=cbind(1,samp$x,samp$z),pi.fn=pi.fn,pi.h=nh/Nh,cutoffs=cutoffs)
        sim.estimators[,r,"ML2x"] <- lm.MLE.2x.poisson.strat(ys=samp$y,x1s=samp$x,x2s=samp$z,pi.h=nh/Nh,cutoffs=cutoffs,pi.s=samp$pi,
                                                             start.beta=sim.estimators[c("b0","b1","b2"),r,"MSL"],
                                                             start.sd=sim.estimators["sd",r,"MSL"])
      } else{
        sim.estimators[c("b0","b1","b2","sd"),r,"MSL"] <- lm.samplik(ys=samp$y,xs=cbind(1,samp$x,samp$z),pi.fn=pi.fn,
                                                                     method="GH")
        sim.estimators[,r,"ML2x"] <- lm.MLE.2x.poisson(ys=samp$y,x1s=samp$x,x2s=samp$z,pi.fn=pi.fn,pi.s=samp$pi,Rquad=150,
                                                       start.beta=sim.estimators[c("b0","b1","b2"),r,"MSL"],
                                                       start.sd=sim.estimators["sd",r,"MSL"])
      }
    }
  }
  E <- as.data.frame(apply(sim.estimators,c(3,1),mean))
  E$est <- rownames(E)
  E$summary <- "E"
  V <- as.data.frame(apply(sim.estimators,c(3,1),var,na.rm=T))
  V$est <- rownames(V)
  V$summary <- "V"
  truevalues <- array(data=c(beta0,beta1,beta2,sd,prob.x2,mu0.x1,mu1.x1,sd.x1),dim=dim(sim.estimators),dimnames=dimnames(sim.estimators))
  sim.errors <- sim.estimators - truevalues
  MSE <- as.data.frame(apply(sim.errors^2,c(3,1),mean))
  MSE$est <- rownames(MSE)
  MSE$summary <- "MSE"
  true <- data.frame(b0=beta0,b1=beta1,b2=beta2,sd=sd,prob.x2=prob.x2,mu0.x1=mu0.x1,mu1.x1=mu1.x1,sd.x1=sd.x1)
  true$est <- rownames(true)
  true$summary <- "true"
  summary.k <- rbind( E , V[,colnames(E)] , MSE[,colnames(E)] , true[,colnames(E)] )
  summary.k$scenario <- k
  rownames(summary.k) <- NULL
  all.simresults <- rbind(all.simresults,summary.k)
}
t2 <- proc.time()[3]
t2-t1

all.simresults$design <- scenarios$design[all.simresults$scenario]
all.simresults$N <- scenarios$N[all.simresults$scenario]
all.simresults$n <- scenarios$n[all.simresults$scenario]
all.simresults$beta1 <- scenarios$beta1[all.simresults$scenario]
all.simresults$longdesign <- paste0(all.simresults$design," with n=",all.simresults$n,
" and N=",all.simresults$N)

save.image("image_cts.Rdata")

################################################################
# Tables of expected values of \hat{\beta}_1 for \beta=0, 1
################################################################

make.table.b1 <- function(simresults,apply.function=I,include.efficiency=FALSE){
  tab <- apply.function(with(simresults,tapply(b1,list(longdesign,est),mean)))
  tab <- tab[,c("OLS","WLS","MSL","ML2x")]
  if(include.efficiency) tab <- cbind(tab,(tab[,"ML2x"]/tab[,"MSL"])^2)
  tab <- format(round(tab,digits=4),nsmall=4)
  tab <- cbind(rownames(tab),tab)
  rownames(tab) <- NULL
  colnames(tab)[1] <- "design"
  tab[grep("NA",tab)] <- " "
  tab
}

make.table.b2 <- function(simresults,apply.function=I,include.efficiency=FALSE){
  tab <- apply.function(with(simresults,tapply(b2,list(longdesign,est),mean)))
  tab <- tab[,c("OLS","WLS","MSL","ML2x")]
  if(include.efficiency) tab <- cbind(tab,(tab[,"ML2x"]/tab[,"MSL"])^2)
  tab <- format(round(tab,digits=4),nsmall=4)
  tab <- cbind(rownames(tab),tab)
  rownames(tab) <- NULL
  colnames(tab)[1] <- "design"
  tab[grep("NA",tab)] <- " "
  tab
}

tab1 <- make.table.b1(all.simresults[with(all.simresults,(beta1==0)&(summary=="E")),])
tab1
tab2 <- make.table.b1(all.simresults[with(all.simresults,(beta1==0.7)&(summary=="E")),])
tab2
tab3 <- make.table.b1(all.simresults[with(all.simresults,(beta1==0.85)&(summary=="E")),])
tab3

################################################################
# Tables of standard deviations of \hat{\beta}_1 for \beta=0, 1
################################################################

tab4 <- make.table.b1(all.simresults[with(all.simresults,(beta1==0)&(summary=="V")),],sqrt)
tab4
tab5 <- make.table.b1(all.simresults[with(all.simresults,(beta1==0.7)&(summary=="V")),],sqrt)
tab5
tab6 <- make.table.b1(all.simresults[with(all.simresults,(beta1==0.85)&(summary=="V")),],sqrt)
tab6

################################################################
# Tables of RMSEs of \hat{\beta}_1 for \beta=0, 1
################################################################

tab7 <- make.table.b1(all.simresults[with(all.simresults,(beta1==0)&(summary=="MSE")),],sqrt,include.efficiency=TRUE)
tab7
tab8 <- make.table.b1(all.simresults[with(all.simresults,(beta1==0.7)&(summary=="MSE")),],sqrt,include.efficiency=TRUE)
tab8
tab9 <- make.table.b1(all.simresults[with(all.simresults,(beta1==0.85)&(summary=="MSE")),],sqrt,include.efficiency=TRUE)
tab9

################################################################
# Table of RMSEs of \hat{\beta}_2 when \beta_0=\beta_1=0
################################################################

tab10 <- make.table.b2(all.simresults[with(all.simresults,(beta1==0)&(summary=="MSE")),],sqrt,include.efficiency=TRUE)
tab10
tab11 <- make.table.b2(all.simresults[with(all.simresults,(beta1==0.7)&(summary=="MSE")),],sqrt,include.efficiency=TRUE)
tab11
tab12 <- make.table.b2(all.simresults[with(all.simresults,(beta1==0.85)&(summary=="MSE")),],sqrt,include.efficiency=TRUE)
tab12


################################################################
# Table 1 for Paper: RMSEs of \hat{\beta}_1 when beta1=0
################################################################

simresults.sub <- all.simresults[with(all.simresults,(beta1==0)&(summary=="MSE")),]
paper.table1 <- sqrt(with(simresults.sub,tapply(b1,list(longdesign,est),mean)))*100
simresults.sub2 <- all.simresults[with(all.simresults,(beta1==0)&(summary=="E")),]
paper.table1.biases <- with(simresults.sub2,tapply(b1,list(longdesign,est),mean))*100
paper.table1 <- as.data.frame(paper.table1)
paper.table1.biases <- as.data.frame(paper.table1.biases)

paper.table1[,c("OLS","WLS","MSL","ML2x")] <- format(round(paper.table1[,c("OLS","WLS","MSL","ML2x")],digits=2),nsmall=2)
paper.table1.biases[,c("OLS","WLS","MSL","ML2x")] <- format(round(paper.table1.biases[,c("OLS","WLS","MSL","ML2x")],digits=2),nsmall=2)

paper.table1$WLS[grep("NA",paper.table1$WLS)] <- paper.table1$OLS[grep("NA",paper.table1$WLS)]
paper.table1$MSL[grep("NA",paper.table1$MSL)] <- paper.table1$OLS[grep("NA",paper.table1$MSL)]
paper.table1$ML2x[grep("NA",paper.table1$ML2x)] <- paper.table1$OLS[grep("NA",paper.table1$ML2x)]

paper.table1.biases$WLS[grep("NA",paper.table1.biases$WLS)] <- paper.table1.biases$OLS[grep("NA",paper.table1.biases$WLS)]
paper.table1.biases$MSL[grep("NA",paper.table1.biases$MSL)] <- paper.table1.biases$OLS[grep("NA",paper.table1.biases$MSL)]
paper.table1.biases$ML2x[grep("NA",paper.table1.biases$ML2x)] <- paper.table1.biases$OLS[grep("NA",paper.table1.biases$ML2x)]

paper.table1$OLS <- paste0(paper.table1$OLS," (",paper.table1.biases$OLS,")")
paper.table1$WLS <- paste0(paper.table1$WLS," (",paper.table1.biases$WLS,")")
paper.table1$MSL <- paste0(paper.table1$MSL," (",paper.table1.biases$MSL,")")
paper.table1$ML2x <- paste0(paper.table1$ML2x," (",paper.table1.biases$ML2x,")")

paper.table1 <- paper.table1[c(2,4,6,8,10,12,14,16,1,3,5,7,9,11,13,15),]
paper.table1$design <- rep( c("EPS","strat1","exclude","strat2",rep("ushaped",4)) ,
                            times=2 )
paper.table1$tuning <- rep( c("","","","","0.1","0.3","0.6","1.0") , times=2)
paper.table1$n <- rep(c(300,1000),each=8)
paper.table1 <- paper.table1[,c("design","tuning","n","OLS","WLS","MSL","ML2x")]
paper.table1
write.table(paper.table1[c(1,2,4:7),-3],file="cts_sim_tab1.txt",sep=" & ",eol=" \\\\ \n" , row.names=FALSE , col.names=FALSE,quote=F)
write.table(paper.table1[8+c(1,2,4:7),-3],file="cts_supp_sim_tab1.txt",sep=" & ",eol=" \\\\ \n" , row.names=FALSE , col.names=FALSE,quote=F)

################################################################
# Table 2 for Paper: RMSEs of \hat{\beta}_1 when beta1=0.7
################################################################

simresults.sub <- all.simresults[with(all.simresults,(beta1==0.7)&(summary=="MSE")),]
paper.table2 <- sqrt(with(simresults.sub,tapply(b1,list(longdesign,est),mean)))*100
simresults.sub2 <- all.simresults[with(all.simresults,(beta1==0.7)&(summary=="E")),]
paper.table2.biases <- (with(simresults.sub2,tapply(b1,list(longdesign,est),mean))-0.7)*100
paper.table2 <- as.data.frame(paper.table2)
paper.table2.biases <- as.data.frame(paper.table2.biases)

paper.table2[,c("OLS","WLS","MSL","ML2x")] <- format(round(paper.table2[,c("OLS","WLS","MSL","ML2x")],digits=2),nsmall=2)
paper.table2.biases[,c("OLS","WLS","MSL","ML2x")] <- format(round(paper.table2.biases[,c("OLS","WLS","MSL","ML2x")],digits=2),nsmall=2)

paper.table2$WLS[grep("NA",paper.table2$WLS)] <- paper.table2$OLS[grep("NA",paper.table2$WLS)]
paper.table2$MSL[grep("NA",paper.table2$MSL)] <- paper.table2$OLS[grep("NA",paper.table2$MSL)]
paper.table2$ML2x[grep("NA",paper.table2$ML2x)] <- paper.table2$OLS[grep("NA",paper.table2$ML2x)]

paper.table2.biases$WLS[grep("NA",paper.table2.biases$WLS)] <- paper.table2.biases$OLS[grep("NA",paper.table2.biases$WLS)]
paper.table2.biases$MSL[grep("NA",paper.table2.biases$MSL)] <- paper.table2.biases$OLS[grep("NA",paper.table2.biases$MSL)]
paper.table2.biases$ML2x[grep("NA",paper.table2.biases$ML2x)] <- paper.table2.biases$OLS[grep("NA",paper.table2.biases$ML2x)]

paper.table2$OLS <- paste0(paper.table2$OLS," (",paper.table2.biases$OLS,")")
paper.table2$WLS <- paste0(paper.table2$WLS," (",paper.table2.biases$WLS,")")
paper.table2$MSL <- paste0(paper.table2$MSL," (",paper.table2.biases$MSL,")")
paper.table2$ML2x <- paste0(paper.table2$ML2x," (",paper.table2.biases$ML2x,")")

paper.table2 <- paper.table2[c(2,4,6,8,10,12,14,16,1,3,5,7,9,11,13,15),]
paper.table2$design <- rep( c("EPS","strat1","exclude","strat2",rep("ushaped",4)) ,
                            times=2 )
paper.table2$tuning <- rep( c("","","","","0.1","0.3","0.6","1.0") , times=2)
paper.table2$n <- rep(c(300,1000),each=8)
paper.table2 <- paper.table2[,c("design","tuning","n","OLS","WLS","MSL","ML2x")]
paper.table2
write.table(paper.table2[c(1,2,4:7),-3],file="cts_sim_tab2.txt",sep=" & ",eol=" \\\\ \n" , row.names=FALSE , col.names=FALSE,quote=F)
write.table(paper.table2[8+c(1,2,4:7),-3],file="cts_supp_sim_tab2.txt",sep=" & ",eol=" \\\\ \n" , row.names=FALSE , col.names=FALSE,quote=F)

################################################################
# Table 3 for Paper: RMSEs (biases) of \hat{\beta}_2 when beta1=0
################################################################

simresults.sub <- all.simresults[with(all.simresults,(beta1==0)&(summary=="MSE")),]
paper.table3 <- sqrt(with(simresults.sub,tapply(b2,list(longdesign,est),mean)))*100
simresults.sub2 <- all.simresults[with(all.simresults,(beta1==0)&(summary=="E")),]
paper.table3.biases <- (with(simresults.sub2,tapply(b2,list(longdesign,est),mean))+0.5)*100
paper.table3 <- as.data.frame(paper.table3)
paper.table3.biases <- as.data.frame(paper.table3.biases)

paper.table3[,c("OLS","WLS","MSL","ML2x")] <- format(round(paper.table3[,c("OLS","WLS","MSL","ML2x")],digits=2),nsmall=2)
paper.table3.biases[,c("OLS","WLS","MSL","ML2x")] <- format(round(paper.table3.biases[,c("OLS","WLS","MSL","ML2x")],digits=2),nsmall=2)

paper.table3$WLS[grep("NA",paper.table3$WLS)] <- paper.table3$OLS[grep("NA",paper.table3$WLS)]
paper.table3$MSL[grep("NA",paper.table3$MSL)] <- paper.table3$OLS[grep("NA",paper.table3$MSL)]
paper.table3$ML2x[grep("NA",paper.table3$ML2x)] <- paper.table3$OLS[grep("NA",paper.table3$ML2x)]

paper.table3.biases$WLS[grep("NA",paper.table3.biases$WLS)] <- paper.table3.biases$OLS[grep("NA",paper.table3.biases$WLS)]
paper.table3.biases$MSL[grep("NA",paper.table3.biases$MSL)] <- paper.table3.biases$OLS[grep("NA",paper.table3.biases$MSL)]
paper.table3.biases$ML2x[grep("NA",paper.table3.biases$ML2x)] <- paper.table3.biases$OLS[grep("NA",paper.table3.biases$ML2x)]

paper.table3$OLS <- paste0(paper.table3$OLS," (",paper.table3.biases$OLS,")")
paper.table3$WLS <- paste0(paper.table3$WLS," (",paper.table3.biases$WLS,")")
paper.table3$MSL <- paste0(paper.table3$MSL," (",paper.table3.biases$MSL,")")
paper.table3$ML2x <- paste0(paper.table3$ML2x," (",paper.table3.biases$ML2x,")")

paper.table3 <- paper.table3[c(2,4,6,8,10,12,14,16,1,3,5,7,9,11,13,15),]
paper.table3$design <- rep( c("EPS","strat1","exclude","strat2",rep("ushaped",4)) ,
                            times=2 )
paper.table3$tuning <- rep( c("","","","","0.1","0.3","0.6","1.0") , times=2)
paper.table3$n <- rep(c(300,1000),each=8)
paper.table3 <- paper.table3[,c("design","tuning","n","OLS","WLS","MSL","ML2x")]
paper.table3
write.table(paper.table3[c(1,2,4:7),-3],file="cts_sim_tab3.txt",sep=" & ",eol=" \\\\ \n" , row.names=FALSE , col.names=FALSE,quote=F)
write.table(paper.table3[8+c(1,2,4:7),-3],file="cts_supp_sim_tab3.txt",sep=" & ",eol=" \\\\ \n" , row.names=FALSE , col.names=FALSE,quote=F)

