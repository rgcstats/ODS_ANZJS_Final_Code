##################################################
# ODS Simulation Study Continuous Response
# Final version for ANZJS Paper
##################################################

#################################
# Load packages
#################################

devtools::install_github("rgcstats/ODS") # note: this is different from the CRAN ODS package (different authors)
library(ODS)

#########################################################################
# First of all, calculate marginal P[Y=1] for different scenarios
#########################################################################

marginal.prob <- function(beta){
  prob.x2 <- 0.5 ; mu0.x1 <- 0 ; mu1.x1 <- 0 ; sd.x1 <- 1
  # Firstly get P[Y=1|x2=0], by integrating f(beta1+beta2*x1+beta2*0) over
  #     x1~N(mu0.x1,sd.x1) conditional on x2
  integrand0 <- function(x1){
    make.link("logit")$linkinv(beta[1]+beta[2]*x1) * dnorm(x1,mu0.x1,sd.x1)
  }
  mu0.y <- integrate(integrand0,-Inf,Inf)$value
  # Similarly get P[Y=1|x2=1], by integrating f(beta1+beta2*x1+beta2*1) over
  #     x1~N(mu1.x1,sd.x1) conditional on x2
  integrand1 <- function(x1){
    make.link("logit")$linkinv(beta[1]+beta[2]*x1+beta[3]) * dnorm(x1,mu1.x1,sd.x1)
  }
  mu1.y <- integrate(integrand1,-Inf,Inf)$value
  # Now we can get mu=E[Y]
  mu <- (1-prob.x2)*mu0.y + prob.x2*mu1.y
  return(mu)
}

ratio.case.control.dist <- function(beta0,beta1,beta2,target.ratio){
  mu <- marginal.prob(c(beta0,beta1,beta2))
  return(((mu/(1-mu))-target.ratio)^2)
}

#################################
# scenarios
#################################

R <- 500
scenarios <- expand.grid( design="CC",
                          n=c(300,1000),fcase=c(1,0.5,0.25),fcontrol=c(0.25,0.1,0.01),
                          beta1=c(0,0.7,0.85) , beta2=-0.5,
                          stringsAsFactors=FALSE)
estimator.names <- c("OLS","WLS","ML2x")
quantity.names <- c("b0","b1","b2","prob.x2","mu0.x1","mu1.x1","sd.x1")
t1 <- proc.time()[3]
all.simresults <- NULL
for(k in c(1:nrow(scenarios))){
  cat("************************ \n")
  cat("STARTING SCENARIO ",k,"\n")
  cat(paste(colnames(scenarios),as.character(scenarios[k,]),sep="=",collapse=" ; "),"\n")
  cat("************************ \n")
  n <- scenarios[k,"n"] ; design <- scenarios[k,"design"]
  beta1 <- scenarios[k,"beta1"] ; beta2 <- scenarios[k,"beta2"]
  fcase <- scenarios[k,"fcase"] ; fcontrol <- scenarios[k,"fcontrol"]
  prob.x2 <- 0.5 ; mu0.x1 <- 0 ; mu1.x1 <- 0 ; sd.x1 <- 1
  sim.estimators <- array(data=NA,dim=c(length(quantity.names),R,length(estimator.names)),
                          dimnames=list(quantity=quantity.names,r=as.character(c(1:R)),est=estimator.names))
  # calculate beta0 and N
  if(fcontrol==0) fcontrol <- fcase
  ncase <- ncontrol <- round(n/2)
  Ncase.expected <- round(ncase/fcase)
  Ncontrol.expected <- round(ncontrol/fcontrol)
  N <- Ncase.expected + Ncontrol.expected
  beta0 <- optimize(f=ratio.case.control.dist,interval=c(-50,50),
                    beta1=beta1,beta2=beta2,
                    target.ratio=Ncase.expected/Ncontrol.expected)$minimum
  pi.case <- fcase ; pi.control <- fcontrol 
  scenarios$N[k] <- N ; scenarios$beta0[k] <- beta0
  set.seed(50155)
  for(r in c(1:R)){
    if(r%%25==0) cat("starting simulation ",r,"\n")
    # generate population
    popn <- data.frame(z=1*(runif(N)<=prob.x2))
    popn$x <- mu0.x1*(1-popn$z) + mu1.x1*popn$z + sd.x1*rnorm(N)
    popn$mu.y <- make.link("logit")$linkinv(beta0 + beta1*popn$x + beta2*(popn$z-prob.x2))
    popn$y <- 1*(runif(N)<=popn$mu.y)
    # select sample according to specified design, and create popn$pi, pi.fn, p.s and specs
    popn$pi <- pmin( 1 , popn$y*pi.case + (1-popn$y)*pi.control )
    s <- c(1:N)[runif(N)<=popn$pi]
    samp <- popn[s,]
    # Calculate linear model coefficient estimators
    ols.fit <- glm(y~x+z,data=samp,family=binomial)
    sim.estimators[c("b0","b1","b2"),r,"OLS"] <- ols.fit$coef
    if(pi.case!=pi.control){
      wls.fit <- glm(y~x+z,data=samp,weights=1/samp$pi,family=binomial)
      sim.estimators[c("b0","b1","b2"),r,"WLS"] <- wls.fit$coef
      sim.estimators[,r,"ML2x"] <- lm.MLE.2x.poisson.logistic(ys=samp$y,x1s=samp$x,x2s=samp$z,pi1=pi.case,
                                                              pi0=pi.control,pi.s=samp$pi)
    }
  }
  E <- as.data.frame(apply(sim.estimators,c(3,1),mean))
  E$est <- rownames(E)
  E$summary <- "E"
  V <- as.data.frame(apply(sim.estimators,c(3,1),var,na.rm=T))
  V$est <- rownames(V)
  V$summary <- "V"
  truevalues <- array(data=c(beta0,beta1,beta2,prob.x2,mu0.x1,mu1.x1,sd.x1),dim=dim(sim.estimators),dimnames=dimnames(sim.estimators))
  sim.errors <- sim.estimators - truevalues
  MSE <- as.data.frame(apply(sim.errors^2,c(3,1),mean))
  MSE$est <- rownames(MSE)
  MSE$summary <- "MSE"
  true <- data.frame(b0=beta0,b1=beta1,b2=beta2,prob.x2=prob.x2,mu0.x1=mu0.x1,mu1.x1=mu1.x1,sd.x1=sd.x1)
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
all.simresults$n <- scenarios$n[all.simresults$scenario]
all.simresults$N <- scenarios$N[all.simresults$scenario]
all.simresults$fcase <- scenarios$fcase[all.simresults$scenario]
all.simresults$fcontrol <- scenarios$fcontrol[all.simresults$scenario]
all.simresults$beta1 <- scenarios$beta1[all.simresults$scenario]
all.simresults$beta0 <- scenarios$beta0[all.simresults$scenario]
all.simresults$longdesign <- paste0("CC with n=",all.simresults$n,
                                    ", fcase=",all.simresults$fcase,
                                         ", fcontrol=",all.simresults$fcontrol,
                                         ", beta1=",all.simresults$beta1)

save.image("image_binary.Rdata")

################################################################
# Tables of expected values of \hat{\beta}_1 for \beta=0, 1
################################################################

make.table.b1 <- function(simresults,apply.function=I){
  tab <- apply.function(with(simresults,tapply(b1,list(longdesign,est),mean)))
  tab <- tab[,c("OLS","WLS","ML2x")]
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

tab7 <- make.table.b1(all.simresults[with(all.simresults,(beta1==0)&(summary=="MSE")),],sqrt)
tab7
tab8 <- make.table.b1(all.simresults[with(all.simresults,(beta1==0.7)&(summary=="MSE")),],sqrt)
tab8
tab9 <- make.table.b1(all.simresults[with(all.simresults,(beta1==0.85)&(summary=="MSE")),],sqrt)
tab9

################################################################
# Table 1 for Paper: RMSEs of \hat{\beta}_1 when beta1=0
################################################################

simresults.sub <- all.simresults[with(all.simresults,(beta1==0)&(summary=="MSE")),]
paper.table1 <- sqrt(with(simresults.sub,tapply(b1,list(longdesign,est),mean)))*100
pi1 <- with(simresults.sub,tapply(fcase,longdesign,mean))
pi0 <- with(simresults.sub,tapply(fcontrol,longdesign,mean))
N <- with(simresults.sub,tapply(N,longdesign,mean))
N1 <- round(N/pi1)
N0 <- round(N/pi0)
simresults.sub2 <- all.simresults[with(all.simresults,(beta1==0)&(summary=="E")),]
paper.table1.biases <- with(simresults.sub2,tapply(b1,list(longdesign,est),mean))*100
paper.table1 <- as.data.frame(paper.table1)
paper.table1.biases <- as.data.frame(paper.table1.biases)

paper.table1[,c("OLS","WLS","ML2x")] <- format(round(paper.table1[,c("OLS","WLS","ML2x")],digits=2),nsmall=2)
paper.table1.biases[,c("OLS","WLS","ML2x")] <- format(round(paper.table1.biases[,c("OLS","WLS","ML2x")],digits=2),nsmall=2)

paper.table1$WLS[grep("NA",paper.table1$WLS)] <- paper.table1$OLS[grep("NA",paper.table1$WLS)]
paper.table1$ML2x[grep("NA",paper.table1$ML2x)] <- paper.table1$OLS[grep("NA",paper.table1$ML2x)]
paper.table1.biases$WLS[grep("NA",paper.table1.biases$WLS)] <- paper.table1.biases$OLS[grep("NA",paper.table1.biases$WLS)]
paper.table1.biases$ML2x[grep("NA",paper.table1.biases$ML2x)] <- paper.table1.biases$OLS[grep("NA",paper.table1.biases$ML2x)]

paper.table1$OLS <- paste0(paper.table1$OLS," (",paper.table1.biases$OLS,")")
paper.table1$WLS <- paste0(paper.table1$WLS," (",paper.table1.biases$WLS,")")
paper.table1$ML2x <- paste0(paper.table1$ML2x," (",paper.table1.biases$ML2x,")")

paper.table1 <- cbind(pi0,pi1,N0/1000,N1/1000,paper.table1)

supp.table1 <- paper.table1[1:9,c("pi0","pi1","N0/1000","N1/1000","OLS","WLS","ML2x")]
paper.table1 <- paper.table1[10:18,c("pi0","pi1","N0/1000","N1/1000","OLS","WLS","ML2x")]
paper.table1
write.table(paper.table1[paper.table1$pi0<0.25,-c(3,4)],file="bin_sim_tab1.txt",sep=" & ",eol=" \\\\ \n" , row.names=FALSE , col.names=FALSE,quote=F)
write.table(supp.table1[supp.table1$pi0<0.25,-c(3,4)],file="bin_supp_sim_tab1.txt",sep=" & ",eol=" \\\\ \n" , row.names=FALSE , col.names=FALSE,quote=F)

################################################################
# Table 2 for Paper: RMSEs of \hat{\beta}_1 when beta1=0.7
################################################################

simresults.sub <- all.simresults[with(all.simresults,(beta1==0.7)&(summary=="MSE")),]
paper.table2 <- sqrt(with(simresults.sub,tapply(b1,list(longdesign,est),mean)))*100
pi1 <- with(simresults.sub,tapply(fcase,longdesign,mean))
pi0 <- with(simresults.sub,tapply(fcontrol,longdesign,mean))
N <- with(simresults.sub,tapply(N,longdesign,mean))
N1 <- round(N/pi1)
N0 <- round(N/pi0)
simresults.sub2 <- all.simresults[with(all.simresults,(beta1==0.7)&(summary=="E")),]
paper.table2.biases <- (with(simresults.sub2,tapply(b1,list(longdesign,est),mean))-0.7)*100
paper.table2 <- as.data.frame(paper.table2)
paper.table2.biases <- as.data.frame(paper.table2.biases)

paper.table2[,c("OLS","WLS","ML2x")] <- format(round(paper.table2[,c("OLS","WLS","ML2x")],digits=2),nsmall=2)
paper.table2.biases[,c("OLS","WLS","ML2x")] <- format(round(paper.table2.biases[,c("OLS","WLS","ML2x")],digits=2),nsmall=2)

paper.table2$WLS[grep("NA",paper.table2$WLS)] <- paper.table2$OLS[grep("NA",paper.table2$WLS)]
paper.table2$ML2x[grep("NA",paper.table2$ML2x)] <- paper.table2$OLS[grep("NA",paper.table2$ML2x)]
paper.table2.biases$WLS[grep("NA",paper.table2.biases$WLS)] <- paper.table2.biases$OLS[grep("NA",paper.table2.biases$WLS)]
paper.table2.biases$ML2x[grep("NA",paper.table2.biases$ML2x)] <- paper.table2.biases$OLS[grep("NA",paper.table2.biases$ML2x)]

paper.table2$OLS <- paste0(paper.table2$OLS," (",paper.table2.biases$OLS,")")
paper.table2$WLS <- paste0(paper.table2$WLS," (",paper.table2.biases$WLS,")")
paper.table2$ML2x <- paste0(paper.table2$ML2x," (",paper.table2.biases$ML2x,")")

paper.table2 <- cbind(pi0,pi1,N0/1000,N1/1000,paper.table2)

supp.table2 <- paper.table2[1:9,c("pi0","pi1","N0/1000","N1/1000","OLS","WLS","ML2x")]
paper.table2 <- paper.table2[10:18,c("pi0","pi1","N0/1000","N1/1000","OLS","WLS","ML2x")]
paper.table2
write.table(paper.table2[paper.table2$pi0<0.25,-c(3,4)],file="bin_sim_tab2.txt",sep=" & ",eol=" \\\\ \n" , row.names=FALSE , col.names=FALSE,quote=F)
write.table(supp.table2[supp.table2$pi0<0.25,-c(3,4)],file="bin_supp_sim_tab2.txt",sep=" & ",eol=" \\\\ \n" , row.names=FALSE , col.names=FALSE,quote=F)

################################################################
# Table 3 for Paper: RMSEs of \hat{\beta}_2 when beta1=0
################################################################

simresults.sub <- all.simresults[with(all.simresults,(beta1==0)&(summary=="MSE")),]
paper.table3 <- sqrt(with(simresults.sub,tapply(b2,list(longdesign,est),mean)))*100
pi1 <- with(simresults.sub,tapply(fcase,longdesign,mean))
pi0 <- with(simresults.sub,tapply(fcontrol,longdesign,mean))
N <- with(simresults.sub,tapply(N,longdesign,mean))
N1 <- round(N/pi1)
N0 <- round(N/pi0)
simresults.sub2 <- all.simresults[with(all.simresults,(beta1==0)&(summary=="E")),]
paper.table3.biases <- (with(simresults.sub2,tapply(b2,list(longdesign,est),mean))+0.5)*100
paper.table3 <- as.data.frame(paper.table3)
paper.table3.biases <- as.data.frame(paper.table3.biases)

paper.table3[,c("OLS","WLS","ML2x")] <- format(round(paper.table3[,c("OLS","WLS","ML2x")],digits=2),nsmall=2)
paper.table3.biases[,c("OLS","WLS","ML2x")] <- format(round(paper.table3.biases[,c("OLS","WLS","ML2x")],digits=2),nsmall=2)

paper.table3$WLS[grep("NA",paper.table3$WLS)] <- paper.table3$OLS[grep("NA",paper.table3$WLS)]
paper.table3$ML2x[grep("NA",paper.table3$ML2x)] <- paper.table3$OLS[grep("NA",paper.table3$ML2x)]
paper.table3.biases$WLS[grep("NA",paper.table3.biases$WLS)] <- paper.table3.biases$OLS[grep("NA",paper.table3.biases$WLS)]
paper.table3.biases$ML2x[grep("NA",paper.table3.biases$ML2x)] <- paper.table3.biases$OLS[grep("NA",paper.table3.biases$ML2x)]

paper.table3$OLS <- paste0(paper.table3$OLS," (",paper.table3.biases$OLS,")")
paper.table3$WLS <- paste0(paper.table3$WLS," (",paper.table3.biases$WLS,")")
paper.table3$ML2x <- paste0(paper.table3$ML2x," (",paper.table3.biases$ML2x,")")

paper.table3 <- cbind(pi0,pi1,N0/1000,N1/1000,paper.table3)

supp.table3 <- paper.table3[1:9,c("pi0","pi1","N0/1000","N1/1000","OLS","WLS","ML2x")]
paper.table3 <- paper.table3[10:18,c("pi0","pi1","N0/1000","N1/1000","OLS","WLS","ML2x")]
paper.table3
write.table(paper.table3[paper.table3$pi0<0.25,-c(3,4)],file="bin_sim_tab3.txt",sep=" & ",eol=" \\\\ \n" , row.names=FALSE , col.names=FALSE,quote=F)
write.table(supp.table3[supp.table3$pi0<0.25,-c(3,4)],file="bin_supp_sim_tab3.txt",sep=" & ",eol=" \\\\ \n" , row.names=FALSE , col.names=FALSE,quote=F)



