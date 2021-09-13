n <- as.integer(commandArgs(TRUE)[1])  ## n = 250, 500, 1000, 2000
set <- as.integer(commandArgs(TRUE)[2])  ## set = 1,2,...

J <- function(hl,hu,rho,x){
  f <- function(z){pnorm(-(rho*z+hu+sqrt(x-max(z-hl,0)^2))/sqrt(1-rho^2),0,1)*dnorm(z,0,1)}
  return(pnorm(hl+sqrt(x),0,1)-integrate(f,lower=-Inf,upper=hl+sqrt(x))$value)
}
J1 <- function(hl,hu,x){return(pnorm(hl+sqrt(x),0,1)-pnorm(-hu-sqrt(x),0,1))}
cv <- function(x) {J(hl,hu,rho,x)-0.95}
cv1 <- function(x) {J1(hl,hu,x)-0.95}


if(set==1){prob <- c(0.13,0.05,0.018,0.002,0.04,0.05,0.01,0.7)}
if(set==2){prob <- c(0.15,0.035,0.035,0.002,0.03,0.01,0.04,0.698)}
if(set==3){prob <- c(0.15,0.008,0.04,0.002,0.06,0.1,0.18,0.46)}


p111 <- prob[1]
p101 <- prob[2]
p011 <- prob[3]
p001 <- prob[4]
p110 <- prob[5]
p100 <- prob[6]
p010 <- prob[7]
p000 <- prob[8]
sen.A <- (p111+p101)/(p111+p101+p011+p001)
sen.B <- (p111+p011)/(p111+p101+p011+p001)
theta.true <- sen.A-sen.B
theta0.true <- (p101-p011)/(p111+p101+p011)
theta1.true <- (p101-p011)/(p111+p101+p011+p001+p000)
thetaL.true <- min(theta0.true,theta1.true)
thetaU.true <- max(theta0.true,theta1.true)
spe.A <- (p010+p000)/(p110+p100+p010+p000)
spe.B <- (p100+p000)/(p110+p100+p010+p000)
psi.true <- (p010+p000)/(p110+p100+p010+p000)-(p100+p000)/(p110+p100+p010+p000)
psi0.true <- (p010-p100)/(p110+p100+p010)
psi1.true <- (p010-p100)/(p110+p100+p010+p001+p000)
psiL.true <- min(psi0.true,psi1.true)
psiU.true <- max(psi0.true,psi1.true)
c(sen.A,sen.B,theta.true,thetaL.true,thetaU.true)
c(spe.A,spe.B,psi.true,psiL.true,psiU.true)

nrep <- 5000
theta0 <- theta1 <- var.theta0 <- var.theta1 <- rep(NA,nrep)
thetaL <- thetaU <- var.thetaL <- var.thetaU <- rep(NA,nrep)
psi0 <- psi1 <- var.psi0 <- var.psi1 <- rep(NA,nrep)
psiL <- psiU <- var.psiL <- var.psiU <- rep(NA,nrep)
LB.theta <- UB.theta <- LB.psi <- UB.psi <- rep(NA,nrep)
corr.theta <- corr.psi <- rep(NA,nrep)
b <- 3
for (k in 1:nrep){
  
  set.seed(k)  
  
  ## Generate Data ##  
  data <- rmultinom(1,n,prob)
  n111 <- data[1]
  n101 <- data[2]
  n011 <- data[3]
  n001 <- data[4]
  n110 <- data[5]
  n100 <- data[6]
  n010 <- data[7]
  n000 <- data[8]
  n00 <- n001+n000
  
  ## Estimate ##
  A0 <- n011/n
  B0 <- n101/n
  a0 <- n010/n
  b0 <- n100/n
  DA <- (n101+n111)/n
  Da <- (n100+n110)/n
  P <- n00/n
  
  # Bounds #
  theta0[k] <- (B0-A0)/(DA+A0)
  theta1[k] <- (B0-A0)/(DA+A0+P)
  psi0[k] <- (a0-b0)/(Da+a0)
  psi1[k] <- (a0-b0)/(Da+a0+P)
  
  # Variance #
  var.theta0[k] <- (B0+A0-(B0-A0)^2/(DA+A0))/(DA+A0)^2
  var.theta1[k] <- (B0+A0-(B0-A0)^2/(DA+A0+P))/(DA+A0+P)^2
  var.psi0[k] <- (a0+b0-(a0-b0)^2/(Da+a0))/(Da+a0)^2
  var.psi1[k] <- (a0+b0-(a0-b0)^2/(Da+a0+P))/(Da+a0+P)^2
  
  
  # Correlation #
  sigma1.theta <- (B0+A0)*P^2+(B0-A0)^2*P*(1-4*P)
  sigma2.theta <- (DA+A0)*(DA+A0+P)*(4*(DA+A0)*(1-DA-A0-P)+P)
  cov12.theta <-  P*(A0-B0)*(A0+DA)*(4*(A0+DA+P)-3)-P^2*(A0-B0)
  vardiff.theta <- sigma1.theta/((DA+A0)^2*(DA+A0+P)^2)-2*(B0-A0)*P*cov12.theta/((DA+A0)^3*(DA+A0+P)^3)+
    (B0-A0)^2*P^2*sigma2.theta/((DA+A0)^4*(DA+A0+P)^4)
  corr.theta[k] <- ((var.theta0[k]+var.theta1[k]-vardiff.theta)/2)/(sqrt(var.theta0[k])*sqrt(var.theta1[k]))
  
  sigma1.psi <- (b0+a0)*P^2+(b0-a0)^2*P*(1-4*P)
  sigma2.psi <- (Da+a0)*(Da+a0+P)*(4*(Da+a0)*(1-Da-a0-P)+P)
  cov12.psi <- P*(a0-b0)*(a0+Da)*(4*(a0+Da+P)-3)-P^2*(a0-b0)
  vardiff.psi <- sigma1.psi/((Da+a0)^2*(Da+a0+P)^2)-2*(b0-a0)*P*cov12.psi/((Da+a0)^3*(Da+a0+P)^3)+
    (b0-a0)^2*P^2*sigma2.psi/((Da+a0)^4*(Da+a0+P)^4)
  corr.psi[k] <- ((var.psi0[k]+var.psi1[k]-vardiff.psi)/2)/(sqrt(var.psi0[k])*sqrt(var.psi1[k]))
  
  
  ## Lower and upper bounds ##
  if(B0>=A0){
    thetaL[k] <- theta1[k]
    thetaU[k] <- theta0[k]
    var.thetaL[k] <- var.theta1[k]
    var.thetaU[k] <- var.theta0[k]
  }
  if(B0<A0){
    thetaL[k] <- theta0[k]
    thetaU[k] <- theta1[k]
    var.thetaL[k] <- var.theta0[k]
    var.thetaU[k] <- var.theta1[k]
  }
  if(a0<=b0){
    psiL[k] <- psi0[k]
    psiU[k] <- psi1[k]
    var.psiL[k] <- var.psi0[k]
    var.psiU[k] <- var.psi1[k]
  }
  if(a0>b0){
    psiL[k] <- psi1[k]
    psiU[k] <- psi0[k]
    var.psiL[k] <- var.psi1[k]
    var.psiU[k] <- var.psi0[k]
  }
  
  
  ## Critical value and CI ##
  hl <- sqrt(n)*(thetaU[k]-thetaL[k])*ifelse(thetaU[k]-thetaL[k]>n^(-1/b),1,0)/max(sqrt(var.thetaL[k]),sqrt(var.thetaU[k]))
  hu <- 0
  rho <- corr.theta[k]
  if((rho-1)>=1e-14){c.theta <- uniroot(cv,c(1,5))$root}
  if((rho-1)<1e-14){c.theta <- uniroot(cv1,c(1,5))$root}
  if(thetaU[k]-thetaL[k]>= (-min(sqrt(var.thetaL[k]),sqrt(var.thetaU[k]))*sqrt(c.theta/n))){
    LB.theta[k] <- thetaL[k]-sqrt(var.thetaL[k])*sqrt(c.theta/n)
    UB.theta[k] <- thetaU[k]+sqrt(var.thetaU[k])*sqrt(c.theta/n)
  }
  
  hl <- sqrt(n)*(psiU[k]-psiL[k])*ifelse(psiU[k]-psiL[k]>n^(-1/b),1,0)/max(sqrt(var.psiL[k]),sqrt(var.psiU[k]))
  hu <- 0
  rho <- corr.psi[k]
  if((rho-1)>=1e-14){c.psi <- uniroot(cv,c(1,5))$root}
  if((rho-1)<1e-14){c.psi <- uniroot(cv1,c(1,5))$root}
  if(psiU[k]-psiL[k]>= (-min(sqrt(var.psiL[k]),sqrt(var.psiU[k]))*sqrt(c.psi/n))){
    LB.psi[k] <- psiL[k]-sqrt(var.psiL[k])*sqrt(c.psi/n)
    UB.psi[k] <- psiU[k]+sqrt(var.psiU[k])*sqrt(c.psi/n)
  }
  
  
}  


c(thetaL.true,mean(thetaL),mean(thetaL)-thetaL.true,var(thetaL),median(var.thetaL/n))
c(thetaU.true,mean(thetaU),mean(thetaU)-thetaU.true,var(thetaU),median(var.thetaU/n))
theta.true
sum(theta.true>=LB.theta&theta.true<=UB.theta)/nrep*100

c(psiL.true,mean(psiL),mean(psiL)-psiL.true,var(psiL),median(var.psiL/n))
c(psiU.true,mean(psiU),mean(psiU)-psiU.true,var(psiU),median(var.psiU/n))
psi.true
sum(psi.true>=LB.psi&psi.true<=UB.psi)/nrep*100




