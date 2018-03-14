#Created by Tao He, Shaoyu Li, Ping-Shou Zhong, and Yuehua Cui
#This R code is designed for the paper
#"An optimal kernel-based method for gene set association analysis"

rm(list=ls())
###################
# load packages ###
###################
library(ttutils)
library(methods)
library(psych)
library(EQL)
library(mvtnorm)
library(Matrix)



##########function to generate kernel matrix######
#Note:for demonstration reasons, only certain commonly used kernels are defined in the function,
#such as Gaussion/Quadratic/Linear/IBS. 
#The users can define the candidate kernel set of their interests.

gen.ker <- function(covx,kernel.index)
#covx:X; 
#kernel.index=c("Gau","Lin","Quad","IBS"): index of kernel function;
{ 
  n <- nrow(covx)
  p <- ncol(covx)
  ker <- matrix(0,n,n)
  for (i in 1:n)
    for (j in i:n)
    {
      x <- covx[i,]
      y <- covx[j,]
      
      if (kernel.index=="Gau")
      {
        ker[i,j] <- exp(-sum((x-y)^2)/p) #gaussian kernel
      }
      
      if (kernel.index=="Lin")
      {
        ker[i,j] <- sum(x*y)  #linear kernel        
      }
      
      if (kernel.index=="Quad")
      {
        ker[i,j]=(sum(x*y)+1)^2 # Quadratic kernel
      }
      if (kernel.index=="IBS")
      {
        ker[i,j] <- 1-sum(abs(x-y))/(2*p) #IBS kernel
      }
    }
    ker <- as.matrix(forceSymmetric(ker))
    ker0 <- ker
    diag(ker0) <- rep(0,n)
    J <- matrix(1,n,n)
    ker.cen <- ker-J%*%ker0/(n-1)-ker0%*%J/(n-1)+J%*%ker0%*%J/n/(n-1)
    v1.cen <- tr(ker.cen)/n
    return(list(ker.cen=ker.cen,v1.cen=v1.cen))
}


##########function to test using a given kernel######
#Note: the rejection region is constructed 
#based on normal asymptotic propoerty

gen.Ustat<-function(y,covx,z,ker)
#y:n by 1 vector
#covx: n by p matrix
#z:n by M matrix
#ker: n by n matrix
{
  n <- length(y)
  p <- ncol(covx)
  h <- diag(n)-z%*%solve(t(z)%*%z)%*%t(z)
  ker1 <- ker
  diag(ker1) <- rep(0,n)
  A <- h%*%ker1%*%h
  hkhk <- A%*%ker1
  hkh.hkh <- A*A
  hk <- h%*%ker1
  sigma.hat <- t(y)%*%h%*%y/(n-ncol(z))
  sd <- as.numeric(sqrt(sigma.hat))
  z0 <- h%*%y/sd
  m4 <- mean(z0^4)
  e1 <- tr(hk)
  e2 <- tr(hkhk)
  mu1 <- 2*e1
  mu2 <- e1*(1+2/(n-1))
  delta <- m4-3
  var <- e1^2*(-2/(n-1))+e2*(2-12/(n-1)) 
  var <- var+delta*(-e1^2/n+6*e2/n+tr(hkh.hkh)) 
  test <- t(h%*%y)%*%ker1%*%(h%*%y)/(sigma.hat)
  test.sd <- test/sqrt(var)
  p.norm <- pnorm(test.sd,lower.tail = F)
  return(list(test.sd=test.sd,p.norm=p.norm, test.sd=test.sd,A.ma=A/sqrt(var))) 
  
}

##########generate correlation matrix################
gen.sig.corr<-function(A.ma.list) 
{
  h <- diag(n)-z%*%solve(t(z)%*%z)%*%t(z)
  sigma.hat <- t(y)%*%h%*%y/(n-ncol(z))
  sd <- as.numeric(sqrt(sigma.hat))
  z0 <- h%*%y/sd
  m4 <- mean(z0^4)
  delta <- m4-3
  
  D <- length(A.ma.list)
  sig.corr <- matrix(0,D,D)
  for (i in  1:D)
    for (j in 1:D)
    {
      Ai <- A.ma.list[[i]]
      Aj <- A.ma.list[[j]]
      a <- -(2+delta)/n
      b <- 2+(6*delta-12)/n
      sig.corr[i,j] <- min(a*tr(Ai)*tr(Aj)+b*tr(Ai%*%Aj)+delta*tr(Ai*Aj),1)
      
    }   
  diag(sig.corr) <- rep(1,D)
  sig.corr
}

#########################
##     main fun   #######
#########################
#Note: for demonstration reasons, consider Gau/Lin/Quad kernels as candidate kernel set
#The users can alter this part accordingly.

opt.ker <- function(x,y,z,alpha=0.05)
{
  z <- cbind(rep(1,n),z)
  ker.G <- gen.ker(x,"Gau")        
  ker.L <- gen.ker(x,"Lin")
  ker.Q <- gen.ker(x,"Quad")
  ######centralized and normalized kernels
  kerM.G <- ker.G$ker.cen/ker.G$v1.cen
  kerM.L <- ker.L$ker.cen/ker.L$v1.cen
  kerM.Q <- ker.Q$ker.cen/ker.Q$v1.cen
  a3 <- rep(1/3,3)
  kerM.a3 <- a3[1]*kerM.G+a3[2]*kerM.L+a3[3]*kerM.Q #average kernel
  
  test.G.ker <- gen.Ustat(y,x,z,kerM.G)
  test.L.ker <- gen.Ustat(y,x,z,kerM.L)
  test.Q.ker <- gen.Ustat(y,x,z,kerM.Q)
  test.a3.ker <- gen.Ustat(y,x,z,kerM.a3)
  
  A.g <- test.G.ker$A.ma
  A.l <- test.L.ker$A.ma
  A.q <- test.Q.ker$A.ma
  sig.corr <- gen.sig.corr(list(A.g,A.l,A.q))
  
  sd.vec <- c(test.G.ker$test.sd,test.L.ker$test.sd,test.Q.ker$test.sd)
  sd.max <- max(sd.vec)
  q.cri <- qmvnorm(1-alpha, c(1.64,2.5),sigma = sig.corr, tail = "lower.tail",algorithm=TVPACK())$quantile
  re.opt <- ifelse(sd.max>q.cri,"reject H0: h(x0)=0","retain H0: h(x0)=0")
  return(list(max.test.result=re.opt, pvalue.averageKer=test.a3.ker$p.norm,pvalue.GauKer=test.G.ker$p.norm,pvalue.LinKer=test.L.ker$p.norm,pvalue.QuadKer=test.Q.ker$p.norm))
}



#load data and test
n <- 100;p <- 100
x <- matrix(rnorm(n*p),n,p)
z <- matrix(rnorm(n*2),n,2)
y <- rnorm(n)
opt.ker(x,y,z,0.05)


