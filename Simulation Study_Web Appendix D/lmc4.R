############################################################################################
########### Coherence project simulation study (LMC coherence) #############################
############################################################################################
library(fields)
library(MASS)
library(doParallel)
#####################################################################
########### Creating required functions #############################
#####################################################################
setwd("/ibex/scratch/qadirga/Simulation_Studies_Old/LMC_coherence_function")

##### Creating function for computing norm #####
norm_vec <- function(x) sqrt(sum(x^2))

#####################################################
########### Parametrising spline coefficients #######
#####################################################
bcoeff<-function(x)
{
  
  temp<-sin(x)
  temp[temp==-1]<--1+1e-10
  temp[temp==1]<-1-1e-10
  return(temp)
}



#####################################
##### Matern spectral density  ######
#####################################

f.matern<-function(w,nu,sigmasq,a,d)
{
  
  ln<-length(w)
  num1<-den1<-temp<-rep(NA,times=ln)
  for(i in 1:ln)
  {
    num1[i]<-gamma(nu+(d/2))*(a^(2*nu))
    den1[i]<-gamma(nu)*(pi^(d/2))*(((a^2)+(w[i]^2))^(nu+(d/2)))
    temp[i]<-num1[i]/den1[i]
  }
  return(sigmasq*temp)
  
}

#######################################################
############### Matern covariance function ############
#######################################################

my.mat2<-function(h,sigmasq,nu,a)
{
  h[h==0]<-1e-10
  num1<-(2^(1-nu))*((a*h)^nu)*besselK(x=a*h,nu=nu)
  den1<-gamma(nu)
  temp<-num1/den1
  return(sigmasq*temp)
}


#############################################
######## lmc matern coherence ###############
#############################################

lmc.coh<-function(w,a1,nu1,sigmasq1,a2,nu2,sigmasq2,b11,b12,b21,b22)
{
  f1<-f.matern(w=w,nu=nu1,sigmasq=sigmasq1,a=a1,d=2)
  f2<-f.matern(w=w,nu=nu2,sigmasq = sigmasq2,a=a2,d=2)
  num<-b11*b21*f1+b12*b22*f2
  den1<-sqrt((b11^2)*f1+(b12^2)*f2)
  den2<-sqrt((b21^2)*f1+(b22^2)*f2)
  den<-den1*den2
  return(list(mar1=den1^2,mar2=den2^2,cross=num,coh=num/den))
}


############################################
##### Function for B spline ################
############################################
Bspline<-function(j,k,delta,x)
{
  fin.val<-numeric(length=length(x))
  for(index in 1:length(x))
  {
    dummy<-rep(NA,times=k+1)
    for(i in 1:(k+1))
    {
      r<-i-1
      num1<-(-1)^(k-r)
      den1<-factorial(k-1)
      kcr<-factorial(k)/(factorial(r)*factorial(k-r))
      f.c<-max(0,(r-x[index]/delta+j)^(k-1))
      dummy[i]<-(num1)*kcr*f.c
    }
    fin.val[index]<-sum(dummy)/den1  
  }
  return(fin.val)
}



# Setting simulation parameter for the LMC model

a1.lmc<-0.5
nu1.lmc<-1
sigmasq1.lmc<-1

a2.lmc<-0.5
nu2.lmc<-2
sigmasq2.lmc<-1

b11.lmc<-1
b12.lmc<-0.4
b21.lmc<-0.9
b22.lmc<-7.5
par(mfrow=c(1,1))
u<-seq(0,9.8,length=300)
par(mfrow=c(2,2))
plot(u,lmc.coh(w=u,a1=a1.lmc,nu1=nu1.lmc,sigmasq1 = sigmasq1.lmc,a2=a2.lmc,sigmasq2 = sigmasq2.lmc,nu2=nu2.lmc,b11=b11.lmc,b12=b12.lmc,b21=b21.lmc,b22=b22.lmc)$coh,ylab="LMC Coherence",main="Coherence for LMC")
plot(u,lmc.coh(w=u,a1=a1.lmc,nu1=nu1.lmc,sigmasq1 = sigmasq1.lmc,a2=a2.lmc,sigmasq2 = sigmasq2.lmc,nu2=nu2.lmc,b11=b11.lmc,b12=b12.lmc,b21=b21.lmc,b22=b22.lmc)$mar1,ylab="f1",main="Marginal spectral density 1")
plot(u,lmc.coh(w=u,a1=a1.lmc,nu1=nu1.lmc,sigmasq1 = sigmasq1.lmc,a2=a2.lmc,sigmasq2 = sigmasq2.lmc,nu2=nu2.lmc,b11=b11.lmc,b12=b12.lmc,b21=b21.lmc,b22=b22.lmc)$mar2,ylab="f2",main="Marginal spectral density 2")
plot(u,lmc.coh(w=u,a1=a1.lmc,nu1=nu1.lmc,sigmasq1 = sigmasq1.lmc,a2=a2.lmc,sigmasq2 = sigmasq2.lmc,nu2=nu2.lmc,b11=b11.lmc,b12=b12.lmc,b21=b21.lmc,b22=b22.lmc)$cross,ylab="f12",main="Cross spectral density")


# n.loc = number of locations
# Creating Irregular grid for simulation on [1,40]^2
n.loc<-500
set.seed(123)
x<-runif(n=n.loc,min=1,max = 40)
set.seed(133)
y<-runif(n=n.loc,min=1,max = 40)
length(unique(x))
length(unique(y))
plot(x,y,xlim=c(0,40),ylim=c(0,40),pch=19,cex=1,main="Simulation grid")



#### Computing distance matrix for the simulation grid

dist.mat<-rdist(cbind(x,y))

COV11<-(b11.lmc^2)*my.mat2(h=dist.mat,sigmasq = sigmasq1.lmc,nu=nu1.lmc,a=a1.lmc)+(b12.lmc^2)*my.mat2(h=dist.mat,sigmasq = sigmasq2.lmc,nu=nu2.lmc,a=a2.lmc)
COV22<-(b21.lmc^2)*my.mat2(h=dist.mat,sigmasq = sigmasq1.lmc,nu=nu1.lmc,a=a1.lmc)+(b22.lmc^2)*my.mat2(h=dist.mat,sigmasq = sigmasq2.lmc,nu=nu2.lmc,a=a2.lmc)
COV12<-(b11.lmc*b21.lmc)*my.mat2(h=dist.mat,sigmasq = sigmasq1.lmc,nu=nu1.lmc,a=a1.lmc)+(b12.lmc*b22.lmc)*my.mat2(h=dist.mat,sigmasq = sigmasq2.lmc,nu=nu2.lmc,a=a2.lmc)

C.gen<-rbind(cbind(COV11,COV12),cbind(t(COV12),COV22))

#plot(c(dist.mat),c(COV22))

#### Simulating 100 realisations of bivariate random field #####

var.bvm<-matrix(NA,nrow=2*n.loc,ncol=100)

sim.time<-system.time(for(i in 1:100)
{
  set.seed(i)
  var.bvm[,i]<-mvrnorm(n=1,mu=rep(0,times=1000),Sigma = C.gen)
})



par(mfrow=c(1,1))
par(cex.axis=2, cex.lab=3, cex.main=1, cex.sub=1,mar=c(2,4,1,2)+.1)

quilt.plot(x,y,var.bvm[1:n.loc,1],xlab="x",ylab="y")
quilt.plot(x,y,var.bvm[(n.loc+1):(2*n.loc),1],xlab="x",ylab="y")
cor(var.bvm[1:n.loc,1],var.bvm[(n.loc+1):(2*n.loc),1])







#### Now we do the optimization ###

ncores<-detectCores()-6
ncores
registerDoParallel(cores = ncores) #run only once

optim.pars<-foreach(i=1:100) %dopar%
{
  library(fields)
  z<-var.bvm[,i]
  
  init.1<-c(0.55,1.1,1.15,0.52,2,56,0.99,0.45,0.40,0.65,0.76,0.81,0.91,0.65)
  
  
  
  # function that returns everything about loglikelihood
  
  
  log_likelihood_all_component<-function(theta)
  {
    
    
    #setting parameters
    a1<-theta[1]
    nu1<-theta[2]
    sigmasq1<-theta[3]
    a2<-theta[4]
    nu2<-theta[5]
    sigmasq2<-theta[6]
    b_3<-theta[7]
    b_2<-theta[8]
    b_1<-theta[9]
    b0<-theta[10]
    b1<-theta[11]
    b2<-theta[12]
    b3<-theta[13]
    b4<-theta[14]
    
    # setting hard constraints on the parameters
    if(a1<=0|nu1<=0|sigmasq1<=0|a2<=0|nu2<=0|sigmasq2<=0|b_3>=1|b_2>=1|b_1>=1|b0>=1|b1>=1|b0>=1|b1>=1|
       b2>=1|b3>=1|b4>=1|b_3<=-1|b_2<=-1|b_1<=-1|b0<=-1|b1<=-1|b2<=-1|b3<=-1|b4<=-1)
    {
      return(10000000000)
      
    } 
    
    else{
      
      u<-seq(0,9.8,length=300) # setting nodes
      f.var1<-f.matern(w=u, nu=nu1, sigmasq = 1, a=a1, d=2)
      f.var2<-f.matern(w=u, nu=nu2, sigmasq = 1, a=a2, d=2)
      delta.s<-2 #setting delta
      coh12<-b_3*Bspline(j=-3,k=4,delta = delta.s,x=u)+b_2*Bspline(j=-2,k=4,delta = delta.s,x=u)+b_1*Bspline(j=-1,k=4,delta = delta.s,x=u)+b0*Bspline(j=0,k=4,delta = delta.s,x=u)+b1*Bspline(j=1,k=4,delta = delta.s,x=u)+b2*Bspline(j=2,k=4,delta = delta.s,x=u)+b3*Bspline(j=3,k=4,delta = delta.s,x=u)++b4*Bspline(j=4,k=4,delta = delta.s,x=u)
      f.var12<-coh12*sqrt(f.var1*f.var2)
      dist.mat<-rdist(cbind(x,y)) ### distance matrix
      
      
      
      #### Computing Block Covariance matrices ####
      
      tempmat1<-tempmat2<-tempmat12<-matrix(0,nrow=nrow(dist.mat),ncol=ncol(dist.mat))
      
      mul.const1<-2*pi*u*f.var1 # multiplier for marginal 1
      mul.const2<-2*pi*u*f.var2 # multiplier for marginal 2
      mul.const12<-2*pi*u*f.var12 #multiplier for cross
      
      for(i in 1:length(u))
      {
        tempbessel<-besselJ(x=dist.mat*u[i],nu=0)
        cmat1<-mul.const1[i]*tempbessel
        cmat2<-mul.const2[i]*tempbessel
        cmat12<-mul.const12[i]*tempbessel
        tempmat1<-tempmat1+cmat1
        tempmat2<-tempmat2+cmat2
        tempmat12<-tempmat12+cmat12
      }
      
      
      ######## Rescaling covariances ########
      scl1<-max(tempmat1)
      scl2<-max(tempmat2)
      COV11.my<-tempmat1*sigmasq1/scl1
      COV22.my<-tempmat2*sigmasq2/scl2
      COV12.my<-(sqrt(sigmasq1*sigmasq2)/sqrt(scl1*scl2))*tempmat12
      
      C<-rbind(cbind(COV11.my,COV12.my),cbind(t(COV12.my),COV22.my))
      
      if(sum(C==Inf)>0||sum(is.nan(C))>0)
      {
        nlog_likelihood <- 1e+12
      }
      else
      { 
        cholS<-chol(C)
        nlog_likelihood <-
          -as.numeric(-0.5 * determinant(C)$modulus -
                        0.5 * t(z) %*% chol2inv(cholS) %*% z -
                        0.5 * length(z)*log(2*pi))
        
        if (abs(nlog_likelihood) == Inf || is.nan(nlog_likelihood)) nlog_likelihood <- 1e+08
      } 
      return(list(cov11=COV11.my,cov22=COV22.my,cov12=COV12.my,coherence=coh12,nloglikelihood=nlog_likelihood))
      
    }
    
  }
  
  
  # function that returns only loglikelihood values     
  n.log.likelihood<-function(p)
  {
    temp<-log_likelihood_all_component(p)
    if(is.recursive(temp))
    {
      return(temp$nloglikelihood)
    }
    else
      return(10000000000)
  }
  
  
  
  
  
  optim(n.log.likelihood,par=init.1,control = list(trace=6,maxit=6000))
}

save.image("estimatedlmcnew4.RData")
est.coh<-foreach(i=1:100) %dopar%
{
  z<-var.bvm[,i] #is of no use here
  
  
  log_likelihood_all_component<-function(theta)
  {
    
    
    #setting parameters
    a1<-theta[1]
    nu1<-theta[2]
    sigmasq1<-theta[3]
    a2<-theta[4]
    nu2<-theta[5]
    sigmasq2<-theta[6]
    b_3<-theta[7]
    b_2<-theta[8]
    b_1<-theta[9]
    b0<-theta[10]
    b1<-theta[11]
    b2<-theta[12]
    b3<-theta[13]
    b4<-theta[14]
    
    # setting hard constraints on the parameters
    if(a1<=0|nu1<=0|sigmasq1<=0|a2<=0|nu2<=0|sigmasq2<=0|b_3>=1|b_2>=1|b_1>=1|b0>=1|b1>=1|b0>=1|b1>=1|
       b2>=1|b3>=1|b4>=1|b_3<=-1|b_2<=-1|b_1<=-1|b0<=-1|b1<=-1|b2<=-1|b3<=-1|b4<=-1)
    {
      return(10000000000)
      
    } 
    
    else{
      
      u<-seq(0,9.8,length=300) # setting nodes
      f.var1<-f.matern(w=u, nu=nu1, sigmasq = 1, a=a1, d=2)
      f.var2<-f.matern(w=u, nu=nu2, sigmasq = 1, a=a2, d=2)
      delta.s<-2 #setting delta
      coh12<-b_3*Bspline(j=-3,k=4,delta = delta.s,x=u)+b_2*Bspline(j=-2,k=4,delta = delta.s,x=u)+b_1*Bspline(j=-1,k=4,delta = delta.s,x=u)+b0*Bspline(j=0,k=4,delta = delta.s,x=u)+b1*Bspline(j=1,k=4,delta = delta.s,x=u)+b2*Bspline(j=2,k=4,delta = delta.s,x=u)+b3*Bspline(j=3,k=4,delta = delta.s,x=u)++b4*Bspline(j=4,k=4,delta = delta.s,x=u)
      f.var12<-coh12*sqrt(f.var1*f.var2)
      dist.mat<-rdist(cbind(x,y)) ### distance matrix
      
      
      
      #### Computing Block Covariance matrices ####
      
      tempmat1<-tempmat2<-tempmat12<-matrix(0,nrow=nrow(dist.mat),ncol=ncol(dist.mat))
      
      mul.const1<-2*pi*u*f.var1 # multiplier for marginal 1
      mul.const2<-2*pi*u*f.var2 # multiplier for marginal 2
      mul.const12<-2*pi*u*f.var12 #multiplier for cross
      
      for(i in 1:length(u))
      {
        tempbessel<-besselJ(x=dist.mat*u[i],nu=0)
        cmat1<-mul.const1[i]*tempbessel
        cmat2<-mul.const2[i]*tempbessel
        cmat12<-mul.const12[i]*tempbessel
        tempmat1<-tempmat1+cmat1
        tempmat2<-tempmat2+cmat2
        tempmat12<-tempmat12+cmat12
      }
      
      
      ######## Rescaling covariances ########
      scl1<-max(tempmat1)
      scl2<-max(tempmat2)
      COV11.my<-tempmat1*sigmasq1/scl1
      COV22.my<-tempmat2*sigmasq2/scl2
      COV12.my<-(sqrt(sigmasq1*sigmasq2)/sqrt(scl1*scl2))*tempmat12
      
      C<-rbind(cbind(COV11.my,COV12.my),cbind(t(COV12.my),COV22.my))
      
      if(sum(C==Inf)>0||sum(is.nan(C))>0)
      {
        nlog_likelihood <- 1e+12
      }
      else
      { 
        cholS<-chol(C)
        nlog_likelihood <-
          -as.numeric(-0.5 * determinant(C)$modulus -
                        0.5 * t(z) %*% chol2inv(cholS) %*% z -
                        0.5 * length(z)*log(2*pi))
        
        if (abs(nlog_likelihood) == Inf || is.nan(nlog_likelihood)) nlog_likelihood <- 1e+08
      } 
      return(list(cov11=COV11.my,cov22=COV22.my,cov12=COV12.my,coherence=coh12,nloglikelihood=nlog_likelihood))
      
    }
    
  }
  
  
  
  log_likelihood_all_component(optim.pars[[i]]$par)$coherence
}

counti<-numeric()
for(i in 1:100)
{
  counti[i]<-optim.pars[[i]]$counts[1]
}
par(mfrow=c(1,2))
plot(u,lmc.coh(w=u,a1=a1.lmc,nu1=nu1.lmc,sigmasq1 = sigmasq1.lmc,a2=a2.lmc,sigmasq2 = sigmasq2.lmc,nu2=nu2.lmc,b11=b11.lmc,b12=b12.lmc,b21=b21.lmc,b22=b22.lmc)$coh,ylab="LMC Coherence",main="Coherence for LMC",ylim=c(-1,1),type="l",lwd=3)

coh.mat<-matrix(NA,nrow=length(u),ncol=100)
for(i in 1:100)
{
  lines(u,est.coh[[i]],col=i+1)
  coh.mat[,i]<-est.coh[[i]]
}

plot(u,lmc.coh(w=u,a1=a1.lmc,nu1=nu1.lmc,sigmasq1 = sigmasq1.lmc,a2=a2.lmc,sigmasq2 = sigmasq2.lmc,nu2=nu2.lmc,b11=b11.lmc,b12=b12.lmc,b21=b21.lmc,b22=b22.lmc)$coh,ylab="LMC Coherence",main="Coherence for LMC",ylim=c(-1,1),type="l",lwd=3)

lines(u,rowMeans(coh.mat),col="red",lty=2,lwd=2)
avg.coh<-rowMeans(coh.mat)
sd.coh<-numeric(length = length(u))
for(i in 1:length(u))
{
  sd.coh[i]<-sd(coh.mat[i,]) 
}

######## Extracting marginal parameters #####
e.a1<-e.a2<-e.sigmasq1<-e.sigmasq2<-e.nu1<-e.nu2<-numeric(length=100)
for(i in 1:100)
{
  e.a1[i]<-optim.pars[[i]]$par[1]
  e.nu1[i]<-optim.pars[[i]]$par[2]
  e.sigmasq1[i]<-optim.pars[[i]]$par[3]
  e.a2[i]<-optim.pars[[i]]$par[4]
  e.nu2[i]<-optim.pars[[i]]$par[5]
  e.sigmasq2[i]<-optim.pars[[i]]$par[6]
}
round(mean(e.sigmasq1),2)
round(sd(e.sigmasq1),2)

round(mean(e.a1),2)
round(sd(e.a1),2)

round(mean(e.nu1),2)
round(sd(e.nu1),2)

round(mean(e.sigmasq2),2)
round(sd(e.sigmasq2),2)

round(mean(e.a2),2)
round(sd(e.a2),2)


round(mean(e.nu2),2)


round(sd(e.nu2),2)


save.image("full_lmcnew4.RData")

library(ggplot2)
up95<-c(rep(NA,times=length(avg.coh)),avg.coh+1.96*sd.coh)
up95[up95>1]<-1

plotdata<-data.frame(
 frequency=c(u,u),
  coh=c(lmc.coh(w=u,a1=a1.lmc,nu1=nu1.lmc,sigmasq1 = sigmasq1.lmc,a2=a2.lmc,sigmasq2 = sigmasq2.lmc,nu2=nu2.lmc,b11=b11.lmc,b12=b12.lmc,b21=b21.lmc,b22=b22.lmc)$coh,avg.coh)
 ,type=c(rep("True coherence function",times=length(avg.coh)),rep("Mean(estimated coherence functions)",times=length(avg.coh))),
up95=up95,low95=c(rep(NA,times=length(avg.coh)),avg.coh-1.96*sd.coh))
par(mfrow=c(1,1))
par(cex.axis=2, cex.lab=3, cex.main=1, cex.sub=1,mar=c(3,3,1,1))
p<-ggplot(data=plotdata, aes(x=frequency, y=coh, colour=type)) + geom_line(aes(linetype=type),lwd=1.0)
p<-p+geom_ribbon(aes(ymin=plotdata$low95, ymax=plotdata$up95), linetype=2, alpha=0.2)+ylim(-1, 1)
p<-p+labs(x= expression(omega),    # works fine
          y= expression(gamma(omega)))+theme(legend.title=element_blank(),legend.position = c(0.65, 0.15), legend.background = element_rect(color = "black", fill = "grey90", size = 0.1, linetype = "solid"))

plot(p)
