############################################################################################
########### Coherence project simulation study (Decreasing coherence) ######################
############################################################################################
library(fields)
library(MASS)
library(doParallel)
#####################################################################
########### Creating required functions #############################
#####################################################################
setwd("/ibex/scratch/qadirga/Simulation_Studies_Old/Decreasing_coherence_function")

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
  ln<-length(h)
  num1<-den1<-temp<-rep(NA,times=ln)
  for(i in 1:ln)
  {
    if(h[i]>0)
    {
      num1[i]<-(2^(1-nu))*((a*h[i])^nu)*besselK(x=a*h[i],nu=nu)
      den1[i]<-gamma(nu)
      temp[i]<-num1[i]/den1[i]
    }
    
    else
    {
      c.dist<-1e-10
      num1[i]<-(2^(1-nu))*((c.dist)^nu)*besselK(x=c.dist,nu=nu)
      den1[i]<-gamma(nu)
      temp[i]<-num1[i]/den1[i]
    }
  }
  return(sigmasq*temp)
}

#############################################
######## Bivariate matern coherence #########
#############################################

biwm_coh<-function(w,a1,a2,v1,v2,a12,v12,d,rho)
{
  temp<-numeric(length=length(w))
  for(i in 1:length(w))
  {
    num<-rho^2*((gamma(v12+(d/2)))^2)*gamma(v1)*gamma(v2)*(a12^(4*v12))*((a1^2+w[i]^2)^(v1+(d/2)))*((a2^2+w[i]^2)^(v2+(d/2)))
    den<-gamma(v1+(d/2))*gamma(v2+(d/2))*(gamma(v12)^2)*(a1^(2*v1))*(a2^(2*v2))*((a12^2+w[i]^2)^(2*v12+d))
    temp[i]<-num/den
  }
  return(((temp)))  
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




#############################################
######## Creating simulation grid ###########
#############################################

n<-30 ##### number of observations = n^2
x<-y<-seq(1,30,length=n)

#######################################################
########## Creating regularly spaced locations ########
#######################################################

x.s<-rep(x,each=length(y))
y.s<-rep(y,times=length(x))

################################################
######### Plotting simulation grid #############
################################################

par(mfrow=c(1,1))
plot(x.s,y.s,main="Simulation grid",xlab="x",ylab="y")

###############################################
##### Computing distance matrix ###############
###############################################

dist.mat<-rdist(cbind(x.s,y.s))


##############################################
###### Computing unique distances ############
##############################################

uniq.dist<-unique(c(dist.mat))
uniq.dist<-sort(uniq.dist)  ##### sorting distances in increasing order #####
d=2
###################################################
####### Computing Max Freq length #################
###################################################
freq.max<-4.5


###########################################################################################
################## Now we simulate from bivariate matern with constant coherence ##########
###########################################################################################

sig1.bvm<-1
sig2.bvm<-1
a1.bvm<-1
a2.bvm<-1
a12.bvm<-1
nu1.bvm<-3
nu2.bvm<-3
nu12.bvm<-4
rho<-0.4



u<-seq(0,freq.max,length=length(uniq.dist)-1)
plot(u,sqrt(biwm_coh(w=u,a1=a1.bvm,a2=a2.bvm,v1=nu1.bvm,v2=nu2.bvm,a12=a12.bvm,v12 = nu12.bvm,d=2,rho = rho)),type="l",ylim=c(0,1),main="Generating Coherence",ylab="Coh",lwd=2)




C11.bvm<-my.mat2(h=uniq.dist,sigmasq = sig1.bvm,nu=nu1.bvm,a=a1.bvm)
C22.bvm<-my.mat2(h=uniq.dist,sigmasq = sig2.bvm,nu=nu2.bvm,a=a2.bvm)
C12.bvm<-rho*sqrt(sig1.bvm*sig2.bvm)*my.mat2(h=uniq.dist,sigmasq = 1,nu=nu12.bvm,a=a12.bvm)

COV11.bvm<-COV22.bvm<-COV12.bvm<-matrix(NA,nrow=n^2,ncol = n^2)

for(i in 1:length(uniq.dist))
{
  COV11.bvm[which(dist.mat==uniq.dist[i],arr.ind = T)]<-C11.bvm[i]
  COV22.bvm[which(dist.mat==uniq.dist[i],arr.ind = T)]<-C22.bvm[i]
  COV12.bvm[which(dist.mat==uniq.dist[i],arr.ind = T)]<-C12.bvm[i]
}


C.bvm<-rbind(cbind(COV11.bvm,COV12.bvm),cbind(t(COV12.bvm),COV22.bvm))
chol(C.bvm)

u<-seq(0,freq.max,length=length(uniq.dist)-1) ####### Setting nodes #########

e.min.lik<-numeric(length=100) ####### For storing minimum of negative likelihood #######
e.a1<-e.a2<-e.sigmasq1<-e.sigmasq2<-e.b_3<-e.b_2<-e.b_1<-e.b0<-e.b1<-e.b2<-e.b3<-e.b4<-numeric(length = 100) ###### For storing the parameters #########
e.coh12<-matrix(NA,nrow=length(u),ncol=100)  ########## For storing estimated coherences #############
var.bvm<-matrix(NA,nrow=1800,ncol=100) #### For storing simulated realisations ######



#############################################
#### Plotting the generating coherence ######
#############################################
par(mfrow=c(1,1))
plot(u,sqrt(biwm_coh(w=u,a1=a1.bvm,a2=a2.bvm,v1=nu1.bvm,v2=nu2.bvm,a12=a12.bvm,v12 = nu12.bvm,d=2,rho = rho)),type="l",ylim=c(0,1),main="Generating Coherence",ylab="Coh",lwd=2)


ncores<-detectCores()-6
registerDoParallel(cores = ncores) #run only once
#registerDoParallel(cores = 28) #run only once
ncores
##################################################################################
###################### Defining negative loglikelihood function ##################
##################################################################################


############### Simulating 50 realisations of zero mean Gaussain process from the specified bivariate matern covariance function ########
sim.time<-system.time(
  
  for(i in 1:100) {
    set.seed(i)
    var.bvm[,i]<-mvrnorm(n=1,mu=rep(0,times=2*(n^2)),Sigma = C.bvm)
    
    
  }
)


cmats.sim<-foreach(i=1:length(u)) %dopar% {
  
  besselJ(x=dist.mat*u[i],nu=0)
  
  
}





loglik_Model.delta.1.allcomponent <- function(theta,locs,data.at.locs,cmats) {
  # theta1: a1
  # theta2: nu1
  # theta3: sigmasq1
  # theta4: a2
  # theta5: nu2
  # theta6: sigmasq2
  # theta7: b_3
  # theta8: b_2
  # theta9: b_1
  # theta10: b0
  # theta11: b1
  # theta12: b2
  # theta13: b3
  # theta14: b4
  # theta15 :nug1
  # theta 16 :nug2
  
  
  
  a1<-theta[1]
  a2<-theta[2]
  sigmasq1<-theta[3]
  sigmasq2<-theta[4]
  b_3<-theta[5]
  b_2<-theta[6]
  b_1<-theta[7]
  b0<-theta[8]
  b1<-theta[9]
  b2<-theta[10]
  b3<-theta[11]
  b4<-theta[12]
  
  
  
  
  
  ## Hard constraints on parameters
  if(a1<=0||sigmasq1<=0||a2<=0||sigmasq2<=0||b_3<=-1||b_2<=-1||b_1<=-1||b0<=-1||b1<=-1||b2<=-1||b3<=-1||b4<=-1||b_3>=1||b_2>=1||b_1>=1||b0>=1||b1>=1||b2>=1||b3>=1||b4>=1)
  {
    return(list(nloglikelihood=1000000000))
  } 
  else {
    ##### Computing discrete frequencies at which we will evaluate the spectral density functions #####
    u<-seq(0,freq.max,length=length(uniq.dist)-1)
    f.var1<-f.matern(w=u, nu=3, sigmasq = 1, a=a1, d=2)  ####### Spectral density for variable 1
    f.var2<-f.matern(w=u, nu=3, sigmasq = 1, a=a2, d=2)  ####### Spectral density for variable 2
    delta.s<-1
    coh12<-b_3*Bspline(j=-3,k=4,delta = delta.s,x=u)+b_2*Bspline(j=-2,k=4,delta = delta.s,x=u)+b_1*Bspline(j=-1,k=4,delta = delta.s,x=u)+b0*Bspline(j=0,k=4,delta = delta.s,x=u)+b1*Bspline(j=1,k=4,delta = delta.s,x=u)+b2*Bspline(j=2,k=4,delta = delta.s,x=u)+b3*Bspline(j=3,k=4,delta = delta.s,x=u)+b4*Bspline(j=4,k=4,delta = delta.s,x=u) ### Coherence for var1 and var2
    f.var12<-coh12*sqrt(f.var1*f.var2) ####### Cross spectral density ########
    dist.mat<-rdist(locs)
    
    ######### Now we compute the covariance matrix #########
    ######### Now we compute the covariance matrix #########
    mul1<-2*pi*u*f.var1
    mul2<-2*pi*u*f.var2
    mul12<-2*pi*u*f.var12
    tempmat1<-tempmat2<-tempmat12<-matrix(0,nrow=nrow(dist.mat),ncol=ncol(dist.mat))
    for(i in 1:length(u))
    {
      tempmat1<-tempmat1+mul1[i]*cmats[[i]]
      tempmat2<-tempmat2+mul2[i]*cmats[[i]]
      tempmat12<-tempmat12+mul12[i]*cmats[[i]]
    }
    
    ######## Rescaling covariances ########
    scl1<-max(tempmat1)
    scl2<-max(tempmat2)
    COV11<-tempmat1*sigmasq1/scl1
    COV22<-tempmat2*sigmasq2/scl2
    COV12<-(sqrt(sigmasq1*sigmasq2)/sqrt(scl1*scl2))*tempmat12
    
    
    
    ##################################################################################
    ############## Creating full covariance matrix ###################################
    ##################################################################################
    C1<-cbind(COV11,COV12)
    C2<-cbind(t(COV12),COV22)
    C<-rbind(C1,C2)
    
    
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
    return(list(cov11=COV11,cov22=COV22,cov12=COV12,Coherence=coh12,FullC=C,nloglikelihood=nlog_likelihood))
  }
}



## loglik_Model.delta.1 returns only the negative loglikelihod value for training data
loglik_Model.delta.1<-function(p)
{
  temp<-loglik_Model.delta.1.allcomponent(theta=p,locs = cbind(x.s,y.s),data.at.locs = z,cmats = cmats.sim)
  return(temp$nloglikelihood)
}


########## initialized coherence ############
tr.sp<-c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2)
#tr.coh<-tr.sp[1]*Bspline(j=-3,k=4,delta = 1,x=u)+tr.sp[2]*Bspline(j=-2,k=4,delta = 1,x=u)+tr.sp[3]*Bspline(j=-1,k=4,delta = 1,x=u)+tr.sp[4]*Bspline(j=0,k=4,delta = 1,x=u)+tr.sp[5]*Bspline(j=1,k=4,delta = 1,x=u)+tr.sp[6]*Bspline(j=2,k=4,delta = 1,x=u)+tr.sp[7]*Bspline(j=3,k=4,delta = 1,x=u)+tr.sp[8]*Bspline(j=4,k=4,delta = 1,x=u)
#lines(u,tr.coh,col="blue")

init.2<-c(0.7,0.7,0.7,0.7,tr.sp)
#z<-var.bvm[,1] rm(z)
#system.time(mval<-loglik_Model.delta.1(init.2))
###### Plotting few realisations #########

par(mfrow=c(1,1))
par(cex.axis=2, cex.lab=3, cex.main=1, cex.sub=1,mar=c(2,4,1,2)+.1)

quilt.plot(x.s,y.s,var.bvm[1:n^2,1],ncol=n,nrow=n,xlab="x",ylab="y")
quilt.plot(x.s,y.s,var.bvm[(n^2+1):(2*n^2),1],ncol=n,nrow=n,xlab="x",ylab="y")

###checkpoint ###
sin(0.4)

optim.pars<-foreach(i=1:100) %dopar%
{
  library(fields)
  z<-var.bvm[,i]
  
  
  loglik_Model.delta.1.allcomponent <- function(theta,locs,data.at.locs,cmats) {
    # theta1: a1
    # theta2: nu1
    # theta3: sigmasq1
    # theta4: a2
    # theta5: nu2
    # theta6: sigmasq2
    # theta7: b_3
    # theta8: b_2
    # theta9: b_1
    # theta10: b0
    # theta11: b1
    # theta12: b2
    # theta13: b3
    # theta14: b4
    # theta15 :nug1
    # theta 16 :nug2
    
    
    
    a1<-theta[1]
    a2<-theta[2]
    sigmasq1<-theta[3]
    sigmasq2<-theta[4]
    b_3<-theta[5]
    b_2<-theta[6]
    b_1<-theta[7]
    b0<-theta[8]
    b1<-theta[9]
    b2<-theta[10]
    b3<-theta[11]
    b4<-theta[12]
    
    
    
    
    
    ## Hard constraints on parameters
    if(a1<=0||sigmasq1<=0||a2<=0||sigmasq2<=0||b_3<=-1||b_2<=-1||b_1<=-1||b0<=-1||b1<=-1||b2<=-1||b3<=-1||b4<=-1||b_3>=1||b_2>=1||b_1>=1||b0>=1||b1>=1||b2>=1||b3>=1||b4>=1)
    {
      return(list(nloglikelihood=1000000000))
    } 
    else {
      ##### Computing discrete frequencies at which we will evaluate the spectral density functions #####
      u<-seq(0,freq.max,length=length(uniq.dist)-1)
      f.var1<-f.matern(w=u, nu=3, sigmasq = 1, a=a1, d=2)  ####### Spectral density for variable 1
      f.var2<-f.matern(w=u, nu=3, sigmasq = 1, a=a2, d=2)  ####### Spectral density for variable 2
      delta.s<-1
      coh12<-b_3*Bspline(j=-3,k=4,delta = delta.s,x=u)+b_2*Bspline(j=-2,k=4,delta = delta.s,x=u)+b_1*Bspline(j=-1,k=4,delta = delta.s,x=u)+b0*Bspline(j=0,k=4,delta = delta.s,x=u)+b1*Bspline(j=1,k=4,delta = delta.s,x=u)+b2*Bspline(j=2,k=4,delta = delta.s,x=u)+b3*Bspline(j=3,k=4,delta = delta.s,x=u)+b4*Bspline(j=4,k=4,delta = delta.s,x=u) ### Coherence for var1 and var2
      f.var12<-coh12*sqrt(f.var1*f.var2) ####### Cross spectral density ########
      dist.mat<-rdist(locs)
      
      ######### Now we compute the covariance matrix #########
      ######### Now we compute the covariance matrix #########
      mul1<-2*pi*u*f.var1
      mul2<-2*pi*u*f.var2
      mul12<-2*pi*u*f.var12
      tempmat1<-tempmat2<-tempmat12<-matrix(0,nrow=nrow(dist.mat),ncol=ncol(dist.mat))
      for(i in 1:length(u))
      {
        tempmat1<-tempmat1+mul1[i]*cmats[[i]]
        tempmat2<-tempmat2+mul2[i]*cmats[[i]]
        tempmat12<-tempmat12+mul12[i]*cmats[[i]]
      }
      
      ######## Rescaling covariances ########
      scl1<-max(tempmat1)
      scl2<-max(tempmat2)
      COV11<-tempmat1*sigmasq1/scl1
      COV22<-tempmat2*sigmasq2/scl2
      COV12<-(sqrt(sigmasq1*sigmasq2)/sqrt(scl1*scl2))*tempmat12
      
      
      
      ##################################################################################
      ############## Creating full covariance matrix ###################################
      ##################################################################################
      C1<-cbind(COV11,COV12)
      C2<-cbind(t(COV12),COV22)
      C<-rbind(C1,C2)
      
      
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
      return(list(cov11=COV11,cov22=COV22,cov12=COV12,Coherence=coh12,FullC=C,nloglikelihood=nlog_likelihood))
    }
  }
  
  
  
  ## loglik_Model.delta.1 returns only the negative loglikelihod value for training data
  loglik_Model.delta.1<-function(p)
  {
    temp<-loglik_Model.delta.1.allcomponent(theta=p,locs = cbind(x.s,y.s),data.at.locs = z,cmats = cmats.sim)
    return(temp$nloglikelihood)
  }
  
  
  ########## initialized coherence ############
  tr.sp<-c(0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2)
  #tr.coh<-tr.sp[1]*Bspline(j=-3,k=4,delta = 1,x=u)+tr.sp[2]*Bspline(j=-2,k=4,delta = 1,x=u)+tr.sp[3]*Bspline(j=-1,k=4,delta = 1,x=u)+tr.sp[4]*Bspline(j=0,k=4,delta = 1,x=u)+tr.sp[5]*Bspline(j=1,k=4,delta = 1,x=u)+tr.sp[6]*Bspline(j=2,k=4,delta = 1,x=u)+tr.sp[7]*Bspline(j=3,k=4,delta = 1,x=u)+tr.sp[8]*Bspline(j=4,k=4,delta = 1,x=u)
  #lines(u,tr.coh,col="blue")
  
  init.2<-c(0.7,0.7,0.7,0.7,tr.sp)
  optim(par=init.2,fn=loglik_Model.delta.1,control=list(maxit=5000))
}

save.image("estimateddec.RData")
par(mfrow=c(1,2))
plot(u,sqrt(biwm_coh(w=u,a1=a1.bvm,a2=a2.bvm,v1=nu1.bvm,v2=nu2.bvm,a12=a12.bvm,v12 = nu12.bvm,d=2,rho = rho)),type="l",ylim=c(-1,1),main="Generating Coherence",ylab="Coh",lwd=2)


###########################################################
############ Extracting parameter values ##################
###########################################################
counti<-numeric()
for(i in 1:100){
  
  e.min.lik[i]<-optim.pars[[i]]$value
  counti[i]<-optim.pars[[i]]$counts[1]
  p<-optim.pars[[i]]$par
  e.a1[i]<-(p[1])
  e.a2[i]<-(p[2])
  e.sigmasq1[i]<-(p[3])
  e.sigmasq2[i]<-(p[4])
  e.b_3[i]<-(p[5])
  e.b_2[i]<-(p[6])
  e.b_1[i]<-(p[7])
  e.b0[i]<-(p[8])
  e.b1[i]<-(p[9])
  e.b2[i]<-(p[10])
  e.b3[i]<-(p[11])
  e.b4[i]<-(p[12])
  
  
  e.coh12[,i]<-e.b_3[i]*Bspline(j=-3,k=4,delta = 1,x=u)+e.b_2[i]*Bspline(j=-2,k=4,delta = 1,x=u)+e.b_1[i]*Bspline(j=-1,k=4,delta = 1,x=u)+e.b0[i]*Bspline(j=0,k=4,delta = 1,x=u)+e.b1[i]*Bspline(j=1,k=4,delta = 1,x=u)+e.b2[i]*Bspline(j=2,k=4,delta = 1,x=u)+e.b3[i]*Bspline(j=3,k=4,delta = 1,x=u)+e.b4[i]*Bspline(j=4,k=4,delta = 1,x=u)
  
  lines(u,e.coh12[,i],col=i+1)
}
lines(u,sqrt(biwm_coh(w=u,a1=a1.bvm,a2=a2.bvm,v1=nu1.bvm,v2=nu2.bvm,a12=a12.bvm,v12 = nu12.bvm,d=2,rho = rho)),type="l",lwd=3,ylim=c(-1,1))

plot(u,sqrt(biwm_coh(w=u,a1=a1.bvm,a2=a2.bvm,v1=nu1.bvm,v2=nu2.bvm,a12=a12.bvm,v12 = nu12.bvm,d=2,rho = rho)),type="l",ylim=c(-1,1),main="Generating Coherence",ylab="Coh",lwd=2)
avg.coh<-rowMeans(e.coh12)
sd.coh<-numeric()
for(i in 1:length(u))
{
  sd.coh[i]<-sd(e.coh12[i,]) 
}

lines(u,avg.coh,col="red",lty=2,lwd=2)
lines(u,avg.coh+3*sd.coh,col="blue",lty=3)
lines(u,avg.coh-3*sd.coh,col="blue",lty=3)
legend("bottomright",lty=c(1,2,3),col=c("black","red","blue"),lwd=c(2,2,1),c("True Coherence","Mean of Estimated Coherences","Mean+-3sigma"))



round(mean(e.a1),2)
round(mean(e.sigmasq1),2)
round(mean(e.a2),2)
round(mean(e.sigmasq2),2)

round(sd(e.a1),2)
round(sd(e.sigmasq1),2)
round(sd(e.a2),2)
round(sd(e.sigmasq2),2)



save.image("Fullanalysisdec.RData")

library(ggplot2)
plotdata<-data.frame(
  frequency=c(u,u),
  coh=c(sqrt(biwm_coh(w=u,a1=a1.bvm,a2=a2.bvm,v1=nu1.bvm,v2=nu2.bvm,a12=a12.bvm,v12 = nu12.bvm,d=2,rho = rho)),avg.coh)
  ,type=c(rep("True coherence function",times=length(avg.coh)),rep("Mean(estimated coherence functions)",times=length(avg.coh))),
  up95=c(rep(NA,times=length(avg.coh)),avg.coh+1.96*sd.coh),low95=c(rep(NA,times=length(avg.coh)),avg.coh-1.96*sd.coh))
par(mfrow=c(1,1))
#par(cex.axis=2, cex.lab=3, cex.main=1, cex.sub=1,mar=c(3,3,1,1))
p<-ggplot(data=plotdata, aes(x=frequency, y=coh, colour=type)) + geom_line(aes(linetype=type),lwd=1.0)
p<-p+geom_ribbon(aes(ymin=plotdata$low95, ymax=plotdata$up95), linetype=2, alpha=0.2)+ylim(-1, 1)
p<-p+labs(x= expression(omega),    # works fine
          y= expression(gamma(omega)))+theme(legend.title=element_blank(),legend.position = c(0.65, 0.15), legend.background = element_rect(color = "black", fill = "grey90", size = 0.1, linetype = "solid"))

plot(p)













































