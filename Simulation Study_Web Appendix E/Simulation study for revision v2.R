#####################################################################
########### Coherence project simulation study ######################
#####################################################################
#### loading libraries #######
library(fields)
library(MASS)
library(doParallel)
library(scoringRules)
library(ggplot2)
library(gridExtra)
#####################################################################
########### Creating required functions #############################
#####################################################################

##### Creating function for computing norm #####
norm_vec <- function(x) sqrt(sum(x^2))
setwd("D:/Project 2/Revision work/Final version")
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


f.matern<-function(w,nu,sigma,a,d)
{
  const<-(sigma^2)*gamma(nu+(d/2))*(a^(2*nu))/(gamma(nu)*(pi^(d/2)))
  varying<-1/((a+w^2)^(nu+(d/2)))
  return(const*varying)
}

######### Plotting Matern spectral density ########
w.plot<-seq(-8,8,length.out = 1000)
plot(w.plot,f.matern(w=w.plot,a=1,nu=1,sigma = 1,d=2),xlab=expression(omega),ylab="f")


#######################################################
############### Matern covariance function ############
#######################################################

my.matern<-function(h,a,sigma,nu)
{
  h[h==0]<-1e-10
  num1<-(sigma^2)*(2^(1-nu))/gamma(nu)
  num2<-(h*a)^nu
  num3<-besselK(x=(h*a), nu=nu)
  return(num1*num2*num3)
}

#######################################################
######## Bivariate matern coherence ###################
#######################################################

biwm_coh<-function(w,a1,a2,v1,v2,a12,v12,d,rho)
{
  temp<-numeric(length=length(w))
  for(i in 1:length(w))
  {
    num<-((gamma(v12+(d/2)))^2)*gamma(v1)*gamma(v2)*(a12^(4*v12))*((a1^2+w[i]^2)^(v1+(d/2)))*((a2^2+w[i]^2)^(v2+(d/2)))
    den<-gamma(v1+(d/2))*gamma(v2+(d/2))*(gamma(v12)^2)*(a1^(2*v1))*(a2^(2*v2))*((a12^2+w[i]^2)^(2*v12+d))
    temp[i]<-rho*sqrt(num/den)
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


freq.max<-7.5 ### (omega_t in the manuscript)

####################################################################################
####### Now we select marginal spectral densities using our proposed method ########
####################################################################################

a1<-(1/5)
nu1<-1
sigma1<-1
a2<-1/10
nu2<-3
sigma2<-1
#a12<-sqrt(0.5*((a1^2)+(a2^2)))
#nu12<-(nu1+nu2)/2
#### Plotting marginal Matern spectral densities #####
w.plot<-seq(0,freq.max,length.out = 1000)
#plot(w.plot,f.matern(w=w.plot,nu=nu1,sigma = sigma1,a=a1,d=2))
#plot(w.plot,f.matern(w=w.plot,nu=nu2,sigma = sigma1,a=a2,d=2))
#plot(w.plot,f.matern(w=w.plot,nu=nu2,sigma = sigma1,a=a2,d=2))
#plot(w.plot,f.matern(w=w.plot,nu=nu12,sigma = sigma1,a=a12,d=2))
#plot(w.plot,0.25*f.matern(w=w.plot,nu=nu12,sigma = sigma1,a=a12,d=2)/sqrt(f.matern(w=w.plot,nu=nu1,sigma = sigma1,a=a1,d=2)*f.matern(w=w.plot,nu=nu2,sigma = sigma1,a=a2,d=2)),ylab="coh")

###################################################
########## Designing coherence functions ##########
###################################################
Delta=2
b_3<-0.95
b_2<-0.95
b_1<--0.30
b0<-0.99
b1<-0.5
b2<--0.3
b3<-0.9

u<-seq(0,freq.max,length.out = 1000)
coh12<-b_3*Bspline(j=-3,k=4,delta = Delta,x=u)+b_2*Bspline(j=-2,k=4,delta = Delta,x=u)+b_1*Bspline(j=-1,k=4,delta = Delta,x=u)+b0*Bspline(j=0,k=4,delta = Delta,x=u)+b1*Bspline(j=1,k=4,delta = Delta,x=u)+b2*Bspline(j=2,k=4,delta = Delta,x=u)+b3*Bspline(j=3,k=4,delta = Delta,x=u)

plot(u,coh12,xlab=expression(omega),ylab = expression(paste(gamma,"(",omega,")")),main="True coherence function")




####### Computing multivariate covariance functions ######
u<-seq(0,freq.max,length=1000)
f.var1<-f.matern(w=u, nu=nu1, sigma = sigma1, a=a1, d=2)
f.var2<-f.matern(w=u, nu=nu2, sigma = sigma2, a=a2, d=2)
par(mfrow=c(1,2))
plot(u,f.var1,main="True Marginal Spectral density (Variable 1)")
plot(u,f.var2,main="True Marginal Spectral density (Variable 2)")

coh12<-b_3*Bspline(j=-3,k=4,delta = Delta,x=u)+b_2*Bspline(j=-2,k=4,delta = Delta,x=u)+b_1*Bspline(j=-1,k=4,delta = Delta,x=u)+b0*Bspline(j=0,k=4,delta = Delta,x=u)+b1*Bspline(j=1,k=4,delta = Delta,x=u)+b2*Bspline(j=2,k=4,delta = Delta,x=u)+b3*Bspline(j=3,k=4,delta = Delta,x=u)
par(mfrow=c(1,1))
plot(u,coh12,main="True Coherence function")
f.var12<-coh12*sqrt(f.var1*f.var2)
plot(u,f.var12,main="True Cross-spectral density")

#### Computing Block Covariance matrices ####


uniq.dist<-unique(c(dist.mat))
uniq.dist<-sort(uniq.dist)  ##### sorting distances in increasing order #####

######### Creating covariance matrix (slow but readable code) ######


C11<-numeric(length = length(uniq.dist))
for(i in 1:length(uniq.dist))
{
  dummy<-numeric(length = length(u))
  {
    for(j in 1:length(u))
    {
      theta<-uniq.dist[i]*u[j]
      J<-besselJ(x=theta,nu=0)
      dummy[j]<-2*pi*u[j]*J*f.var1[j]
    }
    C11[i]<-sum(dummy)
  }
}

####################################
####### Scaling it to 0 to 1 #######
####################################
scl.1<-max(C11)
C11<-C11/scl.1
C11<-(sigma1^2)*C11
plot(uniq.dist,C11,main = "True marginal Matern covariance (Variable 1)")
COV11<-COV22<-COV12<-matrix(NA,nrow = nrow(dist.mat),ncol=ncol(dist.mat))
for(i in 1:length(uniq.dist))
{
  COV11[which(dist.mat==uniq.dist[i],arr.ind = T)]<-C11[i]
  
}




##### C22 #####
C22<-numeric(length = length(uniq.dist))
for(i in 1:length(uniq.dist))
{
  dummy<-numeric(length = length(u))
  {
    for(j in 1:length(u))
    {
      theta<-uniq.dist[i]*u[j]
      J<-besselJ(x=theta,nu=0)
      dummy[j]<-2*pi*u[j]*J*f.var2[j]
    }
    C22[i]<-sum(dummy)
  }
}

####################################
####### Scaling it to 0 to 1 #######
####################################
scl.2<-max(C22)
C22<-C22/scl.2

###############################################
########## Mutiplying marginal variance #######
###############################################

C22<-(sigma2^2)*C22

for(i in 1:length(uniq.dist))
{
  COV22[which(dist.mat==uniq.dist[i],arr.ind = T)]<-C22[i]
  
}



###### C12 #####
C12<-numeric(length = length(uniq.dist))
for(i in 1:length(uniq.dist))
{
  dummy<-numeric(length = length(u))
  {
    for(j in 1:length(u))
    {
      theta<-uniq.dist[i]*u[j]
      J<-besselJ(x=theta,nu=0)
      dummy[j]<-2*pi*u[j]*J*f.var12[j]
    }
    C12[i]<-sum(dummy)
  }
}

#C12<-C12/sqrt(scl.1*scl.2)
#plot(uniq.dist,C12)


####################################
####### Scaling it #################
####################################


###############################################
########## Mutiplying marginal variance #######
###############################################

C12<-(sigma1*sigma2)*C12

for(i in 1:length(uniq.dist))
{
  COV12[which(dist.mat==uniq.dist[i],arr.ind = T)]<-C12[i]
  
}


myC<-rbind(cbind(COV11,COV12),cbind(t(COV12),COV22))

system.time(chol(myC))
set.seed(123)

####### Registering parallel #####
ncores<-detectCores()-4
registerDoParallel(cores = ncores)
getDoParWorkers()



###########################################################################################################################################
#################### Now we again create the full covariance matrix by optimising the code (making it fast computation) ###################
###########################################################################################################################################
theta<-outer(uniq.dist,u,"*")
system.time(Bes_mat<-besselJ(x=theta,nu=0))
index.mat<-NULL
index.val<-NULL
for(i in 1:length(uniq.dist))
{
  index.mat<-rbind(index.mat,which(dist.mat==uniq.dist[i],arr.ind = T))
  index.val<-c(index.val,rep(i,times=nrow(which(dist.mat==uniq.dist[i],arr.ind = T))))
}




full.cov.compute<-function(f1.f,f2.f,f12.f,u.f,index.mat.f,index.val.f,Bes_mat.f,uniq.dist.f,sigma1.f,sigma2.f,dmat.f)
{
  
  mult.func.f1<-function(a)
  {
    return(a*2*pi*u.f*f1.f)
  }
  mult.func.f2<-function(a)
  {
    return(a*2*pi*u.f*f2.f)
  }
  
  mult.func.f12<-function(a)
  {
    return(a*2*pi*u.f*f12.f)
  }
  C11m<-apply(Bes_mat.f,1,"mult.func.f1")
  C11n<-colSums(C11m)
  scl1<-max(C11n)
  C11n<-(C11n/scl1)*(sigma1.f^2)
  
  
  C22m<-apply(Bes_mat.f,1,"mult.func.f2")
  C22n<-colSums(C22m)
  scl2<-max(C22n)
  C22n<-(C22n/scl2)*(sigma2.f^2)
  
  C12m<-apply(Bes_mat.f,1,"mult.func.f12")
  C12n<-colSums(C12m)
  C12n<-(C12n/(sqrt(scl1*scl2)))*(sigma1.f*sigma2.f)
  dist.mat<-dmat.f
  COV11op<-COV22op<-COV12op<-matrix(NA,nrow = nrow(dist.mat),ncol=ncol(dist.mat))
  
  COV11op[index.mat.f]<-C11n[index.val.f]
  COV22op[index.mat.f]<-C22n[index.val.f]
  COV12op[index.mat.f]<-C12n[index.val.f]
  myC4<-rbind(cbind(COV11op,COV12op),cbind(t(COV12op),COV22op))
  return(myC4)
}


system.time(myC2<-full.cov.compute(f1.f=f.var1,f2.f=f.var2,f12.f=f.var12,u.f=u,index.mat.f=index.mat,index.val.f=index.val,Bes_mat.f=Bes_mat,uniq.dist.f=uniq.dist,sigma1.f=sigma1,sigma2.f=sigma2,dmat.f = dist.mat))

####### Checking if myC is equal to myC2 ####
sum(round(myC,13)==round(myC2,13)) ### numerically equivalent with computational time of 0.16 seconds! compared to 6.936 seconds
#############################################

#################### Now we simulate 50 realizations of zero mean Gaussian process with covariance matrix myC2 ################## 
set.seed(123)
mysim<-mvrnorm(n=50,mu=rep(0,times=nrow(myC)),Sigma=myC2)
###################################################
######### Plotting simulated realizations #########
###################################################
par(mfrow=c(1,2))
quilt.plot(x.s,y.s,mysim[1,1:length(x.s)],main="Variable 1",nx=30,ny=30)
quilt.plot(x.s,y.s,mysim[1,(length(x.s)+1):(2*length(x.s))],main="Variable 2",nx=30,ny=30)

par(mar=c(4,4,0.5,4))
quilt.plot(x.s,y.s,mysim[1,1:length(x.s)],nx=30,ny=30,xlab="x",ylab="y")
quilt.plot(x.s,y.s,mysim[1,(length(x.s)+1):(2*length(x.s))],xlab="x",ylab="y",nx=30,ny=30)


cor(mysim[1,1:length(x.s)],mysim[1,(length(x.s)+1):(2*length(x.s))])

#####################################################
###### Estimation of semi-parametric model ##########
#####################################################

######################################################
####### Now we estimate bivariate Matern model #######
######################################################

###checking if parallel is working or not ####
tttt<-foreach(i=1:50)%dopar%{
  i
}


######### Bivariate Matern negative log-likelihood function ##############

mle_bvm<-function(p,z,dmat.ml,index.mat.ml,index.val.ml,uniq.dist.ml)
{
  a1<-p[1]
  nu1<-p[2]
  sigma1<-p[3]
  a2<-p[4]
  nu2<-p[5]
  sigma2<-p[6]
  deltaA<-p[7]
  deltaB<-p[8]
  beta<-sin(p[9])
  if(sum(p[1:6]<0)!=0||p[7]<0||p[8]<0)
  {
    nloglikelihood<-10000000
    return(list(mlv=nloglikelihood,params=NULL))
  }
  else
  {
    
    dist.mat<-dmat.ml
    C11<-my.matern(h=uniq.dist.ml,a=a1,sigma = sigma1,nu=nu1)
    C22<-my.matern(h=uniq.dist.ml,a=a2,sigma = sigma2,nu=nu2)
    a12<-sqrt((a1^2+a2^2)/2+deltaB)
    nu12<-((nu1+nu2)/2)+deltaA
    num1<-beta*(a12^(-2*deltaA-(nu1+nu2)))*gamma(((nu1+nu2)/2)+(2/2))*gamma(nu12)
    den1<-(a1^(-deltaA-(nu1)))*(a2^(-deltaA-(nu2)))*sqrt(gamma(nu1)*gamma(nu2))*gamma(nu12+(2/2))
    rho12<-num1/den1
    C12<-sigma1*sigma2*rho12*my.matern(h=uniq.dist.ml,a=a12,sigma = 1,nu=nu12)
    
    COV11op<-COV22op<-COV12op<-matrix(NA,nrow = nrow(dist.mat),ncol=ncol(dist.mat))
    
    COV11op[index.mat.ml]<-C11[index.val.ml]
    COV22op[index.mat.ml]<-C22[index.val.ml]
    COV12op[index.mat.ml]<-C12[index.val.ml]
    
    
    
    ############## Inverting C11 ##########
    C<-rbind(cbind(COV11op,COV12op),cbind(t(COV12op),COV22op))
    Cchol<-chol(C)
    Cinv<-chol2inv(Cchol)
    logD<-determinant(C)$modulus
    
    nloglikelihood <-(0.5 * logD + 0.5 * t(z) %*% Cinv %*% z+0.5*length(z)*log(2*pi))
    if(abs(nloglikelihood) == Inf || is.nan(nloglikelihood)){ nloglikelihood <- 1e+08}
    return(list(mlv=nloglikelihood,a1=a1,a2=a2,nu1=nu1,nu2=nu2,sigma1=sigma1,sigma2=sigma2,a12=a12,nu12=nu12,rho12=rho12,full.cov=C))
  }
  
}


init.bvm<-c(1/7,1.5,1.2,1/12,2.5,1.2,0,0,0) ### setting initial values close to the true parameter values ####
system.time(mle_bvm(p=init.bvm,z=mysim[1,],dmat.ml=dist.mat,index.mat.ml=index.mat,index.val.ml=index.val,uniq.dist.ml=uniq.dist)$mlv)


######################################################
####### Now we estimate bivariate Matern model #######
######################################################
allbvm.estimates<-foreach(i=1:50)%dopar%{
  library(fields)
  
  mle_bvm_mlv<-function(pars)
  {
    return(mle_bvm(p=pars,z=mysim[i,],dmat.ml=dist.mat,index.mat.ml=index.mat,index.val.ml=index.val,uniq.dist.ml=uniq.dist)$mlv)
  }
  
  
  optim_bvm_loglik <- function(par){
    optim(par=par,
          fn = mle_bvm_mlv,
          hessian=FALSE,
          control=list(trace=6,
                       pgtol=0,
                       parscale=rep(0.1,length(par)),
                       maxit=5000))
  }
  
  init.bvm<-c(1/7,1.5,1.2,1/12,2.5,1.2,0,0,0)
  #system.time(mle_bvm_mlv(init.bvm))
  bvm.estim<-optim_bvm_loglik(init.bvm)
  #bvm.params<-mle_bvm(p=bvm.estim$par,z=mysim[i,],x=x.s,y=y.s)
  returnlist<-list(bvm.estim)
  returnlist
}
save.image("BVMn2.RData")


#######################################
######### LMC estimation ##############
#######################################
lmc.loglikelihood_allcomponents<-function(p,z,dmat.ml,index.mat.ml,index.val.ml,uniq.dist.ml){
  
  #theta1<-a1
  #theta2<-nu1
  #theta3<-sigmasq1
  #theta4<-a2
  #theta5<-nu2
  #theta6<-sigmasq2
  #theta7<-b11
  #theta8<-b12
  #theta9<-b21
  #theta10<-b22
  #theta11<-nug1
  #theta12<-nug2
  theta<-p
  a1<-theta[1]
  nu1<-theta[2]
  sigma1<-1
  a2<-theta[3]
  nu2<-theta[4]
  sigma2<-1
  b11<-theta[5]
  b12<-theta[6]
  b21<-theta[7]
  b22<-theta[8]
  
  ######## Putting hard constraints on the parameters #############
  if( a1<=0 | nu1<=0 | sigma1<=0 | a2<=0 | nu2<=0 | sigma2<=0 )
  {
    return(list(mlv=Inf))
  }
  else
  {
    
    dist.mat<-dmat.ml
    n<-nrow(dist.mat)
    C11<-my.matern(h=uniq.dist.ml,a=a1,nu=nu1,sigma = sigma1)
    C22<-my.matern(h=uniq.dist.ml,a=a2,nu=nu2,sigma = sigma2)
    
    COV11<-(b11^2)*C11+(b12^2)*C22
    COV22<-(b21^2)*C11+(b22^2)*C22
    COV12<-(b11*b21)*C11+(b12*b22)*C22
    
    
    ##################################################################################
    ############## Creating full covariance matrix ###################################
    ##################################################################################
    COV11op<-COV22op<-COV12op<-matrix(NA,nrow = nrow(dist.mat),ncol=ncol(dist.mat))
    
    COV11op[index.mat.ml]<-COV11[index.val.ml]
    COV22op[index.mat.ml]<-COV22[index.val.ml]
    COV12op[index.mat.ml]<-COV12[index.val.ml]
    
    
    
    ############## Inverting C11 ##########
    C<-rbind(cbind(COV11op,COV12op),cbind(t(COV12op),COV22op))
    
    
    
    Cchol<-chol(C)
    Cinv<-chol2inv(Cchol)
    logD<-determinant(C)$modulus
    
    nloglikelihood <-(0.5 * logD + 0.5 * t(z) %*% Cinv %*% z+0.5*length(z)*log(2*pi))
    if(abs(nloglikelihood) == Inf || is.nan(nloglikelihood)){ nloglikelihood <- 1e+08}
    return(list(mlv=nloglikelihood,full.cov=C))
  }
}

init.lmc<-c(1/7,1.5,1/12,2.5,1.2,0,0,1.2)
system.time(lmc.loglikelihood_allcomponents(p=init.lmc,z=mysim[1,],dmat.ml=dist.mat,index.mat.ml=index.mat,index.val.ml=index.val,uniq.dist.ml=uniq.dist)$mlv)
##########################################################################################
##### Now we write the code for log_likelihood of Linear model of coregionalization ######
##########################################################################################

all.lmc.estimates<-foreach(i=1:50)%dopar%{
  library(fields)
  
   lmc.loglikelihood<-function(par)
  {
    
    return(lmc.loglikelihood_allcomponents(p=par,z=mysim[i,],dmat.ml=dist.mat,index.mat.ml=index.mat,index.val.ml=index.val,uniq.dist.ml=uniq.dist)$mlv)
  }
  
  
  ## just to cross check if our lmc.loglikelihood function is without error, checking the log likelihood value for full bivariate mater parameters
  #system.time(lmc.loglikelihood(init.lmc)) ####### Should be nearly equal to fitted bivariate matern -log likelihood #######
  ###value is fine
  
  
  ### Finding mle parametrs for lmc model ####
  optim_lmc.loglikelihood <- function(par){
    optim(par=par,
          fn = lmc.loglikelihood,
          hessian=FALSE,
          control=list(trace=6,
                       pgtol=0,
                       parscale=rep(0.1,length(par)),
                       maxit=5000))
  }
  
  
  #### fit.Model.lmc is the fitted lmc model
  #### lmc.aic is the aic value for the lmc model
  
  fit.Model.lmc <- optim_lmc.loglikelihood(par=init.lmc)
  fit.Model.lmc
}

save.image("LMCn.RData")



###### Reporting ##########
###############################################################
###### Estimation of semi-parametric model (Delta=2) ##########
###############################################################
mle_sp<-function(p,z,dmat.ml,theta.ml,Bes_mat.ml,index.mat.ml,index.val.ml,u.ml,uniq.dist.ml)
{
  a1<-p[1]
  nu1<-p[2]
  sigma1<-p[3]
  a2<-p[4]
  nu2<-p[5]
  sigma2<-p[6]
  b_3<-bcoeff(p[7])
  b_2<-bcoeff(p[8])
  b_1<-bcoeff(p[9])
  b0<-bcoeff(p[10])
  b1<-bcoeff(p[11])
  b2<-bcoeff(p[12])
  b3<-bcoeff(p[13])
  if(sum(p[1:6]<0)!=0)
  {
    nloglikelihood<-10000000
    return(list(mlv=nloglikelihood,params=NULL))
  }
  else
  { 
    f.var1<-f.matern(w=u, nu=nu1, sigma = sigma1, a=a1, d=2)
    f.var2<-f.matern(w=u, nu=nu2, sigma = sigma2, a=a2, d=2)
    Delta=2
    coh12<-b_3*Bspline(j=-3,k=4,delta = Delta,x=u)+b_2*Bspline(j=-2,k=4,delta = Delta,x=u)+b_1*Bspline(j=-1,k=4,delta = Delta,x=u)+b0*Bspline(j=0,k=4,delta = Delta,x=u)+b1*Bspline(j=1,k=4,delta = Delta,x=u)+b2*Bspline(j=2,k=4,delta = Delta,x=u)+b3*Bspline(j=3,k=4,delta = Delta,x=u)
    f.var12<-coh12*sqrt(f.var1*f.var2)
    
    
    
    C<-full.cov.compute(f1.f=f.var1,f2.f=f.var2,f12.f=f.var12,u.f=u.ml,index.mat.f=index.mat.ml,index.val.f=index.val.ml,Bes_mat.f=Bes_mat.ml,uniq.dist.f=uniq.dist.ml,sigma1.f=sigma1,sigma2.f=sigma2,dmat.f = dmat.ml)
    
    
    ############## Inverting C11 ##########
    if(sum(C==Inf)>0||sum(is.nan(C))>0)
    {
      nloglikelihood <- 1e+12
    }
    else
    { 
      #checking due to numerical issues
      
      cholS<-chol(C)
      nloglikelihood <-
        -as.numeric(-0.5 * determinant(C)$modulus -
                      0.5 * t(z) %*% chol2inv(cholS) %*% z -
                      0.5 * length(z)*log(2*pi))
      
    }
    
    if (abs(nloglikelihood) == Inf || is.nan(nloglikelihood)) nloglikelihood <- 1e+08
    return(list(mlv=nloglikelihood,a1=a1,a2=a2,nu1=nu1,nu2=nu2,sigma1=sigma1,sigma2=sigma2,coh12=coh12,u=u,full.cov=C))
  }
  
}



#######################################################################
####### Now we estimate semiparametric model with Delta=2 model #######
#######################################################################
all.sp.estimates<-foreach(i=1:50)%dopar%
{
  library(fields)
    mle_sp_mlv<-function(pars)
  {
    return(mle_sp(p=pars,z=mysim[i,],dmat.ml=dist.mat,theta.ml=theta,Bes_mat.ml=Bes_mat,index.mat.ml=index.mat,index.val.ml=index.val,u.ml=u,uniq.dist.ml=uniq.dist)$mlv)
  }
  
  #system.time(mle_sp_mlv(c(1,1,1,1,1,1,1,1,1,1,1,1,1))) #935.9546
  optim_sp_loglik <- function(par){
    optim(par=par,
          fn = mle_sp_mlv,
          hessian=FALSE,
          control=list(trace=6,
                       pgtol=0,
                       parscale=rep(0.1,length(par)),
                       maxit=7000))
  }
  
  init.sp<-c(1/7,1.5,1.2,1/12,2.5,1.2,0.5,0.5,0.5,0.5,0.5,0.5,0.5)
  sp.estim<-optim_sp_loglik(init.sp)
  sp.estim
}

save.image("SPestnv3.RData")

#######################################################################################################
######################## Now we estimate semiparametric model with Delta=1 ############################
#######################################################################################################
#conv<-feval<-numeric(length=50)
#for(i in 1:50)
#{
#  conv[i]<-all.sp.estimates[[i]]$convergence
#  feval[i]<-all.sp.estimates[[i]]$counts
  
#}
#50-sum(conv)
#sum(feval==5002)

mle_sp_D1<-function(p,z,dmat.ml,theta.ml,Bes_mat.ml,index.mat.ml,index.val.ml,u.ml,uniq.dist.ml)
{
  a1<-p[1]
  nu1<-p[2]
  sigma1<-p[3]
  a2<-p[4]
  nu2<-p[5]
  sigma2<-p[6]
  b_3<-bcoeff(p[7])
  b_2<-bcoeff(p[8])
  b_1<-bcoeff(p[9])
  b0<-bcoeff(p[10])
  b1<-bcoeff(p[11])
  b2<-bcoeff(p[12])
  b3<-bcoeff(p[13])
  b4<-bcoeff(p[14])
  b5<-bcoeff(p[15])
  b6<-bcoeff(p[16])
  b7<-bcoeff(p[17])
  
  if(sum(p[1:6]<0)!=0)
  {
    nloglikelihood<-10000000
    return(list(mlv=nloglikelihood,params=NULL))
  }
  else
  { 
    f.var1<-f.matern(w=u, nu=nu1, sigma = sigma1, a=a1, d=2)
    f.var2<-f.matern(w=u, nu=nu2, sigma = sigma2, a=a2, d=2)
    Delta=1
    coh12<-b_3*Bspline(j=-3,k=4,delta = Delta,x=u)+b_2*Bspline(j=-2,k=4,delta = Delta,x=u)+b_1*Bspline(j=-1,k=4,delta = Delta,x=u)+b0*Bspline(j=0,k=4,delta = Delta,x=u)+b1*Bspline(j=1,k=4,delta = Delta,x=u)+b2*Bspline(j=2,k=4,delta = Delta,x=u)+b3*Bspline(j=3,k=4,delta = Delta,x=u)+b4*Bspline(j=4,k=4,delta = Delta,x=u)+b5*Bspline(j=5,k=4,delta = Delta,x=u)+b6*Bspline(j=6,k=4,delta = Delta,x=u)+b7*Bspline(j=7,k=4,delta = Delta,x=u)
    f.var12<-coh12*sqrt(f.var1*f.var2)
    
    
    
    C<-full.cov.compute(f1.f=f.var1,f2.f=f.var2,f12.f=f.var12,u.f=u.ml,index.mat.f=index.mat.ml,index.val.f=index.val.ml,Bes_mat.f=Bes_mat.ml,uniq.dist.f=uniq.dist.ml,sigma1.f=sigma1,sigma2.f=sigma2,dmat.f = dmat.ml)
    
    
    ############## Inverting C11 ##########
    if(sum(C==Inf)>0||sum(is.nan(C))>0)
    {
      nloglikelihood <- 1e+12
    }
    else
    { 
      #checking due to numerical issues
      
      cholS<-chol(C)
      nloglikelihood <-
        -as.numeric(-0.5 * determinant(C)$modulus -
                      0.5 * t(z) %*% chol2inv(cholS) %*% z -
                      0.5 * length(z)*log(2*pi))
      
    }
    
    if (abs(nloglikelihood) == Inf || is.nan(nloglikelihood)) nloglikelihood <- 1e+08
    return(list(mlv=nloglikelihood,a1=a1,a2=a2,nu1=nu1,nu2=nu2,sigma1=sigma1,sigma2=sigma2,coh12=coh12,u=u,full.cov=C))
  }
  
}

init.spD1<-c(1/7,1.5,1.2,1/12,2.5,1.2,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5)

#system.time(mle_spD1_mlv(c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1))) #935.9546




all.spD1.estimates<-foreach(i=1:50)%dopar%
{
  library(fields)
  mle_spD1_mlv<-function(pars)
  {
    return(mle_sp_D1(p=pars,z=mysim[i,],dmat.ml=dist.mat,theta.ml=theta,Bes_mat.ml=Bes_mat,index.mat.ml=index.mat,index.val.ml=index.val,u.ml=u,uniq.dist.ml=uniq.dist)$mlv)
  }
  
  #system.time(mle_sp_mlv(c(1,1,1,1,1,1,1,1,1,1,1,1,1))) #935.9546
  optim_spD1_loglik <- function(par){
    optim(par=par,
          fn = mle_spD1_mlv,
          hessian=FALSE,
          control=list(trace=6,
                       pgtol=0,
                       parscale=rep(0.1,length(par)),
                       maxit=7000))
  }
  
  init.spD1<-c(1/7,1.5,1.2,1/12,2.5,1.2,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5)
  sp.estimd1<-optim_spD1_loglik(init.spD1)
  sp.estimd1
  #mle_sp_D1(p=init.spD1,z=mysim[i,],dmat.ml=dist.mat,theta.ml=theta,Bes_mat.ml=Bes_mat,index.mat.ml=index.mat,index.val.ml=index.val,u.ml=u,uniq.dist.ml=uniq.dist)$mlv
  #mle_spD1_mlv(init.spD1)
}
system.time(mle_spD1_mlv(init.spD1))
save.image("SPestd1n2.RData")


################################################################################################
################# Now we estimate semiparametric model with Delta=3 ############################
################################################################################################


mle_sp_D3<-function(p,z,dmat.ml,theta.ml,Bes_mat.ml,index.mat.ml,index.val.ml,u.ml,uniq.dist.ml)
{
  a1<-p[1]
  nu1<-p[2]
  sigma1<-p[3]
  a2<-p[4]
  nu2<-p[5]
  sigma2<-p[6]
  b_3<-bcoeff(p[7])
  b_2<-bcoeff(p[8])
  b_1<-bcoeff(p[9])
  b0<-bcoeff(p[10])
  b1<-bcoeff(p[11])
  b2<-bcoeff(p[12])
  #b3<-bcoeff(p[13])
  #b4<-bcoeff(p[14])
  #b5<-bcoeff(p[15])
  #b6<-bcoeff(p[16])
  #b7<-bcoeff(p[17])
  
  if(sum(p[1:6]<0)!=0)
  {
    nloglikelihood<-10000000
    return(list(mlv=nloglikelihood,params=NULL))
  }
  else
  { 
    f.var1<-f.matern(w=u, nu=nu1, sigma = sigma1, a=a1, d=2)
    f.var2<-f.matern(w=u, nu=nu2, sigma = sigma2, a=a2, d=2)
    Delta=3
    coh12<-b_3*Bspline(j=-3,k=4,delta = Delta,x=u)+b_2*Bspline(j=-2,k=4,delta = Delta,x=u)+b_1*Bspline(j=-1,k=4,delta = Delta,x=u)+b0*Bspline(j=0,k=4,delta = Delta,x=u)+b1*Bspline(j=1,k=4,delta = Delta,x=u)+b2*Bspline(j=2,k=4,delta = Delta,x=u)#+b3*Bspline(j=3,k=4,delta = Delta,x=u)+b4*Bspline(j=4,k=4,delta = Delta,x=u)+b5*Bspline(j=5,k=4,delta = Delta,x=u)+b6*Bspline(j=6,k=4,delta = Delta,x=u)+b7*Bspline(j=7,k=4,delta = Delta,x=u)
    f.var12<-coh12*sqrt(f.var1*f.var2)
    
    
    
    C<-full.cov.compute(f1.f=f.var1,f2.f=f.var2,f12.f=f.var12,u.f=u.ml,index.mat.f=index.mat.ml,index.val.f=index.val.ml,Bes_mat.f=Bes_mat.ml,uniq.dist.f=uniq.dist.ml,sigma1.f=sigma1,sigma2.f=sigma2,dmat.f = dmat.ml)
    
    
    ############## Inverting C11 ##########
    if(sum(C==Inf)>0||sum(is.nan(C))>0)
    {
      nloglikelihood <- 1e+12
    }
    else
    { 
        
      
      cholS<-chol(C)
      nloglikelihood <-
        -as.numeric(-0.5 * determinant(C)$modulus -
                      0.5 * t(z) %*% chol2inv(cholS) %*% z -
                      0.5 * length(z)*log(2*pi))
      
    }
    
    if (abs(nloglikelihood) == Inf || is.nan(nloglikelihood)) nloglikelihood <- 1e+08
    return(list(mlv=nloglikelihood,a1=a1,a2=a2,nu1=nu1,nu2=nu2,sigma1=sigma1,sigma2=sigma2,coh12=coh12,u=u,full.cov=C))
  }
  
}

init.spD3<-c(1/7,1.5,1.2,1/12,2.5,1.2,0.5,0.5,0.5,0.5,0.5,0.5)

#system.time(mle_spD3_mlv(c(1,1,1,1,1,1,1,1,1,1,1,1))) #935.9546



all.spD3.estimates<-foreach(i=1:50)%dopar%
{
  library(fields)
  mle_spD3_mlv<-function(pars)
  {
    return(mle_sp_D3(p=pars,z=mysim[i,],dmat.ml=dist.mat,theta.ml=theta,Bes_mat.ml=Bes_mat,index.mat.ml=index.mat,index.val.ml=index.val,u.ml=u,uniq.dist.ml=uniq.dist)$mlv)
  }
  
  #system.time(mle_sp_mlv(c(1,1,1,1,1,1,1,1,1,1,1,1,1))) #935.9546
  optim_spD3_loglik <- function(par){
    optim(par=par,
          fn = mle_spD3_mlv,
          hessian=FALSE,
          control=list(trace=6,
                       pgtol=0,
                       parscale=rep(0.1,length(par)),
                       maxit=7000))
  }
  
  init.spD3<-c(1/7,1.5,1.2,1/12,2.5,1.2,0.5,0.5,0.5,0.5,0.5,0.5)
  sp.estimd3<-optim_spD3_loglik(init.spD3)
  sp.estimd3
  #mle_sp_D1(p=init.spD1,z=mysim[i,],dmat.ml=dist.mat,theta.ml=theta,Bes_mat.ml=Bes_mat,index.mat.ml=index.mat,index.val.ml=index.val,u.ml=u,uniq.dist.ml=uniq.dist)$mlv
  #mle_spD1_mlv(init.spD1)
}
system.time(mle_spD3_mlv(init.spD3))
save.image("SPestd3n2.RData")


################################################################################################
################# Now we estimate semiparametric model with Delta=4 ############################
################################################################################################


mle_sp_D4<-function(p,z,dmat.ml,theta.ml,Bes_mat.ml,index.mat.ml,index.val.ml,u.ml,uniq.dist.ml)
{
  a1<-p[1]
  nu1<-p[2]
  sigma1<-p[3]
  a2<-p[4]
  nu2<-p[5]
  sigma2<-p[6]
  b_3<-bcoeff(p[7])
  b_2<-bcoeff(p[8])
  b_1<-bcoeff(p[9])
  b0<-bcoeff(p[10])
  b1<-bcoeff(p[11])
  #b2<-bcoeff(p[12])
  #b3<-bcoeff(p[13])
  #b4<-bcoeff(p[14])
  #b5<-bcoeff(p[15])
  #b6<-bcoeff(p[16])
  #b7<-bcoeff(p[17])
  
  if(sum(p[1:6]<0)!=0)
  {
    nloglikelihood<-10000000
    return(list(mlv=nloglikelihood,params=NULL))
  }
  else
  { 
    f.var1<-f.matern(w=u, nu=nu1, sigma = sigma1, a=a1, d=2)
    f.var2<-f.matern(w=u, nu=nu2, sigma = sigma2, a=a2, d=2)
    Delta=4
    coh12<-b_3*Bspline(j=-3,k=4,delta = Delta,x=u)+b_2*Bspline(j=-2,k=4,delta = Delta,x=u)+b_1*Bspline(j=-1,k=4,delta = Delta,x=u)+b0*Bspline(j=0,k=4,delta = Delta,x=u)+b1*Bspline(j=1,k=4,delta = Delta,x=u)#+b2*Bspline(j=2,k=4,delta = Delta,x=u)#+b3*Bspline(j=3,k=4,delta = Delta,x=u)+b4*Bspline(j=4,k=4,delta = Delta,x=u)+b5*Bspline(j=5,k=4,delta = Delta,x=u)+b6*Bspline(j=6,k=4,delta = Delta,x=u)+b7*Bspline(j=7,k=4,delta = Delta,x=u)
    f.var12<-coh12*sqrt(f.var1*f.var2)
    
    
    
    C<-full.cov.compute(f1.f=f.var1,f2.f=f.var2,f12.f=f.var12,u.f=u.ml,index.mat.f=index.mat.ml,index.val.f=index.val.ml,Bes_mat.f=Bes_mat.ml,uniq.dist.f=uniq.dist.ml,sigma1.f=sigma1,sigma2.f=sigma2,dmat.f = dmat.ml)
    
    
    ############## Inverting C11 ##########
    if(sum(C==Inf)>0||sum(is.nan(C))>0)
    {
      nloglikelihood <- 1e+12
    }
    else
    { 
      
      
      cholS<-chol(C)
      nloglikelihood <-
        -as.numeric(-0.5 * determinant(C)$modulus -
                      0.5 * t(z) %*% chol2inv(cholS) %*% z -
                      0.5 * length(z)*log(2*pi))
      
    }
    
    if (abs(nloglikelihood) == Inf || is.nan(nloglikelihood)) nloglikelihood <- 1e+08
    return(list(mlv=nloglikelihood,a1=a1,a2=a2,nu1=nu1,nu2=nu2,sigma1=sigma1,sigma2=sigma2,coh12=coh12,u=u,full.cov=C))
  }
  
}

init.spD4<-c(1/7,1.5,1.2,1/12,2.5,1.2,0.5,0.5,0.5,0.5,0.5)

#system.time(mle_spD4_mlv(c(1,1,1,1,1,1,1,1,1,1,1,1))) #935.9546


all.spD4.estimates<-foreach(i=1:50)%dopar%
{
  library(fields)
  mle_spD4_mlv<-function(pars)
  {
    return(mle_sp_D4(p=pars,z=mysim[i,],dmat.ml=dist.mat,theta.ml=theta,Bes_mat.ml=Bes_mat,index.mat.ml=index.mat,index.val.ml=index.val,u.ml=u,uniq.dist.ml=uniq.dist)$mlv)
  }
  
  #system.time(mle_sp_mlv(c(1,1,1,1,1,1,1,1,1,1,1,1,1))) #935.9546
  optim_spD4_loglik <- function(par){
    optim(par=par,
          fn = mle_spD4_mlv,
          hessian=FALSE,
          control=list(trace=6,
                       pgtol=0,
                       parscale=rep(0.1,length(par)),
                       maxit=7000))
  }
  
  init.spD4<-c(1/7,1.5,1.2,1/12,2.5,1.2,0.5,0.5,0.5,0.5,0.5)
  sp.estimd4<-optim_spD4_loglik(init.spD4)
  sp.estimd4
  #mle_sp_D1(p=init.spD1,z=mysim[i,],dmat.ml=dist.mat,theta.ml=theta,Bes_mat.ml=Bes_mat,index.mat.ml=index.mat,index.val.ml=index.val,u.ml=u,uniq.dist.ml=uniq.dist)$mlv
  #mle_spD1_mlv(init.spD1)
}
system.time(mle_spD4_mlv(init.spD4))
save.image("SPestd4n2.RData")



#Sys.sleep(3600)

###############################################################################
######## Now plotting the estimated and true coherence functions ##############
###############################################################################


######### Extracting BVM Parameters ##############
all.bvm.params<-foreach(i=1:50)%dopar%
{
  mle_bvm(p=allbvm.estimates[[i]][[1]]$par,z=mysim[i,],dmat.ml=dist.mat,index.mat.ml=index.mat,index.val.ml=index.val,uniq.dist.ml=uniq.dist)
}

#### Computing the coherence functions and mean coherence function from the estimated Bivariate Matern models ########
bvm.coh.mat<-matrix(NA,nrow = 50,ncol = 1000)
for(i in 1:50)
{
  tempcoh<-(biwm_coh(w=u,a1=all.bvm.params[[i]]$a1,a2=all.bvm.params[[i]]$a2,v1=all.bvm.params[[i]]$nu1,v2=all.bvm.params[[i]]$nu2,a12=all.bvm.params[[i]]$a12,v12 = all.bvm.params[[i]]$nu12,d=2,rho=all.bvm.params[[i]]$rho12 ))
  bvm.coh.mat[i,]<-tempcoh
}

#par(mfrow=c(1,1))
#plot(u,coh12,xlab=expression(omega),ylab = expression(paste(gamma,"(",omega,")")),main="Coherence functions",type = "l",lwd=2,ylim=c(0,1))
#for(i in 1:50)
#{
#  lines(u,bvm.coh.mat[i,],col=i+1)
#}

bvm.coh.mean<-colMeans(bvm.coh.mat)

#lines(u,bvm.coh.mean,col=2,lwd=2)

#legend("bottomright",lty=c(1,1,1,1),col=c(1,2,3,4),lwd=c(2,2,2,2),c("True coherence","BVM coherence","LMC coherence","SP coherence"))




##### Extracting LMC parameters and computing coherence functions and mean coherence functions from the estimate LMC model ####

lmc.coh<-function(w,a1,nu1,sigma1,a2,nu2,sigma2,b11,b12,b21,b22)
{
  num<-b11*b21*f.matern(w=w,nu=nu1,a=a1,sigma = sigma1,d=2)+b12*b22*f.matern(w=w,nu=nu2,sigma = sigma2,a=a2,d=2)
  den<-sqrt((b11^2)*f.matern(w=w,nu=nu1,a=a1,sigma = sigma1,d=2)+(b12^2)*f.matern(w=w,nu=nu2,sigma = sigma2,a=a2,d=2))*sqrt((b21^2)*f.matern(w=w,nu=nu1,a=a1,sigma = sigma1,d=2)+(b22^2)*f.matern(w=w,nu=nu2,sigma = sigma2,a=a2,d=2))
  return(num/den)
}



lmc.param.comp<-function(p)
{
  theta<-p
  a1<-theta[1]
  nu1<-theta[2]
  sigma1<-1
  a2<-theta[3]
  nu2<-theta[4]
  sigma2<-1
  b11<-theta[5]
  b12<-theta[6]
  b21<-theta[7]
  b22<-theta[8]
  
  return(list(a1=a1,nu1=nu1,sigma1=sigma1,a2=a2,nu2=nu2,sigma2=sigma2,b11=b11,b12=b12,b21=b21,b22=b22))
}

#lmc.params<-lmc.param.comp(p=fit.Model.lmc$par)
#lmc.params

all.lmc.params<-foreach(i=1:50)%dopar%
{
  lmc.param.comp(p=all.lmc.estimates[[i]]$par)
}

lmc.coh.mat<-matrix(NA,nrow = 50,ncol = 1000)
for(i in 1:50)
{
  tempcoh<-lmc.coh(w=u,a1=all.lmc.params[[i]]$a1,nu1=all.lmc.params[[i]]$nu1,sigma1=all.lmc.params[[i]]$sigma1,a2=all.lmc.params[[i]]$a2,nu2=all.lmc.params[[i]]$nu2,sigma2=all.lmc.params[[i]]$sigma2,b11=all.lmc.params[[i]]$b11,b12=all.lmc.params[[i]]$b12,b21=all.lmc.params[[i]]$b21,b22=all.lmc.params[[i]]$b22)
  lmc.coh.mat[i,]<-tempcoh
}

lmc.coh.mean<-colMeans(lmc.coh.mat)

#for(i in 1:50)
#{
#  lines(u,lmc.coh.mat[i,],col=2*i+1)
#}

#lines(u,lmc.coh.mean,col="green",lwd=2)



############## Computing coherence functions and mean coherence functions from the estimated Semi-parametric model with different delta #############


##### Delta= 2####

sp.coh.mat<-matrix(NA,nrow = 50,ncol = 1000)

sp.coh.list<-foreach(i=1:50)%dopar%{
  
  mle_sp(p=all.sp.estimates[[i]]$par,z=mysim[i,],dmat.ml=dist.mat,theta.ml=theta,Bes_mat.ml=Bes_mat,index.mat.ml=index.mat,index.val.ml=index.val,u.ml=u,uniq.dist.ml=uniq.dist)$coh12
}
for(i in 1:50)
{
  sp.coh.mat[i,]<-sp.coh.list[[i]]
}

sp.coh.mean<-colMeans(sp.coh.mat)



####### Delta =1 ######

spd1.coh.mat<-matrix(NA,nrow = 50,ncol = 1000)


spd1.coh.list<-foreach(i=1:50)%dopar%{
  
  mle_sp_D1(p=all.spD1.estimates[[i]]$par,z=mysim[i,],dmat.ml=dist.mat,theta.ml=theta,Bes_mat.ml=Bes_mat,index.mat.ml=index.mat,index.val.ml=index.val,u.ml=u,uniq.dist.ml=uniq.dist)$coh12
}
for(i in 1:50)
{
  spd1.coh.mat[i,]<-spd1.coh.list[[i]]
}

spd1.coh.mean<-colMeans(spd1.coh.mat)


######## Delta = 3 ########

spd3.coh.mat<-matrix(NA,nrow = 50,ncol = 1000)


sp.cohd3.list<-foreach(i=1:50)%dopar%{
  
  mle_sp_D3(p=all.spD3.estimates[[i]]$par,z=mysim[i,],dmat.ml=dist.mat,theta.ml=theta,Bes_mat.ml=Bes_mat,index.mat.ml=index.mat,index.val.ml=index.val,u.ml=u,uniq.dist.ml=uniq.dist)$coh12
}
for(i in 1:50)
{
  spd3.coh.mat[i,]<-sp.cohd3.list[[i]]
}


spd3.coh.mean<-colMeans(spd3.coh.mat)

######## Delta = 4 ##########

spd4.coh.mat<-matrix(NA,nrow = 50,ncol = 1000)
sp.cohd4.list<-foreach(i=1:50)%dopar%{
  
  mle_sp_D4(p=all.spD4.estimates[[i]]$par,z=mysim[i,],dmat.ml=dist.mat,theta.ml=theta,Bes_mat.ml=Bes_mat,index.mat.ml=index.mat,index.val.ml=index.val,u.ml=u,uniq.dist.ml=uniq.dist)$coh12
}
for(i in 1:50)
{
  spd4.coh.mat[i,]<-sp.cohd4.list[[i]]
}

spd4.coh.mean<-colMeans(spd4.coh.mat)




##############################################################################################
############################# Now Plotting the coherence functions properly ##################
##############################################################################################
par(mfrow=c(2,3))
plot(u,coh12,xlab=expression(omega),ylab = expression(paste(gamma,"(",omega,")")),main="Coherence functions (BVM)",type = "l",lwd=2,ylim=c(-1,1),col="blue")
for(i in 1:50)
{
  lines(u,bvm.coh.mat[i,],col="grey")
}
lines(u,coh12,lwd=2,col="blue")
lines(u,bvm.coh.mean,col="red",lwd=2)


plot(u,coh12,xlab=expression(omega),ylab = expression(paste(gamma,"(",omega,")")),main="Coherence functions (LMC)",type = "l",lwd=2,ylim=c(-1,1),col="blue")
for(i in 1:50)
{
  lines(u,lmc.coh.mat[i,],col="grey")
}
lines(u,coh12,lwd=2,col="blue")
lines(u,lmc.coh.mean,col="red",lwd=2)


plot(u,coh12,xlab=expression(omega),ylab = expression(paste(gamma,"(",omega,")")),main="Coherence functions (SPD2)",type = "l",lwd=2,ylim=c(-1,1),col="blue")
for(i in 1:50)
{
  lines(u,sp.coh.mat[i,],col="grey")
}
lines(u,coh12,lwd=2,col="blue")
lines(u,sp.coh.mean,col="red",lwd=2)


plot(u,coh12,xlab=expression(omega),ylab = expression(paste(gamma,"(",omega,")")),main="Coherence functions (SPD1)",type = "l",lwd=2,ylim=c(-1,1),col="blue")
for(i in 1:50)
{
  lines(u,spd1.coh.mat[i,],col="grey")
}
lines(u,coh12,lwd=2,col="blue")
lines(u,spd1.coh.mean,col="red",lwd=2)



plot(u,coh12,xlab=expression(omega),ylab = expression(paste(gamma,"(",omega,")")),main="Coherence functions (SPD3)",type = "l",lwd=2,ylim=c(-1,1),col="blue")
for(i in 1:50)
{
  lines(u,spd3.coh.mat[i,],col="grey")
}
lines(u,coh12,lwd=2,col="blue")
lines(u,spd3.coh.mean,col="red",lwd=2)



plot(u,coh12,xlab=expression(omega),ylab = expression(paste(gamma,"(",omega,")")),main="Coherence functions (SPD4)",type = "l",lwd=2,ylim=c(-1,1),col="blue")
for(i in 1:50)
{
  lines(u,spd4.coh.mat[i,],col="grey")
}
lines(u,coh12,lwd=2,col="blue")
lines(u,spd4.coh.mean,col="red",lwd=2)


#for(i in 1:50)
#{
#  lines(u,lmc.coh.mat[i,],col="deepskyblue")
#}

#lines(u,lmc.coh.mean,col="deepskyblue4",lwd=3)


#for(i in 1:50)
#{
#  lines(u,sp.coh.mat[i,],col="chartreuse")
#}

#lines(u,sp.coh.mean,col="chartreuse4",lwd=3)


#for(i in 1:50)
#{
#  lines(u,spd1.coh.mat[i,],col="gray71")
#}

#lines(u,spd1.coh.mean,col="gray28",lwd=3)



#for(i in 1:50)
#{
#  lines(u,spd3.coh.mat[i,],col="lightcyan")
#}

#lines(u,spd3.coh.mean,col="lightcyan4",lwd=3)


#for(i in 1:50)
#{
#  lines(u,spd4.coh.mat[i,],col="lightsalmon")
#}

#lines(u,spd4.coh.mean,col="lightsalmon4",lwd=3)



#lines(u,coh12,lwd=3)

######## First we compute pointwise standard deviation for all the estimated coherence functions ######

bvm.sd<-apply(bvm.coh.mat,2, "sd")
lmc.sd<-apply(lmc.coh.mat,2, "sd")
sp.sd<-apply(sp.coh.mat,2, "sd")
spd1.sd<-apply(spd1.coh.mat,2, "sd")
spd3.sd<-apply(spd3.coh.mat,2, "sd")
spd4.sd<-apply(spd4.coh.mat,2, "sd")


true.coh.d<-data.frame(u=u,coh=coh12,c.min=NA,
                     c.max=NA, Model="True coherence function")
ggplot(data=true.coh.d, aes(x=u, y=coh, colour=Model)) + geom_line(lwd=1.0)+ylim(-1, 1)+
  labs(x= expression(omega),    # works fine
       y= expression(gamma(omega)))+theme(legend.title=element_blank(),legend.position = c(0.3, 0.15), legend.background = element_rect(color = "black", fill = "grey90", size = 0.1, linetype = "solid"))



l.limit.spd2<-sp.coh.mean-1.96*sp.sd
u.limit.spd2<-sp.coh.mean+1.96*sp.sd
u.limit.spd2[u.limit.spd2>1]=1

pl.spd2.data<-data.frame(u=u,coh=sp.coh.mean,
                        c.min=l.limit.spd2,
                        c.max=u.limit.spd2,Model="Mean(estimated coherence functions)")

pl.spd2.data<-rbind(pl.spd2.data,true.coh.d)

ggplot(data=pl.spd2.data, aes(x=u, y=coh, colour=Model)) + geom_line(lwd=1.0)+
  geom_ribbon(aes(ymin=c.min, ymax=c.max), linetype=2, alpha=0.2)+ylim(-1, 1)+
  labs(x= expression(omega),    # works fine
          y= expression(gamma(omega)),title=expression(paste("Semiparametric (",Delta,"=2)")))+theme(legend.title=element_blank(),legend.position = c(0.35, 0.15), legend.background = element_rect(color = "black", fill = "grey90", size = 0.1, linetype = "solid"))




l.limit.spd1<-spd1.coh.mean-1.96*spd1.sd
u.limit.spd1<-spd1.coh.mean+1.96*spd1.sd
u.limit.spd1[u.limit.spd1>1]=1

pl.spd1.data<-data.frame(u=u,coh=spd1.coh.mean,
                         c.min=l.limit.spd1,
                         c.max=u.limit.spd1,Model="Mean(estimated coherence functions)")

pl.spd1.data<-rbind(pl.spd1.data,true.coh.d)

ggplot(data=pl.spd1.data, aes(x=u, y=coh, colour=Model)) + geom_line(lwd=1.0)+
  geom_ribbon(aes(ymin=c.min, ymax=c.max), linetype=2, alpha=0.2)+ylim(-1, 1)+
  labs(x= expression(omega),    # works fine
       y= expression(gamma(omega)),title=expression(paste("Semiparametric (",Delta,"=1)")))+theme(legend.title=element_blank(),legend.position = c(0.35, 0.15), legend.background = element_rect(color = "black", fill = "grey90", size = 0.1, linetype = "solid"))



l.limit.spd3<-spd3.coh.mean-1.96*spd3.sd
u.limit.spd3<-spd3.coh.mean+1.96*spd3.sd
u.limit.spd3[u.limit.spd3>1]=1

pl.spd3.data<-data.frame(u=u,coh=spd3.coh.mean,
                         c.min=l.limit.spd3,
                         c.max=u.limit.spd3,Model="Mean(estimated coherence functions)")

pl.spd3.data<-rbind(pl.spd3.data,true.coh.d)

ggplot(data=pl.spd3.data, aes(x=u, y=coh, colour=Model)) + geom_line(lwd=1.0)+
  geom_ribbon(aes(ymin=c.min, ymax=c.max), linetype=2, alpha=0.2)+ylim(-1, 1)+
  labs(x= expression(omega),    # works fine
       y= expression(gamma(omega)),title=expression(paste("Semiparametric (",Delta,"=3)")))+theme(legend.title=element_blank(),legend.position = c(0.35, 0.15), legend.background = element_rect(color = "black", fill = "grey90", size = 0.1, linetype = "solid"))


l.limit.spd4<-spd4.coh.mean-1.96*spd4.sd
u.limit.spd4<-spd4.coh.mean+1.96*spd4.sd
u.limit.spd4[u.limit.spd4>1]=1

pl.spd4.data<-data.frame(u=u,coh=spd4.coh.mean,
                         c.min=l.limit.spd4,
                         c.max=u.limit.spd4,Model="Mean(estimated coherence functions)")

pl.spd4.data<-rbind(pl.spd4.data,true.coh.d)

ggplot(data=pl.spd4.data, aes(x=u, y=coh, colour=Model)) + geom_line(lwd=1.0)+
  geom_ribbon(aes(ymin=c.min, ymax=c.max), linetype=2, alpha=0.2)+ylim(-1, 1)+
  labs(x= expression(omega),    # works fine
       y= expression(gamma(omega)),title=expression(paste("Semiparametric (",Delta,"=4)")))+theme(legend.title=element_blank(),legend.position = c(0.35, 0.15), legend.background = element_rect(color = "black", fill = "grey90", size = 0.1, linetype = "solid"))

l.limit.lmc<-lmc.coh.mean-1.96*lmc.sd
u.limit.lmc<-lmc.coh.mean+1.96*lmc.sd
u.limit.lmc[u.limit.lmc>1]=1

pl.lmc.data<-data.frame(u=u,coh=lmc.coh.mean,
                         c.min=l.limit.lmc,
                         c.max=u.limit.lmc,Model="Mean(estimated coherence functions)")

pl.lmc.data<-rbind(pl.lmc.data,true.coh.d)

ggplot(data=pl.lmc.data, aes(x=u, y=coh, colour=Model)) + geom_line(lwd=1.0)+
  geom_ribbon(aes(ymin=c.min, ymax=c.max), linetype=2, alpha=0.2)+ylim(-1, 1)+
  labs(x= expression(omega),    # works fine
       y= expression(gamma(omega)),title=expression(paste("LMC")))+theme(legend.title=element_blank(),legend.position = c(0.35, 0.15), legend.background = element_rect(color = "black", fill = "grey90", size = 0.1, linetype = "solid"))



l.limit.bvm<-bvm.coh.mean-1.96*bvm.sd
u.limit.bvm<-bvm.coh.mean+1.96*bvm.sd
u.limit.bvm[u.limit.bvm>1]=1

pl.bvm.data<-data.frame(u=u,coh=bvm.coh.mean,
                        c.min=l.limit.bvm,
                        c.max=u.limit.bvm,Model="Mean(estimated coherence functions)")

pl.bvm.data<-rbind(pl.bvm.data,true.coh.d)

ggplot(data=pl.bvm.data, aes(x=u, y=coh, colour=Model)) + geom_line(lwd=1.0)+
  geom_ribbon(aes(ymin=c.min, ymax=c.max), linetype=2, alpha=0.2)+ylim(-1, 1)+
  labs(x= expression(omega),    # works fine
       y= expression(gamma(omega)),title=expression(paste("Full bivariate Matérn")))+theme(legend.title=element_blank(),legend.position = c(0.35, 0.15), legend.background = element_rect(color = "black", fill = "grey90", size = 0.1, linetype = "solid"))







plot.data<-data.frame(x=rep(u,times=7),y=c(coh12,sp.coh.mean,spd1.coh.mean,spd3.coh.mean,spd4.coh.mean,lmc.coh.mean,bvm.coh.mean),sd=c(rep(NA,times=1000),sp.sd,spd1.sd,spd3.sd,spd4.sd,lmc.sd,bvm.sd),model=c(rep("True",times=1000),rep("SPD2",times=1000),rep("SPD1",times=1000),rep("SPD3",times=1000),rep("SPD4",times=1000),rep("LMC",times=1000),rep("BVM",times=1000)))

plot.d1<-data.frame(x=rep(u,times=2),y=c(coh12,sp.coh.mean),sd=c(rep(NA,times=1000),sp.sd),Model=c(rep("True",times=1000),rep("SPD2",times=1000)))
plot.d2<-data.frame(x=rep(u,times=2),y=c(coh12,spd1.coh.mean),sd=c(rep(NA,times=1000),spd1.sd),Model=c(rep("True",times=1000),rep("SPD1",times=1000)))
plot.d3<-data.frame(x=rep(u,times=2),y=c(coh12,spd3.coh.mean),sd=c(rep(NA,times=1000),spd3.sd),Model=c(rep("True",times=1000),rep("SPD3",times=1000)))
plot.d4<-data.frame(x=rep(u,times=2),y=c(coh12,spd4.coh.mean),sd=c(rep(NA,times=1000),spd4.sd),Model=c(rep("True",times=1000),rep("SPD4",times=1000)))
plot.d5<-data.frame(x=rep(u,times=2),y=c(coh12,lmc.coh.mean),sd=c(rep(NA,times=1000),lmc.sd),Model=c(rep("True",times=1000),rep("LMC",times=1000)))
plot.d6<-data.frame(x=rep(u,times=2),y=c(coh12,bvm.coh.mean),sd=c(rep(NA,times=1000),bvm.sd),Model=c(rep("True",times=1000),rep("BVM",times=1000)))

#myp1 <- ggplot(plot.d1,aes(x=x,y=y,col=Model))+geom_line() + ylim(-1, 1)+
 # geom_line(lwd=1.5) +labs(x= expression(omega),    # works fine
  #                         y= expression(gamma(omega)))+
  #geom_ribbon(aes(ymin=y-1.96*sd,ymax=y+1.96*sd),alpha=0.1,fill=c(rep(1,times=1000),rep(2,times=1000)))# +


#myp2 <- ggplot(plot.d2,aes(x=x,y=y,col=Model))+geom_line() + ylim(-1, 1)+ 
 # geom_line(lwd=1.5) +labs(x= expression(omega),    # works fine
  #                         y= expression(gamma(omega)))+
  #geom_ribbon(aes(ymin=y-1.96*sd,ymax=y+1.96*sd),alpha=0.1,fill=c(rep(1,times=1000),rep(2,times=1000)))# +

#myp3 <- ggplot(plot.d3,aes(x=x,y=y,col=Model))+geom_line() + ylim(-1, 1)+ 
 # geom_line(lwd=1.5) +labs(x= expression(omega),    # works fine
  #                         y= expression(gamma(omega)))+
  #geom_ribbon(aes(ymin=y-1.96*sd,ymax=y+1.96*sd),alpha=0.1,fill=c(rep(1,times=1000),rep(2,times=1000)))# +

#myp4 <- ggplot(plot.d4,aes(x=x,y=y,col=Model))+geom_line() + ylim(-1, 1)+
 # geom_line(lwd=1.5) +labs(x= expression(omega),    # works fine
  #                         y= expression(gamma(omega)))+
  #geom_ribbon(aes(ymin=y-1.96*sd,ymax=y+1.96*sd),alpha=0.1,fill=c(rep(1,times=1000),rep(2,times=1000)))# +

#myp5 <- ggplot(plot.d5,aes(x=x,y=y,col=Model))+geom_line()+ ylim(-1, 1)+
 # geom_line(lwd=1.5) +labs(x= expression(omega),    # works fine
  #                         y= expression(gamma(omega)))+
  #geom_ribbon(aes(ymin=y-1.96*sd,ymax=y+1.96*sd),alpha=0.1,fill=c(rep(1,times=1000),rep(2,times=1000)))# +

#myp6 <- ggplot(plot.d6,aes(x=x,y=y,col=Model))+geom_line()+ ylim(-1, 1)+
 # geom_line(lwd=1.5) +labs(x= expression(omega),    # works fine
  #                         y= expression(gamma(omega)))+
  #geom_ribbon(aes(ymin=y-1.96*sd,ymax=y+1.96*sd),alpha=0.1,fill=c(rep(1,times=1000),rep(2,times=1000)))# +

#grid.arrange(myp1,myp2,myp3,myp4,myp5, myp6, ncol=3, nrow = 2)
################################################################################################
############## Now we compare the prediction performance for all the candidate models ##########
############ Prediction performance ############################################################

library(scoringRules)
pred_performance<-function(C,sample)
{
  #### Using Sherman Morrison updating formula from the article
  ### Relationship between the Inverses of a Matrix and a Submatrix by
  ####E. Ju?rez-Ruiz1
  ####, R. Cort?s-Maldonado2
  ####, F. P?rez-Rodr?guez2
  
  A<-C
  Ainv<-solve(C)
  s.size<-length(sample)
  pred<-pred.var<-numeric(length=s.size)
  for(i in 1:s.size)
  {
    B<-A[-i,-i]
    a<-A[,i]
    e<-rep(0,times=nrow(A))
    e[i]<-1
    u<-a-e
    v<-e
    U<-matrix(u,ncol = 1)
    V<-matrix(v,ncol=1)
    subt<-U%*%t(V)
    num<-(Ainv%*%U)%*%(t(V)%*%Ainv)
    den<-(1-t(V)%*%Ainv%*%U)
    SM_Mat<-Ainv+(c(1/den))*(num)
    
    sigma22inv<-SM_Mat[-i,-i]
    #sigma22<-C[-i,-i]
    z<-matrix(sample[-i],nrow=s.size-1,ncol=1)
    sigma12<-matrix(c(C[i,-i]),nrow=1)
    #sigma22inv<-solve(sigma22)
    weights<-sigma12%*%sigma22inv
    pred[i]<-weights%*%z
    pred.var[i]<-C[i,i]-weights%*%t(sigma12)
    
  }
  
  MSPE<-mean((sample-pred)^2)
  Log_Score<-mean(log(2*pi*pred.var)+(((sample-pred)^2)/pred.var))
  MAE<-mean(abs(sample-pred))
  RMSE<-sqrt(mean((sample-pred)^2))
  NMSE<-1-(sum((sample-pred)^2))/(sum((sample-mean(sample))^2)) #Normalized mean squared error
  CRPS<-mean(crps(sample, "norm", mean = c(pred), sd = c(sqrt(pred.var)))) #Mean CRPS
  
  return(list(kriged.values=pred,kriged.variance=pred.var,Observed=sample,MSPE=MSPE,Log_Score=Log_Score,MAE=MAE,RMSE=RMSE,CRPS=CRPS,NMSE=NMSE))  
}


######  Bivariate Matern ######

bvm.preds<-foreach(i=1:50)%dopar%
{library(scoringRules)
  C_i<-mle_bvm(p=allbvm.estimates[[i]][[1]]$par,z=mysim[i,],dmat.ml=dist.mat,index.mat.ml=index.mat,index.val.ml=index.val,uniq.dist.ml=uniq.dist)$full.cov
  pred_performance(C=C_i,sample =mysim[i,] )
}
save.image("bvmpredsn3.RData")

####### LMC #######

lmc.preds<-foreach(i=1:50)%dopar%
{library(scoringRules)
  C_i<-lmc.loglikelihood_allcomponents(p=all.lmc.estimates[[i]]$par,z=mysim[i,],dmat.ml=dist.mat,index.mat.ml=index.mat,index.val.ml=index.val,uniq.dist.ml=uniq.dist)$full.cov
  pred_performance(C=C_i,sample =mysim[i,] )
}

save.image("lmcpredsn3.RData")

###### Semiparametric Delta=2 ######

sp.preds<-foreach(i=1:50)%dopar%
{library(scoringRules)
  C_i<-mle_sp(p=all.sp.estimates[[i]]$par,z=mysim[i,],dmat.ml=dist.mat,theta.ml=theta,Bes_mat.ml=Bes_mat,index.mat.ml=index.mat,index.val.ml=index.val,u.ml=u,uniq.dist.ml=uniq.dist)$full.cov
  pred_performance(C=C_i,sample =mysim[i,] )
}
save.image("sppredsn3.RData")

###### Semiparametric Delta=1 ######

spd1.preds<-foreach(i=1:50)%dopar%
{library(scoringRules)
  C_i<-mle_sp_D1(p=all.spD1.estimates[[i]]$par,z=mysim[i,],dmat.ml=dist.mat,theta.ml=theta,Bes_mat.ml=Bes_mat,index.mat.ml=index.mat,index.val.ml=index.val,u.ml=u,uniq.dist.ml=uniq.dist)$full.cov
  pred_performance(C=C_i,sample =mysim[i,] )
}
save.image("spd1predsn3.RData")

###### Semiparametric Delta=3 ######

spd3.preds<-foreach(i=1:50)%dopar%
{library(scoringRules)
  C_i<-mle_sp_D3(p=all.spD3.estimates[[i]]$par,z=mysim[i,],dmat.ml=dist.mat,theta.ml=theta,Bes_mat.ml=Bes_mat,index.mat.ml=index.mat,index.val.ml=index.val,u.ml=u,uniq.dist.ml=uniq.dist)$full.cov
  pred_performance(C=C_i,sample =mysim[i,] )
}
save.image("spd3predsn3.RData")


###### Semiparametric Delta=4 ######

spd4.preds<-foreach(i=1:50)%dopar%
{library(scoringRules)
  C_i<-mle_sp_D4(p=all.spD4.estimates[[i]]$par,z=mysim[i,],dmat.ml=dist.mat,theta.ml=theta,Bes_mat.ml=Bes_mat,index.mat.ml=index.mat,index.val.ml=index.val,u.ml=u,uniq.dist.ml=uniq.dist)$full.cov
  pred_performance(C=C_i,sample =mysim[i,] )
}
save.image("spd4predsn3.RData")

######################################################

#### Ploting prediction scores boxplot #####


lmc.rmse<-lmc.mspe<-lmc.mae<-lmc.logs<-lmc.crps<-lmc.nmse<-bvm.rmse<-bvm.mspe<-bvm.mae<-bvm.logs<-bvm.crps<-bvm.nmse<-sp.rmse<-sp.mspe<-sp.mae<-sp.logs<-sp.crps<-sp.nmse<-numeric(length = 50)

spd1.rmse<-spd1.mspe<-spd1.mae<-spd1.logs<-spd1.crps<-spd1.nmse<-spd3.rmse<-spd3.mspe<-spd3.mae<-spd3.logs<-spd3.crps<-spd3.nmse<-spd4.rmse<-spd4.mspe<-spd4.mae<-spd4.logs<-spd4.crps<-spd4.nmse<-numeric(length = 50)
for(i in 1:50)
{
  lmc.rmse[i]<-lmc.preds[[i]]$RMSE
  lmc.mspe[i]<-lmc.preds[[i]]$MSPE
  lmc.mae[i]<-lmc.preds[[i]]$MAE
  lmc.logs[i]<-lmc.preds[[i]]$Log_Score
  lmc.crps[i]<-lmc.preds[[i]]$CRPS
  lmc.nmse[i]<-lmc.preds[[i]]$NMSE
  bvm.rmse[i]<-bvm.preds[[i]]$RMSE
  bvm.mspe[i]<-bvm.preds[[i]]$MSPE
  bvm.mae[i]<-bvm.preds[[i]]$MAE
  bvm.logs[i]<-bvm.preds[[i]]$Log_Score
  bvm.crps[i]<-bvm.preds[[i]]$CRPS
  bvm.nmse[i]<-bvm.preds[[i]]$NMSE
  sp.rmse[i]<-sp.preds[[i]]$RMSE
  sp.mspe[i]<-sp.preds[[i]]$MSPE
  sp.mae[i]<-sp.preds[[i]]$MAE
  sp.logs[i]<-sp.preds[[i]]$Log_Score
  sp.crps[i]<-sp.preds[[i]]$CRPS
  sp.nmse[i]<-sp.preds[[i]]$NMSE
  
  spd1.rmse[i]<-spd1.preds[[i]]$RMSE
  spd1.mspe[i]<-spd1.preds[[i]]$MSPE
  spd1.mae[i]<-spd1.preds[[i]]$MAE
  spd1.logs[i]<-spd1.preds[[i]]$Log_Score
  spd1.crps[i]<-spd1.preds[[i]]$CRPS
  spd1.nmse[i]<-spd1.preds[[i]]$NMSE
  
  
  
  spd3.rmse[i]<-spd3.preds[[i]]$RMSE
  spd3.mspe[i]<-spd3.preds[[i]]$MSPE
  spd3.mae[i]<-spd3.preds[[i]]$MAE
  spd3.logs[i]<-spd3.preds[[i]]$Log_Score
  spd3.crps[i]<-spd3.preds[[i]]$CRPS
  spd3.nmse[i]<-spd3.preds[[i]]$NMSE
  
  
  spd4.rmse[i]<-spd4.preds[[i]]$RMSE
  spd4.mspe[i]<-spd4.preds[[i]]$MSPE
  spd4.mae[i]<-spd4.preds[[i]]$MAE
  spd4.logs[i]<-spd4.preds[[i]]$Log_Score
  spd4.crps[i]<-spd4.preds[[i]]$CRPS
  spd4.nmse[i]<-spd4.preds[[i]]$NMSE
  
}


par(mfrow=c(2,3))

boxplot(sp.rmse,spd1.rmse,spd3.rmse,spd4.rmse,lmc.rmse,bvm.rmse,main="RMSE",names=c("SPD2","SPD1","SPD3","SPD4","LMC","BVM"),col = "blue")
boxplot(sp.mspe,spd1.mspe,spd3.mspe,spd4.mspe,lmc.mspe,bvm.mspe,main="MSPE",names=c("SPD2","SPD1","SPD3","SPD4","LMC","BVM"))
boxplot(sp.mae,spd1.mae,spd3.mae,spd4.mae,lmc.mae,bvm.mae,main="MAE",names=c("SPD2","SPD1","SPD3","SPD4","LMC","BVM"))
boxplot(sp.logs,spd1.logs,spd3.logs,spd4.logs,lmc.logs,bvm.logs,main="LogS",names=c("SPD2","SPD1","SPD3","SPD4","LMC","BVM"))
boxplot(sp.crps,spd1.crps,spd3.crps,spd4.crps,lmc.crps,bvm.crps,main="CRPS",names=c("SPD2","SPD1","SPD3","SPD4","LMC","BVM"))
boxplot(sp.nmse,spd1.nmse,spd3.nmse,spd4.nmse,lmc.nmse,bvm.nmse,main="NMSE",names=c("SPD2","SPD1","SPD3","SPD4","LMC","BVM"))

library(ggplot2)
library(RColorBrewer)
library(viridis)

bp.rmse.data<-data.frame(RMSE=c(sp.rmse,spd1.rmse,spd3.rmse,spd4.rmse,lmc.rmse,bvm.rmse),
                         Model=c(rep("SPD2",times=length(sp.rmse)),
                                 rep("SPD1",times=length(spd1.rmse)),
                                 rep("SPD3",times=length(spd3.rmse)),
                                 rep("SPD4",times=length(spd4.rmse)),
                                 rep("LMC",times=length(lmc.rmse)),
                                 rep("BVM",times=length(bvm.rmse)))
                         )

ggplot(bp.rmse.data, aes(x=Model, y=RMSE, fill=Model)) + 
  geom_boxplot(alpha=0.8,coef=2.6) +
  theme(legend.position="none") +
  scale_fill_brewer(palette="BuPu")


bp.mae.data<-data.frame(MAE=c(sp.mae,spd1.mae,spd3.mae,spd4.mae,lmc.mae,bvm.mae),
                         Model=c(rep("SPD2",times=length(sp.mae)),
                                 rep("SPD1",times=length(spd1.mae)),
                                 rep("SPD3",times=length(spd3.mae)),
                                 rep("SPD4",times=length(spd4.mae)),
                                 rep("LMC",times=length(lmc.mae)),
                                 rep("BVM",times=length(bvm.mae)))
)

ggplot(bp.mae.data, aes(x=Model, y=MAE, fill=Model)) + 
  geom_boxplot(alpha=0.8,coef=2.6) +
  theme(legend.position="none") +
  scale_fill_brewer(palette="BuPu")



bp.nmse.data<-data.frame(NMSE=c(sp.nmse,spd1.nmse,spd3.nmse,spd4.nmse,lmc.nmse,bvm.nmse),
                        Model=c(rep("SPD2",times=length(sp.nmse)),
                                rep("SPD1",times=length(spd1.nmse)),
                                rep("SPD3",times=length(spd3.nmse)),
                                rep("SPD4",times=length(spd4.nmse)),
                                rep("LMC",times=length(lmc.nmse)),
                                rep("BVM",times=length(bvm.nmse)))
)

ggplot(bp.nmse.data, aes(x=Model, y=NMSE, fill=Model)) + 
  geom_boxplot(alpha=0.8,coef=2.4) +
  theme(legend.position="none") +
  scale_fill_brewer(palette="BuPu")



bp.crps.data<-data.frame(mCRPS=c(sp.crps,spd1.crps,spd3.crps,spd4.crps,lmc.crps,bvm.crps),
                         Model=c(rep("SPD2",times=length(sp.crps)),
                                 rep("SPD1",times=length(spd1.crps)),
                                 rep("SPD3",times=length(spd3.crps)),
                                 rep("SPD4",times=length(spd4.crps)),
                                 rep("LMC",times=length(lmc.crps)),
                                 rep("BVM",times=length(bvm.crps)))
)

ggplot(bp.crps.data, aes(x=Model, y=mCRPS, fill=Model)) + 
  geom_boxplot(alpha=0.8,coef=3) +
  theme(legend.position="none") +
  scale_fill_brewer(palette="BuPu")



bp.logs.data<-data.frame(mLogS=c(sp.logs,spd1.logs,spd3.logs,spd4.logs,lmc.logs,bvm.logs),
                         Model=c(rep("SPD2",times=length(sp.logs)),
                                 rep("SPD1",times=length(spd1.logs)),
                                 rep("SPD3",times=length(spd3.logs)),
                                 rep("SPD4",times=length(spd4.logs)),
                                 rep("LMC",times=length(lmc.logs)),
                                 rep("BVM",times=length(bvm.logs)))
)

ggplot(bp.logs.data, aes(x=Model, y=mLogS, fill=Model)) + 
  geom_boxplot(alpha=0.8,coef=1.98) +
  theme(legend.position="none") +
  scale_fill_brewer(palette="BuPu")






save.image("fullsimnew.RData")

deci<-4
pr.summary.dt<-data.frame(
  Model=c("SPD1","SPD2","SPD3","SPD4","LMC","BVM"),
  RSME=round(c(mean(spd1.rmse),mean(sp.rmse),mean(spd3.rmse),mean(spd4.rmse),mean(lmc.rmse),mean(bvm.rmse)),deci),
  MSPE=round(c(mean(spd1.mspe),mean(sp.mspe),mean(spd3.mspe),mean(spd4.mspe),mean(lmc.mspe),mean(bvm.mspe)),deci),
  MAE=round(c(mean(spd1.mae),mean(sp.mae),mean(spd3.mae),mean(spd4.mae),mean(lmc.mae),mean(bvm.mae)),deci),
  NMSE=round(c(mean(spd1.nmse),mean(sp.nmse),mean(spd3.nmse),mean(spd4.nmse),mean(lmc.nmse),mean(bvm.nmse)),deci),
  mCRPS=round(c(mean(spd1.crps),mean(sp.crps),mean(spd3.crps),mean(spd4.crps),mean(lmc.crps),mean(bvm.crps)),deci),
  mLogS=round(c(mean(spd1.logs),mean(sp.logs),mean(spd3.logs),mean(spd4.logs),mean(lmc.logs),mean(bvm.logs)),deci)
)
pr.summary.dt





