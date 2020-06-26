###############################################################################
#################### Code for generating approximation table ##################
###############################################################################

#### Creating  required functions

#####################################################################
########### Creating required functions #############################
#####################################################################

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



###############################################################################################################
########### Here in this section we create a function to compute covariance from our model ####################
###############################################################################################################

semi.p.cov<-function(h,u,a1,nu1,sigmasq1,a2,nu2,a12,nu12,rho12,sigmasq2,coherence)
{
  
  
  
  f.var1<-f.matern(w=u, nu=nu1, sigmasq = 1, a=a1, d=2)
  f.var2<-f.matern(w=u, nu=nu2, sigmasq = 1, a=a2, d=2)
  coh12<-coherence
  f.var12<-coh12*sqrt(f.var1*f.var2)
  
  #### Computing Block Covariance matrices ####
  uniq.dist<-h
  #####C11#####
  
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
  
  ###############################################
  ########## Mutiplying marginal variance #######
  ###############################################
  C11<-sigmasq1*C11
  
  
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
  
  C22<-sigmasq2*C22
  
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
  
  ####################################
  ####### Scaling it #################
  ####################################
  
  C12<-C12/sqrt(scl.1*scl.2)
  
  
  ###############################################
  ########## Mutiplying marginal variance #######
  ###############################################
  
  C12<-sqrt(sigmasq1*sigmasq2)*C12
  
  ###############################################
  ############# Comparison plots ################
  ###############################################
  
  par(mfrow=c(1,3))
  plot(h,C11,xlab="h",ylab="C11",main="Marginal Covariance var1",pch=19,col="red",cex=1)
  points(h,my.mat2(h=h,sigmasq = sigmasq1,nu=nu1,a=a1),pch=19,col="green",cex=0.5)
  
  plot(h,C22,xlab="h",ylab="C22",main="Marginal Covariance var2",pch=19,col="red",cex=1)
  points(h,my.mat2(h=h,sigmasq = sigmasq2,nu=nu2,a=a2),pch=19,col="green",cex=0.5)
  
  plot(h,C12,xlab="h",ylab="C12",main="Cross Covariance",pch=19,col="red",cex=1)
  points(h,rho12*sqrt(sigmasq1*sigmasq2)*my.mat2(h=h,sigmasq = 1,nu=nu12,a=a12),pch=19,col="green",cex=0.5)
  
  c11true<-my.mat2(h=h,sigmasq = sigmasq1,nu=nu1,a=a1)
  c22true<-my.mat2(h=h,sigmasq = sigmasq2,nu=nu2,a=a2)
  c12true<-rho12*sqrt(sigmasq1*sigmasq2)*my.mat2(h=h,sigmasq = 1,nu=nu12,a=a12)
  rmse11<-mean((c11true-C11)^2)
  rmse22<-mean((c22true-C22)^2)
  rmse12<-mean((c12true-C12)^2)
  
  nmse11<-1-(sum((c11true-C11)^2))/(sum((c11true-mean(c11true))^2)) #Normalized mean squared error
  nmse22<-1-(sum((c22true-C22)^2))/(sum((c22true-mean(c22true))^2)) #Normalized mean squared error
  nmse12<-1-(sum((c12true-C12)^2))/(sum((c12true-mean(c12true))^2)) #Normalized mean squared error
  
  
  return(list(C11_RMSE=rmse11,C22_RMSE=rmse22,C12_RMSE=rmse12, C11_NMSE=nmse11,C22_NMSE=nmse22,C12_NMSE=nmse12,cross_covariance=C12))
  
}




#########################################################
######## Full bivariate matern infimum function #########
#########################################################

### Infimum function (to check for upper bound of rho12 such that full bivariate matern is valid in case of v12>0.5(v1+v2))
infi<-function(nu1,nu2,a1,a2,nu12,a12,d)
{
  num<-gamma(nu1+(d/2))*gamma(nu2+(d/2))*((gamma(nu12))^2)*(a1^(2*nu1))*(a2^(2*nu2))*((a12^2)^(2*nu12+d))
  den<-gamma(nu1)*gamma(nu2)*((gamma(nu12+(d/2)))^2)*(a12^(4*nu12))*((a1^2)^(nu1+(d/2)))*((a2^2)^(nu2+(d/2)))
  return(num/den)
}

###############################################################################################################
######### Checking validity of the chosen full bivariate matern parameter (section 2.3 Gneiting 2010) #########
###############################################################################################################

fullbiwm.check<-function(a1,nu1,a2,nu2,a12,nu12,rho12)
{
  nu12.check<-(nu12>=0.5*(nu1+nu2))
  a12.check<-(a12^2>=0.5*((a1^2)+(a2^2)))
  num1<-(a1^(nu1))*(a2^(nu2))*(gamma(nu12))*((exp(1)*((a12^2)-0.5*(a1^2+a2^2)))^(nu12-0.5*(nu1+nu2)))
  den1<-(a12^(2*nu12))*(gamma(nu1)^0.5)*(gamma(nu2)^0.5)*((nu12-0.5*(nu1+nu2))^(nu12-0.5*(nu1+nu2)))
  rhobound<-num1/den1
  rho.check<-(rho12<=rhobound)
  return(list(Validity=(nu12.check & a12.check & rho.check),RhoBound=rhobound))
}
#############################
#### Distance vector ########
#############################

hs<-seq(0,40,length=500)


##############################################################
######## Setting full bivariate matern parameters ############
##############################################################
# case 1
sigmasq1.bvm<-1 
sigmasq2.bvm<-1 
a1.bvm<-0.5   
a2.bvm<-0.5 
a12.bvm<-0.5 
nu1.bvm<-1 
nu2.bvm<-1 
nu12.bvm<-1.5 
rho.bvm<-0.05 
u<-seq(0,9.9,length=200) ### frequencies
fullbiwm.check(a1=a1.bvm,nu1=nu1.bvm,a2=a2.bvm,nu2=nu2.bvm,a12=a12.bvm,nu12=nu12.bvm,rho12=rho.bvm)
par(mfrow=c(1,1))
plot(u,sqrt(biwm_coh(w=u,a1=a1.bvm,a2=a2.bvm,v1=nu1.bvm,v2=nu2.bvm,a12=a12.bvm,v12=nu12.bvm,d=2,rho=rho.bvm)),main="Coherence",ylab="rho")


##### Plot var1 spectral density ######
plot(u,f.matern(w=u,nu=nu1.bvm,sigmasq = sigmasq1.bvm,a=a1.bvm,d=2),xlab="w",ylab="f1",main="Spectral density for var1")
plot(u,f.matern(w=u,nu=nu2.bvm,sigmasq = sigmasq2.bvm,a=a2.bvm,d=2),xlab="w",ylab="f2",main="Spectral density for var2")
plot(u,sqrt(biwm_coh(w=u,a1=a1.bvm,a2=a2.bvm,v1=nu1.bvm,v2=nu2.bvm,a12=a12.bvm,v12=nu12.bvm,d=2,rho=rho.bvm)),xlab="w",ylab="f12",main="Coherence")
plot(u,f.matern(w=u,nu=nu12.bvm,sigmasq = rho.bvm*sqrt(sigmasq2.bvm*sigmasq1.bvm),a=a12.bvm,d=2),xlab="w",ylab="f12",main="Cross-spectral density")





###############################################################################################################
########### Here in this section we decide the node and coherence #############################################
###############################################################################################################



u<-seq(0,9.9,length=200)

truecoh<-sqrt(biwm_coh(w=u,a1=a1.bvm,a2=a2.bvm,v1=nu1.bvm,v2=nu2.bvm,a12=a12.bvm,v12=nu12.bvm,d=2,rho=rho.bvm))
par(mfrow=c(1,1))
plot(u,truecoh)

#########################################################################################################
######## Creating function to obtain least square fit of the coherence function #########################
#########################################################################################################

ls.coh<-function(p,arg_true_coh)
{
  b_3<-bcoeff(p[1])
  b_2<-bcoeff(p[2])
  b_1<-bcoeff(p[3])
  b0<-bcoeff(p[4])
  b1<-bcoeff(p[5])
  b2<-bcoeff(p[6])
  b3<-bcoeff(p[7])
  b4<-bcoeff(p[8])
  #b5<-bcoeff(p[9])
  #b6<-bcoeff(p[10])
  
  delta.s<-2
  mycoh<-b_3*Bspline(j=-3,k=4,delta = delta.s,x=u)+b_2*Bspline(j=-2,k=4,delta = delta.s,x=u)+b_1*Bspline(j=-1,k=4,delta = delta.s,x=u)+b0*Bspline(j=0,k=4,delta = delta.s,x=u)+b1*Bspline(j=1,k=4,delta = delta.s,x=u)+b2*Bspline(j=2,k=4,delta = delta.s,x=u)+b3*Bspline(j=3,k=4,delta = delta.s,x=u)+b4*Bspline(j=4,k=4,delta = delta.s,x=u)
  return(list(Est_coh=mycoh,RSS=sum((mycoh-arg_true_coh)^2)))
}

ls.optim<-function(p)
{
  return(ls.coh(p=p,arg_true_coh = truecoh)$RSS)
}

u<-seq(0,9.9,length=300)

truecoh<-sqrt(biwm_coh(w=u,a1=a1.bvm,a2=a2.bvm,v1=nu1.bvm,v2=nu2.bvm,a12=a12.bvm,v12=nu12.bvm,d=2,rho=rho.bvm))
par(mfrow=c(1,1))
plot(u,truecoh)
init.case1<-c(0.168650436, 0.064034495,  0.006566008,  0.012961652, 0.003252825, 0.007349998, -0.001523553, 0.030244297)
spline.case1<-optim(par=init.case1,ls.optim,control=list(trace=6,maxit=3000))

u<-seq(0,9.9,length=2000)
truecoh<-sqrt(biwm_coh(w=u,a1=a1.bvm,a2=a2.bvm,v1=nu1.bvm,v2=nu2.bvm,a12=a12.bvm,v12=nu12.bvm,d=2,rho=rho.bvm))
mycoh<-ls.coh(p=spline.case1$par,arg_true_coh = truecoh)$Est_coh ### Chosen coherence function for case 1

case1summary<-semi.p.cov(h=hs,u=u,a1=a1.bvm,nu1=nu1.bvm,sigmasq1=sigmasq1.bvm,a2=a2.bvm,nu2=nu2.bvm,a12=a12.bvm,nu12=nu12.bvm,rho12=rho.bvm,sigmasq2=sigmasq2.bvm,coherence=mycoh)
coherence_case1<-data.frame(freq=u,coh=mycoh)
C12_case1<-case1summary$cross_covariance
#############################################
################### Case 2 ##################
#############################################

# case 2
sigmasq1.bvm<-1 #1
sigmasq2.bvm<-1 #1
a1.bvm<-1   #2
a2.bvm<-1 #0.5
a12.bvm<-1.1 #2.5
nu1.bvm<-2 #1
nu2.bvm<-3 #1
nu12.bvm<-5 #1.5
rho.bvm<-0.1 #0.13

fullbiwm.check(a1=a1.bvm,nu1=nu1.bvm,a2=a2.bvm,nu2=nu2.bvm,a12=a12.bvm,nu12=nu12.bvm,rho12=rho.bvm)
par(mfrow=c(1,1))
plot(u,sqrt(biwm_coh(w=u,a1=a1.bvm,a2=a2.bvm,v1=nu1.bvm,v2=nu2.bvm,a12=a12.bvm,v12=nu12.bvm,d=2,rho=rho.bvm)),main="Coherence",ylab="rho",xlim=c(0,4.9))


##### Plot var1 spectral density ######
plot(u,f.matern(w=u,nu=nu1.bvm,sigmasq = sigmasq1.bvm,a=a1.bvm,d=2),xlab="w",ylab="f1",main="Spectral density for var1")
plot(u,f.matern(w=u,nu=nu2.bvm,sigmasq = sigmasq2.bvm,a=a2.bvm,d=2),xlab="w",ylab="f2",main="Spectral density for var2")
plot(u,sqrt(biwm_coh(w=u,a1=a1.bvm,a2=a2.bvm,v1=nu1.bvm,v2=nu2.bvm,a12=a12.bvm,v12=nu12.bvm,d=2,rho=rho.bvm)),xlab="w",ylab="f12",main="Coherence")
plot(u,f.matern(w=u,nu=nu12.bvm,sigmasq = rho.bvm*sqrt(sigmasq2.bvm*sigmasq1.bvm),a=a12.bvm,d=2),xlab="w",ylab="f12",main="Cross-spectral density")





###############################################################################################################
########### Here in this section we decide the node and coherence #############################################
###############################################################################################################



u<-seq(0,6.9,length=300)

truecoh<-sqrt(biwm_coh(w=u,a1=a1.bvm,a2=a2.bvm,v1=nu1.bvm,v2=nu2.bvm,a12=a12.bvm,v12=nu12.bvm,d=2,rho=rho.bvm))
par(mfrow=c(1,1))
plot(u,truecoh)


ls.coh<-function(p,arg_true_coh)
{
  b_3<-bcoeff(p[1])
  b_2<-bcoeff(p[2])
  b_1<-bcoeff(p[3])
  b0<-bcoeff(p[4])
  b1<-bcoeff(p[5])
  b2<-bcoeff(p[6])
  b3<-bcoeff(p[7])
  b4<-bcoeff(p[8])
  b5<-bcoeff(p[9])
  b6<-bcoeff(p[10])
  
  delta.s<-1
  mycoh<-b_3*Bspline(j=-3,k=4,delta = delta.s,x=u)+b_2*Bspline(j=-2,k=4,delta = delta.s,x=u)+b_1*Bspline(j=-1,k=4,delta = delta.s,x=u)+b0*Bspline(j=0,k=4,delta = delta.s,x=u)+b1*Bspline(j=1,k=4,delta = delta.s,x=u)+b2*Bspline(j=2,k=4,delta = delta.s,x=u)+b3*Bspline(j=3,k=4,delta = delta.s,x=u)+b4*Bspline(j=4,k=4,delta = delta.s,x=u)+b5*Bspline(j=5,k=4,delta = delta.s,x=u)+b6*Bspline(j=6,k=4,delta = delta.s,x=u)
  return(list(Est_coh=mycoh,RSS=sum((mycoh-arg_true_coh)^2)))
}

ls.optim<-function(p)
{
  return(ls.coh(p=p,arg_true_coh = truecoh)$RSS)
}

u<-seq(0,6.9,length=300)

truecoh<-sqrt(biwm_coh(w=u,a1=a1.bvm,a2=a2.bvm,v1=nu1.bvm,v2=nu2.bvm,a12=a12.bvm,v12=nu12.bvm,d=2,rho=rho.bvm))
par(mfrow=c(1,1))
plot(u,truecoh)
init.case2<-c(0.1, 0.22,  0.03,  0, 0, 0, 0, 0,0,0)
spline.case2<-optim(par=init.case2,ls.optim,control=list(trace=6,maxit=3000))
points(u,ls.coh(p=spline.case2$par,arg_true_coh = truecoh)$Est_coh,col="blue")

u<-seq(0,6.9,length=2000)
truecoh<-sqrt(biwm_coh(w=u,a1=a1.bvm,a2=a2.bvm,v1=nu1.bvm,v2=nu2.bvm,a12=a12.bvm,v12=nu12.bvm,d=2,rho=rho.bvm))
mycoh<-ls.coh(p=spline.case2$par,arg_true_coh = truecoh)$Est_coh

case2summary<-semi.p.cov(h=hs,u=u,a1=a1.bvm,nu1=nu1.bvm,sigmasq1=sigmasq1.bvm,a2=a2.bvm,nu2=nu2.bvm,a12=a12.bvm,nu12=nu12.bvm,rho12=rho.bvm,sigmasq2=sigmasq2.bvm,coherence=mycoh)

C12_case2<-case2summary$cross_covariance
coherence_case2<-data.frame(freq=u,coh=mycoh)


#############################################
################### Case 3 ##################
#############################################

# case 3
sig1.bvm<-1 #1
sig2.bvm<-1 #1
a1.bvm<-0.6   #2
a2.bvm<-1.4 #0.5
a12.bvm<-1.5 #2.5
nu1.bvm<-3 #1
nu2.bvm<-3 #1
nu12.bvm<-4 #1.5
rho.bvm<-0.1 #0.13

fullbiwm.check(a1=a1.bvm,nu1=nu1.bvm,a2=a2.bvm,nu2=nu2.bvm,a12=a12.bvm,nu12=nu12.bvm,rho12=rho.bvm)
par(mfrow=c(1,1))
plot(u,sqrt(biwm_coh(w=u,a1=a1.bvm,a2=a2.bvm,v1=nu1.bvm,v2=nu2.bvm,a12=a12.bvm,v12=nu12.bvm,d=2,rho=rho.bvm)),main="Coherence",ylab="rho",xlim=c(0,4.9))


##### Plot var1 spectral density ######
plot(u,f.matern(w=u,nu=nu1.bvm,sigmasq = sigmasq1.bvm,a=a1.bvm,d=2),xlab="w",ylab="f1",main="Spectral density for var1")
plot(u,f.matern(w=u,nu=nu2.bvm,sigmasq = sigmasq2.bvm,a=a2.bvm,d=2),xlab="w",ylab="f2",main="Spectral density for var2")
plot(u,sqrt(biwm_coh(w=u,a1=a1.bvm,a2=a2.bvm,v1=nu1.bvm,v2=nu2.bvm,a12=a12.bvm,v12=nu12.bvm,d=2,rho=rho.bvm)),xlab="w",ylab="f12",main="Coherence")
plot(u,f.matern(w=u,nu=nu12.bvm,sigmasq = rho.bvm*sqrt(sigmasq2.bvm*sigmasq1.bvm),a=a12.bvm,d=2),xlab="w",ylab="f12",main="Cross-spectral density")





###############################################################################################################
########### Here in this section we decide the node and coherence #############################################
###############################################################################################################



u<-seq(0,6.9,length=300)

truecoh<-sqrt(biwm_coh(w=u,a1=a1.bvm,a2=a2.bvm,v1=nu1.bvm,v2=nu2.bvm,a12=a12.bvm,v12=nu12.bvm,d=2,rho=rho.bvm))
par(mfrow=c(1,1))
plot(u,truecoh)


ls.coh<-function(p,arg_true_coh)
{
  b_3<-bcoeff(p[1])
  b_2<-bcoeff(p[2])
  b_1<-bcoeff(p[3])
  b0<-bcoeff(p[4])
  b1<-bcoeff(p[5])
  b2<-bcoeff(p[6])
  b3<-bcoeff(p[7])
  b4<-bcoeff(p[8])
  b5<-bcoeff(p[9])
  b6<-bcoeff(p[10])
  b7<-bcoeff(p[11])
  b8<-bcoeff(p[12])
  b9<-bcoeff(p[13])
  b10<-bcoeff(p[14])
  b11<-bcoeff(p[15])
  b12<-bcoeff(p[16])
  b13<-bcoeff(p[17])
  
  delta.s<-0.5
  mycoh<-b_3*Bspline(j=-3,k=4,delta = delta.s,x=u)+b_2*Bspline(j=-2,k=4,delta = delta.s,x=u)+b_1*Bspline(j=-1,k=4,delta = delta.s,x=u)+b0*Bspline(j=0,k=4,delta = delta.s,x=u)+b1*Bspline(j=1,k=4,delta = delta.s,x=u)+b2*Bspline(j=2,k=4,delta = delta.s,x=u)+b3*Bspline(j=3,k=4,delta = delta.s,x=u)+b4*Bspline(j=4,k=4,delta = delta.s,x=u)+b5*Bspline(j=5,k=4,delta = delta.s,x=u)+b6*Bspline(j=6,k=4,delta = delta.s,x=u)+b7*Bspline(j=7,k=4,delta = delta.s,x=u)+b8*Bspline(j=8,k=4,delta = delta.s,x=u)+b9*Bspline(j=9,k=4,delta = delta.s,x=u)+b10*Bspline(j=10,k=4,delta = delta.s,x=u)+b11*Bspline(j=11,k=4,delta = delta.s,x=u)+b12*Bspline(j=12,k=4,delta = delta.s,x=u)+b13*Bspline(j=13,k=4,delta = delta.s,x=u)
  return(list(Est_coh=mycoh,RSS=sum((mycoh-arg_true_coh)^2)))
}

ls.optim<-function(p)
{
  return(ls.coh(p=p,arg_true_coh = truecoh)$RSS)
}

u<-seq(0,6.9,length=300)

truecoh<-sqrt(biwm_coh(w=u,a1=a1.bvm,a2=a2.bvm,v1=nu1.bvm,v2=nu2.bvm,a12=a12.bvm,v12=nu12.bvm,d=2,rho=rho.bvm))
par(mfrow=c(1,1))
plot(u,truecoh)
init.case3<-c(asin(0.14),asin(0.01),asin(0.1),asin(0.25), asin(0.41), asin(0.4),  asin(0.41), asin(0.32), asin(0.31), asin(0.21),asin(0.22),asin(0.19),asin(0.13),asin(0.13),asin(0.12),asin(0.11),asin(0.11))

spline.case3<-optim(par=init.case3,ls.optim,control=list(trace=6,maxit=3000))
points(u,ls.coh(p=spline.case3$par,arg_true_coh = truecoh)$Est_coh,col="blue")
ls.coh(p=spline.case3$par,arg_true_coh = truecoh)$RSS
init.case3<-(spline.case3$par)
u<-seq(0,6.9,length=300)

truecoh<-sqrt(biwm_coh(w=u,a1=a1.bvm,a2=a2.bvm,v1=nu1.bvm,v2=nu2.bvm,a12=a12.bvm,v12=nu12.bvm,d=2,rho=rho.bvm))

spline.case3<-optim(par=init.case3,ls.optim,control=list(trace=6,maxit=3000))
points(u,ls.coh(p=spline.case3$par,arg_true_coh = truecoh)$Est_coh,col="blue")

u<-seq(0,6.9,length=2000)
truecoh<-sqrt(biwm_coh(w=u,a1=a1.bvm,a2=a2.bvm,v1=nu1.bvm,v2=nu2.bvm,a12=a12.bvm,v12=nu12.bvm,d=2,rho=rho.bvm))
mycoh<-ls.coh(p=spline.case3$par,arg_true_coh = truecoh)$Est_coh

case3summary<-semi.p.cov(h=hs,u=u,a1=a1.bvm,nu1=nu1.bvm,sigmasq1=sigmasq1.bvm,a2=a2.bvm,nu2=nu2.bvm,a12=a12.bvm,nu12=nu12.bvm,rho12=rho.bvm,sigmasq2=sigmasq2.bvm,coherence=mycoh)

C12_case3<-case3summary$cross_covariance
coherence_case3<-data.frame(freq=u,coh=mycoh)
###################################################
############### Plots for manuscript ##############
###################################################

par(mfrow=c(1,1))
par(cex.axis=2, cex.lab=3, cex.main=1, cex.sub=1,mar=c(3,3,1,1))


library(ggplot2)
C12_c1<-data.frame(distance=hs,Cross_covariance=C12_case1,var=c(rep("Model 1",times=length(C12_case1))))
C12_c2<-data.frame(distance=hs,Cross_covariance=C12_case2,var=c(rep("Model 2",times=length(C12_case2))))
C12_c3<-data.frame(distance=hs,Cross_covariance=C12_case3,var=c(rep("Model 3",times=length(C12_case3))))

ggplot(data=C12_c1, aes(x=distance, y=Cross_covariance)) +
  geom_line(aes(color=var),show.legend = FALSE) +
  geom_point(aes(color=var),show.legend = FALSE)+labs(x= paste("||",expression(h),"||"),    # works fine
                                  y= expression(C[12]("||h||")),    # works fine
                                  color="")+ theme_grey(base_size = 31)
coherence_case1$var<-rep("gamma12",times=length(coherence_case1$coh))


ggplot(data=coherence_case1, aes(x=freq, y=coh)) +
  geom_line(aes(color=var),show.legend = FALSE) +
  geom_point(aes(color=var),show.legend = FALSE)+labs(x= expression(omega),    # works fine
                                                      y= expression(gamma(omega)),    # works fine
                                                      color="")+ theme_grey(base_size = 31)


ggplot(data=C12_c2, aes(x=distance, y=Cross_covariance)) +
  geom_line(aes(color=var),show.legend = FALSE) +
  geom_point(aes(color=var),show.legend = FALSE)+labs(x= paste("||",expression(h),"||"),    # works fine
                                                      y= expression(C[12]("||h||")),    # works fine
                                                      color="")+ theme_grey(base_size = 31)
coherence_case2$var<-rep("gamma12",times=length(coherence_case2$coh))


ggplot(data=coherence_case2, aes(x=freq, y=coh)) +
  geom_line(aes(color=var),show.legend = FALSE) +
  geom_point(aes(color=var),show.legend = FALSE)+labs(x= expression(omega),    # works fine
                                                      y= expression(gamma(omega)),    # works fine
                                                      color="")+ theme_grey(base_size = 31)



ggplot(data=C12_c3, aes(x=distance, y=Cross_covariance)) +
  geom_line(aes(color=var),show.legend = FALSE) +
  geom_point(aes(color=var),show.legend = FALSE)+labs(x= paste("||",expression(h),"||"),    # works fine
                                                      y= expression(C[12]("||h||")),    # works fine
                                                      color="")+ theme_grey(base_size = 31)
coherence_case3$var<-rep("gamma12",times=length(coherence_case3$coh))


ggplot(data=coherence_case3, aes(x=freq, y=coh)) +
  geom_line(aes(color=var),show.legend = FALSE) +
  geom_point(aes(color=var),show.legend = FALSE)+labs(x= expression(omega),    # works fine
                                                      y= expression(gamma(omega)),    # works fine
                                                      color="")+ theme_grey(base_size = 31)



