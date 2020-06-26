#############################################
##### Oscillatory Cross-covariance ##########
#############################################

#####################################################################
########### Coherence project simulation study ######################
#####################################################################
library(fields)
library(MASS)
#####################################################################
########### Creating required functions #############################
#####################################################################

##### Creating function for computing norm #####
norm_vec <- function(x) sqrt(sum(x^2))



bcoeff<-function(x)
{
  
  temp<-sin(x)
  temp[temp==1]<-1-1e-10
  temp[temp==-1]<-1+1e-10
  return(temp)
}




#####################################################################
############### Spectral density for old parametrisation ############
#####################################################################

f_mat<-function(w,phi,alpha,nu,d)
{ 
  temp<-numeric(length(w))
  for(i in 1:length(w))
  {
    temp[i]<-(phi*((alpha^2+w[i]^2)^(-nu-(d/2))))
  }
  return(temp)
}


###############################################################################
############### Matern covariance function for old parametrisation ############
###############################################################################
my.mat<-function(h,phi,alpha,nu,d)
{
  ln<-length(h)
  num1<-den1<-temp<-rep(NA,times=ln)
  
  for(i in 1:ln)
  {
    
    if(h[i]>0)
    {
      num1[i]<-((pi)^(d/2))*phi*((alpha*h[i])^nu)*besselK(x=alpha*h[i],nu=nu)
      den1[i]<-(2^(nu-1))*(gamma(nu+(d/2)))*(alpha^(2*nu))
      temp[i]<-num1[i]/den1[i]
    }
    
    else
    { 
      d2<-1e-10
      num1[i]<-((pi)^(d/2))*phi*(d2^nu)*besselK(x=d2,nu=nu)
      den1[i]<-(2^(nu-1))*(gamma(nu+(d/2)))*(alpha^(2*nu))
      temp[i]<-num1[i]/den1[i]
    }
  }
  return(temp)
}

#################################################################
##### Matern spectral density fro new parameterisation     ######
#################################################################

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

###############################################################################
############### Matern covariance function for new parametrisation ############
###############################################################################

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

n<-50 ##### number of observations = n^2
x<-y<-seq(1,n,length=n)

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


freq.max<-4.5

###############################################################
####################### Setting nodes #########################
###############################################################
u<-seq(0,freq.max,length=length(uniq.dist)-2)

########### Setting Matern parameters for variable 1  #######
a1<-1
nu1<-1
sigmasq1<-1
########### Setting Matern parameters for variable 2  #######
a2<-1
nu2<-1
sigmasq2<-1



###############################################################
##### Now we compute marginal covariance for variable 1 #######
###############################################################



##################################################################
####### Spectral density for the corresponding correlation #######
##################################################################

f.var1<-f.matern(w=u, nu=nu1, sigmasq = 1, a=a1, d=2)

#################################################################
################# Plotting spectral density #####################
#################################################################

plot(u,f.var1,type="l",main="Spectral density for matern correlation for variable 1")

############################################
######### Computing C11 ####################
############################################

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

#####################################################################################
##### Plotting approximated matern correlation vs true matern correlation  ##########
#####################################################################################

plot(uniq.dist,C11,type="l")
lines(uniq.dist,my.mat2(h=uniq.dist,sigmasq = sigmasq1,nu=nu1,a=a1),col="red",lty=2)
legend("topright",lty=c(1,2),col=c(1,2),c("Approximated Matern covariance","Exact Matern covariance"))


#######################################################
####### Now filling the entries of COV11 ##############
#######################################################

COV11<-matrix(NA,nrow=n^2,ncol = n^2)

for(i in 1:length(uniq.dist))
{
  COV11[which(dist.mat==uniq.dist[i],arr.ind = T)]<-C11[i]
}

########################################################
######### Checking positive definiteness ###############
########################################################

#chol(COV11) ##### Positive definite matrix #####




##################################################################
####### Spectral density for the corresponding correlation #######
##################################################################

f.var2<-f.matern(w=u, nu=nu2, sigmasq = 1, a=a2, d=2)

#################################################################
################# Plotting spectral density #####################
#################################################################

plot(u,f.var2,type="l",main="Spectral density for matern correlation for variable 2")

############################################
######### Computing C22 ####################
############################################

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

#####################################################################################
##### Plotting approximated matern correlation vs true matern correlation  ##########
#####################################################################################

plot(uniq.dist,C22,type="l")
lines(uniq.dist,my.mat2(h=uniq.dist,sigmasq = sigmasq2,nu=nu2,a=a2),col="red",lty=2)
legend("topright",lty=c(1,2),col=c(1,2),c("Approximated Matern covariance","Exact Matern covariance"))


#######################################################
####### Now filling the entries of COV22 ##############
#######################################################

COV22<-matrix(NA,nrow=n^2,ncol = n^2)

for(i in 1:length(uniq.dist))
{
  COV22[which(dist.mat==uniq.dist[i],arr.ind = T)]<-C22[i]
}

########################################################
######### Checking positive definiteness ###############
########################################################

#chol(COV22) ##### Positive definite matrix #####




#############################################################
########## Now computing the cross part #####################
#############################################################
###### Setting spline coefficients #######

b_3<--0.99
b_2<--0.99
b_1<-+0.99
b0<-+0.99
b1<-+0.99
b2<-+0.99
b3<--0.99
b4<--0.99

coh12<-b_3*Bspline(j=-3,k=4,delta = 1,x=u)+b_2*Bspline(j=-2,k=4,delta = 1,x=u)+b_1*Bspline(j=-1,k=4,delta = 1,x=u)+b0*Bspline(j=0,k=4,delta = 1,x=u)+b1*Bspline(j=1,k=4,delta = 1,x=u)+b2*Bspline(j=2,k=4,delta = 1,x=u)+b3*Bspline(j=3,k=4,delta = 1,x=u)+b4*Bspline(j=4,k=4,delta = 1,x=u)
plot(u,coh12,ylim=c(-1,1),type="l",main="Coherence")

f.var12<-coh12*sqrt(f.var1*f.var2)
plot(u,f.var12,main="Cross spectral density")


############################################
######### Computing C12 ####################
############################################

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

plot(uniq.dist,C12,type="l")



###############################################
##### Filling entries of COV12 ################
###############################################

COV12<-matrix(NA,nrow=n^2,ncol = n^2)

for(i in 1:length(uniq.dist))
{
  COV12[which(dist.mat==uniq.dist[i],arr.ind = T)]<-C12[i]
}

##################################################################################
############## Creating full covariance matrix ###################################
##################################################################################


##### Computing the whole covariance matrix #####
C1<-cbind(COV11,COV12)
C2<-cbind(t(COV12),COV22)
C<-rbind(C1,C2)

#chol(C)
##############################################
###### Simulating multivariate normal ########
##############################################

set.seed(123)
var<-mvrnorm(n=1,mu=rep(0,times=2*(n^2)),Sigma = C)
par(mfrow=c(1,2))
quilt.plot(x.s,y.s,var[1:n^2],ncol=n,nrow=n,xlab="x",ylab="y",zlim=c(min(var),max(var)))
quilt.plot(x.s,y.s,var[(n^2+1):(2*n^2)],ncol=n,nrow=n,xlab="x",ylab="y",zlim=c(min(var),max(var)))

####################################################################
############## Plots for manuscript #########################
####################################################################
par(mfrow=c(1,1))
par(cex.axis=2, cex.lab=3, cex.main=1, cex.sub=1,mar=c(3,3,1,1))
marcov12<-data.frame(distance=uniq.dist,cov=C12,var=rep("C12",times=length(uniq.dist)))
marcoh<-data.frame(freq=u,coh=coh12,var=rep("Coh",times=length(u)))

ggplot(data=marcov12, aes(x=distance, y=cov)) +
  geom_line(aes(color=var),show.legend = FALSE) +
  geom_point(aes(color=var),show.legend = FALSE)+labs(x= paste("||",expression(h),"||"),    # works fine
                                                      y= expression(C[12]("||h||")),    # works fine
                                                      color="")+ theme_grey(base_size = 31)


ggplot(data=marcoh, aes(x=freq, y=coh)) +
  geom_line(aes(color=var),show.legend = FALSE) +
  geom_point(aes(color=var),show.legend = FALSE)+labs(x= expression(omega),    # works fine
                                                      y= expression(gamma(omega)),    # works fine
                                                      color="")+ theme_grey(base_size = 31)


quilt.plot(x.s,y.s,var[1:n^2],ncol=n,nrow=n,xlab="x",ylab="y",zlim=c(min(var),max(var)))
quilt.plot(x.s,y.s,var[(n^2+1):(2*n^2)],ncol=n,nrow=n,xlab="x",ylab="y",zlim=c(min(var),max(var)))





################################################
######### Now we do the filtering ##############
################################################
x.s
y.s
nx<-length(unique(x.s))
ny<-length(unique(x.s))

delta1<-x[2]-x[1]
delta2<-y[2]-y[1]
n1<-length(x)
n2<-length(y)
f1<-seq(-n1/2,n1/2,by=1)
f2<-seq(-n2/2,n2/2,by=1)

f11<-f1/(delta1*n1)
f22<-f2/(delta2*n2)
omega1<-2*pi*f11
omega2<-2*pi*f11
w1s<-rep(omega1,times=length(omega2))
w2s<-rep(omega2,each=length(omega1))


############################################
######## Discrete Fourier transform ########
############################################
disc_2d_fft<-function(x,y,Z,w1,w2)
{
  
  x1<-unique(x)
  y1<-unique(y)
  
  N<-length(x1)
  ex1<-complex(length = (N+1)^2)
  ex2<-complex(length = (N+1)^2)
  for(i in 1:((N+1)^2))
  {
    
    
    
    theta<-(x*w1[i])+(y*w2[i])
    temp<-complex(real = cos(theta),imaginary = -sin(theta))
    ex1[i]<-sum(Z*temp)
    ex2[i]<-sum(Z*temp)
    
  }
  B1<-(ex1)
  
  delta1<-x1[2]-x1[1]
  delta2<-y1[2]-y1[1]
  delta<-delta1*delta2
  A<-1#delta/(((2*pi)^2)*N)
  return(A*B1)
  
}

####################################################################
################# Inverse Fourier trnasform ########################
####################################################################


inv_disc_2d_fft<-function(x,y,Z,w1,w2)
{
  
  x1<-unique(x)
  y1<-unique(y)
  
  N<-length(x1)
  ex1<-complex(length = (N-1)^2)
  ex2<-complex(length = (N-1)^2)
  for(i in 1:((N-1)^2))
  {
    
    
    
    theta<-(x*w1[i])+(y*w2[i])
    temp<-complex(real = cos(theta),imaginary = sin(theta))
    ex1[i]<-sum(Z*temp)
    ex2[i]<-sum(Z*temp)
    
  }
  B1<-(ex1)
  delta1<-x1[2]-x1[1]
  delta2<-y1[2]-y1[1]
  delta<-delta1*delta2
  
  A<-(N-1)^2
  return(B1/A)
  
}

########################################################
################# Masking Function #####################
########################################################

masking_fft<-function(spec,w1,w2,mv,mi)
{
  temp<-data.frame(w1=w1,w2=w2,fft=spec)
  img_max<-max(Re(spec))
  img_min<-min(Re(spec))
  par(mfrow=c(1,2))
  quilt.plot(x=temp$w1,y=temp$w2,z=Re(temp$fft),zlim=c(Re(img_min),Re(img_max)),nx=51,ny=51)
  nv<-numeric(length = length(w1))
  for(i in 1:length(w1))
  {
    nv[i]<-norm(c(w1[i],w2[i]))
  }
  log_vec<-(nv>=mv-mi)&(nv<=mv+mi)
  temp[!log_vec,3]<-0
  quilt.plot(x=temp$w1,y=temp$w2,z=Re(temp$fft),zlim=c(Re(img_min),Re(img_max)),nx=51,ny=51)
  par(mfrow=c(1,1))
  
  
  
  return(temp)
}


var1<-var[1:n^2]
var2<-var[(n^2+1):(2*n^2)]

####################################################
########## Variable 1 Fourier transform ############
####################################################

v1fft<-disc_2d_fft(x=x.s,y=y.s,Z=var1,w1=w1s,w2=w2s)
myv1fft<-matrix(v1fft,nrow = 51,ncol = 51,byrow = T)
image.plot(omega1,omega2,Re(myv1fft))


v1_inv_fft<-inv_disc_2d_fft(x=w1s,y=w2s,Z=myv1fft,w1=x.s,w2=y.s)
myv1_inv_fft<-matrix(v1_inv_fft,nrow = 50,ncol = 50,byrow = T)
#image.plot(x,y,Re(myv1_inv_fft),zlim=c(-4,4))
quilt.plot(x.s,y.s,Re(myv1_inv_fft),zlim=c(-4,4),nx=50,ny=50)
quilt.plot(x.s,y.s,var1,zlim=c(-4,4),nx=50,ny=50)


####################################################
########## Variable 2 Fourier transform ############
####################################################

v2fft<-disc_2d_fft(x=x.s,y=y.s,Z=var2,w1=w1s,w2=w2s)
myv2fft<-matrix(v2fft,nrow = 51,ncol = 51,byrow = T)
image.plot(omega1,omega2,Re(myv2fft))


v2_inv_fft<-inv_disc_2d_fft(x=w1s,y=w2s,Z=myv2fft,w1=x.s,w2=y.s)
myv2_inv_fft<-matrix(v2_inv_fft,nrow = 50,ncol = 50,byrow = T)
#image.plot(x,y,Re(myv1_inv_fft),zlim=c(-4,4))
quilt.plot(x.s,y.s,Re(myv2_inv_fft),zlim=c(-4,4),nx=50,ny=50)
quilt.plot(x.s,y.s,var2,zlim=c(-4,4),nx=50,ny=50)





###########################################################
####### Here we set the frequency interval ################
###########################################################
          # set the following band values #
mvs=3.5   #band1=0  band2=3.5
mis=0.5  #band1=0.2 band2=0.5

maskvar1<-masking_fft(spec=v1fft,w1=w1s,w2=w2s,mv=mvs,mi=mis)
maskvar2<-masking_fft(spec=v2fft,w1=w1s,w2=w2s,mv=mvs,mi=mis)



freq_comp_var1<-inv_disc_2d_fft(x=w1s,y=w2s,Z=maskvar1$fft,w1=x.s,w2=y.s)
freq_comp_var2<-inv_disc_2d_fft(x=w1s,y=w2s,Z=maskvar2$fft,w1=x.s,w2=y.s)


####################################################################
######### Plotting each variables frequency component ##############
####################################################################
img_min<-min(c(Re(freq_comp_var1),Re(freq_comp_var2)))
img_max<-max(c(Re(freq_comp_var1),Re(freq_comp_var2)))

par(mfrow=c(1,1))
par(cex.axis=2, cex.lab=3, cex.main=1, cex.sub=1,mar=c(3,3,1,1))
quilt.plot(x.s,y.s,Re(freq_comp_var1),zlim=c(img_min,img_max),nx=50,ny=50,xlab="x",ylab="y")
quilt.plot(x.s,y.s,Re(freq_comp_var2),zlim=c(img_min,img_max),nx=50,ny=50,xlab="x",ylab="y")


rho12<-cor(Re(freq_comp_var1),Re(freq_comp_var2))

rho12  #band1=-0.4972403,band2=0.94266










