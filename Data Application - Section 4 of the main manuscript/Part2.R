#######################################################################
###################### Part 2 of the Data application code ############
#######################################################################

###################################################################################################
######### In this part, we estimate Independent Matern, Bivariate matern, and LMC model ###########
######### for the 100 times randomly splitted training  sample  locations ####################
###################################################################################################


######### Loading libraries #########
library(fields)
library(doParallel)
#####################################


######################################
###### Setting work directory ######## 
######################################
## set the directory where the file "Final dataset to be worked with.RData" is saved from part 1
setwd("/ibex/scratch/qadirga/WNC_final")

##########################################
###### Loading the image from part 1 #####
##########################################

load("Final dataset to be worked with.RData")



######## Loading required functions ########

#####################################
##### Matern spectral density  ######
#####################################


f.matern<-function(w,nu,sigma,a,d)
{
  const<-(sigma^2)*gamma(nu+(d/2))*(a^(2*nu))/(gamma(nu)*(pi^(d/2)))
  varying<-1/((a+w^2)^(nu+(d/2)))
  return(const*varying)
}

#######################################################
############### Matern covariance function ############
#######################################################

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


#############################################
######## Bivariate matern coherence #########
#############################################


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





#######################################
####### Setting number of cores #######
#######################################


ncores<-detectCores()-2
registerDoParallel(cores = ncores)



############## Likelihood Estimation ##############



ind_bvm_lmc_estimations<-foreach(run=1:100)%dopar%{
  library(fields)
  
  
  
  un.grd.train<-un.grd.total[-rand.index[,run],]  ########## Training set ######### 
  un.grd.test<-un.grd.total[rand.index[,run],]    ########## Test set #############
  
  
  dist.mat.train<-rdist(un.grd.train[,-c(3,4)])
  
  
  uniq.dist.train<-unique(c(dist.mat.train))
  uniq.dist.train<-sort(uniq.dist.train)  ##### sorting distances in increasing order #####
  
  index.mat.train<-NULL
  index.val.train<-NULL
  for(i in 1:length(uniq.dist.train))
  {
    index.mat.train<-rbind(index.mat.train,which(dist.mat.train==uniq.dist.train[i],arr.ind = T))
    index.val.train<-c(index.val.train,rep(i,times=nrow(which(dist.mat.train==uniq.dist.train[i],arr.ind = T))))
  }
  
  
  
  mle_ind_mat<-function(p,z,dmat.ml,index.mat.ml,index.val.ml,uniq.dist.ml)
  {
    a1<-p[1]
    nu1<-p[2]
    sigma1<-p[3]
    a2<-p[4]
    nu2<-p[5]
    sigma2<-p[6]
    nug1<-p[7]
    nug2<-p[8]
    if(sum(p[1:6]<0)!=0||nug1<0||nug2<0)
    {
      nloglikelihood<-10000000
      return(list(mlv=nloglikelihood,params=NULL))
    }
    else
    {
      
      dist.mat<-dmat.ml
      C11<-my.matern(h=uniq.dist.ml,a=a1,sigma = sigma1,nu=nu1)
      C22<-my.matern(h=uniq.dist.ml,a=a2,sigma = sigma2,nu=nu2)
      
      COV11op<-COV22op<-COV12op<-matrix(0,nrow = nrow(dist.mat),ncol=ncol(dist.mat))
      
      COV11op[index.mat.ml]<-C11[index.val.ml]
      COV22op[index.mat.ml]<-C22[index.val.ml]
      
      NUG1<-diag(nug1,nrow = nrow(dist.mat),ncol=ncol(dist.mat))
      NUG2<-diag(nug2,nrow = nrow(dist.mat),ncol=ncol(dist.mat))
      
      ############## Inverting C11 ##########
      C<-rbind(cbind(COV11op+NUG1,COV12op),cbind(t(COV12op),COV22op+NUG2))
      Cchol<-chol(C)
      Cinv<-chol2inv(Cchol)
      logD<-determinant(C)$modulus
      
      nloglikelihood <-(0.5 * logD + 0.5 * t(z) %*% Cinv %*% z+0.5*length(z)*log(2*pi))
      if(abs(nloglikelihood) == Inf || is.nan(nloglikelihood)){ nloglikelihood <- 1e+08}
      return(list(mlv=nloglikelihood,a1=a1,a2=a2,nu1=nu1,nu2=nu2,sigma1=sigma1,sigma2=sigma2,full.cov=C))
    }
    
  }
  
  
  init.ind<-c(1,1,1,1,1,1,0,0)
  
  
  mle_ind_mlv<-function(pars)
  {
    return(mle_ind_mat(p=pars,z=c(un.grd.train$PM2_5,un.grd.train$WS),dmat.ml=dist.mat.train,index.mat.ml=index.mat.train,index.val.ml=index.val.train,uniq.dist.ml=uniq.dist.train)$mlv)
  }
  
  
  optim_indmat_loglik <- function(par){
    optim(par=par,
          fn = mle_ind_mlv,
          hessian=FALSE,
          control=list(trace=6,
                       pgtol=0,
                       parscale=rep(0.1,length(par)),
                       maxit=5000))
  }
  
  indmat.estim<-optim_indmat_loglik(init.ind)
  
  
  
  
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
    nug1<-p[10]
    nug2<-p[11]
    if(sum(p[1:6]<0)!=0||p[7]<0||p[8]<0||nug1<0||nug2<0)
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
      
      NUG1<-diag(nug1,nrow = nrow(dist.mat),ncol=ncol(dist.mat))
      NUG2<-diag(nug2,nrow = nrow(dist.mat),ncol=ncol(dist.mat))
      
      ############## Inverting C11 ##########
      C<-rbind(cbind(COV11op+NUG1,COV12op),cbind(t(COV12op),COV22op+NUG2))
      Cchol<-chol(C)
      Cinv<-chol2inv(Cchol)
      logD<-determinant(C)$modulus
      
      nloglikelihood <-(0.5 * logD + 0.5 * t(z) %*% Cinv %*% z+0.5*length(z)*log(2*pi))
      if(abs(nloglikelihood) == Inf || is.nan(nloglikelihood)){ nloglikelihood <- 1e+08}
      return(list(mlv=nloglikelihood,a1=a1,a2=a2,nu1=nu1,nu2=nu2,sigma1=sigma1,sigma2=sigma2,a12=a12,nu12=nu12,rho12=rho12,full.cov=C))
    }
    
  }
  
  
  init.bvm<-c(indmat.estim$par[1:6],0,0,0,indmat.estim$par[7:8])
  
  mle_bvm_mlv<-function(pars)
  {
    return(mle_bvm(p=pars,z=c(un.grd.train$PM2_5,un.grd.train$WS),dmat.ml=dist.mat.train,index.mat.ml=index.mat.train,index.val.ml=index.val.train,uniq.dist.ml=uniq.dist.train)$mlv)
  }
  
  
  optim_bvm_loglik <- function(par){
    optim(par=par,
          fn = mle_bvm_mlv,
          hessian=FALSE,
          control=list(trace=6,
                       pgtol=0,
                       parscale=rep(0.1,length(par)),
                       maxit=10000))
  }
  
  bvm.estim<-optim_bvm_loglik(init.bvm)
  
  
  
  lmc.loglikelihood_allcomponents<-function(p,z,dmat.ml,index.mat.ml,index.val.ml,uniq.dist.ml){
    
    
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
    nug1<-theta[9]
    nug2<-theta[10]
    ######## Putting hard constraints on the parameters #############
    if( a1<=0 | nu1<=0 | sigma1<=0 | a2<=0 | nu2<=0 | sigma2<=0 |nug1<0|nug2<0)
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
      
      NUG1<-diag(nug1,nrow = nrow(dist.mat),ncol=ncol(dist.mat))
      NUG2<-diag(nug2,nrow = nrow(dist.mat),ncol=ncol(dist.mat))
      ##################################################################################
      ############## Creating full covariance matrix ###################################
      ##################################################################################
      COV11op<-COV22op<-COV12op<-matrix(NA,nrow = nrow(dist.mat),ncol=ncol(dist.mat))
      
      COV11op[index.mat.ml]<-COV11[index.val.ml]
      COV22op[index.mat.ml]<-COV22[index.val.ml]
      COV12op[index.mat.ml]<-COV12[index.val.ml]
      
      
      
      ############## Inverting C11 ##########
      C<-rbind(cbind(COV11op+NUG1,COV12op),cbind(t(COV12op),COV22op+NUG2))
      
      
      
      Cchol<-chol(C)
      Cinv<-chol2inv(Cchol)
      logD<-determinant(C)$modulus
      
      nloglikelihood <-(0.5 * logD + 0.5 * t(z) %*% Cinv %*% z+0.5*length(z)*log(2*pi))
      if(abs(nloglikelihood) == Inf || is.nan(nloglikelihood)){ nloglikelihood <- 1e+08}
      return(list(mlv=nloglikelihood,full.cov=C))
    }
  }
  
  init.lmc<-c(indmat.estim$par[c(1,2,4,5)],indmat.estim$par[3],0,0,indmat.estim$par[6],indmat.estim$par[7:8])
  ##########################################################################################
  ##### Now we write the code for log_likelihood of Linear model of coregionalization ######
  ##########################################################################################
  
  
  lmc.loglikelihood<-function(par)
  {
    
    return(lmc.loglikelihood_allcomponents(p=par,z=c(un.grd.train$PM2_5,un.grd.train$WS),dmat.ml=dist.mat.train,index.mat.ml=index.mat.train,index.val.ml=index.val.train,uniq.dist.ml=uniq.dist.train)$mlv)
  }
  
  
  
  
  ### Finding mle parametrs for lmc model ####
  optim_lmc.loglikelihood <- function(par){
    optim(par=par,
          fn = lmc.loglikelihood,
          hessian=FALSE,
          control=list(trace=6,
                       pgtol=0,
                       parscale=rep(0.1,length(par)),
                       maxit=10000))
  }
  
  
  
  fit.Model.lmc <- optim_lmc.loglikelihood(par=init.lmc)
  
  
  rt<-list(ind=indmat.estim,bvm=bvm.estim,lmc=fit.Model.lmc)
  
  
  rt
  
  
  
}


save.image("Estimates from all the other candidate models.RData")



















