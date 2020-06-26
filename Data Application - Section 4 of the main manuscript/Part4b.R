#######################################################################
###################### Part 4b of the Data application code ############
#######################################################################

###################################################################################################
######### In this part, we estimate Semi parametric model with Delta 4 (for i=39:76)     ###########
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
#Sys.sleep(18000)
load("Final dataset to be worked with.RData")
load("aux_var_2_5.RData")
load("Estimates from all the other candidate models.RData")
o_candi.estims<-ind_bvm_lmc_estimations
#######################################
####### Setting number of cores #######
#######################################
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


ncores<-detectCores()-2
registerDoParallel(cores = ncores)



############## Likelihood Estimation ##############



full.cov.compute2<-function(f1.f,f2.f,f12.f,u.f,index.mat.f,index.val.f,Bes_mat.f,uniq.dist.f,sigma1.f,sigma2.f,dmat.f,nug1.f,nug2.f)
{
  
  C11n<-colSums(Bes_mat.f*matrix(rep(2*pi*u*f1.f,length(uniq.dist.f)),
                                 ncol=length(uniq.dist.f),
                                 byrow=F))
  scl1<-max(C11n)
  C11n<-(C11n/scl1)*(sigma1.f^2)
  
 
  C22n<-colSums( Bes_mat.f*matrix(rep(2*pi*u*f2.f,length(uniq.dist.f)),
                                  ncol=length(uniq.dist.f),
                                  byrow=F))
  scl2<-max(C22n)
  C22n<-(C22n/scl2)*(sigma2.f^2)
  
  
  C12n<-colSums(Bes_mat.f*matrix(rep(2*pi*u*f12.f,length(uniq.dist.f)),
                                 ncol=length(uniq.dist.f),
                                 byrow=F))
  C12n<-(C12n/(sqrt(scl1*scl2)))*(sigma1.f*sigma2.f)
  dist.mat<-dmat.f 
  COV11op<-COV22op<-COV12op<-matrix(NA,nrow = nrow(dist.mat),ncol=ncol(dist.mat))
  NUG1<-diag(nug1.f,nrow = nrow(dist.mat),ncol=ncol(dist.mat))
  NUG2<-diag(nug2.f,nrow = nrow(dist.mat),ncol=ncol(dist.mat))
  
  COV11op[index.mat.f]<-C11n[index.val.f]
  COV22op[index.mat.f]<-C22n[index.val.f]
  COV12op[index.mat.f]<-C12n[index.val.f]
  myC4<-rbind(cbind(COV11op+NUG1,COV12op),cbind(t(COV12op),COV22op+NUG2))
  return(myC4)
  
}




#################################################################################
#################################################################################
##### Computing the test and train samples plus all other required matrices #####
##### before we attempt to optimize #############################################
#################################################################################
#################################################################################
bcoeff<-function(x)
{
  
  temp<-sin(x)
  temp[temp==-1]<--1+1e-10
  temp[temp==1]<-1-1e-10
  return(temp)
}

#### Random splitting index ############


freq.max<-9
n.nodes<-500
u<-seq(0,freq.max,length=n.nodes)

full.dist<-rdist(un.grd.total[,-c(3,4)])
full.uniq.dist<-unique(c(full.dist))
full.uniq.dist<-sort(full.uniq.dist)  ##### sorting distances in increasing order #####

full.theta<-outer(u,full.uniq.dist,"*")

full.bessel<-besselJ(x=full.theta,nu=0)

####### Computing indexes to be chosen from columns of full. bessel matrix in the ith run ####
bessel.index<-list()
for(i in 1:100)
{
  tempseq<-1:length(full.uniq.dist)
  temp.index<-tempseq[full.uniq.dist%in%uniq.dist.train.f.list[[i]]]
  bessel.index[[i]]<-temp.index
}


#####################################################
###### Estimation of semi-parametric model ##########
#####################################################

sp2esti<-foreach(i=39:76)%dopar%{
  library(fields)
  freq.max<-9
  n.nodes<-500
  u<-seq(0,freq.max,length=n.nodes)
  
  un.grd.train<-un.grd.total[-rand.index[,i],]  ########## Training set ######### 
  un.grd.test<-un.grd.total[rand.index[,i],]    ########## Test set #############
  dist.mat.train<-dist.mat.train.f.list[[i]]
  uniq.dist.train<-uniq.dist.train.f.list[[i]]  ##### sorting distances in increasing order #####
  
  
  bcoeff<-function(x)
  {
    
    temp<-sin(x)
    temp[temp==-1]<--1+1e-10
    temp[temp==1]<-1-1e-10
    return(temp)
  }
  
  mle_spd4<-function(p,z,dmat.ml,Bes_mat.ml,index.mat.ml,index.val.ml,u.ml,uniq.dist.ml)
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
    nug1<-p[13]
    nug2<-p[14]
    if(sum(p[1:6]<0)!=0||nug1<0||nug2<0)
    {
      nloglikelihood<-10000000
      return(list(mlv=nloglikelihood,params=NULL))
    }
    else
    { 
      f.var1<-f.matern(w=u, nu=nu1, sigma = sigma1, a=a1, d=2)
      f.var2<-f.matern(w=u, nu=nu2, sigma = sigma2, a=a2, d=2)
      Delta=4
      coh12<-b_3*Bspline(j=-3,k=4,delta = Delta,x=u)+b_2*Bspline(j=-2,k=4,delta = Delta,x=u)+b_1*Bspline(j=-1,k=4,delta = Delta,x=u)+b0*Bspline(j=0,k=4,delta = Delta,x=u)+b1*Bspline(j=1,k=4,delta = Delta,x=u)+b2*Bspline(j=2,k=4,delta = Delta,x=u)#+b3*Bspline(j=3,k=4,delta = Delta,x=u)+b4*Bspline(j=4,k=4,delta = Delta,x=u)
      f.var12<-coh12*sqrt(f.var1*f.var2)
      
      
      
      C<-full.cov.compute2(f1.f=f.var1,f2.f=f.var2,f12.f=f.var12,u.f=u.ml,index.mat.f=index.mat.ml,index.val.f=index.val.ml,Bes_mat.f=Bes_mat.ml,uniq.dist.f=uniq.dist.ml,sigma1.f=sigma1,sigma2.f=sigma2,dmat.f = dmat.ml,nug1.f = nug1,nug2.f = nug2)
      
      
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
  mle_spd4_mlv<-function(pars)
  {
    return(mle_spd4(p=pars,z=c(un.grd.train$PM2_5,un.grd.train$WS),dmat.ml=dist.mat.train,Bes_mat.ml=full.bessel[,bessel.index[[i]]],index.mat.ml=index.mat.train.f.list[[i]],index.val.ml=index.val.train.f.list[[i]],u.ml=u,uniq.dist.ml=uniq.dist.train)$mlv)
  }
  
  ###### Finding optimum initial values ########
  
  
  ####### Finding optimized initial values for semiparametric model #########
  bvm.coh.compute<-function(estim.par)
  {
    
    p<-estim.par
    
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
    a12<-sqrt((a1^2+a2^2)/2+deltaB)
    nu12<-((nu1+nu2)/2)+deltaA
    num1<-beta*(a12^(-2*deltaA-(nu1+nu2)))*gamma(((nu1+nu2)/2)+(2/2))*gamma(nu12)
    den1<-(a1^(-deltaA-(nu1)))*(a2^(-deltaA-(nu2)))*sqrt(gamma(nu1)*gamma(nu2))*gamma(nu12+(2/2))
    rho12<-num1/den1
    
    
    coh12<-biwm_coh(w=u,a1=a1,a2=a2,v1=nu1,v2=nu2,a12=a12,v12=nu12,d=2,rho=rho12)
    ############## Inverting C11 ##########
    #C.train<-C22
    #C.test<-C11
    #C.test.train<-C12
    
    return(coh12)
    
    
  }
  
  tvalue<-bvm.coh.compute(o_candi.estims[[i]]$bvm$par)
  ######## Now we find a set of good initial values on the basis of bivariate matern model #######
  l2dist<-function(p)
  {
    b_3<-bcoeff(p[1])
    b_2<-bcoeff(p[2])
    b_1<-bcoeff(p[3])
    b0<-bcoeff(p[4])
    b1<-bcoeff(p[5])
    b2<-bcoeff(p[6])
    #b3<-bcoeff(p[7])
    #b4<-bcoeff(p[8])
    Delta=4
    coh12<-b_3*Bspline(j=-3,k=4,delta = Delta,x=u)+b_2*Bspline(j=-2,k=4,delta = Delta,x=u)+b_1*Bspline(j=-1,k=4,delta = Delta,x=u)+b0*Bspline(j=0,k=4,delta = Delta,x=u)+b1*Bspline(j=1,k=4,delta = Delta,x=u)+b2*Bspline(j=2,k=4,delta = Delta,x=u)#+b3*Bspline(j=3,k=4,delta = Delta,x=u)+b4*Bspline(j=4,k=4,delta = Delta,x=u)
    rv<-sqrt(sum((coh12-tvalue)^2))
    return(rv)
  }
  
  init.value<-optim(par = c(0,0,0,0,0,0),l2dist,
                    hessian=FALSE,
                    control=list(trace=6,
                                 maxit=10000))
  
  init.spd4<-c(o_candi.estims[[i]]$bvm$par[1:6],init.value$par,o_candi.estims[[i]]$bvm$par[10:11])
  
  
  optim_spd4_loglik <- function(par){
    optim(par=par,
          fn = mle_spd4_mlv,
          hessian=FALSE,
          control=list(trace=6,
                       pgtol=0,
                       parscale=rep(0.1,length(par)),
                       maxit=10000))
  }
  spd4.estim<-optim_spd4_loglik(init.spd4)
  
  spd4.estim
  
}
rm(full.bessel,full.theta)
save.image("p4b.RData")







