#########################################################################
###################### Part 2.5 of the Data application code ############
#########################################################################

##############################################################################################################
##### To produce some intermidiate files that will be required for the semi-parametric model estimation ######
##############################################################################################################
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
ncores<-detectCores()-2
registerDoParallel(cores = ncores)
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




###########################################
####### Creating intermediate files #######
###########################################

uniq.dist.train.f.list<-list()
dist.mat.train.f.list<-list()
index.mat.train.f.list<-list()
index.val.train.f.list<-list()


aux_data<-foreach(run=1:100) %dopar% {
  
  
  un.grd.train<-un.grd.total[-rand.index[,run],]  ########## Training set ######### 
  un.grd.test<-un.grd.total[rand.index[,run],]    ########## Test set #############
  dist.mat.train<-rdist(un.grd.train[,-c(3,4)])
  uniq.dist.train<-unique(c(dist.mat.train))
  uniq.dist.train<-sort(uniq.dist.train)  ##### sorting distances in increasing order #####
  
  #dist.mat.train.f.list[[i]]<-dist.mat.train
  #uniq.dist.train.f.list[[i]]<-uniq.dist.train
  index.mat.train<-NULL
  index.val.train<-NULL
  
  for(kk in 1:length(uniq.dist.train))
  {
    index.mat.train<-rbind(index.mat.train,which(dist.mat.train==uniq.dist.train[kk],arr.ind = T))
    index.val.train<-c(index.val.train,rep(kk,times=nrow(which(dist.mat.train==uniq.dist.train[kk],arr.ind = T))))
  }
  #index.mat.train.f.list[[i]]<-index.mat.train
  #index.val.train.f.list[[i]]<-index.val.train
  
  rtvalue<-list(dist.mat.train.f.list=dist.mat.train,
                uniq.dist.train.f.list=uniq.dist.train,
                index.mat.train.f.list=index.mat.train,
                index.val.train.f.list=index.val.train)
  
  rtvalue
  
}

for(i in 1:100)
{
  uniq.dist.train.f.list[[i]]<-aux_data[[i]]$uniq.dist.train.f.list
  dist.mat.train.f.list[[i]]<-aux_data[[i]]$dist.mat.train.f.list
  index.mat.train.f.list[[i]]<-aux_data[[i]]$index.mat.train.f.list
  index.val.train.f.list[[i]]<-aux_data[[i]]$index.val.train.f.list
  
  
}


save.image("aux_var_2_5.RData")






