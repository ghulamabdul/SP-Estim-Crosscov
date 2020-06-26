############################################################################################################################
############ Here we compile the estimated models to evaluate prediction performance and coherence functions ###############
############################################################################################################################

###############################
##### Loading Libraries #######
###############################
library(fields)
library(doParallel)
library(scoringRules)
###############################
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

bcoeff<-function(x)
{
  
  temp<-sin(x)
  temp[temp==-1]<--1+1e-10
  temp[temp==1]<-1-1e-10
  return(temp)
}



###### Setting working directory ######
## set the directory where the file "Final dataset to be worked with.RData" is saved from part 1

setwd("/Users/qadirga/Documents/Project 2/Biometrics submission/Revision work/Revised Data application/Different spatial regions/WNC_Final version/WNC")

###### Loading RData files ########
load("aux_var_2_5.RData")


load("p3a.RData")
part1<-sp2esti
load("p3b.RData")
part2<-sp2esti
load("p3c.RData")
part3<-sp2esti



sp2esti<-c(part1,part2,part3)
spd2esti<-sp2esti

load("p5a.RData")
part1<-sp2esti
load("p5b.RData")
part2<-sp2esti
load("p5c.RData")
part3<-sp2esti


spd5esti<-c(part1,part2,part3)

load("p4a.RData")
part1<-sp2esti
load("p4b.RData")
part2<-sp2esti
load("p4c.RData")
part3<-sp2esti
spd4esti<-c(part1,part2,part3)

load("Estimates from all the other candidate models.RData")

temp1<-temp2<-temp3<-temp4<-temp5<-temp6<-numeric()


############ Now we compute the loglikelihood #######
logl.ind<-logl.bvm<-logl.lmc<-logl.spd2<-logl.spd5<-logl.spd4<-numeric(length = 100)
for(i in 1:100)
{
  logl.ind[i]<--ind_bvm_lmc_estimations[[i]]$ind$value
  logl.bvm[i]<--ind_bvm_lmc_estimations[[i]]$bvm$value
  logl.lmc[i]<--ind_bvm_lmc_estimations[[i]]$lmc$value
  logl.spd2[i]<--spd2esti[[i]]$value
  logl.spd5[i]<--spd5esti[[i]]$value
  logl.spd4[i]<--spd4esti[[i]]$value
  
  temp1[i]<-ind_bvm_lmc_estimations[[i]]$ind$counts[1]
  temp2[i]<-ind_bvm_lmc_estimations[[i]]$bvm$counts[1]
  temp3[i]<-ind_bvm_lmc_estimations[[i]]$lmc$counts[1]
  temp4[i]<-spd2esti[[i]]$counts[1]
  temp5[i]<-spd4esti[[i]]$counts[1]
  temp6[i]<-spd5esti[[i]]$counts[1]
}


##### Plotting loglikelihood values ######
par(mfrow=c(1,1))
boxplot(logl.ind,logl.bvm,logl.lmc,logl.spd2,logl.spd4,logl.spd5,names = c("ind","bvm","lmc","spd2","spd4","spd5"))
mean(logl.ind)
mean(logl.bvm)
mean(logl.lmc)
mean(logl.spd2)
mean(logl.spd4)
mean(logl.spd5)

###### Plotting the AIC values ######
aic.ind<-2*length(ind_bvm_lmc_estimations[[i]]$ind$par)-2*(logl.ind)
aic.bvm<-2*length(ind_bvm_lmc_estimations[[i]]$bvm$par)-2*(logl.bvm)
aic.lmc<-2*length(ind_bvm_lmc_estimations[[i]]$lmc$par)-2*(logl.lmc)
aic.spd2<-2*length(spd2esti[[i]]$par)-2*(logl.spd2)
aic.spd5<-2*length(spd5esti[[i]]$par)-2*(logl.spd5)
aic.spd4<-2*length(spd4esti[[i]]$par)-2*(logl.spd4)


boxplot(aic.ind,aic.bvm,aic.lmc,aic.spd2,aic.spd4,aic.spd5,names = c("ind","bvm","lmc","spd2","spd4","spd5"))

mean(aic.ind)
mean(aic.bvm)
mean(aic.lmc)
mean(aic.spd2)
mean(aic.spd4)
mean(aic.spd5)



model.fit.summary<-data.frame(Model=c("IND","BVM","LMC","SPD2","SPD4","SPD5"),
                              params=c(length(ind_bvm_lmc_estimations[[i]]$ind$par),
                                       length(ind_bvm_lmc_estimations[[i]]$bvm$par),
                                       length(ind_bvm_lmc_estimations[[i]]$lmc$par),
                                       length(spd2esti[[i]]$par),
                                       length(spd4esti[[i]]$par),
                                       length(spd5esti[[i]]$par)),
                              LogL=round(c(mean(logl.ind),
                                           mean(logl.bvm),
                                           mean(logl.lmc),
                                           mean(logl.spd2),
                                           mean(logl.spd4),
                                           mean(logl.spd5)),2),
                              LogL_SE=round(c(sd(logl.ind),
                                              sd(logl.bvm),
                                              sd(logl.lmc),
                                              sd(logl.spd2),
                                              sd(logl.spd4),
                                              sd(logl.spd5)),2),
                              AIC=round(c(mean(aic.ind),
                                          mean(aic.bvm),
                                          mean(aic.lmc),
                                          mean(aic.spd2),
                                          mean(aic.spd4),
                                          mean(aic.spd5)),2),
                              
                              AIC_SE=round(c(sd(aic.ind),
                                             sd(aic.bvm),
                                             sd(aic.lmc),
                                             sd(aic.spd2),
                                             sd(aic.spd4),
                                             sd(aic.spd5)),2)
                              
                              
)

model.fit.summary





########################################################
######### Now we do the prediction performance #########
########################################################

dist.mat<-rdist(un.grd.total[,-c(3,4)])
uniq.dist<-unique(c(dist.mat))
uniq.dist<-sort(uniq.dist)  ##### sorting distances in increasing order #####
pb = txtProgressBar(min = 0, max = length(uniq.dist), initial = 0,style = 3)
index.mat<-NULL
index.val<-NULL
for(i in 1:length(uniq.dist))
{
  index.mat<-rbind(index.mat,which(dist.mat==uniq.dist[i],arr.ind = T))
  index.val<-c(index.val,rep(i,times=nrow(which(dist.mat==uniq.dist[i],arr.ind = T))))
  setTxtProgressBar(pb,i)
}

save.image("before part4.RData")



ind.pred.summary<-function(estim.par,seed)
{
  
  
  train.data<-un.grd.total[-rand.index[,seed],]  ########## Training set ######### 
  test.data<-un.grd.total[rand.index[,seed],]    ########## Test set #############
  full.dist<-rdist(un.grd.total[,-c(3,4)])
  p<-estim.par
  
  a1<-p[1]
  nu1<-p[2]
  sigma1<-p[3]
  a2<-p[4]
  nu2<-p[5]
  sigma2<-p[6]
  nug1<-p[7]
  nug2<-p[8]
  
  
  C11<-my.matern(h=full.dist,a=a1,sigma = sigma1,nu=nu1)
  C22<-my.matern(h=full.dist,a=a2,sigma = sigma2,nu=nu2)
  NUG1<-diag(nug1,nrow = nrow(C11),ncol=ncol(C11))
  NUG2<-diag(nug2,nrow = nrow(C22),ncol=ncol(C22))
  C12<-matrix(0,nrow = nrow(C11),ncol=ncol(C11))
  
  C<-rbind(cbind(C11+NUG1,C12),cbind(t(C12),C22+NUG2))
  
  test.indexes<-c(rand.index[,seed],length(un.grd.total$lon)+rand.index[,seed])
  
  ##################################################################################
  ############## Creating full covariance matrix ###################################
  ##################################################################################
  C.test<-C[test.indexes,test.indexes]
  C.train<-C[-test.indexes,-test.indexes]
  C.test.train<-C[test.indexes,-test.indexes]
  
  ############## Inverting C11 ##########
  #C.train<-C22
  #C.test<-C11
  #C.test.train<-C12
  
  
  
  
  
  
  prediction<-C.test.train%*%solve(C.train)%*%c(train.data$PM2_5,train.data$WS) #Conditional mean of a Multivariate Gaussian
  
  validation.vector<-c(test.data$PM2_5,test.data$WS)
  
  bivar.vaiance<-diag(C.test-C.test.train%*%solve(C.train)%*%t(C.test.train)) #Kriging variance
  
  
  ########### Computing bivariate prediction scores ############
  rmse2<-sqrt(mean((validation.vector-c(prediction))^2)) #Root mean squared error
  mae2<-mean(abs(validation.vector-c(prediction))) #Mean Absolute error
  nmse2<-1-(sum((validation.vector-c(prediction))^2))/(sum((validation.vector-mean(validation.vector))^2)) #Normalized mean squared error
  crps2<-mean(crps(validation.vector, "norm", mean = c(prediction), sd = c(sqrt(bivar.vaiance)))) #Mean CRPS
  logs2<-mean(logs(validation.vector, "norm", mean = c(prediction), sd = c(sqrt(bivar.vaiance)))) #Mean LogS
  return(c(rmse2,mae2,nmse2,crps2,logs2))
}




bvm.pred.summary<-function(estim.par,seed)
{
  
  train.data<-un.grd.total[-rand.index[,seed],]  ########## Training set ######### 
  test.data<-un.grd.total[rand.index[,seed],]    ########## Test set #############
  full.dist<-rdist(un.grd.total[,-c(3,4)])
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
  
  
  C11<-my.matern(h=full.dist,a=a1,sigma = sigma1,nu=nu1)
  C22<-my.matern(h=full.dist,a=a2,sigma = sigma2,nu=nu2)
  a12<-sqrt((a1^2+a2^2)/2+deltaB)
  nu12<-((nu1+nu2)/2)+deltaA
  num1<-beta*(a12^(-2*deltaA-(nu1+nu2)))*gamma(((nu1+nu2)/2)+(2/2))*gamma(nu12)
  den1<-(a1^(-deltaA-(nu1)))*(a2^(-deltaA-(nu2)))*sqrt(gamma(nu1)*gamma(nu2))*gamma(nu12+(2/2))
  rho12<-num1/den1
  C12<-sigma1*sigma2*rho12*my.matern(h=full.dist,a=a12,sigma = 1,nu=nu12)
  
  
  NUG1<-diag(nug1,nrow = nrow(C11),ncol=ncol(C11))
  NUG2<-diag(nug2,nrow = nrow(C22),ncol=ncol(C22))
  
  ############## Inverting C11 ##########
  C<-rbind(cbind(C11+NUG1,C12),cbind(t(C12),C22+NUG2))
  
  
  test.indexes<-c(rand.index[,seed],length(un.grd.total$lon)+rand.index[,seed])
  
  ##################################################################################
  ############## Creating full covariance matrix ###################################
  ##################################################################################
  C.test<-C[test.indexes,test.indexes]
  C.train<-C[-test.indexes,-test.indexes]
  C.test.train<-C[test.indexes,-test.indexes]
  coh12<-biwm_coh(w=u,a1=a1,a2=a2,v1=nu1,v2=nu2,a12=a12,v12=nu12,d=2,rho=rho12)
  ############## Inverting C11 ##########
  #C.train<-C22
  #C.test<-C11
  #C.test.train<-C12
  
  
  
  
  
  
  prediction<-C.test.train%*%solve(C.train)%*%c(train.data$PM2_5,train.data$WS) #Conditional mean of a Multivariate Gaussian
  
  validation.vector<-c(test.data$PM2_5,test.data$WS)
  
  bivar.vaiance<-diag(C.test-C.test.train%*%solve(C.train)%*%t(C.test.train)) #Kriging variance
  
  
  ########### Computing bivariate prediction scores ############
  rmse2<-sqrt(mean((validation.vector-c(prediction))^2)) #Root mean squared error
  mae2<-mean(abs(validation.vector-c(prediction))) #Mean Absolute error
  nmse2<-1-(sum((validation.vector-c(prediction))^2))/(sum((validation.vector-mean(validation.vector))^2)) #Normalized mean squared error
  crps2<-mean(crps(validation.vector, "norm", mean = c(prediction), sd = c(sqrt(bivar.vaiance)))) #Mean CRPS
  logs2<-mean(logs(validation.vector, "norm", mean = c(prediction), sd = c(sqrt(bivar.vaiance)))) #Mean LogS
  return(list(scores=c(rmse2,mae2,nmse2,crps2,logs2),coh=coh12))
  
  
}

# function for the coherence of the lmc model

lmc.coh<-function(w,a1,nu1,a2,nu2,b11,b12,b21,b22)
{
  f1<-f.matern(w=w,nu=nu1,sigma = 1,a=a1,d=2)
  f2<-f.matern(w=w,nu=nu2,sigma = 1,a=a2,d=2)
  num<-b11*b21*f1+b12*b22*f2
  den1<-sqrt((b11^2)*f1+(b12^2)*f2)
  den2<-sqrt((b21^2)*f1+(b22^2)*f2)
  den<-den1*den2
  return(num/den)  
}

lmc.pred.summary<-function(estim.par,seed)
{
  
  train.data<-un.grd.total[-rand.index[,seed],]  ########## Training set ######### 
  test.data<-un.grd.total[rand.index[,seed],]    ########## Test set #############
  full.dist<-rdist(un.grd.total[,-c(3,4)])
  p<-estim.par
  
  
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
  #n<-nrow(dist.mat)
  C11<-my.matern(h=full.dist,a=a1,nu=nu1,sigma = sigma1)
  C22<-my.matern(h=full.dist,a=a2,nu=nu2,sigma = sigma2)
  
  COV11<-(b11^2)*C11+(b12^2)*C22
  COV22<-(b21^2)*C11+(b22^2)*C22
  COV12<-(b11*b21)*C11+(b12*b22)*C22
  
  NUG1<-diag(nug1,nrow = nrow(COV11),ncol=ncol(COV11))
  NUG2<-diag(nug2,nrow = nrow(COV22),ncol=ncol(COV22))
  ##################################################################################
  ############## Creating full covariance matrix ###################################
  ##################################################################################
  
  
  
  ############## Inverting C11 ##########
  C<-rbind(cbind(COV11+NUG1,COV12),cbind(t(COV12),COV22+NUG2))
  
  
  
  test.indexes<-c(rand.index[,seed],length(un.grd.total$lon)+rand.index[,seed])
  
  ##################################################################################
  ############## Creating full covariance matrix ###################################
  ##################################################################################
  C.test<-C[test.indexes,test.indexes]
  C.train<-C[-test.indexes,-test.indexes]
  C.test.train<-C[test.indexes,-test.indexes]
  coh12<-lmc.coh(w=u,a1=a1,nu1=nu1,a2=a2,nu2=nu2,b11=b11,b12=b12,b21=b21,b22=b22)
  ############## Inverting C11 ##########
  #C.train<-C22
  #C.test<-C11
  #C.test.train<-C12
  
  
  
  
  
  
  prediction<-C.test.train%*%solve(C.train)%*%c(train.data$PM2_5,train.data$WS) #Conditional mean of a Multivariate Gaussian
  
  validation.vector<-c(test.data$PM2_5,test.data$WS)
  
  bivar.vaiance<-diag(C.test-C.test.train%*%solve(C.train)%*%t(C.test.train)) #Kriging variance
  
  
  ########### Computing bivariate prediction scores ############
  rmse2<-sqrt(mean((validation.vector-c(prediction))^2)) #Root mean squared error
  mae2<-mean(abs(validation.vector-c(prediction))) #Mean Absolute error
  nmse2<-1-(sum((validation.vector-c(prediction))^2))/(sum((validation.vector-mean(validation.vector))^2)) #Normalized mean squared error
  crps2<-mean(crps(validation.vector, "norm", mean = c(prediction), sd = c(sqrt(bivar.vaiance)))) #Mean CRPS
  logs2<-mean(logs(validation.vector, "norm", mean = c(prediction), sd = c(sqrt(bivar.vaiance)))) #Mean LogS
  return(list(scores=c(rmse2,mae2,nmse2,crps2,logs2),coh=coh12))
  
  
}

full.dist<-rdist(un.grd.total[,-c(3,4)])
full.uniq.dist<-unique(c(full.dist))
full.uniq.dist<-sort(full.uniq.dist)  ##### sorting distances in increasing order #####

full.theta<-outer(u,full.uniq.dist,"*")

full.bessel<-besselJ(x=full.theta,nu=0)

#sp2d.pred.summary(estim.par = spd2esti[[1]]$par,seed = 1)
sp2d.pred.summary<-function(estim.par,seed)
{
  train.data<-un.grd.total[-rand.index[,seed],]  ########## Training set ######### 
  test.data<-un.grd.total[rand.index[,seed],]    ########## Test set #############
  p<-estim.par
  
  
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
  nug1<-p[15]
  nug2<-p[16]
  
  
  f.var1<-f.matern(w=u, nu=nu1, sigma = sigma1, a=a1, d=2)
  f.var2<-f.matern(w=u, nu=nu2, sigma = sigma2, a=a2, d=2)
  Delta=2
  coh12<-b_3*Bspline(j=-3,k=4,delta = Delta,x=u)+b_2*Bspline(j=-2,k=4,delta = Delta,x=u)+b_1*Bspline(j=-1,k=4,delta = Delta,x=u)+b0*Bspline(j=0,k=4,delta = Delta,x=u)+b1*Bspline(j=1,k=4,delta = Delta,x=u)+b2*Bspline(j=2,k=4,delta = Delta,x=u)+b3*Bspline(j=3,k=4,delta = Delta,x=u)+b4*Bspline(j=4,k=4,delta = Delta,x=u)
  f.var12<-coh12*sqrt(f.var1*f.var2)
  
  
  C<-full.cov.compute2(f1.f=f.var1,f2.f=f.var2,f12.f=f.var12,u.f=u,index.mat.f=index.mat,index.val.f=index.val,Bes_mat.f=full.bessel,uniq.dist.f=uniq.dist,sigma1.f=sigma1,sigma2.f=sigma2,dmat.f = dist.mat,nug1.f = nug1,nug2.f = nug2)
  chol(C)
  #C22<-full.cov.compute3(f1.f=f.var1,f2.f=f.var2,f12.f=f.var12,u.f=u,index.mat.f=index.mat.train.f.list[[seed]],index.val.f=index.val.train.f.list[[seed]],Bes_mat.f=full.bessel[,bessel.index[[seed]]],uniq.dist.f=uniq.dist.train.f.list[[seed]],sigma1.f=sigma1,sigma2.f=sigma2,dmat.f = dist.mat.train.f.list[[seed]],nug1.f = nug1,nug2.f = nug2)
  
  #C11<-full.cov.compute3(f1.f=f.var1,f2.f=f.var2,f12.f=f.var12,u.f=u,index.mat.f=index.mat.test.f.list[[seed]],index.val.f=index.val.test.f.list[[seed]],Bes_mat.f=full.bessel[,bessel.test.index[[seed]]],uniq.dist.f=uniq.dist.test.f.list[[seed]],sigma1.f=sigma1,sigma2.f=sigma2,dmat.f = dist.mat.test.f.list[[seed]],nug1.f = nug1,nug2.f = nug2)
  #C12<-full.cov.compute3(f1.f=f.var1,f2.f=f.var2,f12.f=f.var12,u.f=u,index.mat.f=index.mat.test.train.f.list[[seed]][[1]],index.val.f=index.val.test.train.f.list[[seed]][[1]],Bes_mat.f=full.bessel[,bessel.test.train.index[[seed]]],uniq.dist.f=uniq.dist.test.train.f.list[[seed]][[1]],sigma1.f=sigma1,sigma2.f=sigma2,dmat.f = dist.mat.test.train.f.list[[seed]][[1]],nug1.f = 0,nug2.f = 0)
  test.indexes<-c(rand.index[,seed],length(un.grd.total$lon)+rand.index[,seed])
  
  ##################################################################################
  ############## Creating full covariance matrix ###################################
  ##################################################################################
  C.test<-C[test.indexes,test.indexes]
  C.train<-C[-test.indexes,-test.indexes]
  C.test.train<-C[test.indexes,-test.indexes]
  
  ############## Inverting C11 ##########
  #C.train<-C22
  #C.test<-C11
  #C.test.train<-C12
  
  
  
  
  
  
  prediction<-C.test.train%*%solve(C.train)%*%c(train.data$PM2_5,train.data$WS) #Conditional mean of a Multivariate Gaussian
  
  validation.vector<-c(test.data$PM2_5,test.data$WS)
  
  bivar.vaiance<-diag(C.test-C.test.train%*%solve(C.train)%*%t(C.test.train)) #Kriging variance
  
  
  ########### Computing bivariate prediction scores ############
  rmse2<-sqrt(mean((validation.vector-c(prediction))^2)) #Root mean squared error
  mae2<-mean(abs(validation.vector-c(prediction))) #Mean Absolute error
  nmse2<-1-(sum((validation.vector-c(prediction))^2))/(sum((validation.vector-mean(validation.vector))^2)) #Normalized mean squared error
  crps2<-mean(crps(validation.vector, "norm", mean = c(prediction), sd = c(sqrt(bivar.vaiance)))) #Mean CRPS
  logs2<-mean(logs(validation.vector, "norm", mean = c(prediction), sd = c(sqrt(bivar.vaiance)))) #Mean LogS
  return(list(scores=c(rmse2,mae2,nmse2,crps2,logs2),coh=coh12))
  
  
}


sp5d.pred.summary<-function(estim.par,seed)
{
  train.data<-un.grd.total[-rand.index[,seed],]  ########## Training set ######### 
  test.data<-un.grd.total[rand.index[,seed],]    ########## Test set #############
  p<-estim.par
  
  
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
  nug1<-p[12]
  nug2<-p[13]
  
  
  f.var1<-f.matern(w=u, nu=nu1, sigma = sigma1, a=a1, d=2)
  f.var2<-f.matern(w=u, nu=nu2, sigma = sigma2, a=a2, d=2)
  Delta=5
  coh12<-b_3*Bspline(j=-3,k=4,delta = Delta,x=u)+b_2*Bspline(j=-2,k=4,delta = Delta,x=u)+b_1*Bspline(j=-1,k=4,delta = Delta,x=u)+b0*Bspline(j=0,k=4,delta = Delta,x=u)+b1*Bspline(j=1,k=4,delta = Delta,x=u)#+b2*Bspline(j=2,k=4,delta = Delta,x=u)+b3*Bspline(j=3,k=4,delta = Delta,x=u)+b4*Bspline(j=4,k=4,delta = Delta,x=u)
  f.var12<-coh12*sqrt(f.var1*f.var2)
  
  
  C<-full.cov.compute2(f1.f=f.var1,f2.f=f.var2,f12.f=f.var12,u.f=u,index.mat.f=index.mat,index.val.f=index.val,Bes_mat.f=full.bessel,uniq.dist.f=uniq.dist,sigma1.f=sigma1,sigma2.f=sigma2,dmat.f = dist.mat,nug1.f = nug1,nug2.f = nug2)
  #chol(C)
  #C22<-full.cov.compute3(f1.f=f.var1,f2.f=f.var2,f12.f=f.var12,u.f=u,index.mat.f=index.mat.train.f.list[[seed]],index.val.f=index.val.train.f.list[[seed]],Bes_mat.f=full.bessel[,bessel.index[[seed]]],uniq.dist.f=uniq.dist.train.f.list[[seed]],sigma1.f=sigma1,sigma2.f=sigma2,dmat.f = dist.mat.train.f.list[[seed]],nug1.f = nug1,nug2.f = nug2)
  
  #C11<-full.cov.compute3(f1.f=f.var1,f2.f=f.var2,f12.f=f.var12,u.f=u,index.mat.f=index.mat.test.f.list[[seed]],index.val.f=index.val.test.f.list[[seed]],Bes_mat.f=full.bessel[,bessel.test.index[[seed]]],uniq.dist.f=uniq.dist.test.f.list[[seed]],sigma1.f=sigma1,sigma2.f=sigma2,dmat.f = dist.mat.test.f.list[[seed]],nug1.f = nug1,nug2.f = nug2)
  #C12<-full.cov.compute3(f1.f=f.var1,f2.f=f.var2,f12.f=f.var12,u.f=u,index.mat.f=index.mat.test.train.f.list[[seed]][[1]],index.val.f=index.val.test.train.f.list[[seed]][[1]],Bes_mat.f=full.bessel[,bessel.test.train.index[[seed]]],uniq.dist.f=uniq.dist.test.train.f.list[[seed]][[1]],sigma1.f=sigma1,sigma2.f=sigma2,dmat.f = dist.mat.test.train.f.list[[seed]][[1]],nug1.f = 0,nug2.f = 0)
  test.indexes<-c(rand.index[,seed],length(un.grd.total$lon)+rand.index[,seed])
  
  ##################################################################################
  ############## Creating full covariance matrix ###################################
  ##################################################################################
  C.test<-C[test.indexes,test.indexes]
  C.train<-C[-test.indexes,-test.indexes]
  C.test.train<-C[test.indexes,-test.indexes]
  
  ############## Inverting C11 ##########
  #C.train<-C22
  #C.test<-C11
  #C.test.train<-C12
  
  
  
  
  
  
  prediction<-C.test.train%*%solve(C.train)%*%c(train.data$PM2_5,train.data$WS) #Conditional mean of a Multivariate Gaussian
  
  validation.vector<-c(test.data$PM2_5,test.data$WS)
  
  bivar.vaiance<-diag(C.test-C.test.train%*%solve(C.train)%*%t(C.test.train)) #Kriging variance
  
  
  ########### Computing bivariate prediction scores ############
  rmse2<-sqrt(mean((validation.vector-c(prediction))^2)) #Root mean squared error
  mae2<-mean(abs(validation.vector-c(prediction))) #Mean Absolute error
  nmse2<-1-(sum((validation.vector-c(prediction))^2))/(sum((validation.vector-mean(validation.vector))^2)) #Normalized mean squared error
  crps2<-mean(crps(validation.vector, "norm", mean = c(prediction), sd = c(sqrt(bivar.vaiance)))) #Mean CRPS
  logs2<-mean(logs(validation.vector, "norm", mean = c(prediction), sd = c(sqrt(bivar.vaiance)))) #Mean LogS
  return(list(scores=c(rmse2,mae2,nmse2,crps2,logs2),coh=coh12))
  
  
}

sp4d.pred.summary<-function(estim.par,seed)
{
  train.data<-un.grd.total[-rand.index[,seed],]  ########## Training set ######### 
  test.data<-un.grd.total[rand.index[,seed],]    ########## Test set #############
  p<-estim.par
  
  
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
  
  
  f.var1<-f.matern(w=u, nu=nu1, sigma = sigma1, a=a1, d=2)
  f.var2<-f.matern(w=u, nu=nu2, sigma = sigma2, a=a2, d=2)
  Delta=4
  coh12<-b_3*Bspline(j=-3,k=4,delta = Delta,x=u)+b_2*Bspline(j=-2,k=4,delta = Delta,x=u)+b_1*Bspline(j=-1,k=4,delta = Delta,x=u)+b0*Bspline(j=0,k=4,delta = Delta,x=u)+b1*Bspline(j=1,k=4,delta = Delta,x=u)+b2*Bspline(j=2,k=4,delta = Delta,x=u)#+b3*Bspline(j=3,k=4,delta = Delta,x=u)+b4*Bspline(j=4,k=4,delta = Delta,x=u)
  f.var12<-coh12*sqrt(f.var1*f.var2)
  
  
  C<-full.cov.compute2(f1.f=f.var1,f2.f=f.var2,f12.f=f.var12,u.f=u,index.mat.f=index.mat,index.val.f=index.val,Bes_mat.f=full.bessel,uniq.dist.f=uniq.dist,sigma1.f=sigma1,sigma2.f=sigma2,dmat.f = dist.mat,nug1.f = nug1,nug2.f = nug2)
  #chol(C)
  #C22<-full.cov.compute3(f1.f=f.var1,f2.f=f.var2,f12.f=f.var12,u.f=u,index.mat.f=index.mat.train.f.list[[seed]],index.val.f=index.val.train.f.list[[seed]],Bes_mat.f=full.bessel[,bessel.index[[seed]]],uniq.dist.f=uniq.dist.train.f.list[[seed]],sigma1.f=sigma1,sigma2.f=sigma2,dmat.f = dist.mat.train.f.list[[seed]],nug1.f = nug1,nug2.f = nug2)
  
  #C11<-full.cov.compute3(f1.f=f.var1,f2.f=f.var2,f12.f=f.var12,u.f=u,index.mat.f=index.mat.test.f.list[[seed]],index.val.f=index.val.test.f.list[[seed]],Bes_mat.f=full.bessel[,bessel.test.index[[seed]]],uniq.dist.f=uniq.dist.test.f.list[[seed]],sigma1.f=sigma1,sigma2.f=sigma2,dmat.f = dist.mat.test.f.list[[seed]],nug1.f = nug1,nug2.f = nug2)
  #C12<-full.cov.compute3(f1.f=f.var1,f2.f=f.var2,f12.f=f.var12,u.f=u,index.mat.f=index.mat.test.train.f.list[[seed]][[1]],index.val.f=index.val.test.train.f.list[[seed]][[1]],Bes_mat.f=full.bessel[,bessel.test.train.index[[seed]]],uniq.dist.f=uniq.dist.test.train.f.list[[seed]][[1]],sigma1.f=sigma1,sigma2.f=sigma2,dmat.f = dist.mat.test.train.f.list[[seed]][[1]],nug1.f = 0,nug2.f = 0)
  test.indexes<-c(rand.index[,seed],length(un.grd.total$lon)+rand.index[,seed])
  
  ##################################################################################
  ############## Creating full covariance matrix ###################################
  ##################################################################################
  C.test<-C[test.indexes,test.indexes]
  C.train<-C[-test.indexes,-test.indexes]
  C.test.train<-C[test.indexes,-test.indexes]
  
  ############## Inverting C11 ##########
  #C.train<-C22
  #C.test<-C11
  #C.test.train<-C12
  
  
  
  
  
  
  prediction<-C.test.train%*%solve(C.train)%*%c(train.data$PM2_5,train.data$WS) #Conditional mean of a Multivariate Gaussian
  
  validation.vector<-c(test.data$PM2_5,test.data$WS)
  
  bivar.vaiance<-diag(C.test-C.test.train%*%solve(C.train)%*%t(C.test.train)) #Kriging variance
  
  
  ########### Computing bivariate prediction scores ############
  rmse2<-sqrt(mean((validation.vector-c(prediction))^2)) #Root mean squared error
  mae2<-mean(abs(validation.vector-c(prediction))) #Mean Absolute error
  nmse2<-1-(sum((validation.vector-c(prediction))^2))/(sum((validation.vector-mean(validation.vector))^2)) #Normalized mean squared error
  crps2<-mean(crps(validation.vector, "norm", mean = c(prediction), sd = c(sqrt(bivar.vaiance)))) #Mean CRPS
  logs2<-mean(logs(validation.vector, "norm", mean = c(prediction), sd = c(sqrt(bivar.vaiance)))) #Mean LogS
  return(list(scores=c(rmse2,mae2,nmse2,crps2,logs2),coh=coh12))
  
  
}



pb = txtProgressBar(min = 0, max = 100, initial = 0,style = 3)

ind.psummary<-bvm.psummary<-lmc.psummary<-spd2.psummary<-matrix(NA,nrow = 100,ncol = 5)
spd5.psummary<-matrix(NA,nrow = 100,ncol = 5)
spd4.psummary<-matrix(NA,nrow = 100,ncol = 5)

bvm.coh<-lmc.cohv<-spd2.coh<-matrix(NA,nrow = 100,ncol = 500)
spd5.coh<-matrix(NA,nrow = 100,ncol = 500)
spd4.coh<-matrix(NA,nrow = 100,ncol = 500)





for(i in 1:100) 
{
  t0<-ind.pred.summary(estim.par = ind_bvm_lmc_estimations[[i]]$ind$par,seed = i) 
  t1<-bvm.pred.summary(estim.par = ind_bvm_lmc_estimations[[i]]$bvm$par,seed = i) 
  t2<-lmc.pred.summary(estim.par = ind_bvm_lmc_estimations[[i]]$lmc$par,seed = i) 
  t3<-sp2d.pred.summary(estim.par = spd2esti[[i]]$par,seed = i)
  t4<-sp4d.pred.summary(estim.par = spd4esti[[i]]$par,seed = i)
  t5<-sp5d.pred.summary(estim.par = spd5esti[[i]]$par,seed = i)
  ind.psummary[i,]<-t0
  bvm.psummary[i,]<-t1$scores
  bvm.coh[i,]<-t1$coh
  lmc.psummary[i,]<-t2$scores
  lmc.cohv[i,]<-t2$coh
  spd2.psummary[i,]<-t3$scores
  spd2.coh[i,]<-t3$coh
  spd5.psummary[i,]<-t5$scores
  spd5.coh[i,]<-t5$coh
  spd4.psummary[i,]<-t4$scores
  spd4.coh[i,]<-t4$coh
  setTxtProgressBar(pb,i)
  
}

###############################
###### Plotting results #######
###############################
par(mfrow=c(2,3))

boxplot(ind.psummary[,1],bvm.psummary[,1],lmc.psummary[,1],spd2.psummary[,1],spd4.psummary[,1],spd5.psummary[,1],names = c("ind","bvm","lmc","spd2","spd4","spd5"),main="RMSE")
boxplot(ind.psummary[,2],bvm.psummary[,2],lmc.psummary[,2],spd2.psummary[,2],spd4.psummary[,2],spd5.psummary[,2],names = c("ind","bvm","lmc","spd2","spd4","spd5"),main="MAE")
boxplot(ind.psummary[,3],bvm.psummary[,3],lmc.psummary[,3],spd2.psummary[,3],spd4.psummary[,3],spd5.psummary[,3],names = c("ind","bvm","lmc","spd2","spd4","spd5"),main="NMSE")
boxplot(ind.psummary[,4],bvm.psummary[,4],lmc.psummary[,4],spd2.psummary[,4],spd4.psummary[,4],spd5.psummary[,4],names = c("ind","bvm","lmc","spd2","spd4","spd5"),main="CRPS")
boxplot(ind.psummary[,5],bvm.psummary[,5],lmc.psummary[,5],spd4.psummary[,5],spd4.psummary[,5],spd5.psummary[,5],names = c("ind","bvm","lmc","spd2","spd4","spd5"),main="LogS")


######


pred.summary.df<-as.data.frame(rbind(round(colMeans(ind.psummary),4),
                       round(colMeans(bvm.psummary),4),
                       round(colMeans(lmc.psummary),4),
                       round(colMeans(spd2.psummary),4),
                       round(colMeans(spd4.psummary),4),
                       round(colMeans(spd5.psummary),4)))
colnames(pred.summary.df)<-c("RMSE","MAE","NMSE","mCRPS","mLogS")
pred.summary.df$Model<-c("IND","BVM","LMC","SPD2","SPD4","SPD5")
pr.sum.df<-data.frame(pred.summary.df$Model,pred.summary.df$RMSE,pred.summary.df$MAE,
                      pred.summary.df$NMSE,
                      pred.summary.df$mCRPS,
                      pred.summary.df$mLogS
                      )
colnames(pr.sum.df)<-c("Model","RMSE","MAE","NMSE","mCRPS","mLogS")
pr.sum.df

pred.summary.se.df<-as.data.frame(rbind(round(apply(ind.psummary,2,"sd"),3),
                                     round(apply(bvm.psummary,2,"sd"),3),
                                     round(apply(lmc.psummary,2,"sd"),3),
                                     round(apply(spd2.psummary,2,"sd"),3),
                                     round(apply(spd4.psummary,2,"sd"),3),
                                     round(apply(spd5.psummary,2,"sd"),3)))
colnames(pred.summary.se.df)<-c("RMSE_SE","MAE_SE","NMSE_SE","mCRPS_SE","mLogS_SE")
pred.summary.se.df$Model<-c("IND","BVM","LMC","SPD2","SPD4","SPD5")
pr.sum.se.df<-data.frame(pred.summary.se.df$Model,pred.summary.se.df$RMSE_SE,pred.summary.se.df$MAE_SE,
                      pred.summary.se.df$NMSE_SE,
                      pred.summary.se.df$mCRPS_SE,
                      pred.summary.se.df$mLogS_SE
)
colnames(pr.sum.se.df)<-c("Model","RMSE(SE)","MAE(SE)","NMSE(SE)","mCRPS(SE)","mLogS(SE)")

pr.sum.se.df


model.fit.summary
pr.sum.df
pr.sum.se.df



rm(full.bessel,full.theta)
save.image("Fullfinalanalysis.RData")
par(mfrow=c(2,3))
plot(NA,NA,xlim=c(0,freq.max),ylim=c(-1,1),main="LMC coherence")
for(i in 1:100)
{
  lines(u,lmc.cohv[i,],col="grey")
}
lines(u,colMeans(lmc.cohv),lwd=2,col="red")
abline(h=0)


plot(NA,NA,xlim=c(0,freq.max),ylim=c(-1,1),main="BVM coherence")
for(i in 1:100)
{
  lines(u,bvm.coh[i,],col="grey")
}
lines(u,colMeans(bvm.coh),lwd=2,col="red")
abline(h=0)

plot(NA,NA,xlim=c(0,freq.max),ylim=c(-1,1),main="SPD2 coherence")
for(i in 1:100)
{
  lines(u,spd2.coh[i,],col="grey")
}
lines(u,colMeans(spd2.coh),lwd=2,col="red")
abline(h=0)


plot(NA,NA,xlim=c(0,freq.max),ylim=c(-1,1),main="SPD4 coherence")
for(i in 1:100)
{
  lines(u,spd4.coh[i,],col="grey")
}
lines(u,colMeans(spd4.coh),lwd=2,col="red")
abline(h=0)




plot(NA,NA,xlim=c(0,freq.max),ylim=c(-1,1),main="SPD5 coherence")
for(i in 1:100)
{
  lines(u,spd5.coh[i,],col="grey")
}
lines(u,colMeans(spd5.coh),lwd=2,col="red")
abline(h=0)



library(ggplot2)

bp.rmse.data<-data.frame(RMSE=c(ind.psummary[,1],bvm.psummary[,1],lmc.psummary[,1],spd2.psummary[,1],spd4.psummary[,1],spd5.psummary[,1]),
                         Model=c(rep("IND",times=length(ind.psummary[,1])),
                                 rep("BVM",times=length(bvm.psummary[,1])),
                                 rep("LMC",times=length(lmc.psummary[,1])),
                                 rep("SPD2",times=length(spd2.psummary[,1])),
                                 rep("SPD4",times=length(spd4.psummary[,1])),
                                 rep("SPD5",times=length(spd5.psummary[,1])))
)

ggplot(bp.rmse.data, aes(x=Model, y=RMSE, fill=Model)) + 
  geom_boxplot(alpha=0.8) +
  theme(legend.position="none") +
  scale_fill_brewer(palette="BuPu")


bp.mae.data<-data.frame(MAE=c(ind.psummary[,2],bvm.psummary[,2],lmc.psummary[,2],spd2.psummary[,2],spd4.psummary[,2],spd5.psummary[,2]),
                        Model=c(rep("IND",times=length(ind.psummary[,2])),
                                rep("BVM",times=length(bvm.psummary[,2])),
                                rep("LMC",times=length(lmc.psummary[,2])),
                                rep("SPD2",times=length(spd2.psummary[,2])),
                                rep("SPD4",times=length(spd4.psummary[,2])),
                                rep("SPD5",times=length(spd5.psummary[,2])))
)

ggplot(bp.mae.data, aes(x=Model, y=MAE, fill=Model)) + 
  geom_boxplot(alpha=0.8) +
  theme(legend.position="none") +
  scale_fill_brewer(palette="BuPu")



bp.nmse.data<-data.frame(NMSE=c(ind.psummary[,3],bvm.psummary[,3],lmc.psummary[,3],spd2.psummary[,3],spd4.psummary[,3],spd5.psummary[,3]),
                         Model=c(rep("IND",times=length(ind.psummary[,3])),
                                 rep("BVM",times=length(bvm.psummary[,3])),
                                 rep("LMC",times=length(lmc.psummary[,3])),
                                 rep("SPD2",times=length(spd2.psummary[,3])),
                                 rep("SPD4",times=length(spd4.psummary[,3])),
                                 rep("SPD5",times=length(spd5.psummary[,3])))
)

ggplot(bp.nmse.data, aes(x=Model, y=NMSE, fill=Model)) + 
  geom_boxplot(alpha=0.8) +
  theme(legend.position="none") +
  scale_fill_brewer(palette="BuPu")



bp.crps.data<-data.frame(mCRPS=c(ind.psummary[,4],bvm.psummary[,4],lmc.psummary[,4],spd2.psummary[,4],spd4.psummary[,4],spd5.psummary[,4]),
                         Model=c(rep("IND",times=length(ind.psummary[,4])),
                                 rep("BVM",times=length(bvm.psummary[,4])),
                                 rep("LMC",times=length(lmc.psummary[,4])),
                                 rep("SPD2",times=length(spd2.psummary[,4])),
                                 rep("SPD4",times=length(spd4.psummary[,4])),
                                 rep("SPD5",times=length(spd5.psummary[,4])))
)

ggplot(bp.crps.data, aes(x=Model, y=mCRPS, fill=Model)) + 
  geom_boxplot(alpha=0.8) +
  theme(legend.position="none") +
  scale_fill_brewer(palette="BuPu")



bp.logs.data<-data.frame(mLogS=c(ind.psummary[,5],bvm.psummary[,5],lmc.psummary[,5],spd2.psummary[,5],spd4.psummary[,5],spd5.psummary[,5]),
                         Model=c(rep("IND",times=length(ind.psummary[,5])),
                                 rep("BVM",times=length(bvm.psummary[,5])),
                                 rep("LMC",times=length(lmc.psummary[,5])),
                                 rep("SPD2",times=length(spd2.psummary[,5])),
                                 rep("SPD4",times=length(spd4.psummary[,5])),
                                 rep("SPD5",times=length(spd5.psummary[,5])))
)

ggplot(bp.logs.data, aes(x=Model, y=mLogS, fill=Model)) + 
  geom_boxplot(alpha=0.8) +
  theme(legend.position="none") +
  scale_fill_brewer(palette="BuPu")



############ Plotting coherence functions #########


######## First we compute pointwise standard deviation for all the estimated coherence functions ######

bvm.sd<-apply(bvm.coh,2, "sd")
lmc.sd<-apply(lmc.cohv,2, "sd")
spd2.sd<-apply(spd2.coh,2, "sd")
spd4.sd<-apply(spd4.coh,2, "sd")
spd5.sd<-apply(spd5.coh,2, "sd")


######## Pointwise mean ########

bvm.coh.mean<-colMeans(bvm.coh)
lmc.coh.mean<-colMeans(lmc.cohv)
spd2.coh.mean<-colMeans(spd2.coh)
spd4.coh.mean<-colMeans(spd4.coh)
spd5.coh.mean<-colMeans(spd5.coh)




##### Bivariate Matern ######


l.limit.bvm<-bvm.coh.mean-1.96*bvm.sd
u.limit.bvm<-bvm.coh.mean+1.96*bvm.sd


pl.bvm.data<-data.frame(u=u,coh=bvm.coh.mean,
                        c.min=l.limit.bvm,
                        c.max=u.limit.bvm,Model="Mean(estimated coherence functions)")

ggplot(data=pl.bvm.data, aes(x=u, y=coh, colour=Model)) + geom_line(lwd=1.0)+
  geom_ribbon(aes(ymin=c.min, ymax=c.max), linetype=2, alpha=0.2)+ylim(-1, 1)+xlim(0,freq.max)+
  labs(x= expression(omega),    # works fine
       y= expression(gamma(omega)),title=expression(paste("Full bivariate Matérn")))+theme(legend.title=element_blank(),legend.position = c(0.30, 0.85), legend.background = element_rect(color = "black", fill = "grey90", size = 0.1, linetype = "solid"))


###### LMC ######


l.limit.lmc<-lmc.coh.mean-1.96*lmc.sd
u.limit.lmc<-lmc.coh.mean+1.96*lmc.sd


pl.lmc.data<-data.frame(u=u,coh=lmc.coh.mean,
                        c.min=l.limit.lmc,
                        c.max=u.limit.lmc,Model="Mean(estimated coherence functions)")

ggplot(data=pl.lmc.data, aes(x=u, y=coh, colour=Model)) + geom_line(lwd=1.0)+
  geom_ribbon(aes(ymin=c.min, ymax=c.max), linetype=2, alpha=0.2)+ylim(-1, 1)+xlim(0,freq.max)+
  labs(x= expression(omega),    # works fine
       y= expression(gamma(omega)),title=expression(paste("LMC")))+theme(legend.title=element_blank(),legend.position = c(0.30, 0.85), legend.background = element_rect(color = "black", fill = "grey90", size = 0.1, linetype = "solid"))



##### SPD2 #####

l.limit.spd2<-spd2.coh.mean-1.96*spd2.sd
u.limit.spd2<-spd2.coh.mean+1.96*spd2.sd




pl.spd2.data<-data.frame(u=u,coh=spd2.coh.mean,
                         c.min=l.limit.spd2,
                         c.max=u.limit.spd2,Model="Mean(estimated coherence functions)")


ggplot(data=pl.spd2.data, aes(x=u, y=coh, colour=Model)) + geom_line(lwd=1.0)+
  geom_ribbon(aes(ymin=c.min, ymax=c.max), linetype=2, alpha=0.2)+ylim(-1, 1)+
  labs(x= expression(omega),    # works fine
       y= expression(gamma(omega)),title=expression(paste("Semiparametric (",Delta,"=2)")))+theme(legend.title=element_blank(),legend.position = c(0.30, 0.85), legend.background = element_rect(color = "black", fill = "grey90", size = 0.1, linetype = "solid"))




######## SPD4 ##########

l.limit.spd4<-spd4.coh.mean-1.96*spd4.sd
u.limit.spd4<-spd4.coh.mean+1.96*spd4.sd




pl.spd4.data<-data.frame(u=u,coh=spd4.coh.mean,
                         c.min=l.limit.spd4,
                         c.max=u.limit.spd4,Model="Mean(estimated coherence functions)")


ggplot(data=pl.spd4.data, aes(x=u, y=coh, colour=Model)) + geom_line(lwd=1.0)+
  geom_ribbon(aes(ymin=c.min, ymax=c.max), linetype=2, alpha=0.2)+ylim(-1, 1)+
  labs(x= expression(omega),    # works fine
       y= expression(gamma(omega)),title=expression(paste("Semiparametric (",Delta,"=4)")))+theme(legend.title=element_blank(),legend.position = c(0.30, 0.85), legend.background = element_rect(color = "black", fill = "grey90", size = 0.1, linetype = "solid"))



######## SPD5 ##########

l.limit.spd5<-spd5.coh.mean-1.96*spd5.sd
u.limit.spd5<-spd5.coh.mean+1.96*spd5.sd




pl.spd5.data<-data.frame(u=u,coh=spd5.coh.mean,
                         c.min=l.limit.spd5,
                         c.max=u.limit.spd5,Model="Mean(estimated coherence functions)")


ggplot(data=pl.spd5.data, aes(x=u, y=coh, colour=Model)) + geom_line(lwd=1.0)+
  geom_ribbon(aes(ymin=c.min, ymax=c.max), linetype=2, alpha=0.2)+ylim(-1, 1)+
  labs(x= expression(omega),    # works fine
       y= expression(gamma(omega)),title=expression(paste("Semiparametric (",Delta,"=5)")))+theme(legend.title=element_blank(),legend.position = c(0.30, 0.85), legend.background = element_rect(color = "black", fill = "grey90", size = 0.1, linetype = "solid"))


######## IND ##########

l.limit.ind<-rep(0,times=500)-1.96*rep(0,times=500)
u.limit.ind<-rep(0,times=500)+1.96*rep(0,times=500)




pl.ind.data<-data.frame(u=u,coh=rep(0,times=500),
                        c.min=l.limit.ind,
                        c.max=u.limit.ind,Model="Mean(estimated coherence functions)")


ggplot(data=pl.ind.data, aes(x=u, y=coh, colour=Model)) + geom_line(lwd=1.0)+
  geom_ribbon(aes(ymin=c.min, ymax=c.max), linetype=2, alpha=0.2)+ylim(-1, 1)+
  labs(x= expression(omega),    # works fine
       y= expression(gamma(omega)),title=expression(paste("Independent Matérn")))+theme(legend.title=element_blank(),legend.position = c(0.30, 0.85), legend.background = element_rect(color = "black", fill = "grey90", size = 0.1, linetype = "solid"))









save.image("Fulldata_analysiswnc.RData")


