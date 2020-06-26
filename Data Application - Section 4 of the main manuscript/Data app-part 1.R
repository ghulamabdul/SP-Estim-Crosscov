####################################################
############ Loading required libraries ############
####################################################
library(sp)
library(maps)
library(maptools)
library(geosphere)
library(fields)
library(MASS)
library(scoringRules)
library(doParallel)
library(rgdal)
#########################################################
################ Data Creation portion ##################
#########################################################
#### set the working directory to be the one where rawdata files are stored ###
setwd("/Users/qadirga/Documents/Project 2/Biometrics submission/Submission Documents/Rcodes/Data Analysis (Section 4 of the Manuscript)/Raw Data files")
#------ PM25 Data ------#

time_pm25_day_read <- system.time(pm25_day_2013 <- read.table("2013_pm25_daily_average.txt", header = TRUE, sep = ","))
time_pm25_day_read

str(pm25_day_2013); head(pm25_day_2013)
# tmp <- pm25_day_2006[pm25_day_2006$FIPS==1001020100,]

pm25_day_2013$Month <- as.numeric(substr(as.character(pm25_day_2013$Date),6,7))
unique(pm25_day_2013$Month)
time_pm25_mon <- system.time(pm25_mon_2013 <- aggregate(pm25_daily_average ~ Month + FIPS + Longitude + Latitude, data=pm25_day_2013, mean))
time_pm25_mon

str(pm25_mon_2013); head(pm25_mon_2013)
names(pm25_mon_2013)[names(pm25_mon_2013)=="pm25_daily_average"] <- "pm25_monthly_mean"
str(pm25_mon_2013); head(pm25_mon_2013) # 867396 obs

pm25_201301 <- subset(pm25_mon_2013, Month==1)
write.csv(pm25_201301, file = "my_pm25_201301.csv", row.names = FALSE)
pm25_201301<-read.csv("my_pm25_201301.csv")
par(mfrow=c(1,2))


#_______ Data over whole USA _______#
quilt.plot(pm25_201301$Longitude,pm25_201301$Latitude,pm25_201301$pm25_monthly_mean,nx=200,ny=200)


#-------------------------------------  WS DATA--------------------------------------------------------------#

# read raw data ####
grib_GDAL_201301 <- readGDAL("narrmon-a_221_20130101_0000_000.grb")


#==============================================================================================####
# transform LCC coordinate to Long/Lat ####
# str(grib_GDAL_200606)
grib_GDAL_201301@proj4string

# str(grib_GDAL_200606@data)
dim(grib_GDAL_200606@data)

# longlat_200606 <- spTransform(SpatialPoints(coordinates(grib_GDAL_200606), proj4string=grib_GDAL_200606@proj4string), CRS("+proj=longlat +datum=WGS84"))
longlat <- spTransform(SpatialPoints(coordinates(grib_GDAL_201301), proj4string=grib_GDAL_201301@proj4string), CRS("+proj=longlat +datum=WGS84"))
dim(coordinates(grib_GDAL_201301))
dim(coordinates(longlat))

#==============================================================================================####
narr_201301 <- data.frame(Month=201301, Longitude=coordinates(longlat)[,1],	Latitude=coordinates(longlat)[,2],
                          WS=sqrt(grib_GDAL_201301@data[,277]^2 + grib_GDAL_201301@data[,323]^2))

write.csv(narr_201301, file = "my_narr_201301.csv", row.names = FALSE)


##### ______ Merging data files ________#######

pm25_201301 <- read.csv("my_pm25_201301.csv")
narr_201301 <- read.csv("my_narr_201301.csv")

#==============================================================================================####
long_pm25<- unique(pm25_201301$Longitude)
lat_pm25 <- unique(pm25_201301$Latitude)

data_201301 <- subset(narr_201301, Longitude >= min(long_pm25) & Longitude <= max(long_pm25) & Latitude >= min(lat_pm25) & Latitude <= max(lat_pm25))

#==============================================================================================####

LL_pm25 <- cbind(pm25_201301$Longitude, pm25_201301$Latitude)
LL_narr <- cbind(data_201301$Longitude, data_201301$Latitude)
dim(LL_pm25); dim(LL_narr)
quilt.plot(x=pm25_201301$Longitude,y=pm25_201301$Latitude,z=pm25_201301$pm25_monthly_mean)
quilt.plot(x=data_201301$Longitude,y=data_201301$Latitude,z=data_201301$WS)

system.time(DMatrix <- distm(LL_pm25,LL_narr)); dim(DMatrix)    # 1398.215 sec
system.time(tmp <- apply(DMatrix, 1, min))                      # 20.734 sec
system.time(min_ind <- which(DMatrix==tmp, arr.ind=T))          # 3.275 sec
system.time(id_data <- min_ind[order(min_ind[,1]),2])           # 0.002 sec

summary(min_ind[,1]); summary(min_ind[,2])
length(unique(min_ind[,1])); length(unique(min_ind[,2]))

#==============================================================================================####
pm25_201301$id_data <- id_data
data_pm25_201301 <- aggregate(pm25_monthly_mean ~ id_data, data=pm25_201301, mean)

#==============================================================================================####
n_narr <- dim(LL_narr)[1]

data_201301$id_data <- 1:n_narr

#==============================================================================================####
total_201301 <- merge(data_201301, data_pm25_201301, by="id_data"); str(total_201301)
write.csv(total_201301[,-1], file = "my_data_pm25_201301.csv", row.names = FALSE)
total_201301<-read.csv("my_data_pm25_201301.csv")
#####
par(mfrow=c(1,2))
quilt.plot(total_201301$Longitude,total_201301$Latitude,total_201301$WS,main="Wind speed")
quilt.plot(total_201301$Longitude,total_201301$Latitude,total_201301$pm25_monthly_mean,main="PM2.5")

##### -------- Assigning Climatic region --------########
latlong2state <- function(pointsDF) {
  # Prepare SpatialPolygons object with one SpatialPolygon
  # per state (plus DC, minus HI & AK)
  states <- map('state', fill=TRUE, col="transparent", plot=FALSE)
  IDs <- sapply(strsplit(states$names, ":"), function(x) x[1])
  states_sp <- map2SpatialPolygons(states, IDs=IDs,
                                   proj4string=CRS("+proj=longlat +datum=WGS84"))
  
  # Convert pointsDF to a SpatialPoints object 
  pointsSP <- SpatialPoints(pointsDF, 
                            proj4string=CRS("+proj=longlat +datum=WGS84"))
  
  # Use 'over' to get _indices_ of the Polygons object containing each point 
  indices <- over(pointsSP, states_sp)
  
  # Return the state names of the Polygons object containing each point
  stateNames <- sapply(states_sp@polygons, function(x) x@ID)
  stateNames[indices]
}


# Assign States
coords <- as.data.frame(cbind(total_201301$Longitude,total_201301$Latitude))
State <- latlong2state(coords)



#==============================================================================================####
unique(as.factor(State))
summary(as.factor(State))

id_NA <- (1:length(State))[is.na(State)]
id_ST <- (1:length(State))[!is.na(State)]
coords_NA <- coords[id_NA,]
coords_ST <- coords[id_ST,]
dim(coords_NA); dim(coords_ST)

DMatrix <- distm(coords_NA,coords_ST); dim(DMatrix)
tmp <- apply(DMatrix, 1, min)
min_ind <- which(DMatrix==tmp, arr.ind=T)
NA2ST <- id_ST[min_ind[order(min_ind[,1]),2]]

State[id_NA] <- State[NA2ST]
unique(as.factor(State))
summary(as.factor(State))
#==============================================================================================####
# Assing Climatic Regions (CR)
CR_NW <- c("washington", "oregon", "idaho")
CR_W <- c("california", "nevada")
CR_SW <- c("utah", "colorado", "arizona", "new mexico")
CR_WNC <- c("montana", "wyoming", "north dakota", "south dakota", "nebraska")
CR_ENC <- c("minnesota", "iowa", "wisconsin", "michigan")
CR_S <- c("kansas", "oklahoma", "texas", "arkansas", "louisiana", "mississippi")
CR_C <- c("illinois", "indiana", "ohio", "missouri", "kentucky", "west virginia", "tennessee")
CR_NE <- c("maine", "new hampshire", "vermont", "new york", "massachusetts", "connecticut", "rhode island",
           "pennsylvania", "new jersey", "delaware", "maryland")
CR_SE <- c("south carolina", "georgia", "alabama", "florida", "north carolina", "virginia")



CR <- rep(NA,length(State))
CR[State %in% CR_NW] <- "NW"
CR[State %in% CR_W]  <- "W"
CR[State %in% CR_SW] <- "SW"
CR[State %in% CR_WNC]<- "WNC"
CR[State %in% CR_ENC]<- "ENC"
CR[State %in% CR_S]  <- "S"
CR[State %in% CR_C]  <- "C"
CR[State %in% CR_NE] <- "NE"
CR[State %in% CR_SE] <- "SE"



summary(as.factor(CR))
#==============================================================================================####

total_201301$State<-State; total_201301$CR<-CR

coords_NW <- coords[CR == "NW",];  text_NW <- (apply(coords_NW,2,min) + apply(coords_NW,2,max))/2
coords_W  <- coords[CR == "W",];   text_W  <- (apply(coords_W,2,min) + apply(coords_W,2,max))/2
coords_SW <- coords[CR == "SW",];  text_SW <- (apply(coords_SW,2,min) + apply(coords_SW,2,max))/2
coords_WNC<- coords[CR == "WNC",]; text_WNC<- (apply(coords_WNC,2,min) + apply(coords_WNC,2,max))/2
coords_ENC<- coords[CR == "ENC",]; text_ENC<- (apply(coords_ENC,2,min) + apply(coords_ENC,2,max))/2
coords_S  <- coords[CR == "S",];   text_S  <- (apply(coords_S,2,min) + apply(coords_S,2,max))/2
coords_C  <- coords[CR == "C",];   text_C  <- (apply(coords_C,2,min) + apply(coords_C,2,max))/2
coords_NE <- coords[CR == "NE",];  text_NE <- (apply(coords_NE,2,min) + apply(coords_NE,2,max))/2
coords_SE <- coords[CR == "SE",];  text_SE <- (apply(coords_SE,2,min) + apply(coords_SE,2,max))/2

# text_NW <- apply(coords_NW,2,median)
text_W  <- apply(coords_W,2,median)
text_SW <- apply(coords_SW,2,median)
text_WNC<- (apply(coords_WNC,2,median) + text_WNC)/2
text_ENC<- apply(coords_ENC,2,median)
text_S  <- apply(coords_S,2,median)
text_C  <- apply(coords_C,2,median)
text_NE <- apply(coords_NE,2,median)
text_SE <- apply(coords_SE,2,median)


par(mfrow=c(1,1))
map('state'); 
points(total_201301$Longitude, total_201301$Latitude,  col = "red", cex = .1)
map('state', region = CR_NW, lwd=3, interior = FALSE, add = TRUE); text(text_NW[1], text_NW[2], "NW", cex=1.7)
map('state', region = CR_W, lwd=3, interior = FALSE, add = TRUE);  text(text_W[1]+1.5, text_W[2], "W", cex=1.7)
map('state', region = CR_SW, lwd=3, interior = FALSE, add = TRUE); text(text_SW[1], text_SW[2], "SW", cex=1.7)
map('state', region = CR_WNC, lwd=3, interior = FALSE, add = TRUE);text(text_WNC[1], text_WNC[2], "WNC", cex=1.7)
map('state', region = CR_ENC, lwd=3, interior = FALSE, add = TRUE);text(text_ENC[1], text_ENC[2], "ENC", cex=1.7)
map('state', region = CR_S, lwd=3, interior = FALSE, add = TRUE);  text(text_S[1], text_S[2], "S", cex=1.7)
map('state', region = CR_C, lwd=3, interior = FALSE, add = TRUE);  text(text_C[1], text_C[2], "C", cex=1.7)
map('state', region = CR_NE, lwd=3, interior = FALSE, add = TRUE); text(text_NE[1], text_NE[2], "NE", cex=1.7)
map('state', region = CR_SE, lwd=3, interior = FALSE, add = TRUE); text(text_SE[1]-2, text_SE[2], "SE", cex=1.7)
title(main="Climatic Regions", cex.main=2)

c.r.int<-CR_WNC ###Area of interest
a.init<-"WNC"
##########################################################
### Now plotting data of WNC region only ##################
##########################################################
WNC_201301_total<-total_201301[total_201301$CR==a.init,]
par(mfrow=c(1,1))
map('state')
map('state', region = c.r.int, lwd=3, interior = FALSE, add = TRUE); 
points(WNC_201301_total$Longitude,WNC_201301_total$Latitude,col=color.scale(WNC_201301_total$pm25_monthly_mean,tim.colors(),zlim=c(min(WNC_201301_total$pm25_monthly_mean)-0.01,max(WNC_201301_total$pm25_monthly_mean)+0.01)),pch=19,cex=0.5)
image.plot(legend.only = T,horizontal = T,zlim=c(min(WNC_201301_total$pm25_monthly_mean)-0.01,max(WNC_201301_total$pm25_monthly_mean)+0.01))
title("PM25 Data")


map('state')
map('state', region = c.r.int, lwd=3, interior = FALSE, add = TRUE); 
points(WNC_201301_total$Longitude,WNC_201301_total$Latitude,col=color.scale(WNC_201301_total$WS,tim.colors(),zlim=c(min(WNC_201301_total$WS)-0.01,max(WNC_201301_total$WS)+0.01)),pch=19,cex=0.5)
image.plot(legend.only = T,horizontal = T,zlim=c(min(WNC_201301_total$WS)-0.01,max(WNC_201301_total$WS)+0.01))
title("Wind speed")

rm.list<-ls()[(ls()!="WNC_201301_total")&(ls()!="c.r.int")&(ls()!="a.init")]
rm(list=rm.list) ####### removing everything except the data file #######
#------ WNC_201301_total contains the ungridded data of wind speed (unstandardized) and PM2.5 data standardized in the WNC region of the united states ---------------#

######## Data plot of log (PM2.5) and Windspeed #######

#hist(log(WNC_201301_total$pm25_monthly_mean),breaks = 30)

par(mar=c(6.55,4,0.2,0.2))
plot(x=WNC_201301_total$Longitude,y=WNC_201301_total$Latitude,xlab="Longitude",ylab="Latitude",pch=20,col="purple")
points(WNC_201301_total$Longitude,WNC_201301_total$Latitude,col=color.scale(log(WNC_201301_total$pm25_monthly_mean),tim.colors()),pch=20,cex=1)
image.plot(legend.only = T,horizontal = T,zlim=c(min(log(WNC_201301_total$pm25_monthly_mean)),max(log(WNC_201301_total$pm25_monthly_mean))),legend.width = 0.6,legend.mar = 1.9)
map('state', region = c.r.int,add = T)

par(mar=c(6.55,4,0.2,0.2))
plot(x=WNC_201301_total$Longitude,y=WNC_201301_total$Latitude,xlab="Longitude",ylab="Latitude",pch=20,col="purple")
points(WNC_201301_total$Longitude,WNC_201301_total$Latitude,col=color.scale((WNC_201301_total$WS),tim.colors()),pch=20,cex=1)
image.plot(legend.only = T,horizontal = T,zlim=c(min((WNC_201301_total$WS)),max((WNC_201301_total$WS))),legend.width = 0.6,legend.mar = 1.9)
map('state', region = c.r.int,add = T)


##############################################################################################################
########### Now we do standardization and transformation on the ungridded data ################################
##############################################################################################################
# We take log transformation of PM 2.5 data to make it nearly Gaussian
# We then subtract the sample mean and standard deviation to standardize the data

un.grd.total<-data.frame(lon=WNC_201301_total$Longitude,lat=WNC_201301_total$Latitude, PM2_5=WNC_201301_total$pm25_monthly_mean,WS=WNC_201301_total$WS)
hist(un.grd.total$PM2_5)
un.grd.total$PM2_5<-log(un.grd.total$PM2_5)
un.grd.total$PM2_5<-(un.grd.total$PM2_5-mean(un.grd.total$PM2_5))/sd(un.grd.total$PM2_5)
un.grd.total$WS<-(un.grd.total$WS-mean(un.grd.total$WS))/sd(un.grd.total$WS)
hist(un.grd.total$PM2_5) #### PM2_5 is actually log transformed PM2_5 ######
hist(un.grd.total$WS)

#################################################
###### Some required functions ##################
#################################################

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



##### first we divide data in training and test set #####
rm(list=ls()[(ls()!="un.grd.total")&(ls()!="biwm_coh")&(ls()!="Bspline")&(ls()!="f.matern")&(ls()!="my.matern")&(ls()!="c.r.int")&(ls()!="a.init")
             ])
# We removed all the variables from the environment 
# Now un.grd.total is the data file with which we will work



map('state', region = c.r.int)

quilt.plot(x=un.grd.total$lon,y=un.grd.total$lat,z=un.grd.total$PM2_5,add = T,pch=19)


map('state', region = c.r.int)

quilt.plot(x=un.grd.total$lon,y=un.grd.total$lat,z=un.grd.total$WS, add = T,pch=19)


######## Now we generate set of random indices of 20 percent test locations for 100 runs #########
#### Random splitting index ###########

n.test.locs<-round(0.2*length(un.grd.total$PM2_5),0)





rand.index<-matrix(NA,nrow=n.test.locs,ncol=100)

for(i in 1:100)
{
  set.seed(i)
  rand.index[,i]<-sample(size = n.test.locs,1:length(un.grd.total$PM2_5))
  #### ith column is the test index for the ith run of analysis
}


## set the woring directory to be the place where you want to save the final dataset to be worked with.
#setwd("/Users/qadirga/Documents/Project 2/Biometrics submission/Revision work/Revised Data application/Different spatial regions/WNC_Final version/WNC")
save.image("Final dataset to be worked with.RData")
