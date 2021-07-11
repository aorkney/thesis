# This is an example script to illustration how two block PLS is performed
# and its utilities.

# We are going to use the dataset 'SST_and_storms.csv'
# This is a dataset I produced by downloading sea surface temperature data
# from the ARCTIC_REANALYSIS_PHYS_002_003 dataset, 
# and wind-speed from the WIND_GLO_WIND_L3_REP_OBSERVATIONS_012_005 dataset.
# Both of these datasets are available on https://marine.copernicus.eu/

# I averaged data for September, over a series of different years, 
# within 1 degree latitude-longitude boxes, so that the dataset is 
# small and manageable, and so that researchers interested in the original
# data will be motivated to visit https://marine.copernicus.eu/
 
# The csv file contains columns for latitude, longitude, year, 
# sea surface temperature, and storminess. 
# storminess is computed as the average number of days,
# within a grid cell, that wind speeds exceed 10 metres per second. 

setwd('D:/Documents/Github/OCCCI')
# You may store your dataset in a different directory; adjust accordingly

dataset<-read.csv('SST_and_storms.csv')
dataset<-dataset[,-1]
dataset<-dataset[complete.cases(dataset),]
# The dataset has been loaded

library(geomorph)
# 'geomorph' is a dataset intended for evolutionary biologists
# studying variation in shape, but it has a useful function
# for implementing 2BPLS analyses
# If you do not have this package then run 'install.packages('geomorph')'

# Let's explore the covariance structure between 
# geoggraphy (latitude and longitude) 
# and sea surface temperature and storminess

# In the Barents Sea, we expect a south-west to north-east
# decline in sea surface temperature. 
# We do not necessariyl expect stormy days to have a strong relationship
# To geography in a given year. Perhaps there will be more in the warmest waters,
# closest to the Atlantic inflow region, where storms might enter the Barents Sea.

library(ggplot2) # If you have not installed this package run
# install.packages('ggplot2')

sst<-
ggplot(data=dataset[which(dataset$years=='2002'),], aes(x=lon,y=lat))+
geom_raster(aes(x=lon,y=lat,fill=sst))+
labs(fill=expression(paste('SST',,' ('^o,'C)')),y='Latitude',x='Longitude')+
scale_fill_gradientn(colours=rev(rainbow(10))[4:10],na.value='white')
# This is a plot of what SST conditions looked like in September 2002
# A clear geographic gradient along the SW-NE axis is evident. 

storms<-
ggplot(data=dataset[which(dataset$years=='2002'),], aes(x=lon,y=lat))+
geom_raster(aes(x=lon,y=lat,fill=storminess))+
labs(fill='stormy days',y='Latitude',x='Longitude')+
scale_fill_gradientn(colours=rev(rainbow(10))[4:10],na.value='white')
# This is a plot of what storminess conditions looked like in September 2002
# (average number of stormy days within each grid cell)

Geography<-dataset[which(dataset$years=='2002'),c(1,2)]
Weather<-dataset[which(dataset$years=='2002'),c(4,5)]
# Let's define two different blocks of data, one represent geography, another representing 
# the conditions (stormy days and sea surface temperature; we'll just call this 'Weather')

Analysis<-two.b.pls(Geography, Weather)
# We have performed a two block partial least squares analysis, exploring the covariance between these two
# input blocks of data.
# Type 'names(Analysis)' to view what sorts of objects we have produced.
# Type 'summary(Analysis)' to see the summary statistics. Note a P-value well below 0.05, 
# which indicates that the leading axis of covariance between these blocks of data is highly significant.
# A moderate r-PLS value and a high Effect size testify to the strength of the association. 

df<-cbind(Analysis$XScores[,1],Analysis$YScores[,1],dataset[which(dataset$years=='2002'),c(4,5)])
colnames(df)<-c('X','Y','storminess','sst')
# 2BPLS works by decomposing the two blocks of data into left and right orthogonal vectors, 
# over which variance in the data is distributed.
# The pairs are sorted in order from of those which explain the most covariance between the blocks
# to those which explain the least. 

PLS1_sst<-
ggplot(data=df,aes(x=X,y=Y))+
geom_point(aes(x=X,y=Y,col=sst))+
labs(colour=expression(paste('SST',,' ('^o,'C)')),y='PLS1 Block 2',x='PLS1 Block 1')+
scale_colour_gradientn(colours=rev(rainbow(10))[4:10],na.value='white')
# We can see that the first pair of vectors explains variance in sea surface temperature.
# If we review the vectors by typing 'Analysis$left.pls.vectors[,1]' we can see that
# latitude and longitude are increasing. 
# and if we type 'Analysis$right.pls.vectors[,1]' we can see that storminess and sea surface temperature are declining
# the large coefficient for temperature, and the high degree of sorting along these axes, suggest that the first
# explaiend mode of covariance (or 'salience') is a SW-NE gradient in temperature. This conforms to expectations. 

PLS1_storms<-
ggplot(data=df,aes(x=X,y=Y))+
geom_point(aes(x=X,y=Y,col=storminess))+
labs(colour='storminess',y='PLS1 Block 2',x='PLS1 Block 1')+
scale_colour_gradientn(colours=rev(rainbow(10))[4:10],na.value='white')
# If we view variation in storminess, we can see it is poorly sorted along this first axis of covariance. 

# Let us now inspect the second axis of covariance, representing patterns orthogonal to the first axis. 
# This axis explains less of the covariance than the first one. 
Explained_variance<- Analysis$svd[[1]]^2 / sum(Analysis$svd[[1]]^2)
# Indeed, we can see that the first axis explains some 99% of the covariance. 

df<-cbind(Analysis$XScores[,2],Analysis$YScores[,2],dataset[which(dataset$years=='2002'),c(4,5)])
colnames(df)<-c('X','Y','storminess','sst')

PLS2_sst<-
ggplot(data=df,aes(x=X,y=Y))+
geom_point(aes(x=X,y=Y,col=sst))+
labs(colour=expression(paste('SST',,' ('^o,'C)')),y='PLS2 Block 2',x='PLS2 Block 1')+
scale_colour_gradientn(colours=rev(rainbow(10))[4:10],na.value='white')
# Over this second axis, temperature is not very well organised. 

PLS2_storms<-
ggplot(data=df,aes(x=X,y=Y))+
geom_point(aes(x=X,y=Y,col=storminess))+
labs(colour='storminess',y='PLS2 Block 2',x='PLS2 Block 1')+
scale_colour_gradientn(colours=rev(rainbow(10))[4:10],na.value='white')
# Storminess, however, is well organised. 
# If we type Analysis$right.pls.vectors[,2] we can see the coefficient for storminess is much larger along this axis.
# and Analysis$left.pls.vectors[,2] shows that higher storms are associated strongly with a more southern position. 

# We have hence discovered the patterns we expected with this exploratory analysis;
# that there is a very strong SW-NE gradient in temperature, 
# that variation in storms is somewhat independent of this in a given year (because of stochastic processes of weather), 
# and that storms are, overall, more likely to be found in the Southern Barents Sea, because this is the region
# where the atmospheric storm track enters. 

# We have hence seen the utility of using 2BPLS for exploring patterns of variance between different datasets, 
# provided they can be represented as continuously-distributed dense matrices of data, and the 
# fundamental relationships that describe them can be approximated as linear. 

library(ggpubr)# If you have not installed this package run
# install.packages('ggpubr')
dev.new(width=12,height=6,unit='cm')
ggarrange(sst,storms)
dev.new(width=12,height=6,unit='cm')
ggarrange(PLS1_sst,PLS1_storms)
dev.new(width=12,height=6,unit='cm')
ggarrange(PLS2_sst,PLS2_storms)




