# This script will make use of OCCCI ocean colour data
# that can be downloaded with the script 'Download_ocean_colour_data.R'
# We will use Principal Component Analysis (PCA), an exploratory multivariate statistical method,
# to explore variance structure in the data, and identify that 
# the presence of calcifying phytoplankton might bias our capacity to infer
# common place derived-variables from ocean colour, such as chlorophyll-a concentration. 
# Thereafter we will use threshold based on PC1, in order to create a mask to identify pixels
# dominated by calcifying phytoplankton, so that we can identify when derived variables such as 
# chlorophyll-a concentration might be influenced by the presence of calcifiers. 

setwd('D:/Documents/Github/OCCCI/')
# This is the directory I have stored the ocean colour data in
# You may have selected a different directory

latitudes<-rev(seq(70,85,length.out=360))
longitudes<-seq(0,60,length.out=1441)
# Define a system of latitudes and longtiudes

library('reshape2') # If you have not installed this package run
# install.packages('reshape2')

year<-'2016'
month<-'07'
days<-'06'
# This is the time of interest

wavebands<-c('412','443','490','510','560','665')
# These are the spectral channels (nm)

# The following loop will open up the data and arrange it into a single matrix
for(i in 1:length(wavebands)){
	data<-read.csv(paste('Rrs',wavebands[i],'_',year,'_',month,'_',days,'.csv',sep=''), skip = 12, header = F)
	data_store<-as.matrix(data[2:361,2:1442])
	data_store[which(data_store>4e36)]<-NA # Exclude invalid entries
	rownames(data_store)<-latitudes
	colnames(data_store)<-longitudes
	if(i==1){
		complete_data<-melt(data_store)
	} else {
		complete_data<-cbind(complete_data,melt(data_store)[,3])
	}
	rm(data) ; rm(data_store)
}

colnames(complete_data)[3:8]<-wavebands
colnames(complete_data)[1:2]<-c('Lat','Lon')

complete_data<-complete_data[complete.cases(complete_data),]
# Make this sparse matrix into a dense matrix.

PCA<-prcomp(complete_data[,-c(1,2)])
# Perform Principal Compoenent Analysis, to explore variance structure. 
explained_variance<-summary(PCA)$sdev^2/sum(summary(PCA)$sdev^2)
# Compute the variance explained by each Principal Component. 

n<-1
k<-2

plot(PCA$x[,1],PCA$x[,2],xlab=paste('PC',n,'_',round(explained_variance[n]*100),'%',sep=''),
ylab=paste('PC',k,'_',round(explained_variance[k]*100),'%',sep=''),
main='Variance structure of ocean colour data')

# It is evident, from this projection, that almost all of the variance in the ocean colour data
# is summarised along the first Principal Component. An 'arch effect' is also evident, because
# linear trends in PC2 are evident, which are sub-orthogonal to PC1. 

# Let us now produce a Pseudo-True colour image of the data, 
# and compare the appearance with the distribution of the data along PC1. 

# let 560 be the green, let 665 be red, let 443 be blue

red<-complete_data[complete.cases(complete_data),8]
green<-complete_data[complete.cases(complete_data),7]
blue<-complete_data[complete.cases(complete_data),4]

colour_balance<-cbind(red/0.008,green/0.008,blue/0.01)

colour_balance[which(colour_balance[,1]>1),1]<-1
colour_balance[which(colour_balance[,2]>1),2]<-1
colour_balance[which(colour_balance[,3]>1),3]<-1

colour<-rgb(colour_balance[,1],colour_balance[,2],colour_balance[,3])

map_data<-as.data.frame(cbind(complete_data$Lat,complete_data$Lon,PCA$x[,1]))
colnames(map_data)[3]<-'PC1_score'

library(ggplot2) # If you have not installed this package run
# install.packages('ggplot2')

True_colour<-
ggplot(data=map_data[complete.cases(map_data),], aes(x=V2,y=V1))+
labs(y='Latitude',x='Longitude')+
geom_raster(aes(x=V2,y=V1),fill=colour)
# This is a pseudo true colour image. We can see that, within the blue sea, 
# there are regions of green phytoplankton blooms, but also brighter, whiter blooms.
# These bright blooms are probably caused by calcifying phytoplankton, releasing chalky 
# elements of their exoskeletons into the water. 

PC1_val<-
ggplot(data=map_data[complete.cases(map_data),], aes(x=V2,y=V1))+
geom_raster(aes(x=V2,y=V1,fill=PC1_score))+
labs(y='Latitude',x='Longitude')+
scale_fill_gradientn(colours=rev(rainbow(10))[4:10],na.value='white')
# We can easily see that PC1, the main component of variance in the dataset,
# distinguishes the regions with intense blooms of calcifying phytoplankton. 


# We will now use a standard off-the-shelf algorithm to infer 
# chlorophyll-a concentration (a proxy of the biomass of photosynthetic phytoplankton in the 
# surface ocean) from the ocean-colour data
# We will use the OC_6_MERIS algorithm of O'Reilly et al., 2019
# O'Reilly, John E, and P Jeremy Werdell. 2019. Chlorophyll algorithms for ocean color
# sensors-OC4, OC5 & OC6." Remote sensing of environment 229:32-47.
# The resultant concentration is reported in mg per cubic metre

A<-c(0.95087,-3.05489,2.18141,-1.11783,0.15132)


	y<-complete_data[complete.cases(complete_data),-c(1,2,7,8)]
	Max_band<-apply(y,1,function(y) which(y==max(y))[1] )
	i<-cbind(1:nrow(y),Max_band)
	Numerator<-y[i]
	x<-complete_data[complete.cases(complete_data),-c(1,2)]
	Denominator<-apply(x,1,function(x) mean(x[5],x[6]))
	MBR<-Numerator/Denominator
	R<-log10(MBR)
	log_chla<- A[1] + (A[2]*R)+ (A[3]*(R^2))+ (A[4]*(R^3))+
 	(A[5]*(R^4))
	chla<-10^log_chla
	chla[which(chla>20)]<-20 # CAP THE MAXIMUM at 20 mg m-3; values higher than this likely spurious !!!

chla_data<-as.data.frame(cbind(complete_data$Lat,complete_data$Lon,chla))

chla_plot<-
ggplot(data=chla_data, aes(x=V2,y=V1))+
geom_raster(aes(x=V2,y=V1,fill=chla))+
labs(fill=expression(paste('[',italic('chl-a'),']',,' (mg m'^-3,')')),y='Latitude',x='Longitude')+
scale_fill_gradientn(colours=rev(rainbow(10))[4:10],na.value='white')

# Comparison with the pseudo true-colour image shows that the greenest regions
# have the highest chlorophyll-a concentration, but that some regions with calcifying phytoplankton
# are also inferred to have moderate concentrations of chlorophyll-a. 
# This is undesirable, because the ocean colour in these waters is more strongly influenced by 
# chalky exoskeletons than intracellular pigments such as chlorophyll-a. 
# We might hence wish to identify and mask these pixels in any biological analyses. 
# our Principal Component Analysis has shown us that PC1 is highly indicative of coccolithophores, 
# which makes a lot of sense, as they are one of the largest sources of variance in ocean-colour in the image. 

# We can reason that we are therefore not interested in pixels where inferred chlorophyll-a concentration 
# is tightly coupled to variation in PC1. 
# If we plot inferred chlorophyll-a concentration against PC1, we can see that this is true for
# pixels with a PC1 score above 0

plot(PCA$x[,1],chla,xlab=paste('PC',n,'_',round(explained_variance[n]*100),'%',sep=''),
ylab='')
title(ylab=expression(paste('[',italic('chl-a'),']',,' (mg m'^-3,')')),line=1.9, cex.lab=1,)
abline(v=0,col='red',lty=2,lwd='3')

Calcifiers<-rep(0,dim(map_data)[1])
Calcifiers[which(PCA$x[,1]>0)]<-1

Calc_data<-cbind(map_data,Calcifiers)

chla_plot_calc<-
ggplot(data=chla_data, aes(x=V2,y=V1))+
geom_raster(aes(x=V2,y=V1,fill=chla))+
scale_fill_gradientn(colours=rev(rainbow(10))[4:10],na.value='white')+
labs(fill=expression(paste('[',italic('chl-a'),']',,' (mg m'^-3,')')),y='Latitude',x='Longitude')+
geom_contour(data=Calc_data,aes(x=V2,y=V1,z=Calcifiers),col='black',breaks=1,size=1)

# The regions which we should exclude from our analysis, if we wish to infer chlorophyll-a concentration
# from ocean colour, are enclosed in a black contour. 

library(ggpubr) # If you have not installed this package run
# install.packages('ggpubr')

dev.new(width=12,height=6,unit='cm')
ggarrange(True_colour, PC1_val, widths=c(0.8,1))

dev.new(width=12,height=6,unit='cm')
ggarrange(True_colour, chla_plot_calc,common.legend=T)


