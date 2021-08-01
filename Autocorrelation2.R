# This script will explore potential relationships between 
# September Sea Surface Temperature, storminess, Atlantic inflow current 
# speed and Biomass in the Barents Sea between 2002 and 2019
# An attempt will be made to determine whether relationships exhibit
# temporal autocorrelation. 

setwd('D:/Documents/Github/OCCCI')
# You may store your dataset in a different directory; adjust accordingly

dataset<-read.csv('SST_and_storms.csv')
dataset<-dataset[,-1]
# The dataset has been loaded


years<-c(2002:2019)
lat_steps<-seq(70,85,1)
lon_steps<-seq(0,60,1)
# House keeping



owa_ts<-c(212958, 264221, 275861, 217604, 245166, 280260, 206171, 205956, 177643, 227852,
 208346, 240441, 227544, 226351, 277610, 255744, 293658, 262918)
# 'Open Water Area', proxied by the number of unique latitude-longitude combinations observed by remote-sensing ocean-colour 
# satellites each September (a measure of open water area)
# Ocean-colour data available at https://www.oceancolour.org/thredds/catalog-cci.html


deg2rad <- function(deg) {(deg * pi) / (180)}
# We need this to compute weights of pixels as function of latitude

sst_ts<-list()
for(j in years){
	sst_ts[[j]]<-weighted.mean(dataset$sst[which(dataset$years==j)],na.rm=T,
	weights= cos(deg2rad(dataset$lat[which(dataset$years==j)])) )
}
sst_ts<-unlist(sst_ts)
# Sea Surface Temperature series, September 2002-2019, in degrees C


wind_ts<-list()
for(j in years){
	wind_ts[[j]]<-weighted.mean(dataset$storminess[which(dataset$years==j)],na.rm=T,
	weights= cos(deg2rad(dataset$lat[which(dataset$years==j)])) )
}
wind_ts<-unlist(wind_ts)
# Frequency of stormy days (winds over 10 metres per second), September 2002-2019


library(nlme) # We need this package to perform generalised least squares fits
# if you do not have it run 'install.packages('nlme')'

linedata<-as.data.frame(cbind(years,(owa_ts-mean(owa_ts))/sd(owa_ts),
(sst_ts-mean(sst_ts))/sd(sst_ts),(wind_ts-mean(wind_ts))/sd(wind_ts)))
colnames(linedata)[2:4]<-c('OWA','SST','Storms')
# Let's organise the data

library(reshape2) # we need this to rearrange dataframes
# if you do not have it run 'install.packages('reshape2')'
long_line_data<-melt(linedata,id='years') # rearrange data

library(ggplot2) # plotting package
# if you do not have this package run 'install.packages('ggplot2')'
ggplot(data=long_line_data)+
geom_line(aes(x=years,y=value,colour=variable),lwd=1)+
geom_line(aes(x=years,y=value,colour=variable,linetype=variable),lwd=2)+
theme(axis.line=element_line(size=3/2),
      axis.text.x=element_text(size=14,colour='black',angle=-90,vjust=0.5,hjust=1),
      axis.text.y=element_text(size=14,colour='black'),
      axis.ticks=element_line(size=2),
      axis.title.x=element_text(size=10,colour='black'),
      axis.title.y=element_blank(),
      legend.position="bottom",
legend.title=element_text(size=10),
legend.text=element_text(size=10,colour=c('black'),face='bold'),
legend.key.width = unit(3, "line"),
legend.key.height = unit(2, "line"),
      panel.background=element_blank(),
      panel.border=element_blank(),
      panel.grid.major=element_blank(),
      panel.grid.minor=element_blank(),
      plot.background=element_blank())+
  scale_x_continuous("years", labels = as.character(years), breaks = years)
# Visualise the trends in the time series
# Perhaps there are relationships between some of these variables
# Sometimes peaks and troughs in one series may appear to lag another
# so perhaps there is also temporal autocorrelation




# Let's perform geographically resolved regressions. 

# This will be a store for the results 
trends<-cbind(expand.grid(lat_steps,lon_steps),NA)
trends<-as.data.frame(trends)
trends<-cbind(trends,rep(NA,dim(trends)[1]))
trends<-cbind(trends,rep(NA,dim(trends)[1]))
colnames(trends)<-c('lat','lon','trends','odds','auto')


for(a in 1:(length(lat_steps)-1)){ # for each latitude step
	for(b in 1:(length(lon_steps)-1)){ # for each longitude step 
		temp<-dataset[which(dataset$lat==lat_steps[a] & dataset$lon==lon_steps[b]),]$sst # compute the sst time series 
		wind<-dataset[which(dataset$lat==lat_steps[a] & dataset$lon==lon_steps[b]),]$storminess # compute the sotrminess time series 
			if( length(which(is.na(temp)==F))>3 & sd(temp,na.rm=T)>0 & length(which(is.na(wind)==F))>3 & sd(wind,na.rm=T)>0){ 
			# check there are sufficient measurements
			chosen_row<-which(trends$lat==lat_steps[a] & trends$lon==lon_steps[b])
			df1<-as.data.frame(cbind(temp,wind,years))
			df1<-df1[complete.cases(df1),]
			mod1<-gls( temp ~ wind,data=df1)
			 # Make a standard gls model 
			tryCatch( mod2<-gls( temp ~ wind  ,correlation=corAR1(form=~years),data= df1) ,error=function(e) print('test') )
			 # Make an autcorrelated model with lag-1 autocorrelation; if there is an error print 'test'
			pval<-1
			if(exists('mod2')){
					pval<-anova(mod1,mod2)$p[2] # Compare standard and autcorrelated models
			} 

			if(pval <0.05){ # If the autocorrelated model is preferred 
					trends[chosen_row,3]<-summary(mod2)$t[2,1]
					trends[chosen_row,4]<-summary(mod2)$t[2,4]
					trends[chosen_row,5]<-2
					rm('mod2')	
				
			} else {
				trends[chosen_row,3]<-summary(mod1)$t[2,1]
				trends[chosen_row,4]<-summary(mod1)$t[2,4]
				trends[chosen_row,5]<-1
				if(exists('mod2')){
					rm(mod2)
				}
				
			}
		
			
			}
	}
	print(a) # update me on progress; how many of the 15 latitude steps are complete?
}

point_shape<-trends[which(trends$odds<=0.05),5] # Are standard or autcorrelated models preferred? 

library(ggplot2)
df<-as.data.frame(trends)
ggplot()+
geom_tile(data=df,aes(x=lon,y=lat,fill=trends))+
scale_fill_gradient2(high='red',low='blue',mid='white',na.value='white')+
geom_point(data=df[which(trends$odds<=0.05),],aes(x=lon,y=lat),shape=point_shape)
# Circles indicate a standard fit is preferred
# Triangles indicate an autocorrelated fir is preferred
# we can see that SST is modified by wind in different ways in different regions of the Barents Sea
# In the Northern Barents Sea Atlantic waters are overlain by a lid of Polar waters, therefore mixing 
# increases SST.
# In the Southern Barents Sea, the warmest waters are at the top of the water column, so mixing 
# causes cooling of SST. 
# The relationships uncovered indicate that models incorporating autocorrelation tend not to be preferred

fit1<-gls(sst_ts ~ wind_ts)
fit2<-gls(sst_ts ~ wind_ts, corr=corAR1(form=~years))
anova(fit1,fit2)
# An autocorrelated fit does not appear to be preferred
# but the p-value is very close to 0.05 
 


# Let's now investigate any relationship between sst and open water area

trends<-cbind(expand.grid(lat_steps,lon_steps),NA)
trends<-as.data.frame(trends)
trends<-cbind(trends,rep(NA,dim(trends)[1]))
trends<-cbind(trends,rep(NA,dim(trends)[1]))
colnames(trends)<-c('lat','lon','trends','odds','auto')

# as before, but with different variables; I will therefore refrain from commenting
for(a in 1:(length(lat_steps)-1)){
	for(b in 1:(length(lon_steps)-1)){
		temp<-dataset[which(dataset$lat==lat_steps[a] & dataset$lon==lon_steps[b]),]$sst	
			if( length(which(is.na(temp)==F))>3 & sd(temp,na.rm=T)>0 & length(which(is.na(owa_ts)==F))>3 & sd(owa_ts,na.rm=T)>0){
			chosen_row<-which(trends$lat==lat_steps[a] & trends$lon==lon_steps[b])
			df1<-as.data.frame(cbind(owa_ts,temp,years))
			df1<-df1[complete.cases(df1),]
			mod1<-gls( temp~owa_ts ,data=df1)
			tryCatch( mod2<-gls( temp~owa_ts  ,correlation=corAR1(form=~years),data= df1) ,error=function(e) print('test') )
			pval<-1
			if(exists('mod2')){
					pval<-anova(mod1,mod2)$p[2]
			} 

			if(pval <0.05){
					trends[chosen_row,3]<-summary(mod2)$t[2,1]
					trends[chosen_row,4]<-summary(mod2)$t[2,4]
					trends[chosen_row,5]<-2
					rm('mod2')	
				
			} else {
				trends[chosen_row,3]<-summary(mod1)$t[2,1]
				trends[chosen_row,4]<-summary(mod1)$t[2,4]
				trends[chosen_row,5]<-1
				if(exists('mod2')){
					rm(mod2)
				}
				
			}
		
			
			}
	}
	print(a)
}

point_shape<-trends[which(trends$odds<=0.05),5]

library(ggplot2)
df<-as.data.frame(trends)
ggplot()+
geom_tile(data=df,aes(x=lon,y=lat,fill=trends))+
scale_fill_gradient2(high='red',low='blue',mid='white',na.value='white')+
geom_point(data=df[which(trends$odds<=0.05),],aes(x=lon,y=lat),shape=point_shape)
# It looks like there is a relationship between September SST and open water area
# It appears a significant number of grid-cells prefer autocorrelated models. 
# SST and OWA in the Southeast Barents Sea in particular, appear to be associated

fit1<-gls(sst_ts ~ owa_ts)
fit2<-gls(sst_ts ~ owa_ts, corr=corAR1(form=~years))
anova(fit1,fit2)
# An autocorrelated model is preferred for the average pattern across the Barents Sea
# this suggests that there is some kind of inter-annual persistence of periodicity 
# that underlies the relationship between open water area and sea surface temperature



# Let us now explore the relationship between wind and open water area
trends<-cbind(expand.grid(lat_steps,lon_steps),NA)
trends<-as.data.frame(trends)
trends<-cbind(trends,rep(NA,dim(trends)[1]))
trends<-cbind(trends,rep(NA,dim(trends)[1]))
colnames(trends)<-c('lat','lon','trends','odds','auto')



for(a in 1:(length(lat_steps)-1)){
	for(b in 1:(length(lon_steps)-1)){
		wind<-dataset[which(dataset$lat==lat_steps[a] & dataset$lon==lon_steps[b]),]$storminess
			if( length(which(is.na(wind)==F))>3 & sd(wind,na.rm=T)>0 & length(which(is.na(owa_ts)==F))>3 & sd(owa_ts,na.rm=T)>0){
			chosen_row<-which(trends$lat==lat_steps[a] & trends$lon==lon_steps[b])
			df1<-as.data.frame(cbind(wind,owa_ts,years))
			df1<-df1[complete.cases(df1),]
			mod1<-gls(  wind ~ owa_ts  ,data=df1)
			tryCatch( mod2<-gls( wind ~ owa_ts  ,correlation=corAR1(form=~years),data= df1) ,error=function(e) print('test') )
			pval<-1
			if(exists('mod2')){
					pval<-anova(mod1,mod2)$p[2]
			} 

			if(pval <0.05){
					trends[chosen_row,3]<-summary(mod2)$t[2,1]
					trends[chosen_row,4]<-summary(mod2)$t[2,4]
					trends[chosen_row,5]<-2
					rm('mod2')	
				
			} else {
				trends[chosen_row,3]<-summary(mod1)$t[2,1]
				trends[chosen_row,4]<-summary(mod1)$t[2,4]
				trends[chosen_row,5]<-1
				if(exists('mod2')){
					rm(mod2)
				}
				
			}
		
			
			}
	}
	print(a)
}

point_shape<-trends[which(trends$odds<=0.05),5]

library(ggplot2)
df<-as.data.frame(trends)
ggplot()+
geom_tile(data=df,aes(x=lon,y=lat,fill=trends))+
scale_fill_gradient2(high='red',low='blue',mid='white',na.value='white')+
geom_point(data=df[which(trends$odds<=0.05),],aes(x=lon,y=lat),shape=point_shape)
# Years with greater open water area tend to be associated with strong winds in the 
# Southwestern Barents Sea. 
# This could be consistent with greater delivery of Atlantic-derived air-masses
# or winds from the Southwest. These could reduce sea-ice extents. 

fit1<-gls(wind_ts ~ owa_ts)
fit2<-gls(wind_ts ~ owa_ts, corr=corAR1(form=~years))
anova(fit1,fit2)
# An autocorrelated fit does not appear to be preferred

owa_tsf<-owa_ts[-18]
wind_tsf<-wind_ts[-1]
yearsf<-years[-1]
# let us explore a lagged relationship of wind on open water area

trends<-cbind(expand.grid(lat_steps,lon_steps),NA)
trends<-as.data.frame(trends)
trends<-cbind(trends,rep(NA,dim(trends)[1]))
trends<-cbind(trends,rep(NA,dim(trends)[1]))
colnames(trends)<-c('lat','lon','trends','odds','auto')



for(a in 1:(length(lat_steps)-1)){
	for(b in 1:(length(lon_steps)-1)){
		wind<-dataset[which(dataset$lat==lat_steps[a] & dataset$lon==lon_steps[b]),]$storminess
		wind<-wind[-1]
			if( length(which(is.na(wind)==F))>3 & sd(wind,na.rm=T)>0 & length(which(is.na(owa_tsf)==F))>3 & sd(owa_tsf,na.rm=T)>0){
			chosen_row<-which(trends$lat==lat_steps[a] & trends$lon==lon_steps[b])
			df1<-as.data.frame(cbind(wind,owa_tsf,yearsf))
			df1<-df1[complete.cases(df1),]
			mod1<-gls(  wind ~ owa_tsf  ,data=df1)
			tryCatch( mod2<-gls( wind ~ owa_tsf  ,correlation=corAR1(form=~yearsf),data= df1) ,error=function(e) print('test') )
			pval<-1
			if(exists('mod2')){
					pval<-anova(mod1,mod2)$p[2]
			} 

			if(pval <0.05){
					trends[chosen_row,3]<-summary(mod2)$t[2,1]
					trends[chosen_row,4]<-summary(mod2)$t[2,4]
					trends[chosen_row,5]<-2
					rm('mod2')	
				
			} else {
				trends[chosen_row,3]<-summary(mod1)$t[2,1]
				trends[chosen_row,4]<-summary(mod1)$t[2,4]
				trends[chosen_row,5]<-1
				if(exists('mod2')){
					rm(mod2)
				}
				
			}
		
			
			}
	}
	print(a)
}

point_shape<-trends[which(trends$odds<=0.05),5]


library(ggplot2)
df<-as.data.frame(trends)
ggplot()+
geom_tile(data=df,aes(x=lon,y=lat,fill=trends))+
scale_fill_gradient2(high='red',low='blue',mid='white',na.value='white')+
geom_point(data=df[which(trends$odds<=0.05),],aes(x=lon,y=lat),shape=point_shape)
# Years with greater open water area appear to precede years with high winds 

fit1<-gls(wind_tsf ~ owa_tsf)
fit2<-gls(wind_tsf ~ owa_tsf, corr=corAR1(form=~yearsf))
anova(fit1,fit2)
# An autocorrelated fit does not appear to be preferred

