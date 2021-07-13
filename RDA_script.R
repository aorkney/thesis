# This is a script to perform 'Redundancy analysis' 
# We are going to investigate to what extent variation in geography
# constrains and explains variation in the weather conditions,
# sea surface temperature and storminess, in September 2002

setwd('D:/Documents/Github/OCCCI')
# You may store your dataset in a different directory; adjust accordingly

dataset<-read.csv('SST_and_storms.csv')
dataset<-dataset[,-1]
dataset<-dataset[complete.cases(dataset),]
# The dataset has been loaded

Geography<-dataset[which(dataset$years=='2002'),c(1,2)]
Weather<-dataset[which(dataset$years=='2002'),c(4,5)]
# Let's define two different blocks of data, one represent geography, another representing 
# the conditions (stormy days and sea surface temperature; we'll just call this 'Weather')


library('vegan')
subset<-sample(x=c(1:dim(Weather)[1]),size=500)
subset<-1:dim(Weather)[1]

rda<- rda(Weather[subset,] ~ . ,Geography[subset,])
plot(rda) # This is what a default plot looks like
# we can produce a prettier and more intuitive plot though

rda_plot<-cbind(
summary(rda)$sites[,1],
summary(rda)$sites[,2]
)
rda_plot<-as.data.frame(rda_plot)
colnames(rda_plot)<-c('x','y')
rda_plot$x<-rda_plot$x/max(abs(rda_plot$x))
rda_plot$y<-rda_plot$y/max(abs(rda_plot$y))

rda_segments<-as.data.frame(cbind(rep(0,dim(Geography[subset,])[2]),rep(0,dim(Geography[subset,])[2]),
summary(rda)$biplot[,1],summary(rda)$biplot[,2]))/max(abs(summary(rda)$biplot))
rda_segments$V3<-rda_segments$V3/max(abs(rda_plot$x))
rda_segments$V4<-rda_segments$V4/max(abs(rda_plot$y))
rda_segments$V4<- rda_segments$V4 # -
rda_segments$V1<- 0
rda_segments$V2<- 0
magnitudes<-( ((rda_segments$V3)^2)+
((rda_segments$V4)^2) )^.5
rda_segments$V3<-rda_segments$V3/max(magnitudes)
rda_segments$V4<-rda_segments$V4/max(magnitudes)
rda_segments$V3<-rda_segments$V3
rda_segments$V4<-rda_segments$V4 

var_names<-c('Latitude','Longitude')

angles<-(
180*(atan(
(rda_segments$V4-rda_segments$V2)/
(rda_segments$V3-rda_segments$V1)
))/pi)+90
angles[which(angles>90)]<-angles[which(angles>90)]+180

rda_unit<-
rda_segments[,1:2]+ 0.8*(
((rda_segments[,3:4]-rda_segments[,1:2])/sqrt(rowSums((rda_segments[,3:4]-rda_segments[,1:2])^2))) )
colnames(rda_unit)<-c('V3','V4')

rda_size<-sqrt(rowSums((rda_segments[,3:4]-rda_segments[,1:2])^2))*3


library(ggplot2) # This is a plotting package; if you do not have it
# then run 'install.packages('ggplot2')'

sst_plot<-
ggplot() +
geom_point(aes(x= x, y= y,col=Weather$sst[subset]), data=cbind(rda_plot,Weather$sst[subset]),stroke=1,size=4)+
scale_colour_gradientn(colours=rev(rainbow(10))[4:10],na.value='white')+
geom_segment(data=cbind(rda_segments[,1:2],rda_unit),aes(x=V1,y=V2,xend=V3,yend=V4),size=rda_size,
alpha=1,color='darkgrey')+
geom_text(data=cbind(rda_segments[,1:2],rda_unit),aes(x=V3,y=V4,label=paste( var_names ) ),size=6,color='black',parse=T)+
coord_equal()+
labs(colour=expression(paste('SST',,' ('^o,'C)')),y='RDA axis 2',x='RDA axis 1')



#dev.new()
storm_plot<-
ggplot() +
geom_point(aes(x= x, y= y,col=Weather$storminess[subset]), data=cbind(rda_plot,Weather$storminess[subset]),stroke=1,size=4)+
scale_colour_gradientn(colours=rev(rainbow(10))[4:10],na.value='white')+
geom_segment(data=cbind(rda_segments[,1:2],rda_unit),aes(x=V1,y=V2,xend=V3,yend=V4),size=rda_size,
alpha=1,color='darkgrey')+
geom_text(data=cbind(rda_segments[,1:2],rda_unit),aes(x=V3,y=V4,label=paste( var_names ) ),size=6,color='black',parse=T)+
coord_equal()+
labs(colour='storminess',y='RDA axis 2',x='RDA axis 1')






library(ggpubr)# If you have not installed this package run
# 'install.packages('ggpubr')', it allows easy combination of subplots
dev.new(width=12,height=6,unit='cm')
ggarrange(sst_plot,storm_plot)
# The plots show what is rather intuitive; that there is a strong dependency of sea surface temperature
# on latitude. There is also a minor dependency on longitude
# we can observe storms are more common in the southern Barents Sea in September 2002, 
# but that this pattern is somewhat independent of variation in temperature, 
# which is consistent with storms being stochastic weather events that travel into the Barents Sea
# through the southern Atlantic inflow region
# Type 'rda' to see a summary range of statistics, such as the proportion of constrained variance 


