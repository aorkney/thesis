# This script is going to use environmental and geographic variables, 
# from the Barents Sea in September 2002, to discriminate 'Northern'
# and 'Southern' classes. 
# We are going to employ Linear Discriminant Analysis, 
# and perform analyses within the context of a re-sampling procedure, 
# so that we can robustly determine the important variables that
# distinguish the classes. 

library(MASS) # we need this package to perform the analysis;
# if you do not have it run 'install.packages('MASS')'
library(reshape2) # we need this to organise data
# if you do not have it run 'install.packages('reshape2')'


setwd('D:/Documents/Github/OCCCI')
# You may store your dataset in a different directory; adjust accordingly

# You will need to download the dataset 'SST_and_storms.csv' from 
# the github repository if you have not already.

dataset<-read.csv('SST_and_storms.csv')
dataset<-dataset[,-1]
# The dataset has been loaded

Geography_2002<-dataset[which(dataset$years=='2002'),c(1,2)]
Weather_2002<-dataset[which(dataset$years=='2002'),c(4,5)]
# Let's define two different blocks of data, one represent geography, another representing 
# the conditions (stormy days and sea surface temperature; we'll just call this 'Weather')

library(ggplot2) # This is a plotting package; if you do not have it
# then run 'install.packages('ggplot2')'

sst_plot<-
ggplot() +
geom_tile(aes(x= lon, y= lat,fill=Weather_2002$sst), data=cbind(Geography_2002,Weather_2002$sst))+
scale_fill_gradientn(colours=rev(rainbow(10))[4:10],na.value='white')+
labs(fill=expression(paste('SST',,' ('^o,'C)')),y='Latitude',x='Longitude')+
geom_hline(yintercept=76,size=2,linetype='dashed')

# We can see that Northern waters are colder. 

Class<-rep('Southern',length(Geography_2002$lat))
Class[which(Geography_2002$lat>76)]<-'Northern'

# Let's say we want to be able to discriminate data that comes from the Southern and Northern Barents Sea.
# Perhaps SST will turn out to be an effective discriminant. 
# We are therefore going to perform a Linear Discriminant Analysis
# We are going to include storminess, because perhaps this is also predictive of latitude
# We are going to include longitude, because we know this is not predictive of latitude; will our analysis confirm that?

df<-cbind(as.factor(Class),Weather_2002,Geography_2002$lon)[complete.cases(Weather_2002),]
colnames(df)<-c('Class','storminess','sst','lon')

Analysis<-lda(Class ~ ., df)
plot(Analysis)
# It is clear to see that the Linear Discriminant that has been contrived separates Northern and Southern parts of the Barents
Analysis$scaling
# we can see the largest coefficient in this analysis is SST, conforming to our expectations. Neither storminess or longitude
# exhibit a strong capacity to predict weather a grid cell belongs to the Northern or Southern classes. 

# Let's explore the capacity of linear discriminant analysis to successfully predict class membership 
# with a re-sampling procedure

RNorth_PNorth<-list()
RNorth_PSouth<-list()
RSouth_PSouth<-list()
RSouth_PNorth<-list()

storminess_w<-list()
sst_w<-list()
lon_w<-list()

for(i in 1:100){
	train<-sample(1:dim(df)[1],100,replace=T)
	Analysis<-lda(Class ~ ., df, subset=train)
	prediction<-predict(Analysis, df[-train,])$class
	RNorth_PNorth[[i]]<-length(which(prediction==df[-train,]$Class & prediction=='Northern')) / length(which(df[-train,]$Class=='Northern') )
	RSouth_PNorth[[i]]<-length(which(prediction!=df[-train,]$Class & prediction=='Northern')) / length(which(df[-train,]$Class=='Southern') )
	RSouth_PSouth[[i]]<-length(which(prediction==df[-train,]$Class & prediction=='Southern')) / length(which(df[-train,]$Class=='Southern') )
	RNorth_PSouth[[i]]<-length(which(prediction!=df[-train,]$Class & prediction=='Southern')) / length(which(df[-train,]$Class=='Northern') )
	storminess_w[[i]]<-Analysis$scaling[1,]
	sst_w[[i]]<-Analysis$scaling[2,]
	lon_w[[i]]<-Analysis$scaling[3,]

	
}
RNorth_PNorth<-unlist(RNorth_PNorth)
RNorth_PSouth<-unlist(RNorth_PSouth)
RSouth_PSouth<-unlist(RSouth_PSouth)
RSouth_PNorth<-unlist(RSouth_PNorth)

storminess_w<-unlist(storminess_w)
sst_w<-unlist(sst_w)
lon_w<-unlist(lon_w)
	
Confusion_matrix<-matrix(NA,2,2)# actual along left hand side, predicted along top
colnames(Confusion_matrix)<-rownames(Confusion_matrix)<-c('North','South')
Confusion_matrix[1,1]<-mean(RNorth_PNorth)
Confusion_matrix[1,2]<-mean(RNorth_PSouth)
Confusion_matrix[2,1]<-mean(RSouth_PNorth)
Confusion_matrix[2,2]<-mean(RSouth_PSouth)

# Our confusion matrix is a handy was of presenting our results;
# we have found that our Linear Discriminant dominated by SST, on average
# identifies the Northern class correctly almost 90% of the time, and the Southern class about 95% of the time. 

df2<-melt(Confusion_matrix)
ggplot() +
geom_tile(aes(x= Var2, y= Var1 ,fill=value), data=df2 )+
scale_fill_gradientn(colours=rev(rainbow(10))[4:10],na.value='white')+
geom_text(data=df2, aes(x= Var2, y= Var1,label=round(value,digits=2) ),col='white',fontface='bold',size=20)+
labs(colour='Proportion predicted',y='Actual',x='Predicted')+
scale_y_discrete(limits = rev(levels(df2$Var2)))+
scale_x_discrete(expand = c(0, 0),position = 'top')+
theme( axis.text=element_text(size=16),axis.title=element_text(size=20))



# Let's now plot an example of the posterior probabilities that describe the prediction;
# this might provide insight for why prediction of the Northern class membership is less accurate. 


Geo<-Geography_2002[complete.cases(Weather_2002),]

df3<-cbind(Geo[-train,], predict(Analysis, df[-train,])$posterior[,1])
colnames(df3)[3]<-'Northerliness'
PN<-
ggplot() +
geom_tile(aes(x= lon, y= lat,fill=Northerliness), data=df3)+
scale_fill_gradientn(colours=rev(rainbow(10))[4:10],na.value='white')+
labs(fill='Prediction North',y='Latitude',x='Longitude')+
geom_hline(yintercept=76,size=2,linetype='dashed')

df4<-cbind(Geo[-train,], predict(Analysis, df[-train,])$posterior[,2])
colnames(df4)[3]<-'Southerliness'
PS<-
ggplot() +
geom_tile(aes(x= lon, y= lat,fill=Southerliness), data=df4)+
scale_fill_gradientn(colours=rev(rainbow(10))[4:10],na.value='white')+
labs(fill='Prediction South',y='Latitude',x='Longitude')+
geom_hline(yintercept=76,size=2,linetype='dashed')

library(ggpubr) # If you do not already have this package run
# install.packages('ggpubr'), it allows you to combine plots
dev.new(width=12,height=6,unit='cm')
ggarrange(PS,PN)

Analysis$scaling
# we can see the linear discriminant that divides the classes is dominated by SST;
# therefore the region west of Spitsbergen may be less likely to be correctly predicted
# as `Northern' because it has unusually high SST for its latitude

# The other scaling factors are typically small;
mean(lon_w) 
mean(storminess_w)

# Let's visualise the density distribution of the scaling factors from our distribution of re-samplings

pf<-as.data.frame(cbind(c(sst_w,lon_w,storminess_w),c(
rep('SST',length(sst_w)),rep('Longitude',length(lon_w)),rep('Storminess',length(storminess_w)) )) )
pf$V1<-as.numeric(as.character(pf$V1))

ggplot(pf,aes(y=V2,x=V1,color=V2))+
geom_violin(trim=F,size=2)+
geom_vline(xintercept=0,size=2,linetype='dashed')+
labs(x='Scaling coefficient',y='Variable',colour='Variable')+
theme_bw() + theme(axis.text.x=element_text(size=16, angle=90, vjust=0.3,colour='black',face = "bold"),
axis.title.x=element_text(size=24,colour='black',face='bold'),
axis.title.y=element_text(size=24,colour='black',face='bold'),
axis.text.y=element_text(size=16,colour='black',face = "bold"),
plot.title=element_text(size=16,hjust=.5), legend.position="none",panel.background = element_blank())





t.test(sst_w) ; t.test(storminess_w) ; t.test(lon_w)
# storminess plays no role significantly different from zero in the scalings across the analysis
# but longitude does

# This must mean that longitude of the samples is, when combined with temperature, predictive of latitude, 
# even though longitude is by itself agnostic to latitude 
# this could happen if longitude and SST covary significantly

dft<-cbind(Geography_2002$lat,Weather_2002,Geography_2002$lon)[complete.cases(Weather_2002),]
colnames(dft)<-c('lat','storminess','sst','lon')

summary(lm(dft$lat~dft$sst ))
summary(lm(dft$lat~ dft$lon))
summary(lm(dft$lat~dft$sst + dft$lon))

# simple linear regressions appear to demonstrate this; longtiude is not significantly associated with latitude
# by itself, but becomes so in combination with sst

# Let's see what the difference in predictive power is if we exclude longitude:

df<-cbind(as.factor(Class),Weather_2002)[complete.cases(Weather_2002),]
colnames(df)<-c('Class','storminess','sst')


RNorth_PNorth<-list()
RNorth_PSouth<-list()
RSouth_PSouth<-list()
RSouth_PNorth<-list()


for(i in 1:100){
	train<-sample(1:dim(df)[1],100,replace=T)
	Analysis<-lda(Class ~ ., df, subset=train)
	prediction<-predict(Analysis, df[-train,])$class
	RNorth_PNorth[[i]]<-length(which(prediction==df[-train,]$Class & prediction=='Northern')) / length(which(df[-train,]$Class=='Northern') )
	RSouth_PNorth[[i]]<-length(which(prediction!=df[-train,]$Class & prediction=='Northern')) / length(which(df[-train,]$Class=='Southern') )
	RSouth_PSouth[[i]]<-length(which(prediction==df[-train,]$Class & prediction=='Southern')) / length(which(df[-train,]$Class=='Southern') )
	RNorth_PSouth[[i]]<-length(which(prediction!=df[-train,]$Class & prediction=='Southern')) / length(which(df[-train,]$Class=='Northern') )


	
}
RNorth_PNorth<-unlist(RNorth_PNorth)
RNorth_PSouth<-unlist(RNorth_PSouth)
RSouth_PSouth<-unlist(RSouth_PSouth)
RSouth_PNorth<-unlist(RSouth_PNorth)

Confusion_matrix<-matrix(NA,2,2)# actual along left hand side, predicted along top
colnames(Confusion_matrix)<-rownames(Confusion_matrix)<-c('North','South')
Confusion_matrix[1,1]<-mean(RNorth_PNorth)
Confusion_matrix[1,2]<-mean(RNorth_PSouth)
Confusion_matrix[2,1]<-mean(RSouth_PNorth)
Confusion_matrix[2,2]<-mean(RSouth_PSouth)

df2<-melt(Confusion_matrix)
ggplot() +
geom_tile(aes(x= Var2, y= Var1 ,fill=value), data=df2 )+
scale_fill_gradientn(colours=rev(rainbow(10))[4:10],na.value='white')+
geom_text(data=df2, aes(x= Var2, y= Var1,label=round(value,digits=2) ),col='white',fontface='bold',size=20)+
labs(colour='Proportion predicted',y='Actual',x='Predicted')+
scale_y_discrete(limits = rev(levels(df2$Var2)))+
scale_x_discrete(expand = c(0, 0),position = 'top')+
theme( axis.text=element_text(size=16),axis.title=element_text(size=20))

# The predictive power drops slightly




