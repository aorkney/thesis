# This script is going to apply 'Redundancy Analysis', a multivariate statistical method that investigates 
# the variance a set of 'explanatory variables' can constrain in a target dataset of dependent variables. 
# We will also see that it is possible to test which individual explanatory variables constrain variance, 
# and to systematically compare the results of different Redundancy Analyses. 
# This is applied to allow us to test whether geography explains variation in weather conditions in the Barents Sea
# in different years, and whether some years are more similar to one another than others. 

library(vegan) # we need this package to perform the analysis;
# if you do not have it run 'install.packages('vegan')'

setwd('D:/Documents/Github/OCCCI')
# You may store your dataset in a different directory; adjust accordingly

dataset<-read.csv('SST_and_storms.csv')
dataset<-dataset[,-1]
# The dataset has been loaded

Geography_2002<-dataset[which(dataset$years=='2002'),c(1,2)]
Weather_2002<-dataset[which(dataset$years=='2002'),c(4,5)]
# Let's define two different blocks of data, one represent geography, another representing 
# the conditions (stormy days and sea surface temperature; we'll just call this 'Weather')
# We're going to do this for years 2002, 2003 and 2018

Geography_2003<-dataset[which(dataset$years=='2003'),c(1,2)]
Weather_2003<-dataset[which(dataset$years=='2003'),c(4,5)]

Geography_2018<-dataset[which(dataset$years=='2018'),c(1,2)]
Weather_2018<-dataset[which(dataset$years=='2018'),c(4,5)]

# Some of these datasets contain missing values. 
# We are going to investigate the similarity of environmental relationships between the same 
# lat-lon locations in different years, so we need to make sure we 
# only consider samples with complete data. 

check_list<-cbind(Weather_2002,Weather_2003,Weather_2018)
Geography_2002<-Geography_2002[complete.cases(check_list),]
Geography_2003<-Geography_2003[complete.cases(check_list),]
Geography_2018<-Geography_2018[complete.cases(check_list),]

Weather_2002<-Weather_2002[complete.cases(check_list),]
Weather_2003<-Weather_2003[complete.cases(check_list),]
Weather_2018<-Weather_2018[complete.cases(check_list),]

# We will need the following function to explore the relation between
# the individual samples and the geographic vectors explaining variance in environmental
# conditions. 

angle <- function(x,y){
  dot.prod <- x%*%y 
  norm.x <- norm(x,type="2")
  norm.y <- norm(y,type="2")
  theta <- acos(dot.prod / (norm.x * norm.y))
  as.numeric(theta)
}

library(dendextend) # This package will allow us to compare different redundancy analyses
# If you do not have it, then run 'install.packages('dendextend')



# Now that this is finished, let us perform Redundancy Analysis on the weather in each year
# Seeking to explain it with geography
# I am going to use a boot-strapped resampling method, to 
# perform the analyses 100 times, with 300 randomly-drawn locations each time. 
# This will take some time to run. 

recipient<-matrix(NA,100,3)
for(g in 1:100){

	subset<-sample(c(1:dim(Weather_2002)[1]),size=300,replace=T) # Resample

	rda_2002<- rda(Weather_2002[subset,] ~ . ,Geography_2002[subset,]) # Redundancy analysis
	rda_2003<- rda(Weather_2003[subset,] ~ . ,Geography_2003[subset,])
	rda_2018<- rda(Weather_2018[subset,] ~ . ,Geography_2018[subset,])

	x<-scores(choices=c(1:2),rda_2002,scaling=3,display='sites') # Get location scores
	y<-scores(choices=c(1:2),rda_2002,scaling=3,display='bp') # Get explanatory variable vectors

	store<-matrix(NaN,dim(x)[1],dim(y)[1])
	for(i in 1: (dim(x)[1]) ){
		for(j in 1: (dim(y)[1]) ){
			store[i,j]<-(angle(t(as.matrix(x[i,])),as.matrix(y[j,])))
			 # Compute angles between location scores and vectors
		}
	}
	rda_2002_dend <- as.dendrogram(hclust(dist(store)))  # Summarise result as a dendrogram

	x<-scores(choices=c(1:2),rda_2003,scaling=3,display='sites')
	y<-scores(choices=c(1:2),rda_2003,scaling=3,display='bp')

	store<-matrix(NaN,dim(x)[1],dim(y)[1])
	for(i in 1: (dim(x)[1]) ){
		for(j in 1: (dim(y)[1]) ){
			store[i,j]<-(angle(t(as.matrix(x[i,])),as.matrix(y[j,])))
		}
	}
	rda_2003_dend <- as.dendrogram(hclust(dist(store)))


	x<-scores(choices=c(1:2),rda_2018,scaling=3,display='sites')
	y<-scores(choices=c(1:2),rda_2018,scaling=3,display='bp')

	store<-matrix(NaN,dim(x)[1],dim(y)[1])
	for(i in 1: (dim(x)[1]) ){
		for(j in 1: (dim(y)[1]) ){
			store[i,j]<-(angle(t(as.matrix(x[i,])),as.matrix(y[j,])))
		}
	}
	rda_2018_dend <- as.dendrogram(hclust(dist(store)))

	corrs<-cor.dendlist(dendlist(rda_2002_dend,rda_2003_dend,rda_2018_dend))	
	recipient[g,1]<-corrs[1,2] # Store the data 
	recipient[g,2]<-corrs[2,3]
	recipient[g,3]<-corrs[1,3]
	print(g) # Update me on progress 
}

combined_histograms<-
c(recipient[,1],recipient[,2],recipient[,3]) # organise data

combined_histograms<-cbind(
combined_histograms,
c(rep('1 2002 x 2003',length(recipient[,1])),
rep('1 2003 x 2018',length(recipient[,1])),
rep('1 2002 x 2018',length(recipient[,1]))
))
colnames(combined_histograms)<-c('Correlation','Combination')
# Label the organised data 

combined_histograms<-as.data.frame(combined_histograms)
combined_histograms$Combination<-as.factor(combined_histograms$Combination)
combined_histograms$Correlation<-as.numeric(as.character(combined_histograms$Correlation))
# Housekeeping 

library(ggplot2) # This is a plotting package; if you do not have it
# run 'install.packages('ggplot2')

plot_1<-
ggplot(combined_histograms, aes(y=Combination, x=Correlation, color=Combination)) + 
  geom_violin(trim=F,size=1,col='darkgrey',size=2)+
geom_boxplot(width=0.05,size=1,col='black')+ 
lims(x=c(0,1))+ 
theme_bw() + theme(axis.text.x=element_text(size=16, angle=90, vjust=0.3,colour='black',face = "bold"),
			axis.title.x=element_text(size=24,colour='black',face='bold'),
			axis.title.y=element_text(size=24,colour='black',face='bold'),
                     axis.text.y=element_text(size=16,colour='black',face = "bold"),
                     plot.title=element_text(size=16,hjust=.5), legend.position="none",panel.background = element_blank())+
	labs(fill = "Effect Size")
# Violin plot produced




	rda_2002<- rda(Weather_2002 ~ . ,Geography_2002)
	rda_2003<- rda(Weather_2003 ~ . ,Geography_2003)
	rda_2018<- rda(Weather_2018 ~ . ,Geography_2018)
# I have performed redundancy analysis with the full datasets now
# I want to not use permutation analysis to determine the significance of latitude and longitude
# as explanatory variables

rda_2002_pure_effects<-matrix(NA,2,1)
for(i in 1:2){
	temporary_rda<-anova(rda( Weather_2002 ~Geography_2002[,i],Z=Geography_2002[-i]))
	if(temporary_rda$"Pr(>F)" [1] <0.05){
		rda_2002_pure_effects[i,]<-temporary_rda$V[1]/sum(temporary_rda$V)
	} else{
		rda_2002_pure_effects[i,]<-'NA'
	}
}
# This loop calculates 'pure' effects, for each explanatory variable we first
# contrive a model removing all the variance it explains, so that we can 
# see whether the remaining one explains the residual
# this should proxy the effect 'purely' explaiend by the variable of interest

rda_2003_pure_effects<-matrix(NA,2,1)
for(i in 1:2){
	temporary_rda<-anova(rda( Weather_2003 ~Geography_2003[,i],Z=Geography_2003[-i]))
	if(temporary_rda$"Pr(>F)" [1] <0.05){
		rda_2003_pure_effects[i,]<-temporary_rda$V[1]/sum(temporary_rda$V)
	} else{
		rda_2003_pure_effects[i,]<-'NA'
	}
}

rda_2018_pure_effects<-matrix(NA,2,1)
for(i in 1:2){
	temporary_rda<-anova(rda( Weather_2018 ~Geography_2018[,i],Z=Geography_208[-i]))
	if(temporary_rda$"Pr(>F)" [1] <0.05){
		rda_2018_pure_effects[i,]<-temporary_rda$V[1]/sum(temporary_rda$V)
	} else{
		rda_208_pure_effects[i,]<-'NA'
	}
}

pure_effect_matrix<-(cbind(as.numeric(rda_2002_pure_effects),as.numeric(rda_2003_pure_effects),as.numeric(rda_2018_pure_effects)))
rownames(pure_effect_matrix)<-c('Lat','Lon')
# 'pure effect matrix' is just the name of the object where I am storing the results


rda_2002_term_effects<-anova(rda_2002, by="term")
rda_2003_term_effects<-anova(rda_2003, by="term")
rda_2018_term_effects<-anova(rda_2018, by="term")
# Now we can use analysis of variance to look at the term-wise effects
# This explores variation in one term at a time

rda_2002_terms<-rda_2002_term_effects$Variance[1:2]/sum(rda_2002_term_effects$Variance)
rda_2002_terms[which(rda_2002_term_effects$"Pr(>F)"[1:2] >0.05)]<-NA
rda_2003_terms<-rda_2003_term_effects$Variance[1:2]/sum(rda_2003_term_effects$Variance)
rda_2003_terms[which(rda_2003_term_effects$"Pr(>F)"[1:2] >0.05)]<-NA
rda_2018_terms<-rda_2018_term_effects$Variance[1:2]/sum(rda_2018_term_effects$Variance)
rda_2018_terms[which(rda_2018_term_effects$"Pr(>F)"[1:2] >0.05)]<-NA

pure_effect_matrix2<-
cbind(pure_effect_matrix,rep(NA,2),rda_2002_terms,rda_2003_terms,rda_2018_terms)


rda_2002_margin_effects<-anova(rda_2002, by="margin")
rda_2003_margin_effects<-anova(rda_2003, by="margin")
rda_2018_margin_effects<-anova(rda_2018, by="margin")
# Now we're going to look at another alternative; 'marginal' effects
# In this case the effects of the variable of interest are quantified in a multivariate context

rda_2002_margin<-rda_2002_margin_effects$Variance[1:2]/sum(rda_2002_margin_effects$Variance)
rda_2002_margin[which(rda_2002_margin_effects$"Pr(>F)"[1:2] >0.05)]<-NA
rda_2003_margin<-rda_2003_margin_effects$Variance[1:2]/sum(rda_2003_margin_effects$Variance)
rda_2003_margin[which(rda_2003_margin_effects$"Pr(>F)"[1:2] >0.05)]<-NA
rda_2018_margin<-rda_2018_margin_effects$Variance[1:2]/sum(rda_2018_margin_effects$Variance)
rda_2018_margin[which(rda_2018_margin_effects$"Pr(>F)"[1:2] >0.05)]<-NA

pure_effect_matrix3<-
cbind(pure_effect_matrix2,rep(NA,2),rda_2002_margin,rda_2003_margin,rda_2018_margin)

colnames(pure_effect_matrix3)<-
c('2002 P','2003 P','2018 P', '', '2002 T','2003 T','2018 T','','2002 M','2003 M','2018 M' )

library(reshape2) # This allows us to manipulate matrixes of data easily
# If you do not have this package then run 'install.packages('reshape2')'
longData<-melt(t(pure_effect_matrix3)[c(1:3,5:7,9:11),])
longData$Var2<-as.character(longData$Var2)
# organising the data


plot_2<-
ggplot(longData,aes(x=Var2,y=Var1))+
	coord_fixed(ratio=1,xlim=NULL,ylim=NULL,expand=T,clip='on')+
	geom_tile(aes(fill=value))+
	scale_fill_gradientn(colours = c('grey','black'),na.value="white",breaks=c(0,0.25,0.5,0.75)) +
	geom_text(aes(label = round(value, 1)),col='white',fontface='bold',size=4) +
	labs(fill = "Correlation")+
	xlab('')+
	ylab('')+
  	theme_bw() + 
	theme( axis.text=element_text(size=16,colour='black',face = "bold"),
		axis.text.x=element_text(size=16, angle=90, vjust=0.3,colour='black',face = "bold"),
		legend.title=element_text(size=24,colour='black',face='bold'),
		legend.key.width=unit(1,'cm'),
		plot.margin=unit(c(1,1,1,1),'cm'),
		legend.text=element_text(angle=90,size=12,face='bold',colour='black'),
              plot.title=element_text(size=16,hjust=.5), legend.position="bottom",legend.margin=margin(0,100,0,0),panel.background = element_blank())+
	labs(fill = "Explained variance")
# a labelled plot of the explained variance, with pure, term and marginal effects labelled as 'P', 'T','M',



library(ggpubr) # This package makes combining plots easy
# if you do not have it run 'install.packages('ggpubr')'
dev.new(height=8,width=15,unit='cm')
ggarrange(plot_1,plot_2,ncol=2,labels=c('a','b'),font.label=list(size=54))
# The results show us what we expected;
# that the relationships between weather and geography are more similar between 2002 and 2003, 
# than between 2002 and 2018 or 2003 and 2018. 
# we can also see that latitude is the dominant variable explaining variation in the conditions, 
# and that in 2018 lattiudinal gradients in weather conditions must be shallower, because the 
# degree of explained variance has dropped. 
# This is likely consistent with the area occupied by Atlantic inflow waters having expanded
# resulting in a more homogeneous sea on the shelf


