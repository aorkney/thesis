# This script provides an illustrative example of fuzzy-logic clustering (c-means)
# The approach is applied to aid in the categorisation of different phytoplankton assemblages,
# based on their absorption spectra
# absorption spectral data is available through the BODC (British Oceanographic Data Centre)
# [doi:10.5285/97daa7ea-8792-6cff-e053-6c86abc0dd46]
# [doi:10.5285/982b6da2-7e11-060a-e053-6c86abc09389]
# [doi:10.5285/982b6da2-7e12-060a-e053-6c86abc09389]

setwd('D:/Documents/CHEMTAX')
# You will probably store your data in a different directory; adjust accordingly 

data<-read.csv('Ph_abs_w_clusters_and_HPLC_and_microscope.csv')
# You will probably have named your dataset differently, and you will need to stitch this dataset 
# with information from the 'microscope_alone.csv' file available on the Github repository
# An example of how to find matching rows with the BODC data and the microscope_alone.csv file is coded below:

# Convert cruise identities in microscope alone to a vector matching the BODC convention
#cruise_ID<-as.character(microscope_alone$Cruise)
#cruise_ID[which(cruise_ID=='Summer_2018')]<-'JR17006'
#cruise_ID[which(cruise_ID=='Summer_2017')]<-'JR16006'
#cruise_ID[which(cruise_ID=='Spring_2018')]<-'HH180423'

# Find matching rows between 'data1' (some csv file from BODC), and 'microscope_alone.csv'
# The list of matches is stored in 'choices' and represents entries with synonymous event number, depth, lat, julian date, and cruise identity
#choice<-list()
#for(i in 1:(dim(microscope_alone)[1]) ){
#	choice[[i]]<-which( data1$Event==microscope_alone$Event[i] & data1$Depth==microscope_alone$Depth[i] & data1$Lat==microscope_alone$Lat[i] & data1$Julian==microscope_alone$Julian[i] & data1$Cruise==cruise_ID[i])	
#}

# Now that this is assumed to be completed, we will continue
aph<-data[,2:302] # Collect the hyperspectral absorption data (you may need to adjust indices depending on form of your csv file)
# we want the spectra between 400 and 700nm

library(ppclust) # We need this package to perform clustering
# if you do not have it, run the code 'install.packages('ppclust')'


c<-4 # Number of clusters to resolve
m<-2 # Fuzziness

z<-aph 
z<-aph/rowSums(aph) # Normalise the absorption spectra so that we do not explore variation in magnitude


# startign guesses for cluster centres with the k-means++ algorithm
v <- inaparc::kmpp(z,k=c)$v

# starting guesses for membership matrix
u <- inaparc::imembrand(nrow(z),k=c)$u

results<-fcm(z,centers=v,memberships=u,m=m) # Perform the fuzzy-logic clustering 
# This may take a few seconds, depending on the speed of your computer

PCA<-prcomp(z) # Perform a principal component decomposition on the normalised spectra 


df2<-as.data.frame(cbind(PCA$x[,c(1,2,3)],results$u)) # Define a data frame of the results
colnames(df2)[c(4,5,6,7)]<-c('c1','c2','c3','c4')
colnames(df2)[1:3]<-c('V1','V2','V3')

cols<-rgb(df2$c1,df2$c2,df2$c3) # Create a useful colour vector of the memberships (if you have red-green colourblindness you may wish to explore alternatives)





library(ggplot2) # This is a plotting package; if you do not have it
# then run 'install.packages('ggplot2')'

pc12a<-
ggplot(data=df2) +
geom_point(aes(x= V1, y= V2), fill='white',col='black',stroke=1,size=3,shape=21)+
geom_point(aes(x= V1, y= V2), col='black',fill='green',size=5,stroke=1,data=df2[which(data$diatom=='1' & is.na(data$phaeocystis)==T & is.na(data$dinoflagellate)==T),],shape=22)+
geom_point(aes(x= V1, y= V2), col='black',fill='red',size=5,stroke=1,data=df2[which(data$phaeocystis=='1' & is.na(data$diatom)==T & is.na(data$dinoflagellate)==T),],shape=23)+
geom_point(aes(x= V1, y= V2), col='black',fill='black',size=5,stroke=1,data=df2[which(data$dinoflagellate=='1' & is.na(data$diatom)==T & is.na(data$phaeocystis)==T),],shape=24)+
labs(x='PC1',y='PC2')
pc13a<-
ggplot(data=df2) +
geom_point(aes(x= V1, y= V3), fill='white',col='black',stroke=1,size=3,shape=21)+
geom_point(aes(x= V1, y= V3), col='black',fill='green',size=5,stroke=1,data=df2[which(data$diatom=='1' & is.na(data$phaeocystis)==T & is.na(data$dinoflagellate)==T),],shape=22)+
geom_point(aes(x= V1, y= V3), col='black',fill='red',size=5,stroke=1,data=df2[which(data$phaeocystis=='1' & is.na(data$diatom)==T & is.na(data$dinoflagellate)==T),],shape=23)+
geom_point(aes(x= V1, y= V3), col='black',fill='black',size=5,stroke=1,data=df2[which(data$dinoflagellate=='1' & is.na(data$diatom)==T & is.na(data$phaeocystis)==T),],shape=24)+
labs(x='PC1',y='PC3')


pc12b<-
ggplot(data=df2) +
geom_point(aes(x= V1, y= V2), fill=cols,col=cols,stroke=1,size=3,shape=21)+
geom_point(aes(x= V1, y= V2), col='black',fill='green',size=5,stroke=1,data=df2[which(data$diatom=='1' & is.na(data$phaeocystis)==T & is.na(data$dinoflagellate)==T),],shape=0)+
geom_point(aes(x= V1, y= V2), col='black',fill='red',size=5,stroke=1,data=df2[which(data$phaeocystis=='1' & is.na(data$diatom)==T & is.na(data$dinoflagellate)==T),],shape=5)+
geom_point(aes(x= V1, y= V2), col='black',fill='black',size=5,stroke=1,data=df2[which(data$dinoflagellate=='1' & is.na(data$diatom)==T & is.na(data$phaeocystis)==T),],shape=2)+
labs(x='PC1',y='PC2')
pc13b<-
ggplot(data=df2) +
geom_point(aes(x= V1, y= V3), fill=cols,col=cols,stroke=1,size=3,shape=21)+
geom_point(aes(x= V1, y= V3), col='black',fill='green',size=5,stroke=1,data=df2[which(data$diatom=='1' & is.na(data$phaeocystis)==T & is.na(data$dinoflagellate)==T),],shape=0)+
geom_point(aes(x= V1, y= V3), col='black',fill='red',size=5,stroke=1,data=df2[which(data$phaeocystis=='1' & is.na(data$diatom)==T & is.na(data$dinoflagellate)==T),],shape=5)+
geom_point(aes(x= V1, y= V3), col='black',fill='black',size=5,stroke=1,data=df2[which(data$dinoflagellate=='1' & is.na(data$diatom)==T & is.na(data$phaeocystis)==T),],shape=2)+
labs(x='PC1',y='PC3')


library(ggpubr) # this will allows us to view the plots together 
# if you do not have this package run the code 'install.packages('ggpubr')'

ggarrange(pc12a,pc13a,pc12b,pc13b) 
# Note that the colours of the fuzzy clustering are arbitrary

# What are the 'first part the post' results of the memberships to the clusters?
c1_picks<-which(df2$c1 > df2$c2 & df2$c1 > df2$c3 & df2$c1 > df2$c4)
c2_picks<-which(df2$c2 > df2$c1 & df2$c2 > df2$c3 & df2$c2 > df2$c4)
c3_picks<-which(df2$c3 > df2$c2 & df2$c3 > df2$c1 & df2$c3 > df2$c4)
c4_picks<-which(df2$c4 > df2$c2 & df2$c4 > df2$c3 & df2$c4 > df2$c1)


# The following code will generate the mean absorption spectra for each cluster, weighted by cluster memberships 
# The line colours are arbitrary; they may not match the cluster colours in the previous plot

norm<-aph[c1_picks,]/aph$X621[c1_picks]
redact<-which(aph$X621[c1_picks]<0.002) # quality control for extremely low biomass
if(length(redact)>0){
	spectrum<-norm[1,]
	for(i in 1:length(norm)){
		spectrum[i]<-weighted.mean(norm[-redact,i],w=df2$c1[c1_picks][-redact])
	}
} else {
	spectrum<-norm[1,]
	for(i in 1:length(norm)){
		spectrum[i]<-weighted.mean(norm[,i],w=df2$c1[c1_picks])
	}
}
plot(400:700,as.numeric(spectrum),type='l',col='red',ylim=c(0,14),ylab='normalised absorption', xlab='wavelength (nm)',lwd=3)
par(new=T)

norm<-aph[c2_picks,]/aph$X621[c2_picks]
redact<-which(aph$X621[c2_picks]<0.002) # quality control for extremely low biomass
if(length(redact)>0){
	spectrum<-norm[1,]
	for(i in 1:length(norm)){
		spectrum[i]<-weighted.mean(norm[-redact,i],w=df2$c2[c2_picks][-redact])
	}
} else {
	spectrum<-norm[1,]
	for(i in 1:length(norm)){
		spectrum[i]<-weighted.mean(norm[,i],w=df2$c2[c2_picks])
	}
}
plot(400:700,as.numeric(spectrum),type='l',col='green',ylim=c(0,14),ylab='normalised absorption', xlab='wavelength (nm)',lwd=3)
par(new=T)

norm<-aph[c3_picks,]/aph$X621[c3_picks]
redact<-which(aph$X621[c3_picks]<0.002) # quality control for extremely low biomass
if(length(redact)>0){
	spectrum<-norm[1,]
	for(i in 1:length(norm)){
		spectrum[i]<-weighted.mean(norm[-redact,i],w=df2$c3[c3_picks][-redact])
	}
} else {
	spectrum<-norm[1,]
	for(i in 1:length(norm)){
		spectrum[i]<-weighted.mean(norm[,i],w=df2$c3[c3_picks])
	}
}
plot(400:700,as.numeric(spectrum),type='l',col='blue',ylim=c(0,14),ylab='normalised absorption', xlab='wavelength (nm)',lwd=3)
par(new=T)

norm<-aph[c4_picks,]/aph$X621[c4_picks]
redact<-which(aph$X621[c4_picks]<0.002) # quality control for extremely low biomass
if(length(redact)>0){
	spectrum<-norm[1,]
	for(i in 1:length(norm)){
		spectrum[i]<-weighted.mean(norm[-redact,i],w=df2$c4[c4_picks][-redact])
	}
} else {
	spectrum<-norm[1,]
	for(i in 1:length(norm)){
		spectrum[i]<-weighted.mean(norm[,i],w=df2$c4[c4_picks])
	}
}
plot(400:700,as.numeric(spectrum),type='l',col='black',ylim=c(0,14),ylab='normalised absorption', xlab='wavelength (nm)',lwd=3)


# Now, let us explore an example where fuzzy-clustering may not work
# we shall avoid normalising our absorption dataset to its magnitude

PCA<-prcomp(aph)

z<-aph

# startign guesses with the k-means++ algorithm
v <- inaparc::kmpp(z,k=c)$v

# starting guesses for membership matrix
u <- inaparc::imembrand(nrow(z),k=c)$u

results<-fcm(z,centers=v,memberships=u,m=m)


df2<-as.data.frame(cbind(PCA$x[,c(1,2,3)],results$u))
colnames(df2)[c(4,5,6,7)]<-c('c1','c2','c3','c4')
colnames(df2)[1:3]<-c('V1','V2','V3')

cols<-rgb(df2$c1,df2$c2,df2$c3)


library(ggplot2) # This is a plotting package; if you do not have it
# then run 'install.packages('ggplot2')'

pc12a<-
ggplot(data=df2) +
geom_point(aes(x= V1, y= V2), fill='white',col='black',stroke=1,size=3,shape=21)+
geom_point(aes(x= V1, y= V2), col='black',fill='green',size=5,stroke=1,data=df2[which(data$diatom=='1' & is.na(data$phaeocystis)==T & is.na(data$dinoflagellate)==T),],shape=22)+
geom_point(aes(x= V1, y= V2), col='black',fill='red',size=5,stroke=1,data=df2[which(data$phaeocystis=='1' & is.na(data$diatom)==T & is.na(data$dinoflagellate)==T),],shape=23)+
geom_point(aes(x= V1, y= V2), col='black',fill='black',size=5,stroke=1,data=df2[which(data$dinoflagellate=='1' & is.na(data$diatom)==T & is.na(data$phaeocystis)==T),],shape=24)+
labs(x='PC1',y='PC2')
pc13a<-
ggplot(data=df2) +
geom_point(aes(x= V1, y= V3), fill='white',col='black',stroke=1,size=3,shape=21)+
geom_point(aes(x= V1, y= V3), col='black',fill='green',size=5,stroke=1,data=df2[which(data$diatom=='1' & is.na(data$phaeocystis)==T & is.na(data$dinoflagellate)==T),],shape=22)+
geom_point(aes(x= V1, y= V3), col='black',fill='red',size=5,stroke=1,data=df2[which(data$phaeocystis=='1' & is.na(data$diatom)==T & is.na(data$dinoflagellate)==T),],shape=23)+
geom_point(aes(x= V1, y= V3), col='black',fill='black',size=5,stroke=1,data=df2[which(data$dinoflagellate=='1' & is.na(data$diatom)==T & is.na(data$phaeocystis)==T),],shape=24)+
labs(x='PC1',y='PC3')


pc12b<-
ggplot(data=df2) +
geom_point(aes(x= V1, y= V2), fill=cols,col=cols,stroke=1,size=3,shape=21)+
geom_point(aes(x= V1, y= V2), col='black',fill='green',size=5,stroke=1,data=df2[which(data$diatom=='1' & is.na(data$phaeocystis)==T & is.na(data$dinoflagellate)==T),],shape=0)+
geom_point(aes(x= V1, y= V2), col='black',fill='red',size=5,stroke=1,data=df2[which(data$phaeocystis=='1' & is.na(data$diatom)==T & is.na(data$dinoflagellate)==T),],shape=5)+
geom_point(aes(x= V1, y= V2), col='black',fill='black',size=5,stroke=1,data=df2[which(data$dinoflagellate=='1' & is.na(data$diatom)==T & is.na(data$phaeocystis)==T),],shape=2)+
labs(x='PC1',y='PC2')
pc13b<-
ggplot(data=df2) +
geom_point(aes(x= V1, y= V3), fill=cols,col=cols,stroke=1,size=3,shape=21)+
geom_point(aes(x= V1, y= V3), col='black',fill='green',size=5,stroke=1,data=df2[which(data$diatom=='1' & is.na(data$phaeocystis)==T & is.na(data$dinoflagellate)==T),],shape=0)+
geom_point(aes(x= V1, y= V3), col='black',fill='red',size=5,stroke=1,data=df2[which(data$phaeocystis=='1' & is.na(data$diatom)==T & is.na(data$dinoflagellate)==T),],shape=5)+
geom_point(aes(x= V1, y= V3), col='black',fill='black',size=5,stroke=1,data=df2[which(data$dinoflagellate=='1' & is.na(data$diatom)==T & is.na(data$phaeocystis)==T),],shape=2)+
labs(x='PC1',y='PC3')


library(ggpubr)
ggarrange(pc12a,pc13a,pc12b,pc13b)
# It is evident that the clusterng simply partitions variation across PC1, and that this is not clearly related to the phytoplankton
# assemblages we hoped to characterise



