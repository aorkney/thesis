# This code will download OC-CCI data 
# from the thredds server 
# https://www.oceancolour.org/thredds/catalog-cci.html
# This is a multi-channel spectral ocean colour product,
# produced by merging data from several space-borne platforms
# The example will focus on downloading an image of the Barents Sea on the 6th of July, 2016
# You will need to set up a file path 'D:/Documents/Github/OCCCI/' on your computer, or change 'File_name' in the following code.
# The script will download 37.7MB of ocean colour data 
# This script was coded in R version 3.6.3 (2020-02-29) 'Holding the Windsock'

latitudes<-(seq(-89.97916666666666,89.97916666666667,length.out=4319))
longitudes<-(seq(-179.97916666666666,179.97916666666666,length.out=8639))
# I am defining the indices for latitude and longitude

lon_min<-which(abs(longitudes- -0)==min(abs(longitudes- -0)))
lon_max<-which(abs(longitudes- 60)==min(abs(longitudes- 60)))

lat_min<-which(abs(latitudes- -85)==min(abs(latitudes- -85)))
lat_max<-which(abs(latitudes- -70)==min(abs(latitudes- -70)))
# 'lat_max' is actually the southern latitude
# 'lat_min' is actually the northern latitude 
# These lines have defined the region of interest (Barents Sea)


year<-'2016'
wavebands<-c('412','443','490','510','560','665')
# waveband channels (nm)
month<-'07'
day<-'06'
# There is a conspicuous bloom of calcifying phytoplankton
# on this day in the Barents Sea 
# https://earthobservatory.nasa.gov/images/88316/the-barentssea-abloom


# Let's begin downloading the data from the thredds server 
# We will execute this as a for-loop, because there are multiple spectral channels to download

for(i in 1:length(wavebands)){

	URL_name<-
	paste('https://www.oceancolour.org/thredds/dodsC/cci/v5.0-release/geographic/daily/rrs/',year,
	'/ESACCI-OC-L3S-RRS-MERGED-1D_DAILY_4km_GEO_PML_RRS-',year,month,day,'-fv5.0.nc.ascii?Rrs_',wavebands[i],'[0:1:0]',
	'[',lat_min,':1:',lat_max,']','[',lon_min,':1:',lon_max,']',sep='')

	# The URL we will call


	# We wish to store the downloaded data, so we need a file name. You should adjust this according to your preferences.
	File_name<-paste('D:/Documents/Github/OCCCI/','Rrs',wavebands[i],'_',year,'_',month,'_',day,'.csv',sep='')

	download.file(URL_name,
	destfile=File_name)

	# download command
}

# We are now going to view the data we have downloaded
# The bloom of calcifying phytoplankton is conspicuous in the bottom-right of the image
# west of Novaya Zemlya 

for(i in 1:length(wavebands)){
	File_name<-paste('D:/Documents/Github/OCCCI/','Rrs',wavebands[i],'_',year,'_',month,'_',day,'.csv',sep='')

	datafile<-read.csv(File_name, skip = 12, header = F) # Call the data 

	data_store<-as.matrix(datafile[2:(dim(datafile)[1]-6),2:dim(datafile)[2]]) # Make the data into a matrix

	data_store[which(data_store>4e36)]<-NA # Eliminate invalid entries

	image( t(apply(data_store, 2, rev)) )

	# view the data 
	
	Sys.sleep(5)
}




