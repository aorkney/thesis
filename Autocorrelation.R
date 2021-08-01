# This is an example script 
# Illustrating how to investigate autocorrelated time series


library(nlme) # We will need this package for 
# Generalised least squares fits; if you do not have this package
# run 'install.packages('nlme')'



t<-c(1:100) # Time steps 
x<-0.05*t + sample(t)/10 
# x evolves with t, but also has random variation
y<-x
y[-1]<-x[-2]+ x[-length(x)]
y<-y-mean(x)
# y evolves with t too, but it is autocorrelated over a range of
# 1 time step


cor(x[-1],x[-length(x)]) # x has little correlation between time steps
cor(y[-1],y[-length(y)]) # y has a lot more



df<-as.data.frame(cbind(t,x,y)) # make dataframe


fitx<-gls(x~t,data=df) 
# Generalised least squares fit of x as function of t

fity<-gls(y~t,data=df)
# Likewise for y 



library(ggplot2) # this is a plotting package
# if you do not have it installed run 'install.packages('ggplot2')'

df<-cbind(df,predict(fitx),predict(fity))
# make dataframe
colnames(df)[4:5]<-c('px','py')

time_series<-
ggplot()+
geom_line(aes(x=t,y=x),data=df,size=2,col='green')+
geom_line(aes(x=t,y=y),data=df,size=1,col='black')+
geom_line(aes(x=t,y=px),data=df,col='red',linetype='dashed',size=3)+
geom_line(aes(x=t,y=py),data=df,col='red',size=1)+
labs(x='time',y='temperautre')+
theme(axis.line=element_line(size=1),
      axis.title=element_text(size=20)
)
# Make a plot called 'time_series'

p_valuex<-summary(fitx)$tTable[2,4]
p_valuey<-summary(fity)$tTable[2,4]
# These are the p values fo the gls fits, with no autocorrelation represented



df2<-as.data.frame( cbind( residuals(fitx)[-1],residuals(fitx)[1:(length(x)-1)],
residuals(fity)[-1],residuals(fity)[1:(length(y)-1)] ))
# dataframe of residuals and lagged residuals of gls fits for x and y as function t

autocorr<-
ggplot()+
geom_point(aes(x=V2,y=V1),data=df2,size=5,col='green')+
geom_point(aes(x=V4,y=V3),data=df2,size=2)+
labs(x='Residual at time t',y='Residual at time t+1')+
theme(axis.line=element_line(size=1),
      axis.title=element_text(size=20)
)
# plot of the residuals, showing evidence of autocorrelation of y as function t


fit2x<-gls( x ~ t  ,correlation=corAR1(form=~t),data= df) 
fit2y<-gls( y ~ t  ,correlation=corAR1(form=~t),data= df) 
# produce gls fits for x and y as function t with autocorrelation of time step 1


p_valuex2<-summary(fit2x)$tTable[2,4]
p_valuey2<-summary(fit2y)$tTable[2,4]
# p-values of autocorrelated gls fits of x and y as function t 


anova(fitx,fit2x)
anova(fity,fit2y)
# use anova to determine whether standard gls or autocorrelated gls fits are
# superior description of relationships of x as function t and y as function t

library(ggpubr)
dev.new(width=12,height=6,unit='cm')
ggarrange(time_series,autocorr)
# Diagnostic plots




