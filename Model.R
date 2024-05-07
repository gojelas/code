library(StMoMo)
library(demography)
library(urca)
library(tseries)


##########################
### mortality database ###
##########################
#data <- read.demogdata("/media/gouthonabiodunjean-luc/GOJELA'S/Thèse/Correction du modèle/données/Mx_1x1.txt",
#                       "/media/gouthonabiodunjean-luc/GOJELA'S/Thèse/Correction du modèle/données/Exposures_1x1.txt", 
#                      type="mortality", label="France")

#data <- read.demogdata("F:/Thèse/Correction du modèle/données/Mx_1x1.txt",
#                       "F:/Thèse/Correction du modèle/données/Exposures_1x1.txt", 
#                       type="mortality", label="France")
data_fr <- read.demogdata("Mx_1x1.txt",
                          "Exposures_1x1.txt", 
                          type="mortality", label="France")

FraMaleData=StMoMoData(data_fr, series = names(data_fr$rate)[2], type ="central")

FraFeMaleData=StMoMoData(data_fr, series = names(data_fr$rate)[1], type ="central")



#### 1) heat wave indicator
#ind=read.csv("F:/Thèse/Correction du modèle/données/indicateur2.csv",row.names = 1)


#### 2) different temperature scenarios
#temp_=read.csv("F:/Thèse/Correction du modèle/données/Données_globale_fr.csv",row.names = 1)
#temp=temp_[1:70,]
#temp



###################################################
#### 3) données des mois à haute température. #####
###################################################
#write.csv2(df,"F:/Thèse/Correction du modèle/données/heat_month.csv",row.names = FALSE)
df1=read.csv2("heat_month.csv",row.names = 1)

df2=read.csv('rcp2_6_by_year.csv')

df3=read.csv('rcp4_5_by_year.csv')

df4=read.csv('rcp8_5_by_year.csv')

df=c(df1[,6],df2[15:dim(df2)[1],2])
df_=c(df1[,6],(df3[15:dim(df3)[1],2]))
df__=c(df1[,6],(df4[15:dim(df4)[1],2]))

df=data.frame(df,row.names = 1961:(1961+length(df)-1))
colnames(df)=c('RCP2.6')
df['RCP4.5']=df_
df['RCP8.5']=df__


###################################################
#### Fonctions utilisées dans la modélisation #####
###################################################
source(file="models_functions.R")
source(file = "Lee Carter.R")


##############################################
#### Re-ordonner les données de la France ####
##############################################
FraFeMaleData$Dxt=reshape_french_data(FraFeMaleData)$dxt
FraFeMaleData$Ext=reshape_french_data(FraFeMaleData)$ext
female_rate=FraFeMaleData$Dxt/FraFeMaleData$Ext


FraMaleData$Dxt=reshape_french_data(FraMaleData)$dxt
FraMaleData$Ext=reshape_french_data(FraMaleData)$ext
male_rate=FraMaleData$Dxt/FraMaleData$Ext




#######################################################
#### Recherche des meilleurs paramètres du modèle #####
#######################################################
#c=seq(0,1,by=0.01)
#x1=45:53
#x2=54:60
#best_param=get_best_params(FraFeMaleData,c,x1,x2,20:85,1980:2011,Ic)
#best_param$performance

#x1=45 avec 0:85+ on a x1=45, x2=54, c=0 pour les femmes
#x2=60
#c=0.05


### pour les hommes, x1=45, x2=60, c=54

######################################
#### Application avec les données ####
######################################
x1=45
x2=60
ages.fit=20:85
years.fit=1980:2011
Ic=df[(years.fit[1]-1961+1):(tail(years.fit,1)-1961+1),1]
a_m=mean(Ic)
a_c=65
c_x=0.54
I_t=pmax((Ic-mean(Ic)),0)
fr_m=FraMaleData$Dxt[21:86,(years.fit[1]-1816+1):(tail(years.fit,1)-1816+1)]
fr_e=FraMaleData$Ext[21:86,(years.fit[1]-1816+1):(tail(years.fit,1)-1816+1)]
wei=genWeightMat(ages=ages.fit,years=years.fit)
model2=model_function(fr_m,fr_e,wei,Ic,a_m,c_x,a_c,ages.fit,years.fit,x1,x2)


#### Critère de performance du modèle
get_criterion(model2,FraMaleData,ages.fit,years.fit)


#### Back-testing
years.fit=1980:2011
ages.fit=20:85
#mx_test=male_rate[21:86,197:204]
best_model=get_estimation(FraMaleData,ages = ages.fit,
                          years = years.fit,c = 0.54,
                          a_c = 65,x1 = 45,x2 = 60)
mu_pred=get_predict(best_model,years_pred = 2012:2019,3)
mx_test=FraMaleData$Dxt[21:86,(2012-1816+1):(2019-1816+1)]/FraMaleData$Ext[21:86,(2012-1816+1):(2019-1816+1)]

mape_pred=mean(abs(mu_pred$mxt_pred-mx_test)/mx_test)
mape_pred




####################################        
### Comparaison des mxt estimés ####
####################################

ages.fit=20:85
years.fit=1980:2011

mx=data_fr$rate$male[ages.fit+1,]

a=55
mx_a=ts(mx[a,(1980-1816+1):(2019-1816+1)],start=years.fit[1],frequency = 1)
mod_mxt=fit_model(model2,ages.fit,years.fit)$uxt_fit


mu_fit_a=ts(mod_mxt[a,],start = years.fit[1],frequency = 1)
#spo_fit=ts(spo_mxt[a,],start = years.fit[1],frequency = 1)
plot(mx_a,type='l',xlim=c(1979,2019),lwd=1.5,ylim=c((min(mx_a)-0.0009),(max(mx_a)+0.0009)))
lines(mu_fit_a,col="blue")
lines(ts(mu_pred$mxt_pred[a,],start=2012),type='l',col='red',lwd=3)
#lines(spo_fit,type='l',col='blue',lwd=3)





#######################################
######  Intervalle de confiance  ######
#######################################

##############################################
### Intervalle de confiance par simulation des trajectoires futures des kt et de gamma ###
# Cela est basé sur une hypothèse de marche aléatoire multivariée avec drift pour le vecteur 
# (kt1,kt2,kt2,kt3,kt4) et le paramètre kt5 simulé séparément.
##############
sim2=sim_ic_method2(model2,ages.fit = 20:85,years.fit = 1980:2011,
                years.pred = 2012:2040,x1 = 45, x2=60,c=0.54,N=500)

mu_pred=get_predict(best_model,years_pred = 2012:2040,1)
rate=FraMaleData$Dxt[21:86,(1980-1816+1):(2019-1816+1)]/FraMaleData$Ext[21:86,(1980-1816+1):(2019-1816+1)]

par(mfrow=c(2,2))
a=6
probs = c(90, 95)
r=ts(rate[a,],start = 1980)
#text(0.5,0.5,"",cex=2,font=2)
plot(r,col=1,xlim=c(1980,2040),
       ylim=c(min(mu_pred$mxt_pred[a,]),max(rate[a,])),
     xlab="years",ylab="central mortality rate",main="Male-age 25")
fan(sim2$mu_pred[,a,],start = 2012,n.fan = 4,probs = probs,
    fan.col = colorRampPalette(c("gray", "white")), type = "interval")
lines(ts(mu_pred$mxt_pred[a,],start = 2012),col="red")
legend(2009,0.0016,legend =c("observed rate","rate predicted by the model",
                             "central trajectories of simulations"),
       col=c("black","red","orange"),cex=0.7,lty=1,lwd = 2,
       box.lwd = 0.8,box.col="green",box.lty = 2)

a=26
probs = c(90, 95)
r=ts(rate[a,],start = 1980)
plot(r,col=1,xlim=c(1980,2040),
     ylim=c(min(mu_pred$mxt_pred[a,]),max(rate[a,])),
     xlab="years",ylab="central mortality rate",pch=20,main="Male-age 45")
fan(sim2$mu_pred[,a,],start = 2012,n.fan = 4,probs = probs,
    fan.col = colorRampPalette(c("gray", "white")), type = "interval")
lines(ts(mu_pred$mxt_pred[a,],start = 2012),col="red")
legend(2009,0.0052,legend =c("observed rate","rate predicted by the model",
                             "central trajectories of simulations"),
       col=c("black","red","orange"),cex=0.7,lty=1,lwd = 2,
       box.lwd = 0.8,box.col="green",box.lty = 2)


#text(line2user(line=mean(par('mar')[c(2, 4)]), side=2), 
#     line2user(line=2, side=3), 'Central mortality rate for Male', xpd=NA, cex=2, font=2)
a=46
probs = c(90, 95)
r=ts(rate[a,],start = 1980)
plot(r,col=1,xlim=c(1980,2040),
     ylim=c(min(mu_pred$mxt_pred[a,]),max(rate[a,])),
     xlab="years",ylab="central mortality rate",pch=20,main="Male-age 65")
fan(sim2$mu_pred[,a,],start = 2012,n.fan = 4,probs = probs,
    fan.col = colorRampPalette(c("gray", "white")), type = "interval")
lines(ts(mu_pred$mxt_pred[a,],start = 2012),col="red")
legend(2009,0.027,legend =c("observed rate","rate predicted by the model",
                            "central trajectories of simulations"),
       col=c("black","red","orange"),cex=0.7,lty=1,lwd = 2,
       box.lwd = 0.8,box.col="green",box.lty = 2)

a=56
probs = c(90, 95)
r=ts(rate[a,],start = 1980)
plot(r,col=1,xlim=c(1980,2040),
     ylim=c(min(mu_pred$mxt_pred[a,]),max(rate[a,])),
     xlab="years",ylab="central mortality rate",pch=20,main='Male-age 75')
fan(sim2$mu_pred[,a,],start = 2012,n.fan = 4,probs = probs,
    fan.col = colorRampPalette(c("gray", "white")), type = "interval")
lines(ts(mu_pred$mxt_pred[a,],start = 2012),col="red")
legend(2009,0.064,legend =c("observed rate","rate predicted by the model",
                            "central trajectories of simulations"),
       col=c("black","red","orange"),cex=0.7,lty=1,lwd = 2,
       box.lwd = 0.8,box.col="green",box.lty = 2)


############################
#### Intervalle de confiance par boostrap avec 2,5% des taux plus faible
# et plus élévés comme borne.
###########################
IC=int_conf(FraMaleData,20:85,1980:2011,2012:2070,0.54,45,60,B = 5000)
######
model=get_estimation(FraMaleData,20:85,1980:2011,c = 0.54,
                     a_c = 65, x1 = 45, x2 = 60)
mu_pred=get_predict(model,2012:2070,1)

m=IC
rate=FraMaleData$Dxt[21:86,(1980-1816+1):(2019-1816+1)]/FraMaleData$Ext[21:86,(1980-1816+1):(2019-1816+1)]
a=6
plot(ts(rate[a,],start = 1980),col=1,xlim=c(1980,2070),
     ylim=c(min(mu_pred$mxt_pred[a,]),max(rate[a,])))
lines(ts(mu_pred$mxt_pred[a,],start = 2012),col=5)
lines(ts(m$IC_max[a,],start = 2012),col=2)
lines(ts(m$mean_pred[a,],start = 2012),col=3)
lines(ts(m$IC_min[a,],start = 2012),col=4)



###############################################
### Intervalle de confiance par boostrap avec les niveaux de confiance


I2=int_conf(FraMaleData,20:85,1980:2011,2012:2040,0.54,45,60,B = 100)
######
model=get_estimation(FraMaleData,20:85,1980:2011,c = 0.54,
                     a_c = 65, x1 = 45, x2 = 60)
mu_pred_1=get_predict(model,2012:2040,1)
mu_pred_2=get_predict(model,2012:2040,2)
mu_pred_3=get_predict(model,2012:2040,3)

rate=FraMaleData$Dxt[21:86,(1980-1816+1):(2019-1816+1)]/FraMaleData$Ext[21:86,(1980-1816+1):(2019-1816+1)]
a=60
plot(ts(rate[a,],start = 1980),col=1,xlim=c(1980,2040),
     ylim=c(min(mu_pred$mxt_pred[a,]),max(rate[a,])))
plot(ts(mu_pred_1$mxt_pred[a,],start = 2012),col="orange")
lines(ts(mu_pred_2$mxt_pred[a,],start = 2012),col="green")
lines(ts(mu_pred_3$mxt_pred[a,],start = 2012),col="red")






##########################
#### Pour les femmes #####
##########################
x1=45
x2=0
ages.fit=20:85
years.fit=1980:2011
Ic=df[(years.fit[1]-1961+1):(tail(years.fit,1)-1961+1),1]
a_m=mean(Ic)
a_c=65
c_x=0
I_t=pmax((Ic-mean(Ic)),0)
fr_m=FraFeMaleData$Dxt[21:86,(years.fit[1]-1816+1):(tail(years.fit,1)-1816+1)]
fr_e=FraFeMaleData$Ext[21:86,(years.fit[1]-1816+1):(tail(years.fit,1)-1816+1)]
wei=genWeightMat(ages=ages.fit,years=years.fit)
model2_f=model_function(fr_m,fr_e,wei,Ic,a_m,c_x,a_c,ages.fit,years.fit,x1,x2)


get_criterion(model2_f,FraFeMaleData,ages.fit,years.fit)


#### Back-testing
years.fit=1980:2011
ages.fit=20:85
#mx_test=male_rate[21:86,197:204]
best_model_f=get_estimation(FraFeMaleData,ages = ages.fit,
                          years = years.fit,c = 0,
                          a_c = 65,x1 = 45,x2 = 60)
mu_pred_f=get_predict(best_model_f,years_pred = 2012:2040,1)
mx_test_f=FraFeMaleData$Dxt[21:86,(2012-1816+1):(2019-1816+1)]/FraFeMaleData$Ext[21:86,(2012-1816+1):(2019-1816+1)]

mape_pred=mean(abs(mu_pred_f$mxt_pred-mx_test_f)/mx_test_f)
mape_pred


sim2_f=sim_ic_method2(model2_f,ages.fit = 20:85,years.fit = 1980:2011,
                    years.pred = 2012:2040,x1 = 45, x2=60,c=0,N=500)

rate_f=FraFeMaleData$Dxt[21:86,(1980-1816+1):(2019-1816+1)]/FraFeMaleData$Ext[21:86,(1980-1816+1):(2019-1816+1)]

par(mfrow=c(2,2))
a=6
probs = c(90, 95)
r=ts(rate_f[a,],start = 1980)
#text(0.5,0.5,"",cex=2,font=2)
plot(r,col=1,xlim=c(1980,2040),
     ylim=c(min(mu_pred_f$mxt_pred[a,])-0.0001,max(rate_f[a,])),
     xlab="years",ylab="central mortality rate",main="Female-age 25")
fan(sim2_f$mu_pred[,a,],start = 2012,n.fan = 4,probs = probs,
    fan.col = colorRampPalette(c("gray", "white")), type = "interval")
lines(ts(mu_pred_f$mxt_pred[a,],start = 2012),col="red")
legend(2009,0.0007,legend =c("observed rate","rate predicted by the model",
                             "central trajectories of simulations"),
       col=c("black","red","orange"),cex=0.7,lty=1,lwd = 2,
       box.lwd = 0.8,box.col="green",box.lty = 2)

a=26
probs = c(90, 95)
r=ts(rate_f[a,],start = 1980)
plot(r,col=1,xlim=c(1980,2040),
     ylim=c(min(mu_pred_f$mxt_pred[a,])-0.0002,max(rate_f[a,])),
     xlab="years",ylab="central mortality rate",pch=20,main="Female-age 45")
fan(sim2_f$mu_pred[,a,],start = 2012,n.fan = 4,probs = probs,
    fan.col = colorRampPalette(c("gray", "white")), type = "interval")
lines(ts(mu_pred_f$mxt_pred[a,],start = 2012),col="red")
legend(2009,0.0023,legend =c("observed rate","rate predicted by the model",
                             "central trajectories of simulations"),
       col=c("black","red","orange"),cex=0.7,lty=1,lwd = 2,
       box.lwd = 0.8,box.col="green",box.lty = 2)


#text(line2user(line=mean(par('mar')[c(2, 4)]), side=2), 
#     line2user(line=2, side=3), 'Central mortality rate for Male', xpd=NA, cex=2, font=2)
a=46
probs = c(90, 95)
r=ts(rate_f[a,],start = 1980)
plot(r,col=1,xlim=c(1980,2040),
     ylim=c(min(mu_pred_f$mxt_pred[a,])-0.001,max(rate_f[a,])),
     xlab="years",ylab="central mortality rate",pch=20,main="Female-age 65")
fan(sim2_f$mu_pred[,a,],start = 2012,n.fan = 4,probs = probs,
    fan.col = colorRampPalette(c("gray", "white")), type = "interval")
lines(ts(mu_pred_f$mxt_pred[a,],start = 2012),col="red")
legend(2009,0.010,legend =c("observed rate","rate predicted by the model",
                            "central trajectories of simulations"),
       col=c("black","red","orange"),cex=0.7,lty=1,lwd = 2,
       box.lwd = 0.8,box.col="green",box.lty = 2)

a=56
probs = c(90, 95)
r=ts(rate_f[a,],start = 1980)
plot(r,col=1,xlim=c(1980,2040),
     ylim=c(min(mu_pred_f$mxt_pred[a,])-0.008,max(rate_f[a,])),
     xlab="years",ylab="central mortality rate",pch=20,main='Female-age 75')
fan(sim2_f$mu_pred[,a,],start = 2012,n.fan = 4,probs = probs,
    fan.col = colorRampPalette(c("gray", "white")), type = "interval")
lines(ts(mu_pred_f$mxt_pred[a,],start = 2012),col="red")
legend(2009,0.033,legend =c("observed rate","rate predicted by the model",
                            "central trajectories of simulations"),
       col=c("black","red","orange"),cex=0.7,lty=1,lwd = 2,
       box.lwd = 0.8,box.col="green",box.lty = 2)





