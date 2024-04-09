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
df=read.csv2("heat_month.csv",row.names = 1)



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
Ic=df[(years.fit[1]-1961+1):(tail(years.fit,1)-1961+1),6]
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
mx_test=male_rate[21:86,197:204]
best_model=get_estimation(FraMaleData,ages = ages.fit,
                          years = years.fit,c = 0.54,
                          a_c = 65,x1 = 45,x2 = 60)
mu_pred=get_predict(best_model,years_pred = 2012:2019)
mape_pred=mean(abs(mu_pred$mxt_pred-mx_test)/mx_test)
mape_pred



#### A titre comparatif, on considère les résultats 
# avec le modèle de LC 
ages.fit=20:85
years.fit=1980:2011
mx_test=male_rate[21:86,197:204]
wei=genWeightMat(ages=ages.fit,years=years.fit)
LCfit=fit(LC,data=FraMaleData,ages.fit = ages.fit, 
          years.fit = years.fit, wxt = wei)
### Forecasting with PLAT
fore_lc=forecast(LCfit,h=50 )
plot(fore_lc)
mape_lc=mean(abs(fore_lc$rates[,1:8]-mx_test)/mx_test)
mape_lc









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









################################################################
######  Intervalle de confiance par la méthode boostrap   ######
###############################################################


IC=int_conf(FraMaleData,20:85,1980:2011,2012:2019,0.54,45,60)
######







