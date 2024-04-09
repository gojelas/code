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
data_fr <- read.demogdata("/Users/gojelastat/Desktop/Thèse/Thèse/Données/Mx_1x1.txt",
                          "/Users/gojelastat/Desktop/Thèse/Thèse/Données/Exposures_1x1.txt", 
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
df=read.csv2("/Users/gojelastat/Desktop/Thèse/Thèse/Données/heat_month.csv",row.names = 1)



###################################################
#### Fonctions utilisées dans la modélisation #####
###################################################
source(file="/Users/gojelastat/Desktop/Thèse/fonctions_model.R")



##############################################
#### Re-ordonner les données de la France ####
##############################################
FraFeMaleData$Dxt=reshape_french_data(FraFeMaleData)$dxt
FraFeMaleData$Ext=reshape_french_data(FraFeMaleData)$ext
female_rate=FraFeMaleData$Dxt/FraFeMaleData$Ext


FraMaleData$Dxt=reshape_french_data(FraMaleData)$dxt
FraMaleData$Ext=reshape_french_data(FraMaleData)$ext
male_rate=FraMaleData$Dxt/FraMaleData$Ext





### Lee-Carter
## Initialisation du modèle
constLC=function(ax, bx, kt, b0x, gc, wxt, ages)
{
  c1=mean(kt[1,],na.rm=TRUE)
  c2=sum(bx[,1],na.rm=TRUE)
  list(ax=ax+c1*bx,bx=bx/c2,kt=c2*(kt-c1))
}

LC=StMoMo(link="log", staticAgeFun = TRUE, periodAgeFun = "NP", constFun = constLC)

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
plot(fore_lc$rates)
mape_lc=mean(abs(fore_lc$rates[,1:8]-mx_test)/mx_test)
mape_lc


boost_echant=function(data,ages,years,B=1000)
{
  "data représente les données de mortalités"
  ages_tail=tail(ages,1)
  ages_head=ages[1]
  years_head=years[1]
  years_tail=tail(years,1)
  data_ages_head=data$ages[1]
  data_years_head=data$years[1]
  data_ages_tail=tail(data$ages,1)
  data_years_tail=tail(data$years,1)
  
  Dxt=data$Dxt[(ages_head-data_ages_head+1):(ages_tail-data_ages_head+1),
               (years_head-data_years_head+1):(years_tail-data_years_head+1)]
  
  Dxt_b=array(dim = c(B,dim(Dxt)[1],dim(Dxt)[2]))
  for(x in 1:B){
    for(a in 1:dim(Dxt)[1]){
      for(b in 1:dim(Dxt)[2])
      {
        Dxt_b[x,a,b]=rpois(1,Dxt[a,b])
      }
    }
  }
  Dxt_b
}



inf_conf_lc=function(data,ages,years,years_pred,B=50)
{
  dxt=boost_echant(data,0:110,1816:2020,B)
  mu_pred=array(dim=c(B,length(ages),length(years_pred)))
  wei=genWeightMat(ages=ages,years=years)
  
  data$Dxt=dxt[1,,]
  LCfit=fit(LC,data=data,ages.fit = ages, 
            years.fit = years, wxt = wei)
  mu_pred[1,,]=forecast(LCfit,h=length(years_pred))$rates
  
  #plot(ts(mu_pred[1,(a_plot-20+1),],start = years_pred[1]),col='black',
  #    ylim=c(min(mu_pred[1,(a_plot-20+1),]),
  #           max(mu_pred[1,(a_plot-20+1),])))
  
  for(i in 2:B)
  {
    data$Dxt=dxt[i,,]
    LCfit=fit(LC,data=data,ages.fit = ages.fit, 
              years.fit = years.fit, wxt = wei)
    mu_pred[i,,]=fore_lc=forecast(LCfit,h=length(years_pred))$rates
    #lines(ts(mu_pred[i,(a_plot-20+1),],start = years_pred[1]),
    #      col="black",ylim=c(min(mu_pred[i,(a_plot-20+1),]),
    #                        max(mu_pred[i,(a_plot-20+1),])))
  }
  mean_pred=matrix(data=0, nrow = dim(mu_pred)[2],ncol=dim(mu_pred)[3],dimnames = list(ages,years_pred))
  se_pred=matrix(data=0, nrow = dim(mu_pred)[2],ncol=dim(mu_pred)[3],dimnames = list(ages,years_pred))
  IC_min=matrix(data=0, nrow = dim(mu_pred)[2],ncol=dim(mu_pred)[3],dimnames = list(ages,years_pred))
  IC_max=matrix(data=0, nrow = dim(mu_pred)[2],ncol=dim(mu_pred)[3],dimnames = list(ages,years_pred))
  for (x in 1:dim(mu_pred)[2])
  {
    for (y in 1:dim(mu_pred)[3])
    {
      mean_pred[x,y]=mean(mu_pred[,x,y])
      se_pred[x,y]=sd(mu_pred[,x,y])
      IC_min[x,y]=mean_pred[x,y]-1.96*se_pred[x,y]/sqrt(B)
      IC_max[x,y]=mean_pred[x,y]+1.96*se_pred[x,y]/sqrt(B)
    }
  }
  list(mu_pred=mu_pred,
       mean_pred=mean_pred,
       se_pred=se_pred,
       IC_min=IC_min,
       IC_max=IC_max)
  
}




IC_lc=inf_conf_lc(FraMaleData,20:85,1980:2011,2012:2019,B = 10)

ages.fit=20:85
years.fit=1980:2011
LCfit=fit(LC,data=FraMaleData,ages.fit = ages.fit, 
          years.fit = years.fit, wxt = wei)
mu_pred=forecast(LCfit,h=length(2012:2019))$rates

a=65
plot(ts(mu_pred[a,],start = 2012),col=1)
lines(ts(IC_lc$IC_max[a,],start = 2012),col=2)
lines(ts(IC_lc$mean_pred[a,],start = 2012),col=3)
lines(ts(IC_lc$IC_min[a,],start = 2012),col=4)




