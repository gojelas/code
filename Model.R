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
data <- read.demogdata("/Users/gojelastat/Desktop/Thèse/Thèse/Données/Mx_1x1.txt",
                       "/Users/gojelastat/Desktop/Thèse/Thèse/Données/Exposures_1x1.txt", 
                       type="mortality", label="France")

#FraMaleData=StMoMoData(data, series = names(data$rate)[2], type ="central")

FraFeMaleData=StMoMoData(data, series = names(data$rate)[1], type ="central")



#### 1) heat wave indicator
#ind=read.csv("/Volumes/GOJELA'S/Thèse/Correction du modèle/données/indicateur2.csv",row.names = 1)
#ind=read.csv("/media/gouthonabiodunjean-luc/GOJELA'S/Thèse/Correction du modèle/données/indicateur2.csv",row.names = 1)
ind=read.csv("F:/Thèse/Correction du modèle/données/indicateur2.csv",row.names = 1)


#### 2) different temperature scenarios
#temp_=read.csv("/Volumes/GOJELA'S/Thèse/Correction du modèle/données/Données_globale_fr.csv",row.names = 1)
#temp_=read.csv("/media/gouthonabiodunjean-luc/GOJELA'S/Thèse/Correction du modèle/données/Données_globale_fr.csv",row.names = 1)
#temp_=read.csv("F:/Thèse/Correction du modèle/données/Données_globale_fr.csv",row.names = 1)
#temp=temp_[1:70,]
#temp


#### 3) données des mois à haute température
#write.csv2(df,"F:/Thèse/Correction du modèle/données/heat_month.csv",row.names = FALSE)

df=read.csv2("/Users/gojelastat/Desktop/Thèse/Thèse/Données/heat_month.csv",row.names = 1)


ages.fit=20:85
years.fit=1980:2011
fe_m=FraFeMaleData$Dxt[21:86,(1980-1816+1):(2011-1816+1)]
fe_e=FraFeMaleData$Ext[21:86,(1980-1816+1):(2011-1816+1)]
wei=genWeightMat(ages=ages.fit,years=years.fit)
#ind_=ind[(1980-1950+1):(2011-1950+1),1]
#Ic=ind_
Ic=df[(1980-1961+1):(2011-1961+1),6]
m=mean(Ic)
a_c=65
c=0.03
I_t=pmax((Ic-mean(Ic)),0)
x1=40
x2=60
#ind_=as.integer(ind_[,1])
model_f=model_function(fe_m,fe_e,wei,Ic,m,c,x1,x2,a_c,ages.fit,years.fit)


### estimated values function
fit_model=function(model,ages.fit,years.fit)
{
  n.a=length(ages.fit)
  n.t=length(years.fit)
  dxt_fit=matrix(model$model$fitted.values,nrow=n.a,
                 ncol=n.t,dimnames = list(ages.fit,years.fit))
  uxt_fit=dxt_fit/model$Ext
  
  structure(
    list(
      Dxt_fit=dxt_fit,
      uxt_fit=uxt_fit
    )
  )
}


### criterion computation 
ages.fit=20:85
years.fit=1980:2011
wei=genWeightMat(ages=ages.fit,years=years.fit)
loglik_m2=loglik(model_f,FraFeMaleData,ages.fit,years.fit,wei)
aic_m2=aic_mod(model_f,FraFeMaleData,ages.fit,years.fit,wei)
bic_m2=bic_model(model_f,FraFeMaleData,ages.fit,years.fit,wei)

print(c(loglik_m2,model2$model$deviance,aic_m2,bic_m2))
# -10361.117   3024.748  21256.233  22766.223

### Erreur MAPE
years.fit=1980:2011
fit_=fit_model(model_f,ages.fit,years.fit)
mf_uxt=fit_$uxt_fit
mx=model_f$Dxt/model_f$Ext
mf_mape=mean(abs((mx-mf_uxt)/mx))
mf_mape
### mape=0.03047982

#############################################
### Time series forecasts for men's data ####
#############################################

### kt2
kt2=ts(model_f$kt2,start=years.fit[1],frequency = 1)
fit2=auto.arima(kt2)
fore_2=forecast(fit2,h=50)
kt2_pred=fore_2$mean
plot(fore_2)

### kt1
kt1=ts(model_f$kt1,start=years.fit[1],frequency = 1)
fit1=auto.arima(kt1)
fore_1=forecast(fit1, h=50)
kt1_pred=fore_1$mean
plot(fore_1)

### kt3
kt3=ts(model_f$kt3,start=years.fit[1],frequency = 1)
fit3=auto.arima(kt3)
fore_3=forecast(fit3,h=50)
kt3_pred=fore_3$mean
plot(fore_3)

### kt4
kt4=ts(model_f$kt4,start=years.fit[1],frequency = 1)
fit4=auto.arima(kt4)
fore_4=forecast(fit4,h=50)
kt4_pred=fore_4$mean
plot(fore_4)


kt5=ts(model_f$kt5,start=1,frequency = 1)
fit5=auto.arima(kt5)
fore_5=forecast(fit5,h=25)
kt5_pred=fore_5$mean
plot(fore_5)


x = seq(from = 1 , to = length(kt5) , by = 1)
y = kt5
par(mfrow = c(1 , 1))
plot(x , y , type = "l")
points(x , y , col = "red")
l = Reg(x , y)
points(l$x , l$y , col = "blue")

kt5_est=ts(l$y,start=1,frequency = 1)
fit5=auto.arima(kt5_est)
fore_5=forecast(fit5,h=30)
t=seq(3,30,by=3)
kt5_pred=fore_5$mean[t]
plot(fore_5)


###gamma
n.a=length(ages.fit)# age length
n.t=length(years.fit) # time length
cohort=(years.fit[1]-ages.fit[n.a]):(years.fit[n.t]-ages.fit[1])
g=ts(model_f$gc,start=cohort[1],frequency=1)
plot(g)
fit_gc=auto.arima(g)
pc=forecast(fit_gc,h=150)
gc_pred=pc$mean
plot(pc)

gg=cbind(g,gc_pred)
gg

### gc
n.a=length(ages.fit)# age length
n.t=length(years.fit) # time length


indd=df[(2012-1961+1):(2019-1961+1),6]
h=length(indd)
ch_beg=years.fit[n.t]+1-ages.fit[n.a] ### start of cohort forecasts
ch_end=years.fit[n.t]+h-ages.fit[1] ### end of cohort forecasts
gcc=c(gg[(ch_beg-cohort[1]+1):(1950-cohort[1]+1),1],gg[(1950-cohort[1]+2):(ch_end-cohort[1]+1),2])

c1=kt1_pred[1:h]
c2=kt2_pred[1:h]
c3=kt3_pred[1:h]
c4=kt4_pred[1:h]

n.a=tail(ages.fit,1)-ages.fit[1]+1
n.t=length(c1)
n=n.a*n.t
n.c=n.a+n.t-1

x=ages.fit
x.bar=mean(x)
a=mean(Ic)
Ic=indd  
#Ic=ind[(2012-1950+1):(2019-1950+1),1]
I_t=Ic
###Construction de la matrice de kt_5
I_t=pmax((Ic-a),0)
m=diag(I_t)
t=dim(m)[2]
v=c()
u=sapply(1:t,function(x){sum(m[,x])})
for(i in 1:length(u))
{
  if (u[i]!=0)
  {
    v=cbind(v,m[,i])
  }
}
I_t=v

n.t5=dim(I_t)[2]

c5=kt5_pred[1:n.t5]
axx=model_f$ax
teta=c(axx,c1,c2,c3,c4,c5,gcc)




a=65

#### 1) Matrice d'âge
X_a=kronecker(as.matrix(rep(1,n.t)),diag(n.a))

#### 2) Matrice de période kt1
X_t1=kronecker(diag(n.t),as.matrix(rep(1,n.a)))

#### 3) Matrice de période kt2
X_t2=kronecker(diag(n.t),as.matrix(x.bar-x))

#### 4) Matrice de période kt3
X_t3=kronecker(diag(n.t),as.matrix(pmax(x.bar-x,0)))

#X_t5_=kronecker(diag(n.t),matrix(pmax(x.bar-x,0)^2))
#X_t5__=kronecker(diag(I_t),matrix(c_x*pmax(a_c-x,0)^2))
#X_t5=(X_t5_+X_t5__)


#### 4) Matrice de période kt4
X_t4=kronecker(diag(n.t),as.matrix(0.4*pmax(40-x,0)+c_x*pmax(x-50,0))^2)
#X_t4=kronecker(diag(n.t),as.matrix(x-x.bar)^2)

#### 5) Matrice de période kt5
y=(x-a_c)
#y=I(x,a_c)
X_t5=kronecker(I_t,as.matrix(pmax(y,0)))

### Matrice de cohorte
X_c=indic(matrix(0,n,n.c),n.a,n.t)

#### Enfin on obtient
X=cbind(X_a,X_t1,X_t2,X_t3,X_t4,X_t5,X_c)
mu=X%*%teta

log_mu=matrix(mu,nrow = n.a,dimnames =list(ages.fit,2012:2019))
mod_mu=exp(log_mu)

age=55
female_rate=data$rate$female
mx_test=female_rate[21:86,197:204]
plot(ts(mod_mu[age,],start = 2012))
lines(ts(mx_test[age,],start = 2012),col="red")
mape_mod=mean(abs(mod_mu-mx_test)/mx_test)
mape_mod
### mape sur 2012-2019=0.09798836




model_predict=function(model,years.pred)
{
  l=length(years.pred)# l represent the length of prediction
  h=l*5 # to have enough points for kt5 forecast
  ages.fit=model$ages
  years.fit=model$years
  
  
  ### kt1
  kt1=ts(model$kt1,start=years.fit[1],frequency = 1)
  fit1=auto.arima(kt1)
  fore_1=forecast(fit1, h)
  kt1_pred=fore_1$mean
  #plot(fore_1)
  
  ### kt2
  kt2=ts(model$kt2,start=years.fit[1],frequency = 1)
  fit2=auto.arima(kt2)
  fore_2=forecast(fit2,h)
  kt2_pred=fore_2$mean
  #plot(fore_2)
  
  ### kt3
  kt3=ts(model$kt3,start=years.fit[1],frequency = 1)
  fit3=auto.arima(kt3)
  fore_3=forecast(fit3,h)
  kt3_pred=fore_3$mean
  #plot(fore_3)
  
  ### kt4
  kt4=ts(model$kt4,start=years.fit[1],frequency = 1)
  fit4=auto.arima(kt4)
  fore_4=forecast(fit4,h)
  kt4_pred=fore_4$mean
  #plot(fore_4)
  
  
  kt5=ts(model$kt5,start=1,frequency = 1)
  fit5=auto.arima(kt5)
  fore_5=forecast(fit5,h)
  kt5_pred=fore_5$mean
  #plot(fore_5)
  
  
  x = seq(from = 1 , to = length(kt5) , by = 1)
  y = kt5
  par(mfrow = c(1 , 1))
  #plot(x , y , type = "l")
  #points(x , y , col = "red")
  l = Reg(x , y)
  #points(l$x , l$y , col = "blue")
  
  kt5_est=ts(l$y,start=1,frequency = 1)
  fit5=auto.arima(kt5_est)
  fore_5=forecast(fit5,h)
  t=seq(3,30,by=3)
  kt5_pred=fore_5$mean[t]
  #plot(fore_5)
  
  
  ###gamma
  n.a=length(ages.fit)# age length
  n.t=length(years.fit) # time length
  cohort=(years.fit[1]-ages.fit[n.a]):(years.fit[n.t]-ages.fit[1])
  g=ts(model$gc,start=cohort[1],frequency=1)
  #plot(g)
  fit_gc=auto.arima(g)
  pc=forecast(fit_gc,h=150)
  gc_pred=pc$mean
  #plot(pc)
  
  gg=cbind(g,gc_pred)
  gg
  
  ### gc
  n.a=length(ages.fit)# age length
  n.t=length(years.fit) # time length
  
  
  indd=df[(years.pred[1]-1961+1):(tail(years.pred,1)-1961+1),6]
  h=length(indd)
  ch_beg=years.fit[n.t]+1-ages.fit[n.a] ### start of cohort forecasts
  ch_end=years.fit[n.t]+h-ages.fit[1] ### end of cohort forecasts
  gcc=c(gg[(ch_beg-cohort[1]+1):(1950-cohort[1]+1),1],gg[(1950-cohort[1]+2):(ch_end-cohort[1]+1),2])
  
  c1=kt1_pred[1:h]
  c2=kt2_pred[1:h]
  c3=kt3_pred[1:h]
  c4=kt4_pred[1:h]
  
  n.a=tail(ages.fit,1)-ages.fit[1]+1
  n.t=length(c1)
  n=n.a*n.t
  n.c=n.a+n.t-1
  
  x=ages.fit
  x.bar=mean(x)
  Ic=indd  
  #Ic=ind[(2012-1950+1):(2019-1950+1),1]
  ###Construction de la matrice de kt_5
  I_t=pmax((Ic-model$mean_T),0)
  m=diag(I_t)
  t=dim(m)[2]
  v=c()
  u=sapply(1:t,function(x){sum(m[,x])})
  for(i in 1:length(u))
  {
    if (u[i]!=0)
    {
      v=cbind(v,m[,i])
    }
  }
  I_t=v
  
  n.t5=dim(I_t)[2]
  
  c5=kt5_pred[1:n.t5]
  axx=model$ax
  teta=c(axx,c1,c2,c3,c4,c5,gcc)
  
  
  
  #### 1) Matrice d'âge
  X_a=kronecker(as.matrix(rep(1,n.t)),diag(n.a))
  
  #### 2) Matrice de période kt1
  X_t1=kronecker(diag(n.t),as.matrix(rep(1,n.a)))
  
  #### 3) Matrice de période kt2
  X_t2=kronecker(diag(n.t),as.matrix(x.bar-x))
  
  #### 4) Matrice de période kt3
  X_t3=kronecker(diag(n.t),as.matrix(pmax(x.bar-x,0)))
  
  #X_t5_=kronecker(diag(n.t),matrix(pmax(x.bar-x,0)^2))
  #X_t5__=kronecker(diag(I_t),matrix(c_x*pmax(a_c-x,0)^2))
  #X_t5=(X_t5_+X_t5__)
  
  
  #### 4) Matrice de période kt4
  X_t4=kronecker(diag(n.t),as.matrix(pmax((model$x1)-x,0)+(model$c)*pmax(x-(model$x2),0))^2)
  #X_t4=kronecker(diag(n.t),as.matrix(x-x.bar)^2)
  
  #### 5) Matrice de période kt5
  y=(x-(model$a_c))
  #y=I(x,a_c)
  X_t5=kronecker(I_t,as.matrix(pmax(y,0)))
  
  ### Matrice de cohorte
  X_c=indic(matrix(0,n,n.c),n.a,n.t)
  
  #### Enfin on obtient
  X=cbind(X_a,X_t1,X_t2,X_t3,X_t4,X_t5,X_c)
  mu=X%*%teta
  
  log_mu=matrix(mu,nrow = n.a,dimnames =list(ages.fit,years.pred))
  mod_mu=exp(log_mu)
  
  
  structure(list(
    mxt_pred=mod_mu,
    years.pred=years.pred,
    ages=ages.fit
  ))
}

fem_pred=model_predict(model_f,2012:2019)



best_param=function()
  c=seq(0,2,by=0.1)
