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
############################################
### Compute the model selection criteria ###
############################################

# Since APC model is nested in the proposed model the likelihood ratio formula
# can be used to remove the log-likelihood from the proposed model. So, we have:

loglik=function(model,data,ages.fit,years.fit,wxt)
{
  APC=apc(link="log")
  APCfit=fit(APC,data=data,ages.fit = ages.fit, years.fit = years.fit, wxt = wxt)
  dev_apc=APCfit$deviance
  loglik_apc=APCfit$loglik
  dev_mod=model$model$deviance
  
  #### dev_mod-dev_apc=2(loglik_apc-loglik_mod)
  loglik_mod=loglik_apc-(dev_mod-dev_apc)/2
  loglik_mod
}

aic_mod=function(model,data,ages.fit,years.fit,wxt)
{
  loglik_mod=loglik(model,data,ages.fit,years.fit,wxt)
  aic=2*model$npar-2*loglik_mod
  aic
}
bic_model=function(model,data,ages.fit,years.fit,wxt)
{
  loglik_mod=loglik(model,data,ages.fit,years.fit,wxt)
  bic=log(model$nobs)*(model$npar)-2*(loglik_mod)
  bic
}

### MAPE function
mape=function(model)
{
  model_dxt=fitted(model,type = "deaths")
  ext=model$Ext
  dxt=model$Dxt
  mx_fit=model_dxt/ext
  mx=dxt/ext
  mape=mean(abs((mx_fit-mx)/mx))
  mape
}


#################################################
### Function that computes the cohort matrix ####
#################################################
### where df is zero matrix of size n*n.c

indic=function(df,n_a,n_y)
{
  a=c(rep(1:n_a,n_y))
  y=c(rep(1:n_y,each=n_a))
  c=cbind(a,y) ### c représente les indexes de la matrice obtenue après flatten de D_xt.
  
  x=sort(unique(y-a))
  m=dim(df)
  for (i in 1:m[1])
    for(j in 1:m[2])
      if((c[i,2]-c[i,1])==x[j])
        df[i,j]=1
  df
}

I=function(x,a_c)
{
  y=c()
  l=length(x)
  for(t in 1:l){
    if(x[t]-a_c>0)
      y=c(y,1)
    else
      y=c(y,0)
  }
  y
}

non_null=function(I_t)
{
  y=c()
  l=length(I_t)
  for (i in 1:l) 
  {
    if(I_t[i]>0)
      y=c(y,I_t[i])
  }
  y
}

x=5:89
I(x,65)
### Modèle PLAT+TERME de température mais avec le même coeficient par tranche de 5ans
########################
#### Model function ####
########################
model_function=function(dxt,ext,wxt,Ic,m,c,x1,x2,a_c,ages,years)
{
  n.a=dim(dxt)[1]# age length
  n.t=dim(dxt)[2] # time length
  n=n.a*n.t   # lines number of model matrix 
  n.c=n.a+n.t-1 # cohort length
  
  cohort=(years[1]-ages[n.a]):(years[n.t]-ages[1])
  
  
  x=ages
  x.bar=mean(x) # ages mean
  e.log=log(ext) # log of exposure
  
  ###Constructiion of kt_5 matrix
  #I_t=Ic
  I_t=pmax((Ic-m),0)
  #I_t1=pmax((Ic-a),0)
  mi=diag(I_t)
  t=dim(mi)[2]
  v=c()
  u=sapply(1:t,function(x){sum(mi[,x])})
  for(i in 1:length(u))
  {
    if (u[i]!=0)
    {
      v=cbind(v,mi[,i])
    }
  }
  I_t=v
  
  n.t5=dim(I_t)[2]
  
  ### Vectorization of d_xt et e_xt
  v_dxt=as.vector(dxt)
  v_ext=as.vector(e.log)
  v_wxt=as.vector(wxt)
  
  ### Matrix X of model
  
  #### 1) age matrix
  X_a=kronecker(as.matrix(rep(1,n.t)),diag(n.a))
  
  #### 2) period matrix kt1
  X_t1=kronecker(diag(n.t),as.matrix(rep(1,n.a)))
  
  #### 3) period matrix kt2
  X_t2=kronecker(diag(n.t),as.matrix(x.bar-x))
  
  #### 4) period matrix kt3
  X_t3=kronecker(diag(n.t),as.matrix(pmax(x.bar-x,0)))
  
  #### 4) period matrix kt4
  X_t4=kronecker(diag(n.t),as.matrix(pmax(x1-x,0)+c*pmax(x-x2,0))^2)
  #X_t4=kronecker(diag(n.t),as.matrix((x.bar-x)^2))
  
  #### 5) period matrix kt5
  y=(x-a_c)
  
  #X_t5=kronecker(I_t,as.matrix(c_x*I(x,a_c)))
  X_t5=kronecker(I_t,as.matrix(pmax(x-a_c,0)))
  #X_t5_=kronecker(diag(n.t),matrix(pmax(x.bar-x,0)))
  #X_t5__=kronecker(diag(I_t1),matrix(c_x*pmax(a_c-x,0)))
  #X_t5=X_t5_+X_t5__
  ### cohort matrix
  X_c=indic(matrix(0,n,n.c),n.a,n.t)
  
  ### consider cohort that under 1950 years
  n.s=1950-cohort[1]+1 
  X_c=X_c[,1:n.s]
  #n.s=n.c
  #### finally we obtain
  X=cbind(X_a,X_t1,X_t2,X_t3,X_t4,X_t5,X_c)
  #X=cbind(X_a,X_t1,X_t2,X_t3,X_t5,X_c)
  
  #### Constraint matrix
  
  ### Constraint 1
  h1=cbind(t(rep(0,n.a)),t(rep(1,n.t)),t(rep(0,3*n.t+n.t5+n.s)))
  #h1=cbind(t(rep(0,n.a)),t(rep(1,n.t)),t(rep(0,3*n.t+n.s)))
  #h1=cbind(t(rep(0,n.a)),t(rep(1,n.t)),t(rep(0,2*n.t+n.t5+n.s)))
  
  ### Constraint 2
  h2=cbind(t(rep(0,n.a+n.t)),t(rep(1,n.t)),t(rep(0,2*n.t+n.t5+n.s)))
  #h2=cbind(t(rep(0,n.a+n.t)),t(rep(1,n.t)),t(rep(0,2*n.t+n.s)))
  #h2=cbind(t(rep(0,n.a+n.t)),t(rep(1,n.t)),t(rep(0,n.t+n.t5+n.s)))
  
  ### Constraint 3
  h3=cbind(t(rep(0,n.a+2*n.t)),t(rep(1,n.t)),t(rep(0,n.t+n.t5+n.s)))
  #h3=cbind(t(rep(0,n.a+2*n.t)),t(rep(1,n.t)),t(rep(0,n.t5+n.s)))
  #h3=cbind(t(rep(0,n.a+2*n.t)),t(rep(1,n.t)),t(rep(0,n.t+n.s)))
  
  ### Constraint 4                                                                                           ConstraintConstraint 4
  h4=cbind(t(rep(0,n.a+3*n.t)),t(rep(1,n.t)),t(rep(0,n.t5+n.s)))
  #h4=cbind(t(rep(0,n.a+3*n.t)),t(rep(1,n.t)),t(rep(0,n.s)))
  
  ### Constraint 5
  h5=cbind(t(rep(0,n.a+4*n.t+n.t5)),t(rep(1,n.s)))
  #h5=cbind(t(rep(0,n.a+4*n.t)),t(rep(1,n.s)))
  #h5=cbind(t(rep(0,n.a+3*n.t+n.t5)),t(rep(1,n.s)))
  
  ### Constraint 6
  h6=cbind(t(rep(0,n.a+4*n.t+n.t5)),t(1:n.s))
  #h6=cbind(t(rep(0,n.a+4*n.t)),t(1:n.s))
  #h6=cbind(t(rep(0,n.a+3*n.t+n.t5)),t(1:n.s))
  
  ### Constraint 7
  #h7=cbind(t(rep(0,n.a+3*n.t+n.t5)),t((1:n.s)^2))
  
  H=rbind(h1,h2,h3,h4,h5,h6)
  
  
  #### Model 
  model=gnm(v_dxt~-1+offset(v_ext)+X,weights = v_wxt,family = poisson(link = 'log'))
  
  
  
  #### coefficients matrix
  teta=solve(t(X)%*%X+t(H)%*%H,tol=4.29445e-23)%*%t(X)%*%(log(model$fitted.values)-v_ext)
  
  ax=teta[1:n.a]
  kt1=teta[(n.a+1):(n.a+n.t)]
  kt2=teta[(n.a+n.t+1):(n.a+2*n.t)]
  kt3=teta[(n.a+2*n.t+1):(n.a+3*n.t)]
  kt4=teta[(n.a+3*n.t+1):(n.a+4*n.t)]
  kt5=teta[(n.a+4*n.t+1):(n.a+4*n.t+n.t5)]
  gc=teta[(n.a+4*n.t+n.t5+1):dim(teta)[1]]
  
  #kt4=teta[(n.a+3*n.t+1):(n.a+4*n.t)]
  #kt5=teta[(n.a+3*n.t+1):(n.a+3*n.t+n.t5)]
  #gc=teta[(n.a+3*n.t+n.t5+1):dim(teta)[1]]
  
  
  structure(list(
    ax=ts(ax,start = ages.fit[1],frequency = 1),
    kt1=ts(kt1,start=years[1],frequency=1),
    kt2=ts(kt2,start=years[1],frequency=1),
    kt3=ts(kt3,start=years[1],frequency=1),
    kt4=ts(kt4,start=years[1],frequency=1),
    kt5=kt5,
    gc=ts(gc,start = cohort[1],frequency = 1),
    #gc=ts(gc,start = cohort[1],frequency = 1),
    Dxt=dxt,
    Ext=ext,
    wxt=wxt,
    npar=dim(teta)[1],
    nobs=n.a*n.t,
    rank=qr(X)$rank,
    col=dim(X)[2],
    sing=dim(X)[2]-qr(X)$rank,
    #X_t5=X_t5,
    ages=ages.fit,
    years=years.fit,
    mean_T=m,
    x1=x1,
    x2=x2,
    c=c,
    a_c=a_c,
    model=model
  ))
  
}


##############################################
#### Application to male data from France ####
##############################################

### compute model

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



plot(model_f$kt5, type = "l")
plot(model_f$kt1)
plot(model_f$kt2)
plot(model_f$kt3)
plot(model_f$kt4)
plot(model_f$kt5)
plot(model_f$ax)
plot(model_f$gc)



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
