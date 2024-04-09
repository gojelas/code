library(StMoMo)
library(demography)
library(urca)
library(tseries)

##########################
### mortality database ###
##########################
#data <- read.demogdata("/media/gouthonabiodunjean-luc/GOJELA'S/Thèse/Correction du modèle/données/Mx_1x1.txt",
#                       "/media/gouthonabiodunjean-luc/GOJELA'S/Thèse/Correction du modèle/données/Exposures_1x1.txt", 
#                       type="mortality", label="France")

data <- read.demogdata("/Users/gojelastat/Downloads/Thèse/Données/Mx_1x1.txt",
                       "/Users/gojelastat/Downloads/Thèse/Données/Exposures_1x1.txt", 
                       type="mortality", label="France")

FraMaleData=StMoMoData(data, series = names(data$rate)[2], type ="central")


#data <- read.demogdata("/media/gouthonabiodunjean-luc/GOJELA'S/Thèse/Correction du modèle/données/Mx_1x1 (UK).txt",
#                       "/media/gouthonabiodunjean-luc/GOJELA'S/Thèse/Correction du modèle/données/Exposures_1x1 (UK).txt", 
#                       type="mortality", label="France")
#UKMaleData=StMoMoData(data, series = names(data$rate)[2], type ="central")


FraFemaleData=StMoMoData(data, series = names(data$rate)[1], type ="central")

#### 1) heat wave indicator
#ind=read.csv("/Volumes/GOJELA'S/Thèse/Correction du modèle/données/indicateur2.csv",row.names = 1)
#ind=read.csv("/media/gouthonabiodunjean-luc/GOJELA'S/Thèse/Correction du modèle/données/indicateur2.csv",row.names = 1)
ind=read.csv("/Users/gojelastat/Downloads/Thèse/Données/indicateur2.csv",row.names = 1)

source(file="/Users/gojelastat/Desktop/Thèse/fonctions_model.R")


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
mape=function(model,ages.fit,years.fit)
{
  model_dxt=fit_model(model,ages.fit,years.fit)$Dxt_fit
  ext=model$Ext
  dxt=model$Dxt
  mx_fit=model_dxt/ext
  mx=dxt/ext
  mean(abs((mx_fit-mx)/mx))
  
}

mad=function(model,ages.fit,years.fit)
{
  model_dxt=fit_model(model,ages.fit,years.fit)$Dxt_fit
  ext=model$Ext
  dxt=model$Dxt
  mx_fit=model_dxt/ext
  mx=dxt/ext
  mean(abs(mx_fit-mx))
  
}

mse=function(model,ages.fit,years.fit)
{
  model_dxt=fit_model(model,ages.fit,years.fit)$Dxt_fit
  ext=model$Ext
  dxt=model$Dxt
  mx_fit=model_dxt/ext
  mx=dxt/ext
  mean((mx_fit-mx)^2)
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

### Modèle SPO
########################
#### Model function ####
########################
plat_function=function(dxt,ext,wxt,ages,years)
{
  n.a=dim(dxt)[1]# age length
  n.t=dim(dxt)[2] # time length
  n=n.a*n.t   # lines number of model matrix 
  n.c=n.a+n.t-1 # cohort length
  
  cohort=(years[1]-ages[n.a]):(years[n.t]-ages[1])
  
  
  x=ages
  x.bar=mean(x) # ages mean
  e.log=log(ext) # log of exposure
  
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
  
  ### cohort 
  X_c=indic(matrix(0,n,n.c),n.a,n.t)
  
  ### consider cohort that under 1950 years
  n.s=1950-cohort[1]+1 
  X_c=X_c[,1:n.s]
  #n.s=n.c
  X=cbind(X_a,X_t1,X_t2,X_t3,X_c)
  
  #### Constraint matrix
  
  ### Constraint 1
  h1=cbind(t(rep(0,n.a)),t(rep(1,n.t)),t(rep(0,2*n.t+n.s)))
  
  ### Constraint 2
  h2=cbind(t(rep(0,n.a+n.t)),t(rep(1,n.t)),t(rep(0,n.t+n.s)))
  
  ### Constraint 3
  h3=cbind(t(rep(0,n.a+2*n.t)),t(rep(1,n.t)),t(rep(0,n.s)))
  
  ### Constraint 5
  h5=cbind(t(rep(0,n.a+3*n.t)),t(rep(1,n.s)))
  
  ### Constraint 6
  h6=cbind(t(rep(0,n.a+3*n.t)),t(1:n.s))
  
  ### Constraint 7
  #h7=cbind(t(rep(0,n.a+3*n.t+n.t5)),t((1:n.s)^2))
  
  H=rbind(h1,h2,h3,h5,h6)
  
  
  #### Model 
  model=gnm(v_dxt~-1+offset(v_ext)+X,weights = v_wxt,family = poisson(link = 'log'))
  
  
  
  #### coefficients matrix
  teta=solve(t(X)%*%X+t(H)%*%H,tol=4.29445e-23)%*%t(X)%*%(log(model$fitted.values)-v_ext)
  
  ax=teta[1:n.a]
  kt1=teta[(n.a+1):(n.a+n.t)]
  kt2=teta[(n.a+n.t+1):(n.a+2*n.t)]
  kt3=teta[(n.a+2*n.t+1):(n.a+3*n.t)]
  gc=teta[(n.a+3*n.t+1):dim(teta)[1]]
  
  
  structure(list(
    ax=ax,
    kt1=ts(kt1,start=years[1],frequency=1),
    kt2=ts(kt2,start=years[1],frequency=1),
    kt3=ts(kt3,start=years[1],frequency=1),
    gc=ts(gc,start = cohort[1],frequency = 1),
    Dxt=dxt,
    Ext=ext,
    wxt=wxt,
    npar=dim(teta)[1],
    nobs=n.a*n.t,
    rank=qr(X)$rank,
    col=dim(X)[2],
    sing=dim(X)[2]-qr(X)$rank,
    cohort=cohort,
    model=model
  ))
  
}

##############################################
#### Application to male data from France ####
##############################################

### compute model


ages.fit=20:85
years.fit=1980:2011
fr_m=FraMaleData$Dxt[21:86,(1980-1816+1):(2011-1816+1)]
fr_e=FraMaleData$Ext[21:86,(1980-1816+1):(2011-1816+1)]
wei=genWeightMat(ages=ages.fit,years=years.fit)
plat=plat_function(fr_m,fr_e,wei,ages.fit,years.fit)




plat_criteria=get_criterion(plat,FraMaleData,ages.fit,years.fit)
plat_criteria




### fitting PLAT
get_estimation_plat=function(data, ages, years)
{
  age_head=data$ages[1]
  age_tail=tail(data$ages,1)
  year_head=data$years[1]
  year_tail=tail(data$years,1)
  
  fr_m=data$Dxt[(ages[1]-age_head+1):(tail(ages,1)-age_head+1),
                (years[1]-year_head+1):(tail(years,1)-year_head+1)]
  fr_e=data$Ext[(ages[1]-age_head+1):(tail(ages,1)-age_head+1),
                (years[1]-year_head+1):(tail(years,1)-year_head+1)]
  wei=genWeightMat(ages=ages,years=years)
  #Ic=df[(years[1]-1961+1):(tail(years,1)-1961+1),6]
  #a_m=mean(Ic)
  
  plat=plat_function(fr_m,fr_e,wei,ages = ages,years = years)
  
  list(
    years=years,
    ages=ages,
    cohort=plat$cohort,
    kt1=plat$kt1,
    kt2=plat$kt2,
    kt3=plat$kt3,
    #kt4=m2$kt4,
    #kt5=m2$kt5,
    ax=plat$ax,
    gc=plat$gc
    #a_m=a_m,
    #x1=x1,
    #x2=x2,
    #c=c
  )
}



#forecasting fonction for plat model
get_predict_plat=function(estim,years_pred)
{
  h=length(years_pred)
  
  fit1=auto.arima(estim$kt1)
  fore_1=forecast(fit1, h)
  kt1_pred=fore_1$mean
  
  fit2=auto.arima(estim$kt2)
  fore_2=forecast(fit2, h)
  kt2_pred=fore_2$mean
  
  fit3=auto.arima(estim$kt3)
  fore_3=forecast(fit3, h)
  kt3_pred=fore_3$mean
  
  
  fit_gc=auto.arima(estim$g)
  pc=forecast(fit_gc,h=150)
  gc_pred=pc$mean
  
  gg=cbind(estim$g,gc_pred)
  
  ### gc
  n.a=length(estim$ages)# age length
  n.t=length(estim$years) # time length
  cohort=(estim$years[1]-estim$ages[n.a]):(estim$years[n.t]-estim$ages[1])
  
  
  ch_beg=estim$years[n.t]+1-estim$ages[n.a] ### start of cohort forecasts
  ch_end=estim$years[n.t]+h-estim$ages[1] ### end of cohort forecasts
  gcc=c(gg[(ch_beg-estim$cohort[1]+1):(1950-estim$cohort[1]+1),1],
        gg[(1950-estim$cohort[1]+2):(ch_end-estim$cohort[1]+1),2])
  
  c1=kt1_pred
  c2=kt2_pred
  c3=kt3_pred
  
  n.a=length(estim$ages)
  n.t=h
  n=n.a*n.t
  n.c=n.a+n.t-1
  
  axx=estim$ax
  teta=c(axx,c1,c2,c3,gcc)
  x=estim$ages
  x.bar=mean(x)
  
  #### 1) Matrice d'âge
  X_a=kronecker(as.matrix(rep(1,n.t)),diag(n.a))
  #### 2) Matrice de période kt1
  X_t1=kronecker(diag(n.t),as.matrix(rep(1,n.a)))
  #### 3) Matrice de période kt2
  X_t2=kronecker(diag(n.t),as.matrix(x.bar-x))
  #### 4) Matrice de période kt3
  X_t3=kronecker(diag(n.t),as.matrix(pmax(x.bar-x,0)))
  ### Matrice de cohorte
  X_c=indic(matrix(0,n,n.c),n.a,n.t)
  
  #### Enfin on obtient
  X=cbind(X_a,X_t1,X_t2,X_t3,X_c)
  mu=X%*%teta
  
  log_mu=matrix(mu,nrow = n.a,dimnames =list(estim$ages,years_pred))
  mod_mu=exp(log_mu)
  
  list(
    mxt_pred=mod_mu,
    years_pred=years_pred
  )
  
}


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



inf_conf_plat=function(data,ages,years,years_pred,B=50)
{
  dxt=boost_echant(data,0:110,1816:2020,B)
  mu_pred=array(dim=c(B,length(ages),length(years_pred)))
  
  data$Dxt=dxt[1,,]
  plat=get_estimation_plat(data,ages,years)
  mu_pred[1,,]=get_predict_plat(plat,years_pred)$mxt_pred
  
  #plot(ts(mu_pred[1,(a_plot-20+1),],start = years_pred[1]),col='black',
  #    ylim=c(min(mu_pred[1,(a_plot-20+1),]),
  #           max(mu_pred[1,(a_plot-20+1),])))
  
  for(i in 2:B)
  {
    data$Dxt=dxt[i,,]
    plat=get_estimation_plat(data,ages,years)
    mu_pred[i,,]=get_predict_plat(plat,years_pred)$mxt_pred
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




IC_plat=inf_conf_plat(FraMaleData,20:85,1980:2011,2012:2019,B = 1000)


plat=get_estimation_plat(FraMaleData,20:85,1980:2011)
mu_pred=get_predict_plat(plat,2012:2019)$mxt_pred

a=60
plot(ts(mu_pred[a,],start = 2012),col=1)
lines(ts(IC_plat$IC_max[a,],start = 2012),col=2)
lines(ts(IC_plat$mean_pred[a,],start = 2012),col=3)
lines(ts(IC_plat$IC_min[a,],start = 2012),col=4)







#male_rate=data$rate$male
mx_test=male_rate[21:86,197:204]
#plot(ts(plat_mu[20,],start = 2012))
#lines(ts(mx_test[20,],start = 2012),col="red")
mape_plat=mean(abs(plat_mu-mx_test)/mx_test)
mape_plat
### mape sur 2012-2019, 0.2248657

