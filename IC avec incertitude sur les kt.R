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

data_fr <- read.demogdata("Mx_1x1.txt",
                          "Exposures_1x1.txt", 
                          type="mortality", label="France")

FraMaleData=StMoMoData(data_fr, series = names(data_fr$rate)[2], type ="central")

FraFeMaleData=StMoMoData(data_fr, series = names(data_fr$rate)[1], type ="central")

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
  n.s=1946-cohort[1]+1 
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
  #h7=cbind(t(rep(0,n.a+3*n.t)),t((1:n.s)^2))
  
  H=rbind(h1,h2,h3,h5,h6)#,h7)
  
  
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


ages.fit=0:100
years.fit=1950:2017
fr_m=FraMaleData$Dxt[1:101,(1950-1816+1):(2017-1816+1)]
fr_e=FraMaleData$Ext[1:101,(1950-1816+1):(2017-1816+1)]
wei=genWeightMat(ages=ages.fit,years=years.fit)
plat=plat_function(fr_m,fr_e,wei,ages.fit,years.fit)

#fit_model(plat,20:85,1980:2011)$
plot(plat$gc)
  
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
    ax=plat$ax,
    gc=plat$gc
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
  #gcc=c(gg[(ch_beg-estim$cohort[1]+1):(tail(estim$cohort,1)-estim$cohort[1]+1),1],
  #      gg[(tail(estim$cohort,1)-estim$cohort[1]+2):(ch_end-estim$cohort[1]+1),2])
  gcc=c(gg[(ch_beg-estim$cohort[1]+1):(1946-estim$cohort[1]+1),1],
        gg[(1946-estim$cohort[1]+2):(ch_end-estim$cohort[1]+1),2])
  
  
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
  #for(x in 1:B){
  for(a in 1:dim(Dxt)[1]){
    for(b in 1:dim(Dxt)[2])
    {
      Dxt_b[,a,b]=rpois(B,Dxt[a,b])
    }
  }
  Dxt_b
}



inf_conf_plat=function(data,ages,years,years_pred,B=50)
{
  a_min=data$ages[1]
  a_max=tail(data$ages,1)
  y_min=data$years[1]
  y_max=tail(data$years,1)
  dxt=boost_echant(data,a_min:a_max,y_min:y_max,B)
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




IC_plat=inf_conf_plat(FraMaleData,ages.fit,years.fit,2018:2038,B = 100)


plat=get_estimation_plat(FraMaleData,20:85,1980:2011)
mu_pred=get_predict_plat(plat,2012:2019)$mxt_pred

a=50
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











library(StMoMo)
library(demography)
library(urca)
library(tseries)

data_usa <- read.demogdata("/Users/gojelastat/Desktop/Thèse/Thèse/Données/USA/Mx_1x1.txt",
                           "/Users/gojelastat/Desktop/Thèse/Thèse/Données/USA/Exposures_1x1.txt", 
                           type="mortality", label="USA")

USAMaleData=StMoMoData(data_usa, series = names(data_usa$rate)[2], type ="central")

ages.fit=20:84
years.fit=1961:2005
us_m=USAMaleData$Dxt[21:85,(1961-1933+1):(2005-1933+1)]
us_e=USAMaleData$Ext[21:85,(1961-1933+1):(2005-1933+1)]
wei=genWeightMat(ages=ages.fit,years=years.fit)
plat=plat_function(us_m,us_e,wei,ages.fit,years.fit)




dxt_est=fit_model(plat,20:84,1961:2005)$Dxt_fit

#I=information_plat(dxt_est,ages.fit,years.fit)


#IC_plat=inf_conf_plat(USAMaleData,20:84,1961:2005,2006:2019,B = 500)

plat=get_estimation_plat(USAMaleData,20:84,1961:2005)
mu_pred=get_predict_plat(plat,2006:2019)$mxt_pred

a=
  plot(ts(mu_pred[a,],start = 2012),col=1)
lines(ts(IC_plat$IC_max[a,],start = 2012),col=2)
lines(ts(IC_plat$mean_pred[a,],start = 2012),col=3)
lines(ts(IC_plat$IC_min[a,],start = 2012),col=4)

plat$kt1




arima_coef=function(fit1)
{
  "
  Retourne les vecteurs des coefficients ar, ma et drift du fit1
  
  Args:
    fit1 auto.arima(ts): représentant le fit du auto.arima sur une times 
        series
  
  Returns:
    list: Liste contenant les coefficients ar, ma et drift du modèle fit1
  
  "
  
  ### Ordre ar et ma du modèle
  order=arimaorder(fit1)
  ar_order=order[[1]]
  ma_order=order[[3]]
  
  ### Les coefficients du modèle
  coefficient=fit1$coef
  len_coef=length(coefficient)
  
  ### Les coefficients ar, ma et drift
  ar_coef=0
  ma_coef=0
  drift=0
  
  if(len_coef==(ma_order+ar_order))
  {
    if(ar_order!=0)
    {
      if(ma_order!=0)
      {
        ar_coef=coefficient[1:ar_order]
        ma_coef=coefficient[(ar_order+1):(ar_order+ma_order)]
      }
      else
      {
        ar_coef=coefficient[1:ar_order]
        ma_coef=0
      }
        
    }
    else
    {
      ma_coef=coefficient[1:(ar_order+ma_order)]
      ar_coef=0
    }
  }
  
  else if (len_coef==(ma_order+ar_order+1))
  {
    if(ar_order!=0)
    {
      if(ma_order!=0)
      {
        ar_coef=coefficient[1:ar_order]
        ma_coef=coefficient[(ar_order+1):(ar_order+ma_order)]
        drift_coef=coefficient[(ar_order+ma_order+1)]
      }
      else
      {
        ar_coef=coefficient[1:ar_order]
        ma_coef=0
        drift=coefficient[(ar_order+ma_order+1)]
      }
      
    }
    else
    {
      if(ma_order!=0)
      {
        ma_coef=coefficient[1:(ar_order+ma_order)]
        ar_coef=0
        drift=coefficient[(ar_order+ma_order+1)]
      }
      else
      {
        ma_coef=0
        ar_coef=0
        drift=coefficient[(ar_order+ma_order+1)]
      }
      
    }
  }
  
  list(
    ma_coef,
    ar_coef,
    drift
  )
  
}


sim_ic_method=function(model, ages.fit, years.fit, years.pred, N=10)
{
  
  "
  Simuler plusieurs scénarios de kt1,kt2,... et récupérer les taux de mortalité 
  correspondants à chaque scénario pour obtenir le taux moyenne sur les N scénarios
  et le taux moyenne pour les 25% plus bas et 2,5% plus grand.
  
  Args:
    model: modèle plat,spo, modèle proposé,...
    ages.fit (vecteur): Liste des labels d'âges
    years.fit (vecteur): Liste des labels d'année d'estimation
    years.pred (vecteur): Liste des labels d'année de prédiction
    N (int): Entier pour le nombre de stimulation
    
  Return:
    
  "
  n.a=length(ages.fit)
  n.y=length(years.fit)
  n.t=length(years.pred)
  ### Récupérer les valeurs des paramètres du modèle
  kt1=model$kt1
  kt2=model$kt2
  kt3=model$kt3
  gc=model$gc
  ax=model$ax
  
  ### Fiter un auto.arima sur les paramètres 
  fit1=auto.arima(kt1)
  order_kt1=arimaorder(fit1)
  sigma_kt1=sqrt(fit1$sigma2)
  kt1_ma_coef=c(arima_coef(fit1)[[1]])
  kt1_ar_coef=c(arima_coef(fit1)[[2]])
  kt1_drift=c(arima_coef(fit1)[[3]])
  
  fit2=auto.arima(kt2)
  order_kt2=arimaorder(fit2)
  sigma_kt2=sqrt(fit2$sigma2)
  kt2_ma_coef=c(arima_coef(fit2)[[1]])
  kt2_ar_coef=c(arima_coef(fit2)[[2]])
  kt2_drift=c(arima_coef(fit2)[[3]])
  
  fit3=auto.arima(kt3)
  order_kt3=arimaorder(fit3)
  sigma_kt3=sqrt(fit3$sigma2)
  kt3_ma_coef=c(arima_coef(fit3)[[1]])
  kt3_ar_coef=c(arima_coef(fit3)[[2]])
  kt3_drift=c(arima_coef(fit3)[[3]])
  
  ### Dataframe pour enregistrer les simulations
  mu_pred=array(dim=c(N,length(ages.fit),length(years.pred)))
  
  ### Simuler N scénarios
  for(i in range(N))
  {
    sim1=arima.sim(list(order=order_kt1, 
                   ma=kt1_ma_coef),
                   sd=sqrt(sigma_kt1), n=length(kt1),start.innov = kt1)
    kt1_sim=ts(sim1,start = years.fit[1])
    fore_1=forecast(kt1_sim,h=n.t)$mean
    #forecast_kt1=data.frame(forecast_kt1,fore_1)
    
    
    sim2=arima.sim(list(order=order_kt2,
                   ma=kt2_ma_coef), 
                   sd=sqrt(sigma_kt2), n=length(kt2),start.innov = kt2)
    kt2_sim=ts(sim2,start = years.fit[1])
    fore_2=forecast(kt2_sim,h=n.t)$mean
    #forecast_kt2=data.frame(forecast_kt2,fore_2)
    
    sim3=arima.sim(list(order=order_kt3, 
                        ma=kt3_ma_coef,drift=kt3_drift), 
                   sd=sqrt(sigma_kt3), n=length(kt3),start.innov = kt3)
    kt3_sim=ts(sim3,start = years.fit[1])
    fore_3=forecast(kt3_sim,h=n.t)$mean
    #forecast_kt3=data.frame(forecast_kt3,fore_3)
    
    ### Gamma
    fit_gc=auto.arima(model$gc)
    pc=forecast(fit_gc,h=150)
    gc_pred=pc$mean
    gg=cbind(model$gc,gc_pred)
    
    ### gc
    cohort=(years.fit[1]-ages.fit[n.a]):(years.fit[n.y]-ages.fit[1])
    
    
    ch_beg=years.fit[n.y]+1-ages.fit[n.a] ### start of cohort forecasts
    ch_end=years.fit[n.y]+n.t-ages.fit[1] ### end of cohort forecasts
    gcc=c(gg[(ch_beg-cohort[1]+1):(tail(cohort,1)-cohort[1]+1),1],
          gg[(tail(cohort,1)-cohort[1]+2):(ch_end-cohort[1]+1),2])
    
    ## Construire la matrice teta
    ax=model$ax
    teta=c(ax,fore_1,fore_2,fore_3,gcc)
    
    x=ages.fit
    x.bar=mean(x)
    n=n.a*n.t
    n.c=n.a+n.t-1
    
    
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
    
    log_mu=matrix(mu,nrow = n.a,dimnames =list(ages.fit,years.pred))
    mod_mu=exp(log_mu)
    
    mu_pred[i,,]=mod_mu
  }
  
  mu_pred
  
}


sim_ic_method(plat,ages.fit = 20:84,years.fit = 1961:2005,years.pred = 2006:2015)




model=plat
ages.fit= 20:84
years.fit = 1961:2005
years.pred = 2006:2015
n.a=length(ages.fit)
n.y=length(years.fit)
n.t=length(years.pred)
### Récupérer les valeurs des paramètres du modèle
kt1=model$kt1
kt2=model$kt2
kt3=model$kt3
gc=model$gc
ax=model$ax

### Fiter un auto.arima sur les paramètres 
fit1=auto.arima(kt1)
order_kt1=arimaorder(fit1)
sigma_kt1=sqrt(fit1$sigma2)
kt1_ma_coef=c(arima_coef(fit1)[[1]])
kt1_ar_coef=c(arima_coef(fit1)[[2]])
kt1_drift=c(arima_coef(fit1)[[3]])

fit2=auto.arima(kt2)
order_kt2=arimaorder(fit2)
sigma_kt2=sqrt(fit2$sigma2)
kt2_ma_coef=c(arima_coef(fit2)[[1]])
kt2_ar_coef=c(arima_coef(fit2)[[2]])
kt2_drift=c(arima_coef(fit2)[[3]])

fit3=auto.arima(kt3)
order_kt3=arimaorder(fit3)
sigma_kt3=sqrt(fit3$sigma2)
kt3_ma_coef=c(arima_coef(fit3)[[1]])
kt3_ar_coef=c(arima_coef(fit3)[[2]])
kt3_drift=c(arima_coef(fit3)[[3]])
N=100
### Dataframe pour enregistrer les simulations
mu_pred=array(dim=c(N,length(ages.fit),length(years.pred)))

### Simuler N scénarios
for(i in 1:N)
{
  sim1=arima.sim(list(order=order_kt1, 
                      ma=kt1_ma_coef),
                sd=sigma_kt1, n=length(kt1),start.innov = kt1)
  #sim1=simulate(fit1)
  kt1_sim=ts(sim1,start = years.fit[1])
  fore_1=forecast(kt1_sim,h=n.t)$mean
  #forecast_kt1=data.frame(forecast_kt1,fore_1)
  
  
  sim2=arima.sim(list(order=order_kt2,
                      ma=kt2_ma_coef), 
                 sd=sigma_kt2, n=length(kt2),start.innov = kt2)
  #sim2=simulate(fit2)
  kt2_sim=ts(sim2,start = years.fit[1])
  fore_2=forecast(kt2_sim,h=n.t)$mean
  #forecast_kt2=data.frame(forecast_kt2,fore_2)
  
  sim3=arima.sim(list(order=order_kt3), 
                      #ma=kt3_ma_coef,drift=kt3_drift), 
               sd=sigma_kt3, n=length(kt3),start.innov = kt3)
  #sim3=simulate(fit3)
  kt3_sim=ts(sim3,start = years.fit[1])
  fore_3=forecast(kt3_sim,h=n.t)$mean
  #forecast_kt3=data.frame(forecast_kt3,fore_3)
  
  
  ### Gamma
  fit_gc=auto.arima(model$gc)
  pc=forecast(fit_gc,h=150)
  gc_pred=pc$mean
  gg=cbind(model$gc,gc_pred)
  
  ### gc
  cohort=(years.fit[1]-ages.fit[n.a]):(years.fit[n.y]-ages.fit[1])
  
  
  ch_beg=years.fit[n.y]+1-ages.fit[n.a] ### start of cohort forecasts
  ch_end=years.fit[n.y]+n.t-ages.fit[1] ### end of cohort forecasts
  gcc=c(gg[(ch_beg-cohort[1]+1):(1946-cohort[1]+1),1],
        gg[(1946-cohort[1]+2):(ch_end-cohort[1]+1),2])
  
  ## Construire la matrice teta
  ax=model$ax
  teta=c(ax,fore_1,fore_2,fore_3,gcc)
  
  x=ages.fit
  x.bar=mean(x)
  n=n.a*n.t
  n.c=n.a+n.t-1
  
  
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
  
  log_mu=matrix(mu,nrow = n.a,dimnames =list(ages.fit,years.pred))
  mod_mu=exp(log_mu)
  mu_pred[i,,]=mod_mu
}

mu_pred[100,,]


y <- ts(arima.sim(model=list(order = c(1,1,1), ar=.8,ma=.7), 100)) # Simulate data
yseries <- Arima(y,order=c(1,1,1))
simyseries <- ts(replicate(10000, simulate(fit1, nsim=40)),start=end(kt1)+1) # Change the first parameter of replicate() to change the number os simulated paths
matplot(cbind(kt1,simyseries), type='l')

dim(simyseries)
fore_3

X%*%teta
