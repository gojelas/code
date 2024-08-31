library(StMoMo)
library(demography)
library(urca)
library(tseries)



### Importation des données
data_uk <- read.demogdata("Data/Mx_1x1 (UK).txt",
                       "Data/Exposures_1x1 (UK).txt", 
                       type="mortality", label="UK")
UKMaleData=StMoMoData(data_uk, series = names(data_uk$rate)[2], type ="central")


UKMaleData$Dxt=reshape_french_data(UKMaleData)$dxt
UKMaleData$Ext=reshape_french_data(UKMaleData)$ext
male_rate=FraMaleData$Dxt/FraMaleData$Ext


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
  bic=-0.5*log(model$nobs)*(model$npar)+(loglik_mod)
  bic
}


mape=function(model,ages.fit,years.fit)
{
  model_dxt=fit_model(model,ages.fit,years.fit)$Dxt_fit
  ext=model$Ext
  dxt=model$Dxt
  mx_fit=model_dxt/ext
  mx=dxt/ext
  mean(abs((mx_fit-mx)/mx))
  
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
spo_function=function(dxt,ext,wxt,c_x,a_c,ages,years)
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
  
  #### 4) period matrix kt4
  X_t4=kronecker(diag(n.t),as.matrix((pmax(a_c-x,0)+c_x*pmax(x-a_c,0))^2))
  
  ### cohort 
  X_c=indic(matrix(0,n,n.c),n.a,n.t)
  
  ### consider cohort that under 1950 years
  n.s=1945-cohort[1]+1 
  X_c=X_c[,1:n.s]
  n.s=dim(X_c)[2]
  X=cbind(X_a,X_t1,X_t2,X_t3,X_t4,X_c)
  
  #### Constraint matrix
  
  ### Constraint 1
  h1=cbind(t(rep(0,n.a)),t(rep(1,n.t)),t(rep(0,3*n.t+n.s)))
  
  ### Constraint 2
  h2=cbind(t(rep(0,n.a+n.t)),t(rep(1,n.t)),t(rep(0,2*n.t+n.s)))
  
  ### Constraint 3
  h3=cbind(t(rep(0,n.a+2*n.t)),t(rep(1,n.t)),t(rep(0,n.t+n.s)))
  
  ### Constraint 4                                                                                          
  h4=cbind(t(rep(0,n.a+3*n.t)),t(rep(1,n.t)),t(rep(0,n.s)))
  
  ### Constraint 5
  h5=cbind(t(rep(0,n.a+4*n.t)),t(rep(1,n.s)))
  
  ### Constraint 6
  h6=cbind(t(rep(0,n.a+4*n.t)),t(1:n.s))
  
  ### Constraint 7
  #h7=cbind(t(rep(0,n.a+4*n.t)),t((1:n.s)^2))
  
  H=rbind(h1,h2,h3,h4,h5,h6)#,h7)
  
  
  #### Model 
  model=gnm(v_dxt~-1+offset(v_ext)+X,weights = v_wxt,family = poisson(link = 'log'))
  
  
  
  #### coefficients matrix
  teta=solve(t(X)%*%X+t(H)%*%H,tol=5.41168e-27)%*%t(X)%*%(log(model$fitted.values)-v_ext)
  
  ax=teta[1:n.a]
  kt1=teta[(n.a+1):(n.a+n.t)]
  kt2=teta[(n.a+n.t+1):(n.a+2*n.t)]
  kt3=teta[(n.a+2*n.t+1):(n.a+3*n.t)]
  kt4=teta[(n.a+3*n.t+1):(n.a+4*n.t)]
  #kt5=teta[(n.a+4*n.t+1):(n.a+4*n.t+n.t5)]
  gc=teta[(n.a+4*n.t+1):dim(teta)[1]]
  
  
  structure(list(
    ax=ax,
    kt1=ts(kt1,start=years[1],frequency=1),
    kt2=ts(kt2,start=years[1],frequency=1),
    kt3=ts(kt3,start=years[1],frequency=1),
    kt4=ts(kt4,start=years[1],frequency=1),
    gc=ts(gc,start = cohort[1]+3,frequency = 1),
    Dxt=dxt,
    Ext=ext,
    wxt=wxt,
    npar=dim(teta)[1],
    nobs=n.a*n.t,
    rank=qr(X)$rank,
    col=dim(X)[2],
    sing=dim(X)[2]-qr(X)$rank,
    model=model
  ))
  
}

##############################################
#### Application to male data from France ####
##############################################

### compute model

ages.fit=20:85
years.fit=1974:2006
uk_d=UKMaleData$Dxt[(ages.fit+1),(1974-1922+1):(2006-1922+1)]
uk_e=UKMaleData$Ext[(ages.fit+1),(1974-1922+1):(2006-1922+1)]
wei=genWeightMat(ages=ages.fit,years=years.fit)
a_c=50
c_x=-0.5
spo=spo_function(uk_d,uk_e,wei,c_x,a_c,ages.fit,years.fit)

plat=plat_function(uk_d,uk_e,wei,ages.fit,years.fit)



get_criterion=function(model, data, ages.fit, years.fit)
{
  "Compile deviance, bic, mape, mad and mse"
  wei=genWeightMat(ages=ages.fit,years=years.fit)
  loglik_m2=loglik(model,data,ages.fit,years.fit,wei)
  #aic=aic_mod(model,FraMaleData,ages.fit,years.fit,wei)
  bic=bic_model(model,data,ages.fit,years.fit,wei)
  mape=mape(model,ages.fit,years.fit)
  mad=mad(model,ages.fit,years.fit)
  mse=mse(model,ages.fit,years.fit)
  structure(list(
    bic=bic,
    mape=mape,
    mad=mad, 
    mse=mse))
}

get_criterion(spo,UKMaleData,ages.fit = ages.fit,years.fit = years.fit)

plot(spo$kt3)
