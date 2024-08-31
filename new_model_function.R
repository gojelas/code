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

Reg = function(x , y , h = 1){
  n = length(x)
  vect.x = NULL
  vect.y = NULL
  for(i in 1:(n-1)){
    a = (y[i+1] - y[i])/(x[i+1] - x[i])
    b = y[i] - a*x[i]
    dx = (x[i+1] - x[i])/(h+1)
    current.x = x[i]
    vect.x = c(vect.x , x[i])
    vect.y = c(vect.y , y[i])
    for(step in 1:h){
      current.x = current.x + dx
      current.y = a*current.x + b
      vect.x = c(vect.x , current.x)
      vect.y = c(vect.y , current.y)
    }
  }
  vect.x = c(vect.x , x[n])
  vect.y = c(vect.y , y[n])
  return(list("x" = vect.x, 
              "y" = vect.y))
}


reshape_french_data=function(data)
{
  ag=dim(data$Dxt)[1]
  c=c()
  for (i in 1:dim(data$Dxt)[2])
  {
    c=cbind(c,mean(data$Dxt[86:ag,i],na.rm = T))
  }
  c
  
  data_dxt=data$Dxt[1:86,]
  data_dxt[86,]=c
  
  
  ag=dim(data$Ext)[1]
  d=c()
  for (i in 1:dim(data$Ext)[2])
  {
    d=cbind(d,mean(data$Ext[86:ag,i],na.rm = T))
  }
  d
  
  data_ext=data$Ext[1:86,]
  data_ext[86,]=d
  
  data$Dxt[1:86,]=data_dxt
  data$Ext[1:86,]=data_ext
  
  list(dxt=data$Dxt,ext=data$Ext)
  
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



non_null=function(I_t)
{
  "
  Récupère les valeurs non nulles de I_t
  Args: 
    I_t: données de températures ou de d'amplitude thermique
  Return:
    y: vecteur des valeurs non nulles de I_t
  "
  y=c()
  l=length(I_t)
  for (i in 1:l) 
  {
    if(I_t[i]>0)
      y=c(y,I_t[i])
  }
  y
}

indicatrice_It=function(I.t,mean_I.t)
{
  "
  Indicateur d'année chaude, qui vaut 1 lorsque l'écart à la moyenne
  de température est positive et zéro sinon.
  
  Args: 
    I_t: données de température
    mean_It: moyenne de la température sur la période d'estimation du modèle
  
  Return:
    Ind: vecteur d'indicatrice d'année chaude.
  "
  
  Ind = c()
  ecart_moy = pmax((I.t-mean_I.t),0)
  for (i in 1:length(I.t))
  {
    if (ecart_moy[i]>0)
      Ind = cbind(Ind,1)
    else
      Ind = cbind(Ind,0)
  }
  
  c(Ind)
}



new_model_function=function(dxt,ext,wxt,Ic,a,c_x,a_c,ages,years,x1,x2)
{
  "
  C'est une variante du premier modèle, ici on prend l'indicatrice d'année chaude 
  plutôt que de prendre l'amplitude pour cette année.
  "
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
  #I_t=pmax((Ic-a),0)
  I_t=indicatrice_It(Ic,a)
  #I_t1=pmax((Ic-a),0)
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
  X_t4=kronecker(diag(n.t),as.matrix((pmax(x1-x,0)+c_x*pmax(x-x2,0))^2))
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
  n.s=dim(X_c)[2]
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
  
  ### Constraint 4                                                                                           
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
  #h7=cbind(t(rep(0,n.a+4*n.t+n.t5)),t((1:n.s)^2))
  
  H=rbind(h1,h2,h3,h4,h5,h6)#,h7)
  
  
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
    cohort=cohort,
    a_m=a,
    years.fit = years,
    model=model
  ))
  
}




### compute model
get_estimation_new_model=function(data, ages, years, c, a_c=65, x1, x2)
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
  Ic=df[(years[1]-1961+1):(tail(years,1)-1961+1),1]
  a_m=mean(Ic)
  
  m2=new_model_function(fr_m,fr_e,wei,Ic,a_m,c,a_c,ages,years,x1,x2)
  
  list(
    years=years,
    ages=ages,
    cohort=m2$cohort,
    kt1=m2$kt1,
    kt2=m2$kt2,
    kt3=m2$kt3,
    kt4=m2$kt4,
    kt5=m2$kt5,
    ax=m2$ax,
    gc=m2$gc,
    a_m=a_m,
    x1=x1,
    x2=x2,
    c=c
  )
}


get_predict_new_model=function(estim,years_pred,i,method="ARIMAX")
{
  "
  Prédire la taux de mortalité sur la période years_pred.
  Args:
    estim: modèle donné par la fonction get_estimation
    years_pred: la période de prévison
    i indique le choix du scénario RCP, 1 pour RCP2.6, 2 pour RCP4.5 et 3 pour RCP8.5
  
  Returns:
    mxt_pred: matrice des prévisions de taux de mortalité
    years_pred: période de prévision.
  "
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
  
  fit4=auto.arima(estim$kt4)
  fore_4=forecast(fit4, h)
  kt4_pred=fore_4$mean
  
  x = seq(from = 1 , to = length(estim$kt5) , by = 1)
  y = estim$kt5
  l = Reg(x, y, 1)# Data augmentation by inserting two points.
  kt5_est=ts(l$y,start=1,frequency = 1)
  
  years.fit = estim$years
  I_t = df[(years.fit[1]-1961+1):(tail(years.fit,1)-1961+1),1] # Température sur le train
  
  Ic=df[(years_pred[1]-1961+1):(tail(years_pred,1)-1961+1),i]
  a_m=estim$a_m
  
  if (method== "ARIMAX")
  {
    I_t_pred=pmax((Ic-a_m),0)
    I_t_pred=I_t_pred[I_t_pred>0]
    
    I_t = pmax((I_t-a_m),0)
    I_t = I_t[I_t>0] # Récupérer les points non nuls dans I_t
    
    y2 = I_t
    x2 = seq(1, length(I_t))
    l2 = Reg(x2 , y2, h=1) # Augmentation des points de I_t
    I_t_est=ts(l2$y, start = 1, frequency = 1)
    
    y3 = I_t_pred
    x3 = seq(1, length(I_t_pred))
    l3 = Reg(x3, y3, 1) # Augmentation des points de I_t_pred
    I_t_pred=ts(l3$y, start = 1, frequency = 1)
    
    
    data3=ts(cbind(kt5_est,I_t_est), start = 1)
    fit1=auto.arima(data3[,"kt5_est"], xreg = data3[,"I_t_est"])
    
    forecast_values=forecast(fit1,xreg = I_t_pred)
    t=seq(1,(length(I_t_pred)),by=2)
    kt5_pred= forecast_values$mean[t]
  } else if(method=="ARIMA"){
    fit5=auto.arima(kt5_est)
    #fit5=Arima(kt5_est,order = kt5_order,include.mean = TRUE)
    fore_5=forecast(fit5,20*h)
    t=seq(2,(20*h),by=2)
    kt5_pred=fore_5$mean[t]
    #plot(fore_5)
    #print(kt5_pred)
  }
  
  
  
  fit_gc=Arima(estim$g,order = c(1,1,1))
  #fit_gc = auto.arima(estim$gc)
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
  c4=kt4_pred
  
  n.a=length(estim$ages)
  n.t=h
  n=n.a*n.t
  n.c=n.a+n.t-1
  
  x=estim$ages
  x.bar=mean(x)
  a_m=estim$a_m
  Ic=df[(years_pred[1]-1961+1):(tail(years_pred,1)-1961+1),i]
  
  
  
  ###Construction de la matrice de kt_5
  #I_t=pmax((Ic-a_m),0)
  I_t= indicatrice_It(Ic, a_m)
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
  axx=estim$ax
  teta=c(axx,c1,c2,c3,c4,c5,gcc)
  
  a_c=65
  
  #### 1) Matrice d'âge
  X_a=kronecker(as.matrix(rep(1,n.t)),diag(n.a))
  #### 2) Matrice de période kt1
  X_t1=kronecker(diag(n.t),as.matrix(rep(1,n.a)))
  #### 3) Matrice de période kt2
  X_t2=kronecker(diag(n.t),as.matrix(x.bar-x))
  #### 4) Matrice de période kt3
  X_t3=kronecker(diag(n.t),as.matrix(pmax(x.bar-x,0)))
  #### 4) Matrice de période kt4
  X_t4=kronecker(diag(n.t),as.matrix(pmax((estim$x1)-x,0)+(estim$c)*pmax(x-(estim$x2),0))^2)
  #### 5) Matrice de période kt5
  y=(x-a_c)
  #y=I(x,a_c)
  X_t5=kronecker(I_t,as.matrix(pmax(y,0)))
  ### Matrice de cohorte
  X_c=indic(matrix(0,n,n.c),n.a,n.t)
  
  #### Enfin on obtient
  X=cbind(X_a,X_t1,X_t2,X_t3,X_t4,X_t5,X_c)
  mu=X%*%teta
  
  log_mu=matrix(mu,nrow = n.a,dimnames =list(estim$ages,years_pred))
  mod_mu=exp(log_mu)
  
  list(
    mxt_pred=mod_mu,
    years_pred=years_pred
  )
  
}




get_best_params_new_model=function(data,c,x1,x2,ages,years,Ic,method="ARIMAX")
{
  #initialize parameters
  c_0=1
  x1_0=0
  x2_0=0
  mape1=1
  mape2=1
  
  y1=as.integer(0.85*length(years))
  years1=years.fit[1:y1]
  years2=years[(y1+1):length(years)]
  age_head=data$ages[1]
  age_tail=tail(data$ages,1)
  year_head=data$years[1]
  year_tail=tail(data$years,1)
  
  dx_obs=data$Dxt[(ages[1]-age_head+1):(tail(ages,1)-age_head+1),
                  (years2[1]-year_head+1):(tail(years2,1)-year_head+1)]
  ex_obs=data$Ext[(ages[1]-age_head+1):(tail(ages,1)-age_head+1),
                  (years2[1]-year_head+1):(tail(years2,1)-year_head+1)]
  mu_obs=dx_obs/ex_obs
  
  for (i in 1:length(x1))
  {
    for (j in 1:length(x2))
    {
      for (x in 1:length(c)){
        est=get_estimation_new_model(data,ages,years1,c[x],a_c = 65,x1[i],x2[j])
        a_m=est$a_m
        Ic=df[(years1[1]-1961+1):(tail(years1,1)-1961+1),1]
        fr_m=data$Dxt[21:86,(years1[1]-1816+1):(tail(years1,1)-1816+1)]
        fr_e=data$Ext[21:86,(years1[1]-1816+1):(tail(years1,1)-1816+1)]
        wei=genWeightMat(ages=ages,years=years1)
        model=new_model_function(fr_m,fr_e,wei,Ic,a_m,c[x],a_c,ages,years,x1[i],x2[j])
        mu_pred=get_predict_new_model(est,years_pred = years2,1,method)
        m2=mean(abs((mu_pred$mxt_pred-mu_obs)/mu_obs))
        m1=mape(model,ages,years1)
        if(m2 < mape2 & m1<mape1)
        {
          #mape2=rbind(mape2,m2)
          #c_0=rbind(c_0,c[x])
          #mape1=rbind(mape1,m1)
          mape2=m2
          mape1=m1
          c_0=c[x]
          x1_0=x1[i]
          x2_0=x2[j]
          print(m1)
          print(m2)
          print(c_0)
          print(x1_0)
          print(x2_0)
        }
      }
    }
  }
  #  })
  perf=data.frame(c_0,x1_0,x2_0,mape1,mape2,row.names = 1)
  colnames(perf)=c('x1','x2','mape_train','mape_test')
  list(
    c=c_0,
    x1=x1_0,
    x2=x2_0,
    mape1=mape1,
    mape2=mape2,
    performance=perf
  )
}




sim_ic_method_new_model=function(model, ages.fit, years.fit, years.pred, x1, x2, c, 
                                 a_c=65, N=100,i=1,method="ARIMAX")
{
  
  "
  Simuler plusieurs scénarios de kt1,kt2,... et récupérer les taux de mortalité 
  correspondants à chaque scénario pour obtenir le taux moyenne sur les N scénarios
  et le taux moyenne pour les 25% plus bas et 2,5% plus grand.
  (kt1,kt2,kt3,kt4) est considéré comme série temporelle multivariée et est
  estimée par MRWD (Multivariate Random Walk with Drift)
  
  Args:
    model: modèle plat,spo, modèle proposé,...
    ages.fit (vecteur): Liste des labels d'âges
    years.fit (vecteur): Liste des labels d'année d'estimation
    years.pred (vecteur): Liste des labels d'année de prédiction
    N (int): Entier pour le nombre de stimulation
    
  Return:
    
  "
  
  #model=model_m_new
  #years.fit=1980:2011
  #ages.fit = 20:85
  #years.pred = 2012:2030
  
  #x1=45
  #x2=60
  #c=.58
  #a_c=65
  #N=50
  #i=1
  #method="ARIMAX"
  n.a=length(ages.fit)
  n.y=length(years.fit)
  n.t=length(years.pred)
  ### Récupérer les valeurs des paramètres du modèle
  kt1=model$kt1
  kt2=model$kt2
  kt3=model$kt3
  kt4=model$kt4
  kt5=model$kt5
  gc=model$gc
  ax=model$ax
  
  ### Fiter un auto.arima sur les paramètres
  
  #fit1=auto.arima(kt1)
  #sim_kt1=ts(replicate(N,simulate(fit1,nsim = n.t)),start = end(kt1))
  
  #fit2=auto.arima(kt2)
  #sim_kt2=ts(replicate(N,simulate(fit2,nsim = n.t)),start = end(kt2))
  
  #fit3=auto.arima(kt3)
  #sim_kt3=ts(replicate(N,simulate(fit3,nsim = n.t)),start = end(kt3))
  
  #fit4=auto.arima(kt4)
  #sim_kt4=ts(replicate(N,simulate(fit4,nsim = n.t)),start = end(kt4))
  
  mrwd_vect=rbind(kt1,kt2,kt3,kt4)
  colnames(mrwd_vect)=years.fit
  
  mrwd_kt=mrwd(mrwd_vect)
  
  sim_kt=replicate(N,simulate(mrwd_kt,nsim=n.t))
  
  #matplot(sim_kt)
  
  x = seq(from = 1 , to = length(kt5))
  y = kt5
  l = Reg(x , y, 1)# Data augmentation by inserting two points.
  
  kt5_est=ts(l$y,start=1,frequency = 1)
  #fit5=Arima(kt5_est,order = kt5_order,include.mean = val)
  #fit5=auto.arima(kt5_est)
  #sim_kt5=ts(replicate(N,simulate(fit5,nsim = 8*n.t)),start = end(kt5_est))
  #t=seq(3,(8*n.t),by=3)
  #sim_kt5=sim_kt5[t,]
  #fore_5=forecast(fit5,20*h)
  #t=seq(3,(20*h),by=3)
  #kt5_pred=fore_5$mean[t]
  
  #years.pred=2012:2040
  #i=1
  Ic = df[(years.pred[1]-1961+1):(tail(years.pred,1)-1961+1),i]
  a_m = model$a_m
  
  I_t_pred=pmax((Ic-a_m),0)
  I_t_pred=I_t_pred[I_t_pred>0]
  
  I_c = df[(years.fit[1]-1961+1):(tail(years.fit,1)-1961+1),i]
  
  I_t_est = pmax((I_c-a_m),0)
  I_t_est = I_t_est[I_t_est>0]
  
  
  if (method=="ARIMAX")
  {
    amy2 = I_t_est
    amx2 = seq(1, length(I_t_est))
    log2 = Reg(amx2, amy2, h=1)
    I_t_est=ts(log2$y, start = 1, frequency = 1)
    
    amy3 = c(tail(I_t_est,1),I_t_pred)
    amx3 = seq(1, length(amy3))
    log3 = Reg(amx3, amy3, h=1) # Augmentation des points de I_t_pred
    I_t_pred=ts(log3$y[2:length(log3$y)], start = 1, frequency = 1)
    
    
    dat3=ts(cbind(kt5_est,I_t_est), start = 1)
    fit_1=auto.arima(dat3[,"kt5_est"], xreg = dat3[,"I_t_est"])
    
    sim_kt5=ts(replicate(N,simulate(fit_1,nsim = length(I_t_pred),
                                    xreg = I_t_pred)),start = end(kt5_est))
    #forecast_values=forecast(fit1,xreg = I_t_pred)
    tu=seq(2,length(I_t_pred),by=2)
    sim_kt5=sim_kt5[tu,]
    #colMeans(sim_kt5)
    #matplot(sim_kt5)
    #kt5_pred= forecast_values$mean[t]
  } else if(method=="ARIMA"){
    fit_1=auto.arima(kt5_est)
    sim_kt5=ts(replicate(N,simulate(fit_1,nsim = length(I_t_pred)*2)),start = end(kt5_est))
    t=seq(1,(length(I_t_pred)*2),by=2)
    sim_kt5=sim_kt5[t,]
    
    #fore_5=forecast(fit5,20*h)
    #t=seq(2,(20*h),by=2)
    #kt5_pred=fore_5$mean[t]
    #plot(fore_5)
    #print(kt5_pred)
  }
  
  
  
  
  #Ic=df[(years.pred[1]-1961+1):(tail(years.pred,1)-1961+1),i]
  
  ###Construction de la matrice de kt_5
  #I_t=pmax((Ic-a_m),0)
  Ic = df[(years.pred[1]-1961+1):(tail(years.pred,1)-1961+1),i]
  a_m = model$a_m
  
  I_t = indicatrice_It(Ic, a_m)
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
  sim_kt5=sim_kt5[1:n.t5,]
  #matplot(sim_kt5)
  #plot(sim_kt5[,7],type="l")
  #plot(colMeans(sim_kt5),type="l")
  #sim_kt5=kt5_pred[1:n.t5]
  
  ### Gamma
  #fit_gc=Arima(model$gc,order = c(1,0,0),include.mean = FALSE)
  fit_gc=Arima(model$gc,order = c(1,1,1))
  pc=forecast(fit_gc,h=150)
  sim_gc=ts(replicate(N,simulate(fit_gc,nsim = 150)),start = 1951)
  gc_pred=pc$mean
  #gc_pred=sim_gc
  #gg=cbind(model$gc,gc_pred)
  
  ### gc
  #cohort=(years.fit[1]-ages.fit[n.a]):(years.fit[n.y]-ages.fit[1])
  
  
  #ch_beg=years.fit[n.y]+1-ages.fit[n.a] ### start of cohort forecasts
  #ch_end=years.fit[n.y]+n.t-ages.fit[1] ### end of cohort forecasts
  #gcc=c(gg[(ch_beg-cohort[1]+1):(1950-cohort[1]+1),1],
  #     gg[(1950-cohort[1]+2):(ch_end-cohort[1]+1),2])
  
  
  ### Dataframe pour enregistrer les simulations
  mu_pred=array(dim=c(N,length(ages.fit),length(years.pred)))
  
  ### Simuler N scénarios
  for(i in 1:N)
  {
    
    #gg=cbind(model$gc,gc_pred[,i])
    gg=cbind(model$gc,gc_pred)
    ### gc
    cohort=(years.fit[1]-ages.fit[n.a]):(years.fit[n.y]-ages.fit[1])
    
    
    ch_beg=years.fit[n.y]+1-ages.fit[n.a] ### start of cohort forecasts
    ch_end=years.fit[n.y]+n.t-ages.fit[1] ### end of cohort forecasts
    gcc=c(gg[(ch_beg-cohort[1]+1):(1950-cohort[1]+1),1],
          gg[(1950-cohort[1]+2):(ch_end-cohort[1]+1),2])
    
    ## Construire la matrice teta
    ax=model$ax
    #teta=c(ax,sim_kt1[,i],sim_kt2[,i],sim_kt3[,i],sim_kt4[,i],sim_kt5,gcc)
    
    teta=c(ax,sim_kt[1,,i],sim_kt[2,,i],sim_kt[3,,i],sim_kt[4,,i],sim_kt5[,i],gcc)
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
    #### 4) Matrice de période kt4
    X_t4=kronecker(diag(n.t),as.matrix(pmax(x1-x,0)+c*pmax(x-x2,0))^2)
    #### 5) Matrice de période kt5
    y=(x-a_c)
    #y=I(x,a_c)
    X_t5=kronecker(I_t,as.matrix(pmax(y,0)))
    ### Matrice de cohorte
    X_c=indic(matrix(0,n,n.c),n.a,n.t)
    
    #### Enfin on obtient
    X=cbind(X_a,X_t1,X_t2,X_t3,X_t4,X_t5,X_c)
    mu=X%*%teta
    
    log_mu=matrix(mu,nrow = n.a,dimnames =list(ages.fit,years.pred))
    mod_mu=exp(log_mu)
    
    mu_pred[i,,]=mod_mu
  }
  
  
  mean_pred=matrix(data=0, nrow = dim(mu_pred)[2],ncol=dim(mu_pred)[3],dimnames = list(ages.fit,years.pred))
  se_pred=matrix(data=0, nrow = dim(mu_pred)[2],ncol=dim(mu_pred)[3],dimnames = list(ages.fit,years.pred))
  IC_min=matrix(data=0, nrow = dim(mu_pred)[2],ncol=dim(mu_pred)[3],dimnames = list(ages.fit,years.pred))
  IC_max=matrix(data=0, nrow = dim(mu_pred)[2],ncol=dim(mu_pred)[3],dimnames = list(ages.fit,years.pred))
  
  B=N
  for (x in 1:dim(mu_pred)[2])
  {
    for (y in 1:dim(mu_pred)[3])
    {
      mean_pred[x,y]=mean(mu_pred[,x,y])
      #se_pred[x,y]=sd(mu_pred[,x,y])
      #IC_min[x,y]=mean(sort(mu_pred[,x,y])[1:500])
      #IC_max[x,y]=mean(sort(mu_pred[,x,y],decreasing = TRUE)[1:500])
      se_pred[x,y]=sd(mu_pred[,x,y])
      IC_min[x,y]=mean_pred[x,y]-1.96*se_pred[x,y]/sqrt(B)
      IC_max[x,y]=mean_pred[x,y]+1.96*se_pred[x,y]/sqrt(B)
      
    }
  }
  list(mu_pred=mu_pred,
       mean_pred=mean_pred,
       #se_pred=se_pred,
       IC_min=IC_min,
       IC_max=IC_max)
  
}
