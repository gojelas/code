library(StMoMo)
library(demography)
library(urca)
library(tseries)


source(file="/Users/gojelastat/Desktop/Thèse/fonctions_model.R")


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
  n.s=dim(X_c)[2]
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
    gc=ts(gc,start = (cohort[1]),frequency = 1),
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
  
  #fit1=auto.arima(estim$kt1)
  fit1=Arima(estim$kt1,order = c(0,1,0), include.mean = FALSE)
  fore_1=forecast(fit1, h)
  kt1_pred=fore_1$mean
  
  #fit2=auto.arima(estim$kt2)
  fit2=Arima(estim$kt2,order = c(1,0,0), include.mean = FALSE)
  fore_2=forecast(fit2, h)
  kt2_pred=fore_2$mean
  
  #fit3=auto.arima(estim$kt3)
  fit3=Arima(estim$kt3,order = c(1,0,0), include.mean = FALSE)
  fore_3=forecast(fit3, h)
  kt3_pred=fore_3$mean
  
  
  #fit_gc=auto.arima(estim$g)
  fit_gc=Arima(estim$g,order = c(1,0,0), include.mean = FALSE)
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


plot(plat$gc)
fit1=Arima(plat$gc,order = c(1,0,0),include.mean = FALSE)
fit1


sim_ic_method=function(model, ages.fit, years.fit, years.pred, N=100)
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
  fit1=Arima(kt1, order=c(0,1,0), include.drift = TRUE)
  sim_kt1=ts(replicate(10000,simulate(fit1,nsim = n.t)),start = end(kt1)+1)
  
  fit2=Arima(kt2, order=c(1,0,0),include.mean = FALSE)
  sim_kt2=ts(replicate(10000,simulate(fit2,nsim = n.t)),start = end(kt2)+1)
  
  
  fit3=Arima(kt3, order=c(1,0,0),include.mean = FALSE)
  sim_kt3=ts(replicate(10000,simulate(fit3,nsim = n.t)),start = end(kt3)+1)
  
  ### Gamma
  fit_gc=Arima(model$gc, order = c(1,0,0),include.mean = FALSE)
  pc=forecast(fit_gc,h=150)
  gc_pred=pc$mean
  gg=cbind(model$gc,gc_pred)
  
  ### gc
  cohort=(years.fit[1]-ages.fit[n.a]):(years.fit[n.y]-ages.fit[1])
  
  
  ch_beg=years.fit[n.y]+1-ages.fit[n.a] ### start of cohort forecasts
  ch_end=years.fit[n.y]+n.t-ages.fit[1] ### end of cohort forecasts
  #gcc=c(gg[(ch_beg-cohort[1]+1):(tail(cohort,1)-cohort[1]+1),1],
   #     gg[(tail(cohort,1)-cohort[1]+2):(ch_end-cohort[1]+1),2])
  gcc=c(gg[(ch_beg-cohort[1]+1):(1946-cohort[1]+1),1],
        gg[(1946-cohort[1]+2):(ch_end-cohort[1]+1),2])
  
  
  ### Dataframe pour enregistrer les simulations
  mu_pred=array(dim=c(N,length(ages.fit),length(years.pred)))
  
  ### Simuler N scénarios
  for(i in 1:N)
  {
    ## Construire la matrice teta
    ax=model$ax
    teta=c(ax,sim_kt1[,i],sim_kt2[,i],sim_kt3[,i],gcc)
    
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
      se_pred[x,y]=sd(mu_pred[,x,y])
      #IC_min[x,y]=mean_pred[x,y]-1.96*se_pred[x,y]/sqrt(B)
      #IC_max[x,y]=mean_pred[x,y]+1.96*se_pred[x,y]/sqrt(B)
      IC_min[x,y]=mean(sort(mu_pred[,x,y])[1:500])
      IC_max[x,y]=mean(sort(mu_pred[,x,y],decreasing = TRUE)[1:500])
    }
  }
  list(mu_pred=mu_pred,
       mean_pred=mean_pred,
       #se_pred=se_pred,
       IC_min=IC_min,
       IC_max=IC_max)
  
}


m = sim_ic_method(plat, ages.fit = 20:84, years.fit = 1961:2005, years.pred = 2006:2050, N = 10000)

IC_plat = m

plat=get_estimation_plat(USAMaleData,20:84,1961:2005)
mu_pred=get_predict_plat(plat,2006:2050)$mxt_pred
rate=data_usa$rate$male[21:85,(1961-1933+1):(2011-1933+1)]
a=65
plot(ts(rate[a,],start = 1961),col=1,xlim=c(1961,2050),
     ylim=c(min(IC_plat$mean_pred[a,]),max(rate[a,])))
lines(ts(IC_plat$IC_max[a,],start = 2006),col=2)
lines(ts(IC_plat$mean_pred[a,],start = 2006),col=3)
lines(ts(IC_plat$IC_min[a,],start = 2006),col=4)

matplot(IC_plat$mu_pred)

sort(IC_plat$mu_pred[,2,3])[1:25]

### En utilisant la méthode dans PLAT avec des ordres de ARIMA appropriés, j'arrive à reproduire approximativement les r
# résultats de PLAT, mais cette méthode appliquée à notre modèle fait n'importe quoi dans les prévisions (en gros les prévisions
#sont très volatiles tout âge confondu.

### En reprenant la méthode boostrap (avec la technique utilisée dans PLAT c'est-à-dire le choix de 2,5% des données 
# en sus et en sous), on tombe sur des résultats un peu mieux avec notre modèle, mais seulement que avec 
# cette approche le modèle risque d'être moins bon.

### Avec un auto.arima sur les kt on arrive à quelque chose mais à partir d'un certain temps, le
### taux de mortalité se redresse sur les intervalles de confiance.
### Par contre lorsqu'on essaie PLAT avec gamma simulé normalement sur toute la période de cohorte, on a de meilleur résultat
### pour l'IC.


### J'ai aussi essayé de voir le cas où je fite les paramètres temporelles avec des Arima comme dans PLAT, mais là également 
## on remarque une divergence de l'IC

##@ Une 


###########################################################################
### En utilisant le boostrap pour B petit (100) et en utilisant les quantiles
### d'ordre 1-a/2 pour calculer l'IC, on se retrouve finalement avec IC qui contient 
### nos prévisions



