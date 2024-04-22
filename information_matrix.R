#### Matrix d'information de Fisher

i_alpha=function(dxt_est,w_xt,ages.fit,years.fit)
{
  """dxt_est represents the matrix of number of deaths estimated by the model"""
  m=dim(dxt_est)[1]
  d=sapply(1:m, function(x){sum(d[x,])})
  diag(d)
}

i_kt1=function(dxt_est,w_xt,ages.fit,years.fit)
{
  """dxt_est represents the matrix of number of deaths estimated by the model"""
  n=dim(dxt_est)[2]
  d=sapply(1:n, function(t){sum(d[,t])})
  diag(d)
}


i_kt2=function(dxt_est,w_xt,ages.fit,years.fit)
{
  """dxt_est represents the matrix of number of deaths estimated by the model"""
  n=dim(dxt_est)[2] ### length of years
  x_bar=mean(ages)
  a=matrix((x_bar-ages)^2,nrow=1) ### the main idea is to recover sum_x{(bar{x}-x)*hat{dxt_est}
  d=a%*%dxt_est 
  diag(d)
}

i_kt3=function(dxt_est,w_xt,ages.fit,years.fit)
{
  """dxt_est represents the matrix of number of deaths estimated by the model"""
  n=dim(dxt_est)[2] ### length of years
  x_bar=mean(ages)
  a=matrix((pmax(x_bar-ages,0))^2,nrow=1) ### the main idea is to recover sum_x{(bar{x}-x)*hat{dxt_est}
  d=a%*%dxt_est 
  diag(d)
}  
  
  
i_kt4=function(dxt_est,w_xt,ages.fit,years.fit,x1,x2,c)
{
  """dxt_est represents the matrix of number of deaths estimated by the model"""
  n=dim(dxt_est)[2] ### length of years
  #x_bar=mean(ages)
  a=pmax(x1-ages,0)+c*pmax(ages-x2,0)
  a=matrix(a^2,nrow=1) ### the main idea is to recover sum_x{(bar{x}-x)*hat{dxt_est}
  d=a%*%dxt_est 
  diag(d)
}  

i_kt5=function(dxt_est,w_xt,ages.fit,years.fit,a_c,Ic)
{
  """dxt_est represents the matrix of number of deaths estimated by the model"""
  n=dim(dxt_est)[2] ### length of years
  #x_bar=mean(ages)
  
  ### reshape a and t to have (x-a_c)^{+}(Ic-Ic_bar)^+
  a=matrix(pmax(ages-a_c,0),ncol=1) ## a in X*1 
  Ic_bar=mean(Ic)
  t=matrix(pmax(Ic-Ic_bar,0),nrow=1) ## t in 1*T
  
  b=a%*%t ## b in X*1*1*T
  d=(b^2)*dxt_est ### b^2*dxt_est pas un produit matriciel
  e=sapply(1:n,function(t){sum(d[,t])}) ### sum_x d
  diag(e)
}  



i_gamma=function(dxt_est, ages.fit, years.fit)
{
  " 
  L'idée c'est de sommer toutes les valeurs de dxt_est où t-x=c, c étant la cohorte
  puis retourner une matrice de taille (nbr_coh X nbr_coh) avec les sommes trouvées en diagonale.
  "
  n=length(ages.fit)
  k=length(years.fit)
  c=seq(1-n,k-1)
  
  sum_d=matrix(0,ncol=length(c),nrow=1)
  for (t in 1:length(c))
  {
    a=c()
    for(j in 1:k)
      for (i in 1:n)
        if(j-i==c[t]) 
        {
          a=cbind(a,dxt_est[i,j])
        }
    sum_d[t]=sum(a)
  }
  
  
  diag(sum_d)
}


i_alpha_kt1=function(dxt_est,w_xt,ages.fit,years.fit)
{
  dxt_est
}
  

i_alpha_kt2=function(dxt_est,w_xt,ages.fit,years.fit)
{
  m=dim(dxt_est)[1]
  x_bar=mean(ages.fit)
  x=(x_bar-ages)
  sapply(1:m,function(a){x[a]*dxt_est[a,]}) ### (x_bar-x)*dxt_est without sum 
}

 
i_alpha_kt3=function(dxt_est,w_xt,ages.fit,years.fit)
{
  m=dim(dxt_est)[1]
  x_bar=mean(ages.fit)
  x=pmax(x_bar-ages,0)
  sapply(1:m,function(a){x[a]*dxt_est[a,]}) ### (x_bar-x)*dxt_est without sum 
}


i_alpha_kt4=function(dxt_est,w_xt,ages.fit,years.fit,x1,x2,c)
{
  m=dim(dxt_est)[1]
  #x_bar=mean(ages.fit)
  x=pmax(x1-ages,0)+c*pmax(ages-x2,0)
  sapply(1:m,function(a){x[a]*dxt_est[a,]}) ### (x_bar-x)*dxt_est without sum 
}
  

i_alpha_kt5=function(dxt_est,w_xt,ages.fit,years.fit,a_c,Ic)
{
  """dxt_est represents the matrix of number of deaths estimated by the model"""
  n=dim(dxt_est)[2] ### length of years
  #x_bar=mean(ages)
  
  ### reshape a and t to have (x-a_c)^{+}(Ic-Ic_bar)^+
  a=matrix(pmax(ages-a_c,0),ncol=1) ## a in X*1 
  Ic_bar=mean(Ic)
  t=matrix(pmax(Ic-Ic_bar,0),nrow=1) ## t in 1*T
  
  b=a%*%t ## b in X*1*1*T
  d=b*dxt_est ### b*dxt_est pas un produit matriciel
  d
}


i_alpha_gamma=function(dxt_est, ages.fit, years.fit)
{
  "
  L'idée ici c'est qu'on fixe l'age et la cohorte et on parcourt le temps pour voir le dxt_est correspond.
  "
  n=length(ages.fit)
  k=length(years.fit)
  c=seq(years.fit[1]-tail(ages.fit,1),tail(years.fit,1)-ages.fit[1])
  
  alpha_gamma=matrix(nrow = n,ncol = length(c),dimnames = list(ages.fit,c))
  
  ### i pour l'âge, j pour l'année et t pour la cohorte
  for (i in 1:n)
    for(t in 1:length(c))
    {
      a=c()
      for(j in 1:k)
        if(j-i==c[t])
          a=cbind(a,dxt_est[i,j])
      alpha_gamma[i,t]=sum(a)
    }
  alpha_gamma
      
}




i_kt1_kt2=function(dxt_est,w_xt,ages.fit,years.fit)
{
  n=dim(dxt_est)[2] ### length of years
  x_bar=mean(ages)
  a=matrix((x_bar-ages),nrow=1) ### the main idea is to recover sum_x{(bar{x}-x)*hat{dxt_est}
  d=a%*%dxt_est 
  diag(d)
}


i_kt1_kt3=function(dxt_est,w_xt,ages.fit,years.fit)
{
  n=dim(dxt_est)[2] ### length of years
  x_bar=mean(ages)
  a=matrix(pmax(x_bar-ages,0),nrow=1) ### the main idea is to recover sum_x{(bar{x}-x)*hat{dxt_est}
  d=a%*%dxt_est 
  diag(d)
}


i_kt1_kt4=function(dxt_est,w_xt,ages.fit,years.fit,x1,x2,c)
{
  n=dim(dxt_est)[2] ### length of years
  #x_bar=mean(ages)
  a=pmax(x1-ages,0)+c*pmax(ages-x2,0)
  a=matrix(a,nrow=1) ### the main idea is to recover sum_x{(bar{x}-x)*hat{dxt_est}
  d=a%*%dxt_est 
  diag(d)
}

i_kt1_kt5=function(dxt_est,w_xt,ages.fit,years.fit,a_c,Ic)
{
  "
  dxt_est represents the matrix of number of deaths estimated by the model
  "
  n=dim(dxt_est)[2] ### length of years
  #x_bar=mean(ages)
  
  ### reshape a and t to have (x-a_c)^{+}(Ic-Ic_bar)^+
  a=matrix(pmax(ages-a_c,0),ncol=1) ## a in X*1 
  Ic_bar=mean(Ic)
  t=matrix(pmax(Ic-Ic_bar,0),nrow=1) ## t in 1*T
  
  b=a%*%t ## b in X*1*1*T
  d=b*dxt_est ### b^2*dxt_est pas un produit matriciel
  e=sapply(1:n,function(t){sum(d[,t])}) ### sum_x d
  diag(e)
}

i_kt1_gamma=function(dxt_est,ages.fit,years.fit)
{
  "
  C'est la même idée que alpha_gamma sauf que ici c'est l'age qui varie
  "
  n=length(ages.fit)
  k=length(years.fit)
  c=seq(years.fit[1]-tail(ages.fit,1),tail(years.fit,1)-ages.fit[1])
  
  kt1_gamma=matrix(nrow = k,ncol = length(c),dimnames = list(years.fit,c))
  
  ### i pour l'âge, j pour l'année et t pour la cohorte
  for (j in 1:k)
    for(t in 1:length(c))
    {
      a=c()
      for(i in 1:n)
        if(j-i==c[t])
          a=cbind(a,dxt_est[i,j])
      kt1_gamma[j,t]=sum(a)
    }
  kt1_gamma
}


i_kt2_kt3=function(dxt_est,w_xt,ages.fit,years.fit)
{
  n=dim(dxt_est)[2] ### length of years
  x_bar=mean(ages)
  a=matrix((pmax(x_bar-ages,0)*(x_bar-ages)),nrow=1) ### the main idea is to recover sum_x{(bar{x}-x)*hat{dxt_est}
  d=a%*%dxt_est 
  diag(d)
}


i_kt2_kt4=function(dxt_est,w_xt,ages.fit,years.fit,x1,x2,c)
{
  n=dim(dxt_est)[2] ### length of years
  x_bar=mean(ages)
  a=pmax(x1-ages,0)+c*pmax(ages-x2,0)
  a=matrix((x_bar-ages)*a,nrow=1) ### the main idea is to recover sum_x{(bar{x}-x)*hat{dxt_est}
  d=a%*%dxt_est 
  diag(d)
}

i_kt2_kt5=function()
{
  "
  dxt_est represents the matrix of number of deaths estimated by the model
  "
  n=dim(dxt_est)[2] ### length of years
  #x_bar=mean(ages)
  
  ### reshape a and t to have (x-a_c)^{+}(Ic-Ic_bar)^+
  a=matrix(pmax(ages-a_c,0),ncol=1) ## a in X*1 
  Ic_bar=mean(Ic)
  t=matrix(pmax(Ic-Ic_bar,0),nrow=1) ## t in 1*T
  
  b=a%*%t ## b in X*1*1*T
  d=b*dxt_est ### b^2*dxt_est pas un produit matriciel
  a=matrix((x_bar-ages),nrow=1) ### the main idea is to recover sum_x{(bar{x}-x)*hat{b*dxt_est}
  diag(a%*%d)
}



i_kt2_gamma=function(dxt_est,ages.fit,years.fit)
{
  "
  C'est la même idée que alpha_gamma sauf que ici c'est l'age qui varie et on 
  met le coefficient (x_bar-x) devant chaque fois.
  "
  n=length(ages.fit)
  k=length(years.fit)
  coh=seq(years.fit[1]-tail(ages.fit,1),tail(years.fit,1)-ages.fit[1])
  
  x=ages.fit
  x_bar=mean(ages.fit)
  kt2_gamma=matrix(nrow = k,ncol = length(coh),dimnames = list(years.fit,coh))
  
  ### i pour l'âge, j pour l'année et t pour la cohorte
  for (j in 1:k)
    for(t in 1:length(coh))
    {
      a=c()
      for(i in 1:n)
        if(j-i==coh[t])
          a=cbind(a,((x_bar-x[i])*dxt_est[i,j]))
      kt2_gamma[j,t]=sum(a)
    }
  kt2_gamma
}

i_kt3_kt4=function()
{
  n=dim(dxt_est)[2] ### length of years
  x_bar=mean(ages)
  a=pmax(x1-ages,0)+c*pmax(ages-x2,0)
  a=matrix(pmax(x_bar-ages,0)*a,nrow=1) ### the main idea is to recover sum_x{(bar{x}-x)*hat{dxt_est}
  d=a%*%dxt_est 
  diag(d)
}


i_kt3_kt5=function()
{
  "
  dxt_est represents the matrix of number of deaths estimated by the model
  "
  n=dim(dxt_est)[2] ### length of years
  x_bar=mean(ages)
  
  ### reshape a and t to have (x-a_c)^{+}(Ic-Ic_bar)^+
  a=matrix(pmax(ages-a_c,0),ncol=1) ## a in X*1 
  Ic_bar=mean(Ic)
  t=matrix(pmax(Ic-Ic_bar,0),nrow=1) ## t in 1*T
  
  b=a%*%t ## b in X*1*1*T
  d=b*dxt_est ### b*dxt_est pas un produit matriciel
  
  a=matrix(pmax(x_bar-ages,0),nrow=1) ### the main idea is to recover sum_x{(bar{x}-x)*hat{dxt_est}
  diag(a%*%d)
}


i_kt3_gamma=function(dxt_est,ages.fit,years.fit)
{
  "
  C'est la même idée que alpha_gamma sauf que ici c'est l'age qui varie et on 
  met le coefficient (x_bar-x)^+ devant chaque fois.
  "
  n=length(ages.fit)
  k=length(years.fit)
  coh=seq(years.fit[1]-tail(ages.fit,1),tail(years.fit,1)-ages.fit[1])
  
  x=ages.fit
  x_bar=mean(ages.fit)
  kt3_gamma=matrix(nrow = k,ncol = length(coh),dimnames = list(years.fit,coh))
  
  ### i pour l'âge, j pour l'année et t pour la cohorte
  for (j in 1:k)
    for(t in 1:length(coh))
    {
      a=c()
      for(i in 1:n)
        if(j-i==coh[t])
          a=cbind(a,(pmax(x_bar-x[i],0)*dxt_est[i,j]))
      kt3_gamma[j,t]=sum(a)
    }
  kt3_gamma
}


i_kt4_kt5=function(dxt_est,w_xt,ages.fit,years.fit,x1,x2,c,a_c,Ic)
{
  "
  dxt_est represents the matrix of number of deaths estimated by the model
  "
  n=dim(dxt_est)[2] ### length of years
  x_bar=mean(ages)
  
  ### reshape a and t to have (x-a_c)^{+}(Ic-Ic_bar)^+
  a=matrix(pmax(ages-a_c,0),ncol=1) ## a in X*1 
  Ic_bar=mean(Ic)
  t=matrix(pmax(Ic-Ic_bar,0),nrow=1) ## t in 1*T
  
  coef2=a%*%t ## b in X*1*1*T
  d=b*dxt_est ### b*dxt_est pas un produit matriciel
  
  coef1=matrix((pmax(x1-ages,0)+c*pmax(ages-x2,0)),nrow=1) ### the main idea is to recover sum_x{(bar{x}-x)*hat{dxt_est}
  diag(a%*%d)
}

i_kt4_gamma=function(dxt_est,ages.fit,years.fit,x1,x2,c)
{
  "
  C'est la même idée que alpha_gamma sauf que ici c'est l'age qui varie et on 
  met le coefficient (x1-x)^+c*(x-x2)^+ devant chaque fois.
  "
  n=length(ages.fit)
  k=length(years.fit)
  coh=seq(years.fit[1]-tail(ages.fit,1),tail(years.fit,1)-ages.fit[1])
  
  x=ages.fit
  coef=(pmax(x1-x,0)+c*pmax(x-x2,0))^2
  kt4_gamma=matrix(nrow = k,ncol = length(coh),dimnames = list(years.fit,coh))
  
  ### i pour l'âge, j pour l'année et t pour la cohorte
  for (j in 1:k)
    for(t in 1:length(coh))
    {
      a=c()
      for(i in 1:n)
        if(j-i==coh[t])
          a=cbind(a,(coef[i]*dxt_est[i,j]))
      kt4_gamma[j,t]=sum(a)
    }
  kt4_gamma
}

i_kt5_gamma=function(dxt_est,w_xt,ages.fit,years.fit,a_c,Ic)
{
  "
  C'est la même idée que alpha_gamma sauf que ici c'est l'age qui varie et on 
  met le coefficient de kt5 devant chaque fois.
  "
  n=dim(dxt_est)[2] ### length of years
  x_bar=mean(ages)
  
  ### reshape a and t to have (x-a_c)^{+}(Ic-Ic_bar)^+
  a=matrix(pmax(ages-a_c,0),ncol=1) ## a in X*1 
  Ic_bar=mean(Ic)
  t=matrix(pmax(Ic-Ic_bar,0),nrow=1) ## t in 1*T
  
  coef=a%*%t ## b in X*1*1*T
  
  kt5_gamma=matrix(nrow = k,ncol = length(coh),dimnames = list(years.fit,coh))
  
  ### i pour l'âge, j pour l'année et t pour la cohorte
  for (j in 1:k)
    for(t in 1:length(coh))
    {
      a=c()
      for(i in 1:n)
        if(j-i==coh[t])
          a=cbind(a,(coef[i]*dxt_est[i,j]))
      kt5_gamma[j,t]=sum(a)
    }
  kt5_gamma
}




=matrix(seq(1,10),nrow=5)
d
d=sapply(1:2, function(x){sum(d[,x])})
