#### Matrix d'information de Fisher

i_alpha=function(dxt_est,ages.fit,years.fit)
{
  "
  dxt_est represents the matrix of number of deaths estimated by the model
  "
  m=dim(dxt_est)[1]
  d=sapply(1:m, function(x){sum(dxt_est[x,])})
  diag(d)
}

i_kt1=function(dxt_est,ages.fit,years.fit)
{
  "
  dxt_est represents the matrix of number of deaths estimated by the model
  "
  n=dim(dxt_est)[2]
  d=sapply(1:n, function(t){sum(dxt_est[,t])})
  diag(d)
}


i_kt2=function(dxt_est,ages.fit,years.fit)
{
  "
  dxt_est represents the matrix of number of deaths estimated by the model
  "
  n=dim(dxt_est)[2] ### length of years
  x_bar=mean(ages.fit)
  a=matrix((x_bar-ages.fit)^2,nrow=1) ### the main idea is to recover sum_x{(bar{x}-x)*hat{dxt_est}
  d=a%*%dxt_est 
  diag(c(d))
}

i_kt3=function(dxt_est,ages.fit,years.fit)
{
  "
  dxt_est represents the matrix of number of deaths estimated by the model
  "
  n=dim(dxt_est)[2] ### length of years
  x_bar=mean(ages.fit)
  a=matrix((pmax(x_bar-ages.fit,0))^2,nrow=1) ### the main idea is to recover sum_x{(bar{x}-x)*hat{dxt_est}
  d=a%*%dxt_est 
  diag(c(d))
}  
  
  
i_kt4=function(dxt_est,ages.fit,years.fit,x1,x2,c)
{
  "
  dxt_est represents the matrix of number of deaths estimated by the model
  "
  n=dim(dxt_est)[2] ### length of years
  #x_bar=mean(ages)
  a=pmax(x1-ages.fit,0)+c*pmax(ages.fit-x2,0)
  a=matrix(a^2,nrow=1) ### the main idea is to recover sum_x{(bar{x}-x)*hat{dxt_est}
  d=a%*%dxt_est 
  diag(c(d))
}  

i_kt5=function(dxt_est,ages.fit,years.fit,a_c,Ic)
{
  "
  dxt_est represents the matrix of number of deaths estimated by the model
  Il faut penser à récupérer après juste les colonnes non entièrement pour faire 
  correspondre à kt5. Parce que ce sont ces colonnes qui permettent d'avoir le kt5
  "
  n=dim(dxt_est)[2] ### length of years
  #x_bar=mean(ages)
  
  ### reshape a and t to have (x-a_c)^{+}(Ic-Ic_bar)^+
  a=matrix(pmax(ages.fit-a_c,0),ncol=1) ## a in X*1 
  Ic_bar=mean(Ic)
  t=matrix(pmax(Ic-Ic_bar,0),nrow=1) ## t in 1*T
  
  b=a%*%t ## b in X*1*1*T
  d=(b^2)*dxt_est ### b^2*dxt_est pas un produit matriciel
  e=sapply(1:n,function(t){sum(d[,t])}) ### sum_x d
  kt5_with_0=diag(e) ### correspond au kt5 pour toute la période de temps considéré. Sauf 
  #que le kt5 concerne seulement les années où pmax(Ic-mean(Ic),0)!=0.
  
  kt5_without_0=c()
  for (i in 1:n)
  {
    if (sum(kt5_with_0[,i])!=0)
    {
      kt5_without_0=cbind(kt5_without_0,kt5_with_0[,i])
    }
  }
  kt5_without_0
  
}  



i_gamma=function(dxt_est, ages.fit, years.fit)
{
  " 
  L'idée c'est de sommer toutes les valeurs de dxt_est où t-x=c, c étant la cohorte
  puis retourner une matrice de taille (nbr_coh X nbr_coh) avec les sommes trouvées en diagonale.
  "
  n=length(ages.fit)
  k=length(years.fit)
  coh=seq(years.fit[1]-tail(ages.fit,1),tail(years.fit,1)-ages.fit[1])
  
  n.s=1950-coh[1]+1
  
  gamma=matrix(0,ncol=length(coh),nrow=1)
  for (t in 1:length(coh))
  {
    a=c()
    for(j in 1:k)
      for (i in 1:n)
        if(j-i==coh[t]) 
        {
          a=cbind(a,dxt_est[i,j])
        }
    gamma=sum(a)
  }
  
  
  diag(c(gamma[1:n.s]))
}


i_alpha_kt1=function(dxt_est,ages.fit,years.fit)
{
  dxt_est
}
  

i_alpha_kt2=function(dxt_est,ages.fit,years.fit)
{
  "
  t(sapply) in R^(a)x(t) a=dim(ages) et t=dim(years)
  "
  m=dim(dxt_est)[1]
  x_bar=mean(ages.fit)
  x=(x_bar-ages.fit)
  t(sapply(1:m,function(a){x[a]*dxt_est[a,]})) ### (x_bar-x)*dxt_est without sum 
}

 
i_alpha_kt3=function(dxt_est,ages.fit,years.fit)
{
  "
  t(sapply) in R^(a)x(t) a=dim(ages) et t=dim(years)
  "
  m=dim(dxt_est)[1]
  x_bar=mean(ages.fit)
  x=pmax(x_bar-ages.fit,0)
  t(sapply(1:m,function(a){x[a]*dxt_est[a,]})) ### (x_bar-x)*dxt_est without sum 
}


i_alpha_kt4=function(dxt_est,ages.fit,years.fit,x1,x2,c)
{
  "
  t(sapply) in R^(a)x(t) a=dim(ages) et t=dim(years)
  "
  m=dim(dxt_est)[1]
  #x_bar=mean(ages.fit)
  x=pmax(x1-ages.fit,0)+c*pmax(ages.fit-x2,0)
  t(sapply(1:m,function(a){x[a]*dxt_est[a,]}))  ### (x_bar-x)*dxt_est without sum 
}
  

i_alpha_kt5=function(dxt_est,ages.fit,years.fit,a_c,Ic)
{
  "
  dxt_est represents the matrix of number of deaths estimated by the model
  "
  n=dim(dxt_est)[2] ### length of years
  #x_bar=mean(ages)
  
  ### reshape a and t to have (x-a_c)^{+}(Ic-Ic_bar)^+
  a=matrix(pmax(ages.fit-a_c,0),ncol=1) ## a in X*1 
  Ic_bar=mean(Ic)
  t=matrix(pmax(Ic-Ic_bar,0),nrow=1) ## t in 1*T
  
  b=a%*%t ## b in X*1*1*T
  d_with_0=b*dxt_est ### b*dxt_est pas un produit matriciel
  
  d_without_0=c()
  for (i in 1:n)
  {
    if (sum(d_with_0[,i])!=0)
    {
      d_without_0=cbind(d_without_0,d_with_0[,i])
    }
  }
  d_without_0
  
}


i_alpha_gamma=function(dxt_est, ages.fit, years.fit)
{
  "
  L'idée ici c'est qu'on fixe l'age et la cohorte et on parcourt le temps pour voir le dxt_est correspond.
  "
  n=length(ages.fit)
  k=length(years.fit)
  coh=seq(years.fit[1]-tail(ages.fit,1),tail(years.fit,1)-ages.fit[1])
  
  n.s=1950-coh[1]+1
  alpha_gamma=matrix(nrow = n,ncol = length(coh),dimnames = list(ages.fit,coh))
  
  ### i pour l'âge, j pour l'année et t pour la cohorte
  for (i in 1:n)
    for(t in 1:length(coh))
    {
      a=c()
      for(j in 1:k)
        if(j-i==coh[t])
          a=cbind(a,dxt_est[i,j])
      alpha_gamma[i,t]=sum(a)
    }
  alpha_gamma[,1:n.s]
      
}




i_kt1_kt2=function(dxt_est,ages.fit,years.fit)
{
  n=dim(dxt_est)[2] ### length of years
  x_bar=mean(ages.fit)
  a=matrix((x_bar-ages.fit),nrow=1) ### the main idea is to recover sum_x{(bar{x}-x)*hat{dxt_est}
  d=a%*%dxt_est 
  diag(c(d))
}


i_kt1_kt3=function(dxt_est,ages.fit,years.fit)
{
  n=dim(dxt_est)[2] ### length of years
  x_bar=mean(ages.fit)
  a=matrix(pmax(x_bar-ages.fit,0),nrow=1) ### the main idea is to recover sum_x{(bar{x}-x)*hat{dxt_est}
  d=a%*%dxt_est 
  diag(c(d))
}


i_kt1_kt4=function(dxt_est,ages.fit,years.fit,x1,x2,c)
{
  n=dim(dxt_est)[2] ### length of years
  #x_bar=mean(ages)
  a=pmax(x1-ages.fit,0)+c*pmax(ages.fit-x2,0)
  a=matrix(a,nrow=1) ### the main idea is to recover sum_x{(bar{x}-x)*hat{dxt_est}
  d=a%*%dxt_est 
  diag(c(d))
}

i_kt1_kt5=function(dxt_est,ages.fit,years.fit,a_c,Ic)
{
  "
  dxt_est represents the matrix of number of deaths estimated by the model
  "
  n=dim(dxt_est)[2] ### length of years
  #x_bar=mean(ages)
  
  ### reshape a and t to have (x-a_c)^{+}(Ic-Ic_bar)^+
  a=matrix(pmax(ages.fit-a_c,0),ncol=1) ## a in X*1 
  Ic_bar=mean(Ic)
  t=matrix(pmax(Ic-Ic_bar,0),nrow=1) ## t in 1*T
  
  b=a%*%t ## b in X*1*1*T
  d=b*dxt_est ### b^2*dxt_est pas un produit matriciel
  e=sapply(1:n,function(t){sum(d[,t])}) ### sum_x d
  kt1_kt5_with0=diag(e)
  
  kt1_kt5_without0=c()
  for (i in 1:n)
  {
    if (sum(kt1_kt5_with0[,i])!=0)
    {
      kt1_kt5_without0=cbind(kt1_kt5_without0,kt1_kt5_with0[,i])
    }
  }
  kt1_kt5_without0
  
}


i_kt1_gamma=function(dxt_est,ages.fit,years.fit)
{
  "
  C'est la même idée que alpha_gamma sauf que ici c'est l'age qui varie
  "
  n=length(ages.fit)
  k=length(years.fit)
  coh=seq(years.fit[1]-tail(ages.fit,1),tail(years.fit,1)-ages.fit[1])
  n.s=1950-coh[1]+1
  
  kt1_gamma=matrix(nrow = k,ncol = length(coh),dimnames = list(years.fit,coh))
  
  ### i pour l'âge, j pour l'année et t pour la cohorte
  for (j in 1:k)
    for(t in 1:length(coh))
    {
      a=c()
      for(i in 1:n)
        if(j-i==coh[t])
          a=cbind(a,dxt_est[i,j])
      kt1_gamma[j,t]=sum(a)
    }
  kt1_gamma[,1:n.s]
}


i_kt2_kt3=function(dxt_est,ages.fit,years.fit)
{
  n=dim(dxt_est)[2] ### length of years
  x_bar=mean(ages.fit)
  a=matrix((pmax(x_bar-ages.fit,0)*(x_bar-ages.fit)),nrow=1) ### the main idea is to recover sum_x{(bar{x}-x)*hat{dxt_est}
  d=a%*%dxt_est 
  diag(c(d))
}


i_kt2_kt4=function(dxt_est,ages.fit,years.fit,x1,x2,c)
{
  n=dim(dxt_est)[2] ### length of years
  x_bar=mean(ages.fit)
  a=pmax(x1-ages.fit,0)+c*pmax(ages.fit-x2,0)
  a=matrix((x_bar-ages.fit)*a,nrow=1) ### the main idea is to recover sum_x{(bar{x}-x)*hat{dxt_est}
  d=a%*%dxt_est 
  diag(c(d))
}

i_kt2_kt5=function(dxt_est,ages.fit,years.fit,a_c,Ic)
{
  "
  dxt_est represents the matrix of number of deaths estimated by the model
  "
  n=dim(dxt_est)[2] ### length of years
  #x_bar=mean(ages)
  
  ### reshape a and t to have (x-a_c)^{+}(Ic-Ic_bar)^+
  a=matrix(pmax(ages.fit-a_c,0),ncol=1) ## a in X*1 
  Ic_bar=mean(Ic)
  t=matrix(pmax(Ic-Ic_bar,0),nrow=1) ## t in 1*T
  
  b=a%*%t ## b in X*1*1*T
  d=b*dxt_est ### b^2*dxt_est pas un produit matriciel
  a=matrix((x_bar-ages.fit),nrow=1) ### the main idea is to recover sum_x{(bar{x}-x)*hat{b*dxt_est}
  kt2_kt5_with0=diag(c(a%*%d))
  
  kt2_kt5_without0=c()
  for (i in 1:n)
  {
    if (sum(kt2_kt5_with0[,i])!=0)
    {
      kt2_kt5_without0=cbind(kt2_kt5_without0,kt2_kt5_with0[,i])
    }
  }
  kt2_kt5_without0
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
  n.s=1950-coh[1]+1
  
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
  
  kt2_gamma[,1:n.s]
  
}

i_kt3_kt4=function(dxt_est,ages.fit,years.fit,x1,x2,c)
{
  n=dim(dxt_est)[2] ### length of years
  x_bar=mean(ages.fit)
  a=pmax(x1-ages.fit,0)+c*pmax(ages.fit-x2,0)
  a=matrix(pmax(x_bar-ages.fit,0)*a,nrow=1) ### the main idea is to recover sum_x{(bar{x}-x)*hat{dxt_est}
  d=a%*%dxt_est 
  diag(c(d))
}


i_kt3_kt5=function(dxt_est,ages.fit,years.fit,a_c,Ic)
{
  "
  dxt_est represents the matrix of number of deaths estimated by the model
  "
  n=dim(dxt_est)[2] ### length of years
  x_bar=mean(ages.fit)
  
  ### reshape a and t to have (x-a_c)^{+}(Ic-Ic_bar)^+
  a=matrix(pmax(ages.fit-a_c,0),ncol=1) ## a in X*1 
  Ic_bar=mean(Ic)
  t=matrix(pmax(Ic-Ic_bar,0),nrow=1) ## t in 1*T
  
  b=a%*%t ## b in X*1*1*T
  d=b*dxt_est ### b*dxt_est pas un produit matriciel
  
  a=matrix(pmax(x_bar-ages.fit,0),nrow=1) ### the main idea is to recover sum_x{(bar{x}-x)*hat{dxt_est}
  diag(c(a%*%d))
  #a%*%d
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
  n.s=1950-coh[1]+1
  
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
  kt3_gamma[,1:n.s]
}


i_kt4_kt5=function(dxt_est,ages.fit,years.fit,x1,x2,c,a_c,Ic)
{
  "
  dxt_est represents the matrix of number of deaths estimated by the model
  "
  n=length(ages.fit)
  k=length(years.fit)
  coh=seq(years.fit[1]-tail(ages.fit,1),tail(years.fit,1)-ages.fit[1])
  ### length of years
  x_bar=mean(ages.fit)
  
  ### reshape a and t to have (x-a_c)^{+}(Ic-Ic_bar)^+
  a=matrix(pmax(ages.fit-a_c,0),ncol=1) ## a in X*1 
  Ic_bar=mean(Ic)
  t=matrix(pmax(Ic-Ic_bar,0),nrow=1) ## t in 1*T
  
  coef2=a%*%t ## b in X*1*1*T
  d=coef2*dxt_est ### b*dxt_est pas un produit matriciel
  
  coef1=matrix((pmax(x1-ages.fit,0)+c*pmax(ages.fit-x2,0)),nrow=1) ### the main idea is to recover sum_x{(bar{x}-x)*hat{dxt_est}
  kt4_kt5_with0=diag(c(coef1%*%d))
  
  kt4_kt5_without0=c()
  for (i in 1:n)
  {
    if (sum(kt4_kt5_with0[,i])!=0)
    {
      kt4_kt5_without0=cbind(kt4_kt5_without0,kt4_kt5_with0[,i])
    }
  }
  kt4_kt5_without0
  
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
  n.s=1950-coh[1]+1
  
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
  kt4_gamma[,1:n.s]
}

i_kt5_gamma=function(dxt_est,ages.fit,years.fit,a_c,Ic)
{
  "
  C'est la même idée que alpha_gamma sauf que ici c'est l'age qui varie et on 
  met le coefficient de kt5 devant chaque fois.
  "
  n=length(ages.fit)
  k=length(years.fit)
  coh=seq(years.fit[1]-tail(ages.fit,1),tail(years.fit,1)-ages.fit[1])
  n.s=1950-coh[1]+1
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
  kt5_gamma[,1:n.s]
}



information_plat=function(dxt_est,ages.fit,years.fit)
{
  ### First line of the matrix
  alpha=i_alpha(dxt_est,ages.fit,years.fit)
  alpha_kt1=i_alpha_kt1(dxt_est,ages.fit,years.fit)
  alpha_kt2=i_alpha_kt2(dxt_est,ages.fit,years.fit)
  alpha_kt3=i_alpha_kt3(dxt_est,ages.fit,years.fit)
  alpha_gamma=i_alpha_gamma(dxt_est,ages.fit,years.fit)
  
  l1=cbind(alpha,alpha_kt1,alpha_kt2,alpha_kt3,alpha_gamma)
  
  ### Second line 
  kt1_alpha=t(alpha_kt1)
  kt1=i_kt1(dxt_est,ages.fit,years.fit)
  kt1_kt2=i_kt1_kt2(dxt_est,ages.fit,years.fit)
  kt1_kt3=i_kt1_kt3(dxt_est,ages.fit,years.fit)
  kt1_gamma=i_kt1_gamma(dxt_est,ages.fit,years.fit)
  
  l2=cbind(kt1_alpha,kt1,kt1_kt2,kt1_kt3,kt1_gamma)
  
  ### Third line
  kt2_alpha=t(alpha_kt2)
  kt2_kt1=t(kt1_kt2)
  kt2=i_kt2(dxt_est,ages.fit,years.fit)
  kt2_kt3=i_kt2_kt3(dxt_est,ages.fit,years.fit)
  kt2_gamma=i_kt2_gamma(dxt_est,ages.fit,years.fit)
  
  l3=cbind(kt2_alpha,kt2_kt1,kt2,kt2_kt3,kt2_gamma)
  
  ### Fourth line
  kt3_alpha=t(alpha_kt3)
  kt3_kt1=t(kt1_kt3)
  kt3_kt2=t(kt2_kt3)
  kt3=i_kt3(dxt_est,ages.fit,years.fit)
  kt3_gamma=i_kt3_gamma(dxt_est,ages.fit,years.fit)
  
  l4=cbind(kt3_alpha,kt3_kt1,kt3_kt2,kt3,kt3_gamma)
  
  ### Fifth line
  gamma_alpha=t(alpha_gamma)
  gamma_kt1=t(kt1_gamma)
  gamma_kt2=t(kt2_gamma)
  gamma_kt3=t(kt3_gamma)
  gamma=i_gamma(dxt_est,ages.fit,years.fit)
  
  l5=cbind(gamma_alpha,gamma_kt1,gamma_kt2,gamma_kt3,gamma)
  
  inf_matrix_plat=rbind(l1,l2,l3,l4,l5)
  
  
  inf_matrix_plat
}


information_model=function(dxt_est,ages.fit,years.fit,x1,x2,a_c,c,Ic)
{
  ### First line of the matrix
  alpha=i_alpha(dxt_est,ages.fit,years.fit)
  alpha_kt1=i_alpha_kt1(dxt_est,ages.fit,years.fit)
  alpha_kt2=i_alpha_kt2(dxt_est,ages.fit,years.fit)
  alpha_kt3=i_alpha_kt3(dxt_est,ages.fit,years.fit)
  alpha_kt4=i_alpha_kt4(dxt_est,ages.fit,years.fit,x1,x2,c)
  alpha_kt5=i_alpha_kt5(dxt_est,ages.fit,years.fit,a_c,Ic)
  alpha_gamma=i_alpha_gamma(dxt_est,ages.fit,years.fit)
  
  l1=cbind(alpha,alpha_kt1,alpha_kt2,alpha_kt3,alpha_kt4,alpha_kt5,alpha_gamma)
  
  ### Second line 
  kt1_alpha=t(alpha_kt1)
  kt1=i_kt1(dxt_est,ages.fit,years.fit)
  kt1_kt2=i_kt1_kt2(dxt_est,ages.fit,years.fit)
  kt1_kt3=i_kt1_kt3(dxt_est,ages.fit,years.fit)
  kt1_kt4=i_kt1_kt4(dxt_est,ages.fit,years.fit,x1,x2,c)
  kt1_kt5=i_kt1_kt5(dxt_est,ages.fit,years.fit,a_c,Ic)
  kt1_gamma=i_kt1_gamma(dxt_est,ages.fit,years.fit)
  
  l2=cbind(kt1_alpha,kt1,kt1_kt2,kt1_kt3,kt1_kt4,kt1_kt5,kt1_gamma)
  
  ### Third line
  kt2_alpha=t(alpha_kt2)
  kt2_kt1=t(kt1_kt2)
  kt2=i_kt2(dxt_est,ages.fit,years.fit)
  kt2_kt3=i_kt2_kt3(dxt_est,ages.fit,years.fit)
  kt2_kt4=i_kt2_kt4(dxt_est,ages.fit,years.fit,x1,x2,c)
  kt2_kt5=i_kt2_kt5(dxt_est,ages.fit,years.fit,a_c,Ic)
  kt2_gamma=i_kt2_gamma(dxt_est,ages.fit,years.fit)
  
  l3=cbind(kt2_alpha,kt2_kt1,kt2,kt2_kt3,kt2_kt4,kt2_kt5,kt2_gamma)
  
  ### Fourth line
  kt3_alpha=t(alpha_kt3)
  kt3_kt1=t(kt1_kt3)
  kt3_kt2=t(kt2_kt3)
  kt3=i_kt3(dxt_est,ages.fit,years.fit)
  kt3_kt4=i_kt3_kt4(dxt_est,ages.fit,years.fit,x1,x2,c)
  kt3_kt5=i_kt3_kt5(dxt_est,ages.fit,years.fit,a_c,Ic)
  kt3_gamma=i_kt3_gamma(dxt_est,ages.fit,years.fit)
  
  l4=cbind(kt3_alpha,kt3_kt1,kt3_kt2,kt3,kt3_kt4,kt3_kt5,kt3_gamma)
  
  ###Fifth line
  kt4_alpha=t(alpha_kt4)
  kt4_kt1=t(kt1_kt4)
  kt4_kt2=t(kt2_kt4)
  kt4_kt3=t(kt3_kt4)
  kt4=i_kt4(dxt_est,ages.fit,years.fit,x1,x2,c)
  kt4_kt5=i_kt4_kt5(dxt_est,ages.fit,years.fit,x1,x2,c,a_c,Ic)
  kt4_gamma=i_kt4_gamma(dxt_est,ages.fit,years.fit,x1,x2,c)
  
  l5=cbind(kt4_alpha,kt4_kt1,kt4_kt2,kt4_kt3,kt4,kt4_kt5,kt4_gamma)
  
  ### Sixth line
  kt5_alpha=t(alpha_kt5)
  kt5_kt1=t(kt1_kt5)
  kt5_kt2=t(kt2_kt5)
  kt5_kt3=t(kt3_kt5)
  kt5_kt4=t(kt4_kt5)
  kt5=i_kt5(dxt_est,ages.fit,years.fit,a_c,Ic)
  kt5_gamma=i_kt5_gamma(dxt_est,ages.fit,years.fit,a_c,Ic)
  
  l6=cbind(kt5_alpha,kt5_kt1,kt5_kt2,kt5_kt3,kt5_kt4,kt5,kt5_gamma)
  
  ### seventh
  gamma_alpha=t(alpha_gamma)
  gamma_kt1=t(kt1_gamma)
  gamma_kt2=t(kt2_gamma)
  gamma_kt3=t(kt3_gamma)
  gamma_kt4=t(kt4_gamma)
  gamma_kt5=t(kt5_gamma)
  gamma=i_gamma(dxt_est,ages.fit,years.fit)
  
  l7=cbind(gamma_alpha,gamma_kt1,gamma_kt2,gamma_kt3,gamma_kt4,gamma_kt5,gamma)
  
  inf_matrix_model=rbind(l1,l2,l3,l4,l5,l6,l7)
  
  
  inf_matrix_model
}



d=matrix(1:35,nrow = 5)
ages=1:5
years=1:7

i_kt1(dxt_est = d,ages.fit = ages,years.fit = years)

