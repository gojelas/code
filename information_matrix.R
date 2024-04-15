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



i_gamma=function()
{
  
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


i_alpha_gamma=function()
{
  
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
  """dxt_est represents the matrix of number of deaths estimated by the model"""
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

i_kt1_gamma=function(dxt_est,w_xt,ages.fit,years.fit,a_c,Ic)
{
  
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
  """dxt_est represents the matrix of number of deaths estimated by the model"""
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



i_kt2_gamma=function()
{
  
}

i_kt3_kt4=function()
{
  n=dim(dxt_est)[2] ### length of years
  x_bar=mean(ages)
  a=pmax(x1-ages,0)+c*pmax(ages-x2,0)
  a=matrix((pmax(x_bar-ages,0)*a,nrow=1) ### the main idea is to recover sum_x{(bar{x}-x)*hat{dxt_est}
  d=a%*%dxt_est 
  diag(d)
}


i_kt3_kt5=function()
{
  """dxt_est represents the matrix of number of deaths estimated by the model"""
  n=dim(dxt_est)[2] ### length of years
  x_bar=mean(ages)
  
  ### reshape a and t to have (x-a_c)^{+}(Ic-Ic_bar)^+
  a=matrix(pmax(ages-a_c,0),ncol=1) ## a in X*1 
  Ic_bar=mean(Ic)
  t=matrix(pmax(Ic-Ic_bar,0),nrow=1) ## t in 1*T
  
  b=a%*%t ## b in X*1*1*T
  d=b*dxt_est ### b*dxt_est pas un produit matriciel
  
  a=matrix((pmax(x_bar-ages,0),nrow=1) ### the main idea is to recover sum_x{(bar{x}-x)*hat{dxt_est}
  diag(a%*%d)
}


i_kt3_gamma=function()
{
  
}


i_kt4_kt5=function(dxt_est,w_xt,ages.fit,years.fit,x1,x2,c,a_c,Ic)
{
  """dxt_est represents the matrix of number of deaths estimated by the model"""
  n=dim(dxt_est)[2] ### length of years
  x_bar=mean(ages)
  
  ### reshape a and t to have (x-a_c)^{+}(Ic-Ic_bar)^+
  a=matrix(pmax(ages-a_c,0),ncol=1) ## a in X*1 
  Ic_bar=mean(Ic)
  t=matrix(pmax(Ic-Ic_bar,0),nrow=1) ## t in 1*T
  
  b=a%*%t ## b in X*1*1*T
  d=b*dxt_est ### b*dxt_est pas un produit matriciel
  
  a=matrix((pmax(x1-ages,0)+c*pmax(ages-x2,0)),nrow=1) ### the main idea is to recover sum_x{(bar{x}-x)*hat{dxt_est}
  diag(a%*%d)
}

i_kt4_gamma=function()
{
  
}

i_kt5_gamma=function()
{
  
}




=matrix(seq(1,10),nrow=5)
d
d=sapply(1:2, function(x){sum(d[,x])})
