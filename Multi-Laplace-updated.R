# Simulation Study Using Multivariate Laplace Distribution
# for the paper
# Robust Mixture of Linear Mixed Modeling via Multivariate Laplace Distribution


packages = c("MASS", "mvtnorm")
## Now load or install&load all
package.check = lapply(
  packages,
  FUN = function(x) 
  {
    if (!require(x, character.only = TRUE)) 
    {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  })

y.fun=function(I,r,rs,e,x,u,beta,b)
{
  ###########################################################################
  # Model: yi|(zi=j)=xi*betaj+ui*bij+e
  #     I: index number of i
  #     r: repeated measure in balanced case. 
  #    rs: repeated measure in unbalanced case. 
  #        Radomly sampled in 2:r. rs is a vector with I=100 terms. 
  #     e: the error term in the model above. IXr
  #     x: x matrix in model above. (pXr)XI
  #     u: u matrix in model above. (qXr)XI
  #  beta: fixed parameters in model above. (pXm)
  #     b: random parameters in model above. (IXq)
  #
  #  Usage: 
  #     z_yub=y.fun(I=I,r=r,rs=n,e=er,x=z_xub,u=z_uub,beta=beta12,b=b12)
  ###########################################################################
  
  z_yub=array(0,c(r,1,I))  #  empty 3 way array to contain unbalanced y.
  ee=matrix(0,I,r)      #  empty Ixr matrix to contain unbalanced errors.
  
  for (i in 1:I)
  {
    ee[i,1:rs[i]]=e[i,1:rs[i]]
    if (i<=rn)
    {
      j=1
      z_yub[,,i]=t(x[,,i])%*%beta[,j] + t(u[,,i])%*%b[i,] + ee[i,] 
      # 1st component of y
    } else
    { 
      j=2
      z_yub[,,i]=t(x[,,i])%*%beta[,j] + t(u[,,i])%*%b[i,] + ee[i,] 
      # 2nd component of y.
    } #end of 'else'
  }#end of ith loop
  z_yub
}#end of function


dmvl=function(x,mu,Sig)
{
  nd=length(x);
  Q=t(x-mu)%*%solve(Sig)%*%(x-mu);
  dL=besselK(sqrt(2*Q),nd/2-1)*(Q/2)^((1-nd/2)/2)*2/
    ((2*pi)^(nd/2)*sqrt(det(Sig)))
  return(dL)
}

##################################################################
# density of  function of yi's. It is 3 way array (1xmxI) 
# note: for each yi ,we don't know which component it comes from
##################################################################

deny.fun=function(y,I,m,r,beta,x,u,si,sigma)
{
  ##############################################################################
  # model: (yi|zi=j) = xi*betaj + ui*bij + eij
  #     y:  unbalanced y in model above.
  #     I:  index number of i. I=100
  #     m:  number of compomemts. m=2
  #     r:  repeated measure in balanced case. r=4.
  #    rn:  Radomly sampled index number from I=100 index by binom distribution. 
  #     x:  x matrix in model above
  #     u:  u matrix in model above
  #  beta:  fixed parameters in model above.
  #    si:  covariance of b
  # sigma:  variance of error term in model above. 
  #          it has 2 terms relatively to 2 components.
  # Usage: dy=deny.fun(y=z_yub,I=I,m=m,r=r,beta=beta12,x=z_xub,u=z_uub,
  #                    si=si.theta,sigma=sigma.theta)
  ##############################################################################
  dy=array(0,c(1,m,I))
  mu_m=array(0,c(r,m,I))
  Gamma=array(0,c(r,r,m,I))#4 way array to contain variance of each yi
  
  for (i in 1:I)
  {
    for (j in 1:m)
    {
      #mu_m[,j,i]=t(x[,,i])%*%beta[,j]
      mu_m[1:n[i],j,i]=t(x[,1:n[i],i])%*%beta[,j]
      #Gamma[,,j,i]=t(u[,,i])%*%si[,,j]%*%u[,,i]+diag(sigma[j])#variance of yi.
      Gamma[1:n[i],1:n[i],j,i]=t(u[,1:n[i],i])%*%si[,,j]%*%u[,1:n[i],i]+
        diag(sigma[j],n[i])
      # variance of yi.
      dy[,j,i]=dmvl(y[1:n[i],,i],mu_m[(1:n[i]),j,i],Gamma[(1:n[i]),(1:n[i]),j,i])
    }#end of jth function
  }#end of ith function
  #list(dy=dy,mu=mu_m,gamma=Gamma)
  dy=pmax(10^(-50),dy)
  dy=array(dy,c(1,m,I))
  dy
}#end of function

##############################################################################
#function of prop  in e-step. Prop is a Ixm(100x2) matrix.    
#
#That is the probability that each yi is in 1st comp and 2nd comp             #
##############################################################################

pr.fun=function(y,I,m,pr,dy)
{
  ########################################################################
  # definition in model: (yi|zi=j) = xi*betaj + ui*bij + eij
  # y:  unbalanced y in model above.
  # I:  index number of i. I=100
  # m:  number of compomemts. m=2
  # pr: propability of each i in jth component. i.e pij.
  # dy: densities of y.
  # Usage: pr.E=pr.fun(y=z_yub,I=I,m=m,pr=pr.theta,dy=dy )
  #######################################################################
  
  p=matrix(0,I,m)
  num=matrix(0,I,m);
  den=vector();
  
  for (j in 1:m){
      #find the proportion of each yi in jth component
      num[,j]=pr[j]*dy[,j,];
    }#end of jth component

    den=apply(num, 1, sum);
  p=num/den
  #list(num=num,den=den,p=p)
  p
}#end of function


###########################################################
#function of bij. It is a 3 way array with dim q x m*q x I#
#since we don't know which comp each yi comes from.       #
###########################################################

b.fun=function(y,I,q,m,beta,x,u,si,sigma)
{
  #############################################################################
  # model: (yi|zi=j) = xi*betaj + ui*bij + eij
  #     y:  unbalanced y in model above.
  #     I:  index number of i. I=100
  #     m:  number of compomemts. m=2
  #     q:  length of each bi's. q=5 
  #     x:  x matrix in model above
  #     u:  u matrix in model above
  #  beta:  fixed parameters in model above.
  #    si:  covariance of b
  # sigma:  variance of error term in model above. 
  #         it has 2 terms relatively to 2 components.
  # Usage: b.E=b.fun(y=z_yub,I=I,q=q,m=m,x=z_xub,u=z_uub,si=si.theta,
  #  sigma=sigma.theta,beta=beta12)
  ######################################################################
  
  bij=matrix(0,q,I*m)
  z_bij=array(t(bij),c(q,m,I))  
  # empty 3 way array to contain estimate of bij's.
  
  for (i in 1:I)
  {
    for (j in 1:m)
    {
      at=solve(solve(si[,,j])+(1/sigma[j])*u[,1:n[i],i]%*%
                 solve(diag(1,n[i]))%*%t(u[,1:n[i],i]) )
      bt=(1/sigma[j])*u[,1:n[i],i]%*%solve(diag(1,n[i]))%*%
        (y[1:n[i],,i]-t(x[,1:n[i],i])%*%beta[,j])
      z_bij[,j,i]=at%*%bt
    }#end of jth loop
  }#end of ith loop
  z_bij
}#end of function



###############################################################################
# Function of deltasquare: from the result of Fang and Zhang, 1990,p.4.
# For each i, we get m delta squares since we don't know which comp the i 
# comes from.
###############################################################################

delt.fun<-function(y,I,m,beta,x,u,si,sigma,b)
{
  #############################################################################
  #definition in model: (yi|zi=j) = xi*betaj + ui*bij + eij
  #y:  unbalanced y in model above.
  #I:  index number of i. I=100
  #m:  number of compomemts. m=2
  #x:  x matrix in model above
  #u:  u matrix in model above
  #beta: fixed parameters in model above.
  #si:   covariance of b
  #sigma:variance of error term in model above. 
  #  it has 2 terms relatively to 2 components.
  #b: random parameters in model above.
  #Usage: deltsq=delt.fun(y=z_yub,I=I,m=m,beta=beta12,
  #          x=z_xub,u=z_uub,si=si.theta,sigma=sigma.theta, b=b.E) 
  #  here the  random effect b, we use true value of b=b12
  #######################################################################
  
  deltsq=matrix(0,I,m) #empty Ixm (Ixm)matrix to contain estimate of delta squares.
  for (i in 1:I)
  {
    for (j in 1:m)
    {
      del1=t(b[,j,i])%*%solve(si[,,j])%*%b[,j,i];
      del2=solve(sigma[j])%*%
        t(y[1:n[i],,i]-t(x[,1:n[i],i])%*%beta[,j]-t(u[,1:n[i],i])%*%b[,j,i])%*%
        solve(diag(1,n[i]))%*%
        (y[1:n[i],,i]-t(x[,1:n[i],i])%*%beta[,j]-t(u[,1:n[i],i])%*%b[,j,i])
      deltsq[i,j]=del1+del2+10^(-6)
    }#end of jth loop
  }#end of ith loop
  deltsq
  #list(del1=del1,del2=del2,deltsq=deltsq)
}#end of function


tau.fun=function(y,I,m,rn,rs,beta,x,u,si,sigma,b,delt)
{
  ## tau.E=tau.fun(y=z_yub,I=I,m=m,rn=rn,rs=n,x=z_xub,u=z_uub,beta=beta12,          
  ##           b=b12,si=si.theta,sigma=sigma.theta,df=df.theta,delt=deltsq)
  ## for random effect b, we use the true values
  
  tau =matrix(0,I,m) # empty Ixm matrix  to contain estimates of tau.
  for (i in 1:I)
  {
    for (j in 1:m)
    {
      tau[i,j]=sqrt(2)*besselK(sqrt(2*delt[i,j]),n[i]/2-1)/
        (sqrt(delt[i,j])*besselK(sqrt(2*delt[i,j]),1-n[i]/2))
    }
  }
  tau
}


omega.fun=function(y,I,m,q,si,sigma,x,u)
{
  #############################################################################
  # model: (yi|zi=j) = xi*betaj + ui*bij + eij
  #y:  unbalanced y in model above.
  #I:  index number of i. I=100
  #m:  number of compomemts. m=2
  #q: length of each bi's. q=5
  #x:  x matrix in model above
  #u:  u matrix in model above
  #si:   covariance of b
  #sigma:variance of error term in model above.
  #      it has 2 terms relatively to 2 components.
  # Usage: omega.E=omega.fun(y=z_yub,I=I,m=m,q=q,si=si.theta,sigma=sigma.theta, 
  #                   x=z_xub,u=z_uub) 
  ######################################################################
  
  omega=array(0,c(q,q,m,I)) #empty 4 way array to contain estimate of omega
  
  for (i in 1:I)
  {
    for (j in 1:m)
    {
      omega[,,j,i]=solve(solve(si[,,j])+(1/sigma[j])*u[,1:n[i],i]%*%
                           solve(diag(1,n[i]))%*%t(u[,1:n[i],i]))
    }
  }
  omega
}

beta.fun<-function(y,pr,I,m,p,tau,sigma,x,u,b)
{
  #############################################################################
  # model: (yi|zi=j) = xi*betaj + ui*bij + eij
  #y:  unbalanced y in model above.
  #pr: current proportion for each i obs in jth component.(ie: pr.E in E step)
  #I:  index number of i. I=100
  #m:  number of compomemts. m=2
  #q: length of each bi's. q=5
  #tau: current tau values (ie, tau.E)
  #x:  x matrix in model above
  #u:  u matrix in model above
  #b:  current random parameters b (ie. b.E)
  #sigma:variance of error term in model above. it has 2 terms relatively to 2 components.
  # Usage: beta.M=beta.fun(y=z_yub,I=I,p=p,pr=pr.E,m=m,tau=tau.E, 
  #                   sigma=sigma.theta,x=z_xub,u=z_uub,b=b.E) 
  ######################################################################
  
  beta.new=matrix(0,p,m)  #empty matrix to contain beta estimates
  t_num=array(0,c(p,m*p,I))
  #empty 3 dim matrix to contain numerator of beta estimates
  t_den=array(0,c(p,m,I) )
  #empty 3 dim matrix to contain denominator of beta estimates
  for (i in 1:I)
  {
    t_num[,,i]= t(pr[i,]*tau[i,]*(1/sigma))%x%
      (x[,1:n[i],i]%*%solve(diag(1,n[i]))%*%t(x[,1:n[i],i]))
    t_den[,,i]=rep(1,p)%x%t(pr[i,]*tau[i,]*(1/sigma))*(x[,1:n[i],i]%*%                                solve(diag(1,n[i]))%*%
                                                         (cbind(y[1:n[i],,i],y[1:n[i],,i])-t(u[,1:n[i],i])%*%b[,,i]))
  }
  
  num_m=apply(t_num,c(1,2),sum) # get the sum of numerator by index I.
  den_m=apply(t_den,c(1,2),sum) # get the sum of denominater by index I.
  
  for (j in 1:m)
   {
    l=((j-1)*p+1):(j*p)   
    beta.new[,j] = solve(num_m[,l])%*%den_m[,j]
   }
  ##list(t_num=t_num,t_den=t_den,num_m=num_m,den_m=den_m,beta.new=beta.new)
  beta.new
}#end of function

sigma.fun<-function(y,m,I,r,p,q,pr,tau,omega,u,x,beta,b,rs)
{
  #############################################################################
  #definition in model: (yi|zi=j) = xi*betaj + ui*bij + eij
  #y:  unbalanced y in model above.
  #pr: current proportion for each i obs in jth component.(ie: pr.E in E step)
  #I:  index number of i. I=100
  #r:  repeated measure in balanced case. r=4
  #m:  number of compomemts. m=2
  #p: length of betai's
  #q: length of each bi's. q=5
  #tau: current tau values (ie, tau.E)
  #omega: current estimate of omega.
  #x:  x matrix in model above
  #u:  u matrix in model above
  #b:  current random parameters b (ie. b.E)
  #rs: repeated measure in unbalanced case.
  #    Radomly sampled in 2:4. rs is a vector with I=100 terms.
  #beta: current estimate of beta (i.e, beta.M)
  # usage: sigma.M=sigma.fun(y=z_yub,m=m,I=I,r=r,p=p,q=q,pr=pr.E,tau=tau.E,
  #                    omega=omega.E,u=z_uub,x=z_xub,beta=beta.M,b=b.E,rs=n)
  ######################################################################
  
  sigma.new=matrix(0,m,1)
  trace=matrix(0,1,m)
  num=matrix(0,I,m)
  den=matrix(0,1,m)
  res=array(matrix(0,m*r,I),c(r,m,I))
  
  
  for (j in 1:m){
    den[,j]=sum(pr[,j]*rs)
    for (i in 1:I){     
      #trace[,j]=sum(diag(omega[,,j,i]%*%u[,,i]%*%diag(dim(u[,,i])[2])%*%t(u[,,i])))
      trace[,j]=sum(diag(omega[,,j,i]%*%u[,1:n[i],i]%*%diag(n[i])%*%t(u[,1:n[i],i])))
      
      #res[,j,i]=y[,,i]-t(x[,,i])%*%beta[,j]-t(u[,,i])%*%b[,j,i]#res_ij is a rx1 vectro
      res[1:n[i],j,i]=y[1:n[i],,i]-t(x[,1:n[i],i])%*%beta[,j]-t(u[,1:n[i],i])%*%b[,j,i]  #res_ij is a rx1 vector
      
      
      #num[i,j]=pr[i,j]*(tau[i,j]*(t(res[,j,i])%*%diag(1,dim(t(res[,j,i]))[2])%*%res[,j,i])+trace[,j])
      num[i,j]=pr[i,j]*(tau[i,j]*(t(res[1:n[i],j,i])%*%diag(1,n[i])%*%res[1:n[i],j,i])+trace[,j])   
      
    }#end of ith loop
    #sigma.new[j]=sum(num[,j])/den[,j]
    
  }#end of jth loop
  
  #list(num=num,sigma.new=sigma.new)
  sigma.new=matrix(sum(sum(num))/sum(sum(den)),m,1)
  sigma.new
} #end of function


si.fun=function(q,I,m,pr,tau,b,omega)
{
  #############################################################################
  #definition in model: (yi|zi=j) = xi*betaj + ui*bij + eij
  #pr: current proportion for each i obs in jth component.(ie: pr.E in E step)
  #I:  index number of i. I=100
  #m:  number of compomemts. m=2
  #q: length of each bi's. q=5
  #tau: current tau values (ie, tau.E)
  #omega: current estimate of omega.
  #b:  current random parameters b (ie. b.E)
  # Usage: si.M=si.fun(q=q,I=I,m=m,pr=pr.E,tau=tau.E,b=b.E,omega=omega.E) 
  ######################################################################
  
  si.new=array(matrix(0,q,m*q),c(q,q,m))#empty 3 array to contain si matrix.
  num=array(matrix(0,I,m*q),c(q,q,m,I)) #empty 4 array to cintain numerators of 'si' estimate. 
  
  den=apply(pr,2,sum) #denominators of 'si' estimate,get the sum proportion in E-step by columns.
  
  
  for (j in 1:m){
    for (i in 1:I){
      #den=apply(pr,2,sum)
      #denominators of 'si' estimate,get the sum proportion in E-step by columns.
      num[,,j,i]=pr[i,j]*( tau[i,j]*b[,j,i]%*%t(b[,j,i])+omega[,,j,i] ) #numerator of 'si' estimate before got the sum of each component.
      #si.new[,,j]=apply(num,c(1,2,3),sum)[,,j]/den[j]
    }#end of ith loop
    #si.new[,,j]=apply(num[,,j,],c(1,2),sum)/den[j]
  }#end of jth loop
  #list(num=num,si.new=si.new,den=den)
  si.new1=apply(num,c(1,2),sum)/sum(den)##assume same variance of bij
  si.new=array(si.new1,c(q,q,m))      
  si.new
}

loglik.fun<-function(y,x,u,beta,si,sigma,pai,I,m,r)
{
  mu=array(0,c(r,m,I))
  lambda=array(matrix(0,m*r,r*I),c(r,r,m,I)) 
  d_y=array(0,c(1,m,I))
  loglik2=vector()
  for (i in 1:I)
  {
    for (j in 1:m)
    {
      ## mu[,j,i]=t(x[,,i])%*%beta[,j]# balance case
      mu[1:n[i],j,i]=t(x[,1:n[i],i])%*%beta[,j]
      ##lambda[,,j,i]=t(u[,,i])%*%si[,,j]%*%u[,,i]+sigma[j]*diag(1,r)##balance case
      lambda[1:n[i],1:n[i],j,i]=t(u[,1:n[i],i])%*%si[,,j]%*%
        u[,1:n[i],i]+sigma[j]*diag(1,n[i])
      ##d_y[,j,i]=dmvt(y[,,i],mu[,j,i],lambda[,,j,i],df[j],log=FALSE)##balance case
      d_y[,j,i]=dmvl(y[,,i][1:n[i]],mu[(1:n[i]),j,i],lambda[(1:n[i]),(1:n[i]),j,i])
    }#end of jth loop
  }#end of ith loop
  d_y=pmax(10^(-50),d_y)
  d_y=array(d_y,c(1,m,I))
  
  loglik2=log(apply((pai*d_y),3,sum))
  loglik=sum(loglik2)
  #list(d_y=d_y,loglik2=loglik2,loglik=loglik,lam=lambda,mu=mu) 
  return(loglik)
}#end of function

## loglik=loglik.fun(y=z_yub,I=I,m=m,r=r,x=z_xub,u=z_uub,beta=beta.M,
##                   sigma=sigma.M,si=si.M,df=df.M,pai=pr.M)

LMIX_LMM<-function(z_yub,z_xub,z_uub,pr.init,beta.init,sigma.init,si.init)
{
  #######################
  #definition of model  #
  #######################
  # y_{ir_m}=x_{i}beta_{m}+u_{i}b_{i_m}+epsilon_{ir_m}
  # i=1,,,100 is the index.
  # r=4 is the balanced repeated measure.
  # n[i] is the unbalanced repeated measure.
  # m=2 is the 2 component.
  # x,y and u are 3 way arrays:
  # dim(y)=(r,1,I);#response
  # dim(x)=(p,r,I);#matrix in fixed effect
  # dim(u)=(q,r,I);#design matrix in random effect.
  
  I <- dim(z_yub)[3]
  r <- dim(z_yub)[1]
  q <- dim(z_uub)[1]
  p <- dim(z_xub)[1]
  
  ##set initial value
  pr.new=pr.init;
  beta.new=beta.init;
  sigma.new=sigma.init;
  si.new=si.init
  
  loglik.theta <- -1e10
  loglik.new <- -1e10
  
  ######################
  #start EM Algorithm  #
  ######################
  
  loglik.all <- vector()
  sigma.all<-matrix(0,m,iter)
  beta.all<-array(0,c(p,m,iter))
  
  for (g in 1:iter)
  {
    pr.theta=pr.new
    beta.theta=beta.new
    sigma.theta=sigma.new
    si.theta=si.new
    loglik.theta=loglik.new
    
    ###########
    # E-step #
    ###########
    
    ##################
    #get density of y#
    ##################
    
    dy=deny.fun(y=z_yub,I=I,m=m,r=r,beta=beta.theta,x=z_xub,u=z_uub,si=si.theta,
                sigma=sigma.theta)
    
    ####################################
    #conditional probability  in E step#
    ####################################
    
    pr.E=pr.fun(y=z_yub,I=I,m=m,pr=pr.theta,dy=dy)
    b.E=b.fun(y=z_yub,I=I,q=q,m=m,x=z_xub,u=z_uub,si=si.theta,
              sigma=sigma.theta,beta=beta.theta)
    deltsq=delt.fun(y=z_yub,I=I,m=m,beta=beta.theta,
                    x=z_xub,u=z_uub,si=si.theta,sigma=sigma.theta,b=b.E)
    # here the  random effect b, we use true value of b=b12
    tau.E=tau.fun(y=z_yub,I=I,m=m,rn=rn,rs=n,x=z_xub,u=z_uub,beta=beta.theta,
                  b=b.E,si=si.theta,sigma=sigma.theta,delt=deltsq)
    #for random effect b, we use the true values
    
    omega.E=omega.fun(y=z_yub,I=I,m=m,q=q,si=si.theta,sigma=sigma.theta,x=z_xub,u=z_uub)
    
    ##########################
    #M step for EM algorithm #
    ##########################
    
    ##############################
    #get new values of parameters#
    ##############################
    
    pr.new=apply(pr.E,2,sum)/I
    #I is the index of y
    #pr.new (Ixm=100x2) is the proportions for total I index in E-step
    beta.new=beta.fun(y=z_yub,I=I,p=p,pr=pr.E,m=m,tau=tau.E,
                      sigma=sigma.theta,x=z_xub,u=z_uub,b=b.E)
    
    sigma.new=sigma.fun(y=z_yub,m=m,I=I,r=r,p=p,q=q,pr=pr.E,tau=tau.E,
                        omega=omega.E,u=z_uub,x=z_xub,beta=beta.new,b=b.E,rs=n)
    
    ##for(j in 1:m) sigma.new[j] <- ifelse(sigma.new[j]<100,1,sigma.new[j])
    si.new=si.fun(q=q,I=I,m=m,pr=pr.E,tau=tau.E,b=b.E,omega=omega.E)
    
    
    #################################################
    ##compute the loglikelihood for the iteration   #
    #################################################
    
    loglik.new=loglik.fun(y=z_yub,I=I,m=m,r=r,x=z_xub,u=z_uub,beta=beta.new,
                          sigma=sigma.new,si=si.new,pai=pr.new)
    
    loglik.all[g] <- loglik.new
    sigma.all[,g]<-sigma.new
    beta.all[,,g]<-beta.new
    
    dif=abs(loglik.new-loglik.theta)
    if(dif<eps | g >99) break
  }#end of EM iteration
  # invisible(list(loglik=loglik.all))##just for loglike
  # invisible (list(pr.theta=pr.theta,sigma.theta=sigma.theta,
  #  beta.theta=beta.theta,df.theta=df.th      
  #  eta,si.theta=si.theta,loglik=loglik.all,sigma=sigma.all,
  #  beta=beta.all,df.all=df.all,dif=dif))
  invisible (list(pr.theta=pr.new,sigma.theta=sigma.new,
                  beta.theta=beta.new,si.theta=si.new,
                  loglik=loglik.all,dif=dif))
}#end of function
##est_t=LMIX_LMM(z_yub,z_xub, z_uub,pr.init,beta.init,df.init,sigma.init,si.init)


time0=Sys.time()
  
r=4 # repeated measures for each subject in balanced case
I=400   # total subject

for(datatype in c("N","T3","T5","L","CN"))
 {   
p=4 # number of parameters in beta.
q=2  # dimension of random effect
m=2  # total 2 classes
n=rep(r,I)  # repeated measure for each subject in balanced case
nb=r*I  #total sample size of 100 groups with balanced 4 repeated measure.
#n=sample(2:r,I,replace=T) # random select repeated measure for each index.
#nunb=sum(n)         # total number in unbalanced case.

rp=200; # number of replications


#######################
#    repetition;     #
#######################

pm=10 ## number of parameters estimated
## beta(2p),sigma(m),pi(m)

## empty matrix for t_data parameters ;
est_Lpar=array(0,c(pm,m,rp))

##empty matrix for initial parameters;
est_init=array(0,c(pm,m,rp))


set.seed(99)

for(ii in 1:rp)
  
{ # start the repetition loop;
  
  # Generating X matrix
  
  mu_x=rep(0,p)                  # mean matrix of x.
  sigma_x=diag(1,p)              # variance matrix of x, identity diagnal matrix.
  xb=mvrnorm(nb,mu_x,sigma_x)    # big x matrix r*I x p.
  
  ########################################################################################
  #creat 3 dimension diagonal matrix of x
  #here z matrix is a 3 dimension matrix: I matrix, each has different repeated measure
  ########################################################################################
  
  z_x=array(t(xb),c(p,r,I))
  # 3 way array . 'Array' function  reads the matrix with column, so tranpose the x matrix.
  # empty 3 way array used to contain the balanced x matrix.
  z_xub=array(0,c(p,r,I))
  # empty 3 way array used to contain the unbalanced x matrix.
  
  for (i in 1:I)
  {
    for (l in 1:p)
    {
      z_xub[l,,i][1:n[i]]=z_x[l,,i][1:n[i]]
      # put the terms in z_x to the 3 way array unbalanced x,
      # and missing values are all 0's.
    } # end of lth loop.
  } # end of ith loop.
  
  ######################################
  #generate U matrix with q=2.
  ######################################
  
  mu_u=rep(0,q);             #mean vector of u.
  sigma_u=diag(rep(1,q))
  u=mvrnorm(nb,mu_u,sigma_u) #big balanced u matrix with dim nxq (r*I x q)
  z_u=array(t(u),c(q,r,I))   #3 way array u.
  z_uub=array(0,c(q,r,I))    #empty 3 way array used to contain unbalanced u matrix.
  
  for (i in 1:I)
  {
    for (l in 1:q)
    {
      z_uub[l,,i][1:n[i]]=z_u[l,,i][1:n[i]]   
      #put the terms in z_u to the 3 dim unbalanced u matrix, 
      #and missing values are all 0's.
    }                                        #end of lth loop.
  }                                           #end of ith loop
  
  #z_uub  #3 way array u with unbalance repeated measure, all missins are 0's.
  
  ########################################
  #generate parameters of fixed part beta#
  ########################################
  
  beta1=c(1,1,0,0) #fixed part parameters for component1
  beta2=c(0,0,1,1)    #fixed part parameters fro component2
  
  beta12=cbind(beta1,beta2)       #combine the beta parameters for 2 components.
  
  pr=c(0.4,0.6) # mixture propotion vector. 
  # It means the probability of each ith index vector (i.e 'yi') 
  # belonging to component 1 or 2  is 0.5 at the same time. 
  
  #############################################################################
  #use random binomial method to separate all (I=100)index into  2 components. #
  #############################################################################
  
  rn=rbinom(1,I,0.4) # randomly generate index number for component1. 
  # I-rn are the index number for component2.
  
  # generate random effect b with Ar(1) covariance matrix:
  
  #qseq=1:q
  #rho=0.5
  #H=abs(outer(qseq, qseq, "-"))   #times=dim(x),times=dim(y),"-"=x-y
  #sigma1=1*0.5^H                     #covariance matrix for bi in 1st component.
  #sigma2=1*0.5^H                     #covariance of bi in 2nd component.
  
  sigma1=sigma2=matrix(0.5,nrow=q,ncol=q)
  diag(sigma1)=diag(sigma2)=1
  
  if(datatype=="N") # 1. Normal 
  {
    mu=rep(0,r);
    sigma=diag(1,r)
    e11=mvrnorm(rn,mu,sigma)
    e21=mvrnorm(I-rn,mu,sigma)
    er=rbind(e11,e21)
    
    ## conditional b
    
    mu1=mu2=rep(0,q)
    b11=mvrnorm(rn,mu1,sigma1)
    b21=mvrnorm(I-rn,mu2,sigma2)
    b12=rbind(b11,b21)
  } else if(datatype=="L")
  {  
    # 2. Multivariate  Laplace
    
    ##############################
    #   conditional error and b  #
    ##############################
    
    #############################
    #generate exponential distribution#
    #############################
    
    tauij1=rexp(rn)
    tauij2=rexp(I-rn)
    
    ## conditional er;
    mu=rep(0,r);
    sigma=diag(1,r)
    
    e11=matrix(0,rn,r)
    for (i in 1:rn)
    {
      sigma=diag(1,r)*tauij1[i]
      e11[i,]=mvrnorm(1,mu,sigma)
    }
    
    e21=matrix(0,I-rn,r)
    for (i in 1:(I-rn))
    {
      sigma=diag(1,r)*tauij2[i]
      e21[i,]=mvrnorm(1,mu,sigma)
    }
    
    er=rbind(e11,e21)
    
    ## conditional b
    
    mu1=mu2=rep(0,q)
    b11=matrix(0,rn,q);
    for(i in 1:rn)
    {
      b11[i,]=mvrnorm(1,mu1,sigma1*tauij1[i])
    }
    
    b21=matrix(0,I-rn,q);
    for(i in 1:(I-rn))
    {
      b21[i,]=mvrnorm(1,mu2,sigma2*tauij2[i])
    }
    
    b12=rbind(b11,b21)
  } else if((datatype=="T1")|(datatype=="T3")|(datatype=="T5")) # Multivariate t
  {  
    #############################
    #generate Gamma distribution#
    #############################
    
    dfv=1*(datatype=="T1")+3*(datatype=="T3")+5*(datatype=="T5")
    df=c(dfv,dfv)  
    tauij1=1/(rgamma(rn, df[1]/2,df[1]/2))
    tauij2=1/(rgamma(I-rn, df[1]/2,df[2]/2))
    
    ## conditional er;
    mu=rep(0,r);
    sigma=diag(1,r)
    
    e11=matrix(0,rn,r)
    for (i in 1:rn)
    {
      sigma=diag(1,r)*tauij1[i]
      e11[i,]=mvrnorm(1,mu,sigma)
    }
    
    e21=matrix(0,I-rn,r)
    for (i in 1:(I-rn))
    {
      sigma=diag(1,r)*tauij2[i]
      e21[i,]=mvrnorm(1,mu,sigma)
    }
    
    er=rbind(e11,e21)
    
    ## conditional b
    
    mu1=mu2=rep(0,q)
    b11=matrix(0,rn,q);
    for(i in 1:rn)
    {
      b11[i,]=mvrnorm(1,mu1,sigma1*tauij1[i])
    }
    
    b21=matrix(0,I-rn,q);
    for(i in 1:(I-rn))
    {
      b21[i,]=mvrnorm(1,mu2,sigma2*tauij2[i])
    }
    
    b12=rbind(b11,b21)
  } else if(datatype=="CN")   # 4. Contaminated Normal
  {
    crate=0.95
    mu=rep(0,r);
    sigma=diag(1,r)
    u1=runif(I,0,1)
    er=(u1<crate)*mvrnorm(I,mu,sigma)+(u1>crate)*mvrnorm(I,mu,25*sigma)
    
    ## conditional b
    
    mu1=rep(0,q)
    mu2=rep(0,q)
    u2=runif(I,0,1)
    b12=(u2<0.95)*mvrnorm(I,mu1,sigma1)+(u2>0.95)*mvrnorm(I,mu2,25*sigma2)
  }
  ####################################################### 
  
  si.init=array(cbind(sigma1,sigma2),c(q,q,2));
  
  #############################
  # 3 way matrix of response  #
  #############################
  
  z_yub=y.fun(I=I,r=r,rs=n,e=er,x=z_xub,u=z_uub,beta=beta12,b=b12)
  
  pr.theta=c(0.4,0.6);           #probability of ith obs in jth component.
  beta.theta=beta12;             #fixed parameters
  sigma.theta=matrix(c(1,1),2);  #variance about errors
  si.theta=array(cbind(sigma1,sigma2),c(q,q,2));#variance about bij's
  
  pr.init=pr.theta;     #probability of ith obs in jth component.
  beta.init=beta.theta; #fixed parameters
  sigma.init=sigma.theta;
  si.init=si.theta;
  iter=200
  eps=0.0001
  
  est_L=try(
  LMIX_LMM(z_yub=z_yub,z_xub=z_xub, z_uub=z_uub,pr.init=pr.init,
                 beta.init=beta.init,sigma.init=sigma.init,
                 si.init=si.init))
  if(inherits(est_L,"try-error"))
  {next}
  
  cat(est_L$dif,"\n")
  
  ###########################################################
  # end the estimate part by the functions
  ###########################################################
  
  ## Dealing with Label Switching
  
  for (jj in 1:m)
  { ## start loop for group;
    ## construct the parameter matrix with 2 columns
    est_Lpar[,jj,ii]=rbind(est_L$pr.theta[jj],
                           as.matrix(est_L$beta.theta[,jj],ncol=1),
                           t(est_L$sigma.theta)[,jj],
                           as.matrix(as.vector(est_L$si.theta[,,jj]),ncol=1))
    ###find the initial parameter value matrix that can match the estimate matrix
    est_init[,jj,ii]=rbind(pr.init[jj],as.matrix(beta.init[,jj],ncol=1),
                           t(sigma.init)[,jj],
                           as.matrix(as.vector( si.init[,,jj]),ncol=1))
  }##end the loop for group
}####end the loop of rp:repeatition


est_Llbpar= est_Lpar


#######################################################################
#######################################################################
## Check the label and then get the Mse and Bias

est_init_tr=array(0,c(pm,m,rp))
est_Ldfsq=array(0,c(pm,m,rp))
est_Llb=array(0,c(pm,m,rp))


for (ii in 1:rp)
{
  est_init_tr[,,ii]=cbind(est_init[,1,ii],est_init[,1,ii])
  ## get same column in initial value array;
  est_Ldfsq[,,ii]=(est_Llbpar[,,ii]-est_init_tr[,,ii])^2
  
  if(sum(est_Ldfsq[,1,ii])>sum(est_Ldfsq[,2,ii]))
  {
    est_Llb[,1,ii]=est_Llbpar[,2,ii]
    est_Llb[,2,ii]=est_Llbpar[,1,ii]
  } else
  {
    est_Llb[,1,ii]=est_Llbpar[,1,ii]
    est_Llb[,2,ii]=est_Llbpar[,2,ii]
  }
} ## end of rp loop

time1=Sys.time()
time1-time0

Lname=paste0("L-",datatype,"-n",as.character(r),"-I",as.character(I),".Rdata")
fLname=paste0("C:\\Users\\weixing\\Dropbox\\LmmSimul\\",Lname)
save(est_Llb,file=fLname)

mdese_L=apply((est_Llb[1:5,1:2,]-est_init[1:5,1:2,])^2,c(1,2),median)  
# Multivariate Laplace
mdese_L
#medeff=mdese_n/mdese_L

mse_L=apply((est_Llb[1:5,1:2,]-est_init[1:5,1:2,])^2,c(1,2),mean)  # Multivariate T
mse_L
#mseeff=mse_n/mse_L

# library(latex2exp)
# xlabs=c(TeX("$\\hat{\\pi}_{1}$"),
#         TeX("$\\hat{\\beta}_{11}$"),
#         TeX("$\\hat{\\beta}_{12}$"),
#         TeX("$\\hat{\\beta}_{13}$"),
#         TeX("$\\hat{\\beta}_{14}$"),
#         TeX("$\\hat{\\pi}_{2}$"),
#         TeX("$\\hat{\\beta}_{21}$"),
#         TeX("$\\hat{\\beta}_{22}$"),
#         TeX("$\\hat{\\beta}_{23}$"),
#         TeX("$\\hat{\\beta}_{24}$"))

#par(mfrow=c(2,5))
#for(j in seq(2))
#{
#  for(i in seq(5))
#  {
    #boxplot.matrix(cbind(est_Llb[i,j,],est_tlb[i,j,],est_zlb[i,j,]),
    #               names=c('L','T','N'),xlab=xlabs[i+5*(j-1)])
#    boxplot(est_Llb[i,j,], names=c('L'),xlab=xlabs[i+5*(j-1)])
#  }
#}

}