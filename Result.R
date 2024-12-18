# To-Boxplot

library(latex2exp)
xlabs=c(TeX("$\\hat{p}_{1}$"),
        TeX("$\\hat{\\beta}_{11}$"),
        TeX("$\\hat{\\beta}_{12}$"),
        TeX("$\\hat{\\beta}_{13}$"),
        TeX("$\\hat{\\beta}_{14}$"),
        TeX("$\\hat{p}_{2}$"),
        TeX("$\\hat{\\beta}_{21}$"),
        TeX("$\\hat{\\beta}_{22}$"),
        TeX("$\\hat{\\beta}_{23}$"),
        TeX("$\\hat{\\beta}_{24}$"))


est_init=array(matrix(c(0.4,1.0,1.0,0.0,0.0,1.0,1.0,0.5,0.5,1.0,
         0.6,0.0,0.0,1.0,1.0,1.0,1.0,0.5,0.5,1.0)),c(10,2,200))

resTable1=resTable2=resTable3=array(0,c(50,5,4))


k=1
# for(r in c(4,8))
#  {
#   if(r==4) 
#    {Iseq=c(200,400)}
#   else if(r==8)
#    {Iseq=c(100,200)}
#   for(I in Iseq)for(r in c(4,8))
for(r in c(8))
{
  if(r==4) 
  {Iseq=c(200,400)}
  else if(r==8)
  {Iseq=c(200)}
  for(I in Iseq)
    
    
   {
    j=1
    for(datatype in c("N","T3","T5","L","CN"))
    {
    zname=paste0("N-",datatype,"-n",as.character(r),"-I",as.character(I),".Rdata")
    fzname=paste0("C:\\Users\\weixing\\Dropbox\\LmmSimul\\",zname)
    load(fzname)
    tname=paste0("T-",datatype,"-n",as.character(r),"-I",as.character(I),".Rdata")
    ftname=paste0("C:\\Users\\weixing\\Dropbox\\LmmSimul\\",tname)
    load(ftname)
    Lname=paste0("L-",datatype,"-n",as.character(r),"-I",as.character(I),".Rdata")
    fLname=paste0("C:\\Users\\weixing\\Dropbox\\LmmSimul\\",Lname)
    load(fLname)
    
    ## -------- Mean ------------------------------------------------
    
    abias_t=abs(apply((est_tlb[1:5,1:2,]-est_init[1:5,1:2,]),c(1,2),mean))  
    # Multivariate T
    abias_L=abs(apply((est_Llb[1:5,1:2,]-est_init[1:5,1:2,]),c(1,2),mean))  
    # Multivariate Laplace
    abias_n=abs(apply((est_zlb[1:5,1:2,]-est_init[1:5,1:2,]),c(1,2),mean))  
    
    ## ----medse-----------------------------------------------------------------------
    mdese_t=apply((est_tlb[1:5,1:2,]-est_init[1:5,1:2,])^2,c(1,2),median)  
        # Multivariate T
    mdese_L=apply((est_Llb[1:5,1:2,]-est_init[1:5,1:2,])^2,c(1,2),median)  
        # Multivariate Laplace
    mdese_n=apply((est_zlb[1:5,1:2,]-est_init[1:5,1:2,])^2,c(1,2),median)  
        # Multivariate N
    medeffnt=mdese_n/mdese_t
    medeffnL=mdese_n/mdese_L
    
    #temp1=cbind(c(mdese_n),c(mdese_t),c(medeffnt),c(mdese_L),c(medeffnL))
    temp1=cbind(c(abias_n),c(abias_t),c(medeffnt),c(abias_L),c(medeffnL))
    resTable1[,j,k]=c(t(temp1))

    ## ----mse-------------------------------------------------------------------------
    mse_t=apply((est_tlb[1:5,1:2,]-est_init[1:5,1:2,])^2,c(1,2),mean)  
      # Multivariate T
    mse_L=apply((est_Llb[1:5,1:2,]-est_init[1:5,1:2,])^2,c(1,2),mean)  
      # Multivariate Laplace
    mse_n=apply((est_zlb[1:5,1:2,]-est_init[1:5,1:2,])^2,c(1,2),mean)  
      # Multivariate N
    mseeffnt=mse_n/mse_t
    mseeffnL=mse_n/mse_L
    
    temp2=cbind(c(abias_n),c(abias_t),c(mseeffnt),c(abias_L),c(mseeffnL))
    resTable2[,j,k]=c(t(temp2))
    
    temp3=cbind(c(mse_n),c(mse_t),c(mseeffnt),c(mse_L),c(mseeffnL))
    resTable3[,j,k]=c(t(temp3))
    
    j=j+1
    
    x11()
    par(mfrow=c(2,5))
    for(jj in seq(2))
    {
      for(ii in seq(5))
      {
        boxplot.matrix(cbind(est_tlb[ii,jj,],est_Llb[ii,jj,],
                             est_zlb[ii,jj,]),
                       names=c('T','L','N'),xlab=xlabs[ii+5*(jj-1)])
        if(jj==1)
        {
          hv=0.4*(ii==1)+1*(ii==2)+1*(ii==3)+0*(ii==4)+0*(ii==5)
        } else
        {
          hv=0.6*(ii==1)+0*(ii==2)+0*(ii==3)+1*(ii==4)+1*(ii==5)
        }
        abline(h=hv,lty=1)
      }
    }
    #mtext(paste0("Boxplots: ",datatype,'-n',as.character(r),'-m',as.character(I)), side = 3, line = - 2, outer = TRUE)
   }
    k=k+1
   } 
 }

namelist=list(c('pi(N)','pi(T)','pi(T/N)','pi(L)','pi(L/N)',
    'beta11(N)','beta11(T)','beta11(T/N)','beta11(L)','beta11(L/N)',
    'beta12(N)','beta12(T)','beta12(T/N)','beta12(L)','beta12(L/N)',
    'beta13(N)','beta13(T)','beta13(T/N)','beta13(L)','beta13(L/N)',
    'beta14(N)','beta14(T)','beta14(T/N)','beta14(L)','beta14(L/N)',
    'beta21(N)','beta21(T)','beta21(T/N)','beta21(L)','beta21(L/N)',
    'beta22(N)','beta22(T)','beta22(T/N)','beta22(L)','beta22(L/N)',
    'beta23(N)','beta23(T)','beta23(T/N)','beta23(L)','beta23(L/N)',
    'beta24(N)','beta24(T)','beta24(T/N)','beta24(L)','beta24(L/N)'),
  c("N","T3","T5","L","CN"),
  c("r=4,I=200","r=4,I=400","r=8,I=100","r=8,I=200"))

resTable1Final=resTable1[-c(26:30),,]
dimnames(resTable1Final)=namelist
round(resTable1Final,4)

resTable2Final=resTable2[-c(26:30),,]
dimnames(resTable2Final)=namelist
round(resTable2Final,4)

resTable3Final=resTable3[-c(26:30),,]
dimnames(resTable3Final)=namelist
round(resTable3Final,4)

resTable4Final=round(resTable1Final[-c(3,5,8,10,13,15,18,20,23,25,28,30,33,35,38,40,43,45),,],4)
resTable5Final=round(resTable3Final[-c(3,5,8,10,13,15,18,20,23,25,28,30,33,35,38,40,43,45),,],4)

# MSE Comparison

# n=4, I=200 and 400


resTable5Final[,,1][]>=resTable5Final[,,2]

# n=8, I=100 and 200
resTable5Final[,,3]>=resTable5Final[,,4]

# I=200, n=4 and 8

resTable5Final[,,1]>=resTable5Final[,,4]