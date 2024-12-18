# Line Plot of Lung Function Growth Data 
# for the paper
# Robust Mixture of Linear Mixed Modeling via Multivariate Laplace Distribution

# Import bootstrapped data from the R-Workspace "bigy" by Bai et al. (2016) 

lung=read.table("C:\\Users\\weixi\\Dropbox\\LmmSimul\\codes\\topeka.dat")

xmin=min(lung[,3])
xmax=max(lung[,3])
ymin=-0.04082
ymax=max(lung[,6])

plot(lung[,3][lung[,1]==1],lung[,6][lung[,1]==1],type="l",xlim=c(xmin,xmax),
      ylim=c(ymin,ymax),xlab="Age",ylab="Log(FEV1)",col="darkgray")
for(j in seq(1,300,by=1))
{
  if(length(lung[,1]==j)>1)
  {
   lines(lung[,3][lung[,1]==j],lung[,6][lung[,1]==j],col="darkgray")
  }  
}

lines(ksmooth(lung[,3],lung[,6],bandwidth=3),lwd=2)
