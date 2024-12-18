# Stability Study of Lung Function Growth Data 
# for the paper
# Robust Mixture of Linear Mixed Modeling via Multivariate Laplace Distribution

# Outliers are 2,4,6,8,10.

outliers=c(2,4,6,8,10)

x2=matrix( c(0.189, 0.164,  0.192,  0.165,
           0.644, 0.169,  0.645,  0.171,
           0.118, 0.160,  0.120,  0.161,
           0.060, 0.013,  0.060,  0.013,
           -0.365, 0.087, -0.369,  0.090,
           0.096, 0.007,  0.096,  0.007,
           -0.610, 0.088, -0.609,  0.093,
           0.103, 0.008,  0.103,  0.009),byrow=T,nrow=8)

x4=matrix(c(0.191, 0.167,  0.197, 0.173, 
0.644, 0.172,  0.643, 0.179, 
0.118, 0.160,  0.125, 0.160, 
0.060, 0.013,  0.059, 0.013, 
-0.370, 0.090, -0.378, 0.094, 
0.096, 0.007,  0.097, 0.007, 
-0.608, 0.095, -0.606, 0.103, 
0.103, 0.009,  0.103, 0.009),byrow=T,nrow=8)


x6=matrix(c(0.193, 0.170,  0.197, 0.176,
              0.643, 0.175,  0.642, 0.182,
              0.121, 0.159,  0.130, 0.162,
              0.060, 0.013,  0.059, 0.013,
              -0.373, 0.091, -0.383, 0.097,
              0.096, 0.007,  0.097, 0.007,
              -0.609, 0.096, -0.606, 0.107,
              0.103, 0.009,  0.103, 0.009)
          , byrow=T,nrow=8)

x8=matrix(c( 0.194, 0.172,  0.198, 0.178,
             0.642, 0.177,  0.639, 0.184,
             0.123, 0.160,  0.135, 0.160,
             0.060, 0.013,  0.059, 0.013,
             -0.375, 0.093, -0.387, 0.097,
             0.096, 0.007,  0.097, 0.007,
             -0.609, 0.096, -0.605, 0.110,
             0.103, 0.009,  0.103, 0.010)
          , byrow=T,nrow=8)

x10=matrix(c( 0.239, 0.161, 0.253, 0.162,
              0.577, 0.169, 0.561, 0.165,
              0.130, 0.236, 0.148, 0.269,
              0.059, 0.018, 0.057, 0.021,
              -0.311, 0.247,-0.283, 0.268,
              0.092, 0.021, 0.089, 0.022,
              -0.607, 0.124,-0.615, 0.158,
              0.103, 0.011, 0.104, 0.016),byrow=T,nrow=8)


X=array(0,c(8,4,5))

X[,,1]=x2
X[,,2]=x4
X[,,3]=x6
X[,,4]=x8
X[,,5]=x10

y1min=min(X[,1,1],X[,1,2],X[,1,3],X[,1,4],X[,1,5])
y1max=max(X[,1,1],X[,1,2],X[,1,3],X[,1,4],X[,1,5])

y2min=min(X[,3,1],X[,3,2],X[,3,3],X[,3,4],X[,3,5])
y2max=max(X[,3,1],X[,3,2],X[,3,3],X[,3,4],X[,3,5])

par(mfrow=c(1,2))

plot(1, type = "n", axes = TRUE, xlab = "Outliers in first subject", ylab = "Estimates", 
     xlim = c(2,10), ylim = c(y1min, y1max),)

# Add axes manually (optional, ensures clear display)
axis(1)  # X-axis
axis(2)  # Y-axis
box(lty = "blank")  # Remove the surrounding box

  for(k in seq(8))
  {
    lines(outliers,X[k,1,])
  }


plot(1, type = "n", axes = TRUE, xlab = "Outliers in the first two subjects", ylab = "Estimates", 
     xlim = c(2,10), ylim = c(y2min, y2max))

# Add axes manually (optional, ensures clear display)
axis(1)  # X-axis
axis(2)  # Y-axis
box(lty = "blank")  # Remove the surrounding box

  for(k in seq(8))
  {
    lines(outliers,X[k,3,])
  }



