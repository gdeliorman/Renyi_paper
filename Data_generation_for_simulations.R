##to find metrics from random parameters for simulations

library(MASS) ##TO GENERATE MULTIVARIATE NORMAL DATA

##schizo data parameters
sigma1<-matrix(c(544.3285, 0, 300.9816, 0, 
                 0,550.6597,0, 304.4222, 
                 300.9816,0, 180.6831, 0, 
                 0,304.4222, 0,180.9433), nrow=4, ncol=4)

mu<-c(-11.47191,-16.0151,-6.646067,-8.938326)


##generate random parameters
##for data 1 mu
mu_1<-rnorm(750, mean = 1.5, sd = 0.5)
mu_2<-rnorm(750, mean = 4, sd = 2)
mu_3<-rnorm(750, mean = 7, sd = 0.75)
mu_4<-rnorm(750, mean = 10, sd = 0.25)


##for data 2 mu
mu_1<-rnorm(750, mean = 45, sd = 2)
mu_2<-rnorm(750, mean = 4, sd = 5)
mu_3<-rnorm(750, mean = 1, sd = 0.95)
mu_4<-rnorm(750, mean = 6, sd = 3.25)


##generate random sigmas
n <- 100  # Number of matrices
p <- 4  # Dimension
df <- 10  # Degrees of freedom
Sigma <- toeplitz((p:1)/p)  # the matrix parameter of the distribution
# Draw n Wishart distributed matrices
Sigmas <- drop(rWishart(1000, df, Sigma))


for(i in 1:1) {
  mu<- c(mu_1[i], mu_2[i], mu_3[i], mu_4[i])
  mu<- abs(mu)
  Sigma<- Sigmas[, ,i:i]
  Sigma<- abs(Sigma)
  print(is.positive.definite(Sigma)) }

##n=sample size
data_normal <- as.data.frame(mvrnorm(n=500, mu=mu, Sigma=sigma1), empirical = TRUE)

T0<-data_normal[,1]
T1<-data_normal[,2]
S0<-data_normal[,3]
S1<-data_normal[,4]
treatment<-rbinom(500,1,0.5)
unobserved_data<-data.frame(T0,T1,S0,S1, treatment)


##observed data
observed_data<-matrix(NA,nrow= 500, ncol=3)
for (i in 1:500) {
  if(unobserved_data$treatment[i]==0)
  {
    observed_data[i,]<-c(unobserved_data$T0[i], unobserved_data$S0[i], unobserved_data$treatment[i])
  }
  else
    observed_data[i,]<-cbind(unobserved_data$T1[i], unobserved_data$S1[i], unobserved_data$treatment[i])
}

##observed data
placebo<-observed_data[observed_data[,3] == "0",]
experimental<-observed_data[observed_data[,3] == "1",]

True_end<-observed_data[,1]
Surr_end<-observed_data[,2]

T0<- placebo[,1]
T1<- experimental[,1]
S0<- placebo[,2]
S1<- experimental[,2]


##mean
mT0<-mean(T0)
mT1<-mean(T1)
mS0<-mean(S0)
mS1<-mean(S1)

beta<-mT1-mT0
alpha<- mS1-mS0
muu<- c(mT0, mT1, mS0, mS1)

mudelta<-c(beta,alpha)

##variances
vT0T0<-var(T0,T0)
vT1T1<-var(T1,T1)
vS1S1<-var(S1,S1)
vS0S0<-var(S0,S0)
vT1S1<-var(T1,S1)
vT0S0<-var(T0,S0)


##correlations
cT0S0<-cor(T0,S0)
cT1S1<-cor(T1,S1)
cTS<-cor(True_end,Surr_end)

sigma<- matrix( c(vT0T0, 0, vT0S0,0,
                  0,vT1T1,0,vT1S1,
                  vT0S0,0,vS0S0,0,
                  0,vT1S1,0,vS1S1), nrow = 4, ncol = 4)
A<-matrix(c(-1,0,1,0,0,-1,0,1),  ncol=4,nrow = 2)
sigmadelta<- A%*%sigma%*% t(A)
det(sigmadelta)


ICA_armd<-ICA.ContCont(0.7692674, 0.711794,  264.797, 231.7709, 188.9261, 132.6380, T0T1=seq(-1, 1, by=.05), 
                       T0S1=seq(-1, 1, by=.05), T1S0=seq(-1, 1, by=.1), S0S1=seq(-1, 1, by=.1))

#ICAs_armd<- cbind(ICA_armd$Pos.Def,ICA_armd$GoodSurr, ICA_armd$GoodSurr$ICA^2)

ICA_schizo<-ICA.ContCont(0.961795, 0.946760, 544.3285, 550.6597, 180.6831, 180.9433, T0T1=seq(-1, 1, by=.05), 
                  T0S1=seq(-1, 1, by=.05), T1S0=seq(-1, 1, by=.05), S0S1=seq(-1, 1, by=.05))


ICA<-ICA.ContCont(cT0S0, cT1S1, vT0T0, vT1T1, vS0S0, vS1S1, T0T1=seq(-1, 1, by=.05), 
                  T0S1=seq(-1, 1, by=.05), T1S0=seq(-1, 1, by=.05), S0S1=seq(-1, 1, by=.05))


##ICA from package
ICA<-ICA.ContCont(0.3937989, 0.3872117, 14.526506, 16.355106, 10.058608, 6.414937, T0T1=seq(-1, 1, by=.05), 
                  T0S1=seq(-1, 1, by=.1), T1S0=seq(-1, 1, by=.1), S0S1=seq(-1, 1, by=.1))

ICAs_2<- cbind(ICA$Pos.Def,ICA$GoodSurr, ICA$GoodSurr$ICA^2)

##parameters 1
##to genarte the same data save the parameters
mu<- c(1.150252,  8.030536,  6.705058, 10.194465)
sigma1<- matrix(c(14.526506, 12.213278,  4.763175, 2.078263,
                  12.213278, 16.355106,  8.573634, 3.992938,
                  4.763175,  8.573634, 10.058608, 5.898407,
                  2.078263,  3.992938,  5.898407, 6.414937), nrow = 4, ncol = 4)

is.positive.definite(sigma)

##choose parameters
high_ICA_n100<- matrix(c(0.4, 0.3937989, -0.5, -0.5, 0.3872117, -0.4, 0.9686586, 0.9382996), ncol = 8)
high_ICA_n500<- matrix(c(-0.6, 0.3937989, -0.9, -0.9, 0.3872117, -0.1, 0.9691966, 0.939342), ncol = 8)
high_ICA_n1000<- matrix(c(-0.5, 0.3937989, -0.6, -0.5, 0.3872117, 0.4, 0.9704561, 0.941785), ncol = 8)


medium_ICA_n1000<-  matrix(c(-0.9, 0.3937989, -0.6, -0.7, 0.3872117, 0.1, 0.7959334,0.63351), ncol = 8)
medium_ICA_n500<-  matrix(c(0.1, 0.3937989, -0.7, -0.5, 0.3872117, -0.9, 0.7466279, 0.5574532), ncol = 8)
medium_ICA_n100<- matrix(c(-0.1,0.3937989, -0.7, -0.6, 0.3872117, -0.9, 0.7150889, 0.5113521), ncol = 8)

low_ICA_n1000<- matrix(c(-0.6, 0.3937989, -0.7, -0.2, 0.3872117, -0.9,0.4616699,0.2131391), ncol = 8)
low_ICA_n500<- matrix(c(-0.9,0.3937989,-0.7, -0.1, 0.3872117, -0.9, 0.3935466, 0.1548789), ncol = 8)
low_ICA_n100<- matrix(c(-0.8,0.3937989, -0.2, -0.4, 0.3872117,-0.9,0.3809287, 0.1451067), ncol = 8)


##parameters 2
##to genarte the same data save the parameters
mu<- c(42.441383,  2.008065,  0.672878,  7.393521)
sigma<- matrix(c(21.37292, 17.27286, 14.93269, 11.27921,
                 17.27286, 18.72757, 15.63501, 12.17874,
                 14.93269, 15.63501, 17.86141, 15.42664,
                 11.27921, 12.17874, 15.42664, 15.61918), nrow = 4, ncol = 4)

is.positive.definite(sigma)

##choose parameters
high_ICA_n1000<- matrix(c(0.8, 0.759604,  0.4,  0.3, 0.7111961,  0.2, 0.96536978, 0.931938804), ncol = 8)
high_ICA_n500<- matrix(c(-0.3, 0.759604, -0.7, -0.8, 0.7111961, -0.9, 0.9458685,  0.8946672), ncol = 8)
high_ICA_n100<- matrix(c(-0.7, 0.759604, -0.9, -0.9, 0.7111961, -0.9, 0.9106026,  0.8291971), ncol = 8)


medium_ICA_n1000<-  matrix(c(0.5,0.759604,-0.2, 0.7, 0.7111961, 0.1, 0.7247139, 0.52521021), ncol = 8)
medium_ICA_n500<-  matrix(c(-0.6,0.759604,-0.7, -0.4,0.7111961,-0.7, 0.7800917, 0.6085430), ncol = 8)
medium_ICA_n100<- matrix(c(0.8, 0.759604, 0.5,  0.4, 0.7111961, 0.2, 0.71676942,0.513758398), ncol = 8)

low_ICA_n1000<- matrix(c(0.7,0.759604, 0.9, 0.7, 0.7111961, 0.9, -0.3511859, 0.1233315), ncol = 8)
low_ICA_n500<- matrix(c(0.7, 0.759604, 0.6, 0.7, 0.7111961, 0.9,  0.50263889, 0.252645856), ncol = 8)
low_ICA_n100<- matrix(c(0.4, 0.759604, 0.8, 0.2, 0.7111961, 0.5, 0.4343410,0.1886521), ncol = 8)

write.csv(ICAs, "ICAS_schizo_data")


##armd data parameters
sigma1<- matrix( c(264.7975 , 0, 172.0601,0,
                   0,231.7709,0,124.8009,
                   172.0601,0,188.9261,0,
                   0,124.8009,0,132.6380), nrow = 4, ncol = 4)

mu<- c(-11.742268 -14.654762  -6.030928  -7.809524)
