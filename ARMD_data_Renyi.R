##ARMDphrenia data Renyi measures alpha=0.1, 0.5, 1, 1.25, 1.5, 2, 2.5, 3

library(pracma) ##integral2
library(extraDistr) ##bivariate normal
library(Surrogate) ##to get data and to calculate ICA normal
library(LaplacesDemon) ##to check positive definite matrix
library(dplyr) ##for merge

data("ARMD")
placebo<-ARMD[ARMD$Treat == "-1",]
experimental<-ARMD[ARMD$Treat == "1",]

Tr<-ARMD$Diff52
S<-ARMD$Diff24

T0<- placebo$Diff52
T1<- experimental$Diff52
S0<- placebo$Diff24
S1<- experimental$Diff24


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
cTS<-cor(Tr,S)

sigma<- matrix( c(vT0T0, 0, vT0S0,0,
                  0,vT1T1,0,vT1S1,
                  vT0S0,0,vS0S0,0,
                  0,vT1S1,0,vS1S1), nrow = 4, ncol = 4)
A<-matrix(c(-1,0,1,0,0,-1,0,1),  ncol=4,nrow = 2)
sigmadelta<- A%*%sigma%*% t(A)
det(sigmadelta)

##ICA from package
ICA<-ICA.ContCont(cT0S0, cT1S1, vT0T0, vT1T1, vS0S0, vS1S1, T0T1=seq(-1, 1, by=.1), 
                  T0S1=seq(-1, 1, by=.1), T1S0=seq(-1, 1, by=.1), S0S1=seq(-1, 1, by=.1))

ICA_sur<- cbind(ICA$GoodSurr$ICA, ICA$GoodSurr$Sigma.Delta.T, (ICA$GoodSurr$ICA)^2)
colnames(ICA_sur)<- c("rh", "sigmadelta_t","ICA_formul")
ICA_sur<- as.data.frame(ICA_sur)

Sigma_matrix<- matrix(NA, nrow=4, ncol=4)
Sigma_matrix[1,1]<- vT0T0
Sigma_matrix[2,2]<- vT1T1
Sigma_matrix[3,3]<- vS0S0
Sigma_matrix[4,4]<- vS1S1
Sigma_matrix[1,3]<-Sigma_matrix[3,1]<- vT0S0
Sigma_matrix[2,4]<-Sigma_matrix[4,2]<- vT1S1

##possible correlations:

T0T1=seq(-1, 1, by=.1)
T0S1=seq(-1, 1, by=.1)
T1S0=seq(-1, 1, by=.1)
S0S1=seq(-1, 1, by=.1)
combins<- expand.grid(T0T1, T0S1, T1S0, S0S1)

ralpha_01<-0.1
ralpha_05<-0.5
ralpha_125<-1.25
ralpha_15<-1.5
ralpha_2<-2
ralpha_25<-2.5
ralpha_3<-3

df<- data.frame()
for (i in 1:nrow(combins)) {
  
  T0T1 <- combins[i, 1]  
  T0S1 <- combins[i, 2]  
  T1S0 <- combins[i, 3]  
  S0S1 <- combins[i, 4]    
  
  Sigma_matrix[1,2]<-Sigma_matrix[2,1]<- T0T1*sd(T0)*sd(T1)
  Sigma_matrix[1,4]<-Sigma_matrix[4,1]<- T0S1*sd(T0)*sd(S1)
  Sigma_matrix[2,3]<-Sigma_matrix[3,2]<- T1S0*sd(T1)*sd(S0)
  Sigma_matrix[3,4]<-Sigma_matrix[4,3]<- S0S1*sd(S0)*sd(S1)
  
  
  if(is.positive.definite(Sigma_matrix)==TRUE){
    print(i)
    pay<-(sqrt(vT0T0*vS0S0)*cT0S0)+(sqrt(vT1T1*vS1S1)*cT1S1)-(sqrt(vT1T1*vS0S0)*T1S0)-(sqrt(vT0T0*vS1S1)*T0S1)
    payda<-sqrt((vT0T0+vT1T1-2*sqrt(vT0T0*vT1T1)*T0T1)*(vS0S0+vS1S1-2*sqrt(vS0S0*vS1S1)*S0S1))
    rh<-pay/payda
    ICA_formul<-rh^2
    
    sigmadelta_s<-var(S0,S0)+var(S1,S1)- (2*S0S1*sqrt(var(S0,S0)*var(S1,S1)))
    sigmadelta_t<-var(T0,T0)+var(T1,T1)- (2*T0T1*sqrt(var(T0,T0)*var(T1,T1)))
    vv<- rh*sqrt(sigmadelta_t)*sqrt(sigmadelta_s)
    sigmadelta_ts<- matrix(c(sigmadelta_t, vv, sigmadelta_s, vv), nrow = 2, ncol=2)
    det(sigmadelta_ts)
    
    fx<-function(x) dnorm(x,mean = beta,sd=sqrt(sigmadelta_t))
    fy<-function(y) dnorm(y,mean=alpha,sd=sqrt( sigmadelta_s))
    fxy<- function(x,y) dbvnorm(x,y, mean1 = beta, mean2 = alpha, sd1 =sqrt(sigmadelta_t) , sd2 = sqrt(sigmadelta_s), cor =rh)
    
    KLfun<-  function(x,y) fxy(x,y)* ( log(  (fxy(x,y)) / (fx(x)*fy(y)) ))
    RY_01<-  function(x,y) fxy(x,y)^ralpha_01 * ( (fx(x)*fy(y)) ^(1-ralpha_01))
    RY_05<-  function(x,y) fxy(x,y)^ralpha_05 * ( (fx(x)*fy(y)) ^(1-ralpha_05))
    RY_125<- function(x,y) fxy(x,y)^ralpha_125 * ( (fx(x)*fy(y))^(1-ralpha_125))
    RY_15<-  function(x,y) fxy(x,y)^ralpha_15 * ( (fx(x)*fy(y)) ^(1-ralpha_15))
    RY_2<-   function(x,y) fxy(x,y)^ralpha_2 * ( (fx(x)*fy(y)) ^(1-ralpha_2))
    RY_25<-  function(x,y) fxy(x,y)^ralpha_25 * ( (fx(x)*fy(y)) ^(1-ralpha_25))
    RY_3<-   function(x,y) fxy(x,y)^ralpha_3 * ( (fx(x)*fy(y)) ^(1-ralpha_3))
    
    
    ##find min and max
    fT_bound<-matrix(NA, ncol=2)
    for (j in -300:300) {
      
      if (fx(j)>=0.00001)
      {
        
        c<-cbind(j, fx(j))
        fT_bound<-rbind(c,fT_bound)
      }
    }
    
    ###finding bound for fT and fS
    fS_bound<-matrix(NA, ncol=2)
    for (k in -300:300) {
      if (fy(k)>=0.00001)
      {
        
        c<-cbind(k,fy(k))
        fS_bound<-rbind(c,fS_bound)
        
      }
    }
    
    
    xmax<- max(na.omit(fT_bound[,1]))
    xmin<- min(na.omit(fT_bound[,1]))
    ymax<- max(na.omit(fS_bound[,1]))
    ymin<- min(na.omit(fS_bound[,1]))
    
    ##integral 2
    #KL<-integral2(KLfun, xmin,xmax, ymin, ymax, reltol = 1e-10)
    #ICA_1<- 1-exp(-2*KL$Q)
    
    
    KLfun2 <- function(X) {
      x <- X[1]
      y <- X[2]
      fxy(X[1],X[2])* ( log(  (fxy(X[1],X[2])) / (fx(X[1])*fy(X[2])) ))
    }
    
    
    
    KL<-cubature::cubintegrate(f = KLfun2, lower = c(xmin,ymin), upper = c(xmax,ymax),  method = "hcubature")
    ICA_1<- 1-exp(-2*KL$integral)
    
    ##same as previous one as expected
    new_formul_ica_1<- ICA_formul
    
    R_01<-integral2(RY_01, xmin,xmax, ymin, ymax, reltol = 1e-10)
    R01<-(1/ (ralpha_01-1))* log(R_01$Q)
    ICA_01<- 1-exp(-2*R01)
    
    new_formul01<- 1-(1-rh^2)*(1-(1-ralpha_01)^2*rh^2)^(-1/(1-ralpha_01))
    
    R_05<-integral2(RY_05, xmin,xmax, ymin, ymax, reltol = 1e-10)
    R05<-(1/ (ralpha_05-1))* log(R_05$Q)
    ICA_05<- 1-exp(-2*R05)
    
    new_formul05<- 1-(1-rh^2)*(1-(1-ralpha_05)^2*rh^2)^(-1/(1-ralpha_05))
    
    R_125<-integral2(RY_125, xmin,xmax, ymin, ymax, reltol = 1e-10)
    R125<-( 1/(ralpha_125-1))* log(R_125$Q)
    ICA_125<-1-exp(-2*R125)
    
    new_formul125<- 1-(1-rh^2)*(1-(1-ralpha_125)^2*rh^2)^(-1/(1-ralpha_125))
    
    
    R_15<-integral2(RY_15, xmin,xmax, ymin, ymax, reltol = 1e-10)
    R15<- log(R_15$Q)* (1/ (ralpha_15-1))
    ICA_15<-1-exp(-2*R15)
    
    new_formul15<- 1-(1-rh^2)*(1-(1-ralpha_15)^2*rh^2)^(-1/(1-ralpha_15))
    
    R_2<-integral2(RY_2,xmin,xmax, ymin, ymax,reltol = 1e-10)
    R2<- log(R_2$Q)* (1/ (ralpha_2-1))
    ICA_2<-1-exp(-2*R2)
    
    new_formul2<- 1-(1-rh^2)*(1-(1-ralpha_2)^2*rh^2)^(-1/(1-ralpha_2))
    
    R_25<-integral2(RY_25, xmin,xmax, ymin, ymax,, reltol = 1e-10)
    R25<- log(R_25$Q)* (1/ (ralpha_25-1))
    ICA_25<- 1-exp(-2*R25)
    
    new_formul25<- 1-(1-rh^2)* ((1-(1-ralpha_25)^2*rh^2))^((-1/(1-ralpha_25)))
    
    R_3<-integral2(RY_3,xmin,xmax, ymin, ymax, reltol = 1e-10)
    R3<- log(R_3$Q)* (1/ (ralpha_3-1))
    ICA_3<- 1-exp(-2*R3)
    
    new_formul3<- 1-(1-rh^2)*(1-(1-ralpha_3)^2*rh^2)^(-1/(1-ralpha_3))
    
    print(i)
    output<- as.data.frame(cbind(ICA_01, new_formul01, ICA_05, new_formul05,ICA_formul, ICA_1,  ICA_125,new_formul125,
                                 ICA_15, new_formul15, ICA_2,new_formul2, ICA_25, new_formul25, ICA_3,new_formul3, sigmadelta_t))
    df = rbind(df, output) }
  
  
}

##to check values 
df_full<- cbind(df, ICA_sur)

library(utils)
write.csv(df_full, 'ARMD_data_renyis.csv')

##box plot ICA_alphas and ICA
boxplot(df_full[,3], df_full[,4], df_full[,2], df_full[,1], df_full[,5], df_full[,6], df_full[,7], df_full[,8], df_full[,9], 
        names = c(expression(ICA[0.1]), expression(ICA[0.5]), expression(ICA[1]), expression(ICA[]), expression(ICA[1.25]), expression(ICA[1.5]), expression(ICA[2])
                  , expression(ICA[2.5]), expression(ICA[3]))
        #, border = "brown", col="orange")
)


boxplot(df_full[,1], df_full[,2],
        names = c("using integral", "using formula"), main=bquote('Box plot of'~ICA[0.1]), ylim=c(0,1) )


boxplot(df_full[,9], df_full[,10],
        names = c("using integral", "using formula"), main=bquote('Box plot of'~ICA[1.5]), ylim=c(0,1) )

##box plot ICA_alphas 
boxplot(df_full[,3], df_full[,4], df_full[,2], df_full[,5], df_full[,6], df_full[,7], df_full[,8], df_full[,9], 
        names = c(expression(ICA[0.1]), expression(ICA[0.5]), expression(ICA[1]), expression(ICA[1.25]), expression(ICA[1.5]), expression(ICA[2])
                  , expression(ICA[2.5]), expression(ICA[3]))
        #, border = "brown", col="orange")
)


hist(df_full[,5], main=bquote('Histogram of'~ICA[1]), xlab = "ICA_1")
hist(df_full[,18], main=bquote('Histogram of'~rho[Delta]), xlab = "rho")
hist(df_full[,5], main=bquote('Histogram of'~ICA[1]), xlab = "ICA_1")


##summarize
summary(df_full[,1])
sd(df_full[,1])
mm_ICA<-max(df_full[,1])-min(df_full[,1])
hist(df_full[,1], main=bquote('Histogram of'~ICA[01]), xlab = "ICA_01")


summary(as.numeric(na.omit(df_full[,2])))
sd(as.numeric(na.omit(df_full[,2])))
mm_ICA1<-max(as.numeric(na.omit(df_full[,2])))-min(as.numeric(na.omit(df_full[,2])))

#ICA_1
summary(df_full[,5])
sd(df_full[,5])
mm_ICA1<-max(df_full[,5])-min(df_full[,5])

##ICA_01
summary(df_full[,2])
sd(df_full[,2])
mm_ICA01<-max(df_full[,2])-min(df_full[,2])

##ICA_05
summary(df_full[,4])
sd(df_full[,4])
mm_ICA05<-max(df_full[,4])-min(df_full[,4])
hist(df_full[,4], main=bquote('Histogram of'~ICA[05]), xlab = "ICA_05")


##ICA_125
summary(df_full[,8])
sd(df_full[,8])
mm_ICA125<-max(df_full[,8])-min(df_full[,8])
hist(df_full[,8], main=bquote('Histogram of'~ICA[1.25]), xlab = "ICA_1.25")


##ICA_15
summary(df_full[,10])
sd(df_full[,10])
mm_ICA15<-max(df_full[10])-min(df_full[,10])
hist(df_full[,10], main=bquote('Histogram of'~ICA[1.5]), xlab = "ICA_1.5")


##ICA_2
summary(df_full[,12])
sd(df_full[,12])
mm_ICA2<-max(df_full[12])-min(df_full[,12])
hist(df_full[,12], main=bquote('Histogram of'~ICA[2]), xlab = "ICA_2")


##ICA_25  ##INTEGRAL
summary(df_full[,13])
sd(df_full[,13])
mm_ICA25<-max(df_full[13])-min(df_full[,13])
hist(df_full[,13], main=bquote('Histogram of'~ICA[2.5]), xlab = "ICA_2.5")


##ICA_3 ##INTEGRAL
summary(df_full[,15])
sd(df_full[,15])
mm_ICA3<-max(df_full[15])-min(df_full[,15])
hist(df_full[,15], main=bquote('Histogram of'~ICA[3]), xlab = "ICA_3")


##max-min vs alpha
mm_ICAS<- c(mm_ICA01, mm_ICA05, mm_ICA1, mm_ICA125, mm_ICA15, mm_ICA2, mm_ICA25, mm_ICA3)
alphas<- c(0.1, 0.5, 1, 1.25, 1.5, 2, 2.5, 3)
plot(alphas,mm_ICAS, ylab="max-min", xlab = expression(ICA[2.5]))

