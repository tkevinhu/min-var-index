###############################################################
# Project: Minimum Variance Portfolio
# Editors: Xinyuan Liu, Kevin Hu, Quinn Song, Tina Behroozi
# Date: November 30th 2019
###############################################################

###########################################
# Table of Contents
# 1. Data Preparation
# 2. Log-Return Calculation
# 3. Minimum Weight Optimization
# 4. ARIMA Model
# 5. Copula Model
# 6. Garch Model
###########################################

###########################################
## 1: DATA PREPARATION
###########################################

###########################################
## 1.1: Install Packages
###########################################
#install.packages("zoo")
#install.packages("devtools")
#install_github("cvxgrp/CVXR")
#install.packages("CVXR")
#install.packages("tidyverse")
#install.packages("tseries")
#install.packages("forecast")
#install.packages("Metrics")
#install.packages("copula")
#install.packages("readxl")
#install.packages("MASS")
#install.packages("forecast")
#install.packages('rugarch')
#install.packages('Metrics')


###########################################
## 1.2: Call Libraries
###########################################
library(lubridate)
library(zoo)
library(devtools)
library(CVXR)
library(tidyverse)
library(tseries)
library(forecast)
library(Metrics)
library(copula)
library("readxl")
library("MASS")
library("forecast")
library('rugarch')
library('Metrics')

###########################################
## 1.3: Read Raw Dataset
###########################################
rawdata <- read.csv('C:/Users/hukevin/Desktop/Apps/energy_data.csv')

###########################################
## 1.4: Imputation
###########################################
#  Carry the lastest observation forward to missing values
data<-na.locf(rawdata)
#  Use the mean of first six months to impute the first three days data
data[,3:14]<-as.numeric(unlist(data[,3:14]))
data$gas[1]<-mean(data$gas[4:120])
data$gas[2]<-mean(data$gas[4:120])
data$gas[3]<-mean(data$gas[4:120])
#  Check the missing consective data in xdr (our base asset)
missingxdr<-rawdata[which(is.na(rawdata$xdr)),]
data$Date<-as.Date(data$Date,format="%Y-%m-%d")

###########################################
## 1.5: Plot a few data for observation
###########################################
plot(x=data$Date,y=data$xdr)
plot(x=data$Date,y=data$cny)

###########################################
## 1.6: Convert all data from USD to XDR
###########################################
multiply <- function(x) 
{ 
  data[,x]*data[,3]
}
data[,4]<-data[,4]*data[,3]
data[,5]<-multiply(5)
data[,6]<-multiply(6)
data[,7]<-multiply(7)
data[,8]<-multiply(8)
data[,9]<-multiply(9)
data[,10]<-multiply(10)
data[,11]<-multiply(11)
data[,12]<-multiply(12)
data[,13]<-multiply(13)
data[,14]<-multiply(14)

###########################################
## 2: LOG-RETURN CALCULATION
###########################################

###########################################
## 2.1: Calculate Log Return
###########################################
log_returns<-data[-1,2]
# take a look at the log_returns2 which is the xdr
log_returns2<- diff(log(data[,3]), lag=1)
plot(log_returns2)
log_returns3<- diff(log(data[,4]), lag=1)
log_returns4<- diff(log(data[,5]), lag=1)
log_returns5<- diff(log(data[,6]), lag=1)
log_returns6<- diff(log(data[,7]), lag=1)
log_returns7<- diff(log(data[,8]), lag=1)
log_returns8<- diff(log(data[,9]), lag=1)
log_returns9<- diff(log(data[,10]), lag=1)
log_returns10<- diff(log(data[,11]), lag=1)
log_returns11<- diff(log(data[,12]), lag=1)
log_returns12<- diff(log(data[,13]), lag=1)
log_returns13<- diff(log(data[,14]), lag=1)
# Combine log returns data
newdata<-cbind(log_returns,log_returns2,log_returns3,log_returns4,
               log_returns5,log_returns6,log_returns7,log_returns8,
               log_returns9,log_returns10,log_returns11,log_returns12,
               log_returns13)
newdata<-as.data.frame(newdata)
newdata$log_returns<-data[-1,2]
# Assign columns name
colnames(newdata,do.NULL=FALSE)
colnames(newdata)<-c("DATE","XDR","EUR","JPY","CNY","BTC","ETH",
                     "RIP","MON","OIL","GAS","NAT","URA")

###########################################
## 2.2: Categorize data to each quarters
###########################################
newdata$Quarter <- ifelse(
  newdata$DATE< as.Date('2016-12-31'),0,
  ifelse(
    newdata$DATE < as.Date('2017-03-31'),1,
    ifelse(newdata$DATE< as.Date('2017-06-30'),2, 
           ifelse(newdata$DATE < as.Date('2017-09-30'),3,
                  ifelse(newdata$DATE < as.Date('2017-12-31'),4,
                         ifelse(newdata$DATE < as.Date('2018-03-31'),5,
                                ifelse(newdata$DATE< as.Date('2018-06-30'),6,
                                       ifelse(newdata$DATE < as.Date('2018-09-30'),7,
                                              ifelse(newdata$DATE < as.Date('2018-12-31'),8,
                                                     ifelse(newdata$DATE< as.Date('2019-03-31'),9,
                                                            ifelse(newdata$DATE < as.Date('2019-06-30'),10,11)))))))))))

###########################################
## 3: Minimum Weight Optimization
###########################################

###########################################
## 3.1: Minimum Variance Weight Calculation
###########################################
# Initialize empty variance, covariance and weights
LogRetCov <- matrix(,12,1212)
LogRetVar <- matrix(,12,11)
optimw    <- matrix(,12,11)
# Loop the optimization model
for (i in 1:11)
{
  z <- 12*i+1
  k <- 12*(i+1)
  # calculate covariance matrix
  LogRetCov[,z:k]<- cov(subset(newdata[,2:13],newdata$Quarter==i|newdata$Quarter==i-1))
  # calculate the variance
  LogRetVar_1<-sapply(subset(newdata[,2:13],newdata$Quarter==i|newdata$Quarter==i-1),var)
  dim(LogRetVar_1)<-c(12,1)
  LogRetVar[,i]<-LogRetVar_1
  # quadration programming for weight optimization
  w <- Variable(12)
  objective <- Minimize(LogRetVar[,i] %*% w + quad_form(w,LogRetCov[,z:k]))
  constraints <- list(sum(w) == 1, w >= 0.01, w <= 0.25)
  prob <- Problem(objective, constraints)
  # output the weight
  result <- solve(prob)
  optimw[,i]<- result$getValue(w)
}

###########################################
## 3.2: Tidying Up the Results
###########################################
# Tidy-up the results
colnames(optimw,do.NULL=FALSE)
colnames(optimw)<-c("Q1","Q2","Q3","Q4","Q5","Q6","Q7",
                    "Q8","Q9","Q10","Q11")
# Label quarters on the results
colnames(data)[colnames(data) == "Date"] <- "DATE"
data$Quarter <- ifelse(data$DATE < as.Date('2017-03-31'),0,
                       ifelse(data$DATE< as.Date('2017-06-30'),1, 
                              ifelse(data$DATE < as.Date('2017-09-30'),2,
                                     ifelse(data$DATE < as.Date('2017-12-31'),3,
                                            ifelse(data$DATE < as.Date('2018-03-31'),4,
                                                   ifelse(data$DATE< as.Date('2018-06-30'),5,
                                                          ifelse(data$DATE < as.Date('2018-09-30'),6,
                                                                 ifelse(data$DATE < as.Date('2018-12-31'),7,
                                                                        ifelse(data$DATE< as.Date('2019-03-31'),8,
                                                                               ifelse(data$DATE < as.Date('2019-06-30'),9,10))))))))))

###########################################
## 3.3: Calculating Quarter-End Results
###########################################
data$Quarter <- ifelse(data$DATE < as.Date('2017-03-31'),0,
                       ifelse(data$DATE< as.Date('2017-06-30'),1, 
                              ifelse(data$DATE < as.Date('2017-09-30'),2,
                                     ifelse(data$DATE < as.Date('2017-12-31'),3,
                                            ifelse(data$DATE < as.Date('2018-03-31'),4,
                                                   ifelse(data$DATE< as.Date('2018-06-30'),5,
                                                          ifelse(data$DATE < as.Date('2018-09-30'),6,
                                                                 ifelse(data$DATE < as.Date('2018-12-31'),7,
                                                                        ifelse(data$DATE< as.Date('2019-03-31'),8,
                                                                               ifelse(data$DATE < as.Date('2019-06-30'),9,
                                                                                      ifelse(data$DATE==as.Date('2019-09-30'),10,9)))))))))))
data <- data[order(data$Quarter,data$DATE, decreasing=FALSE),]
data1<-subset(data,data$Quarter!=0)

#calculate the weighted index 
data$Quarter <- ifelse(
  data$DATE< as.Date('2016-12-31'),0,
  ifelse(
    data$DATE < as.Date('2017-03-31'),1,
    ifelse(data$DATE< as.Date('2017-06-30'),2, 
           ifelse(data$DATE < as.Date('2017-09-30'),3,
                  ifelse(data$DATE < as.Date('2017-12-31'),4,
                         ifelse(data$DATE < as.Date('2018-03-31'),5,
                                ifelse(data$DATE< as.Date('2018-06-30'),6,
                                       ifelse(data$DATE < as.Date('2018-09-30'),7,
                                              ifelse(data$DATE < as.Date('2018-12-31'),8,
                                                     ifelse(data$DATE< as.Date('2019-03-31'),9,
                                                            ifelse(data$DATE < as.Date('2019-06-30'),10,
                                                                   ifelse(data$DATE==as.Date('2019-09-30'),11,10))))))))))))

#quarter_end price                                                                               
data <- data[order(data$Quarter,data$DATE, decreasing=TRUE),]
qe_price<-data[!duplicated(data$Quarter),]
#quarter_end_index
qe_price<-qe_price[order(qe_price$Quarter,decreasing=FALSE),]
matrix1<-qe_price[-1,-c(1,2,15)]
matrix2<-optimw
matrix3<-as.matrix(matrix1)%*%as.matrix(matrix2)
qe_index<-cbind(matrix3[1,1],matrix3[2,2],matrix3[3,3],matrix3[4,4],
                matrix3[5,5],matrix3[6,6],matrix3[7,7],matrix3[8,8],
                matrix3[9,9],matrix3[10,10],matrix3[11,11])

###########################################
## 3.4: Calculate Realized Index
###########################################
qe_index2 <- as.data.frame(lapply(qe_index, rep,12))
colnames(qe_index2) <- c("Q1","Q2","Q3","Q4","Q5","Q6","Q7","Q8","Q9","Q10","Q11")
matrix1t <- data.frame(t(matrix1))
c <- qe_index2/matrix1t
# Append weights for each quater
append<-cbind(data1,t(c[,data1[,15]]))
index<-rowSums(append[,3:14]*append[,16:27])
index2<-c(rep(index[1],length(index)))
# Index Calculation
index<-index/index2*100
cbind(append[,c(2,15)],index)
# Plot the index
plot(index)
# Label the index with dates
index3 <- cbind(append[,c(2,15)],index)
plot(index3[,1],index3[,3])
#correlation matrix of assets value
corr<-newdata[,c(-1,-14)]
source("http://www.sthda.com/upload/rquery_cormat.r")
corrplot<-rquery.cormat(corr)

###########################################
## 4: ARIMA MODEL
###########################################

###########################################
## 4.1: Fitting ARIMA Model
###########################################
cat("\014")  
rm(list = ls())
graphics.off()

# in this case, we have 9 quarters
index3 <- subset(index3,index3$Quarter<=10)
ARIMARMSE <- vector()  
ARIMAMAPE <-vector()
ARIMAMSE <-vector()
for (i in 1:9)
{
  IndexRealized <- index3 # Rename the realized index dataframe
  TrainIndexi   <- subset(IndexRealized, IndexRealized$Quarter <= i) # Our current training set
  TrainIndexi   <- TrainIndexi$index # Strip out date columns
  TestIndexi    <- subset(IndexRealized, IndexRealized$Quarter == i+1) # Our current testing set
  TestIndexi    <- TestIndexi$index # Strip out date columns
  ARIMAi        <- auto.arima(TrainIndexi
                              ,max.p      = 10
                              ,max.q      = 5
                              ,ic         = "aic") # Fit ARIMA Model
  show(ARIMAi)
  Forecasti <- forecast(ARIMAi, h = length(TestIndexi),bootstrap = T)
  meanboot  = Forecasti$mean
  RMSEi         <- rmse(meanboot, TestIndexi)
  ARIMARMSE     <- cbind(ARIMARMSE, RMSEi)
  MAPEi         <- mape(meanboot, TestIndexi)
  ARIMAMAPE     <- cbind(ARIMAMAPE, MAPEi)
  MSEi          <- mse(meanboot, TestIndexi)
  ARIMAMSE     <- cbind(ARIMAMSE, MSEi)
}
ARIMARMSE
ARIMAMAPE
ARIMAMSE
measure <- rbind(ARIMARMSE,ARIMAMAPE,ARIMAMSE)
rownames(measure)<-c('rmse','mape','mse')
measuret <- t(measure)


###########################################
## 4.2: Backtesting ARIMA Model
###########################################
TrainIndex   <- subset(index3, index3$Quarter <= 7) # Our current training set
TrainIndex   <- TrainIndex$index # Strip out date columns
TestIndex    <- subset(index3, index3$Quarter > 7) # Our current testing set
TestIndex    <- TestIndex$index # Strip out date columns
ARIMA        <- auto.arima(TrainIndex
                           ,max.p      = 10
                           ,max.q      = 5
                           ,ic         = "aic") # Fit ARIMA Model
Forecasti <- forecast(ARIMA, h = length(TestIndex)+10,bootstrap = T)
meanboot  = Forecasti$mean
upperboot = Forecasti$upper[,1]
lowerboot = Forecasti$lower[,1]
xdr = Index2$xdr
btc = Index2$btc
gas = Index2$gas
q = as.character(index3$Quarter)
q[1] = '17/1'
q[277] = '18/4'
q[642] = '19/8'
a =  c(1,which.max(index3$Quarter==2),which.max(index3$Quarter==3),277,which.max(index3$Quarter==5),
       which.max(index3$Quarter==6),which.max(index3$Quarter==7),642,which.max(index3$Quarter==9),which.max(index3$Quarter==10))
plot(index3$index, xaxt = "n",xlab='year/quarter',ylab='Index',main='ARIMA model')
axis(1, at=a, labels=q[a])
lines(c(fitted(Forecasti),upperboot),col='blue')
lines(c(fitted(Forecasti),lowerboot),col='green')
lines(c(fitted(Forecasti),meanboot),col='red')

plot(c(fitted(Forecasti),meanboot), xaxt = "n",xlab='year/quarter',ylab='Index',main='ARIMA model',col='black')
axis(1, at=a, labels=q[a])
lines(gas,col='red', lty=2)
lines(xdr,col='blue', lty=3)
lines(btc,col='green', lty=4)
legend("topleft", legend=c('index',"gas", "xdr",'btc'),
       col=c('black',"red", "blue",'green'), lty=1:4, cex=0.8)

##ACF and PACF
kpss.test(TrainIndex)
kpss.test(diff(TrainIndex))
kpss.test(diff(diff(TrainIndex)))
Acf(diff(diff(TrainIndex)))
Pacf(diff(diff(TrainIndex)))
MAresult = arma(diff(diff(TrainIndex)),order = c(0, 1))
MAresult = as.numeric(na.omit(MAresult$fitted.values))
Acf(MAresult)
Pacf(MAresult)


###########################################
## 5: Copula Simulation of Index
###########################################

# Simulates index using copula:
# We use previous month's data to create empirical marginal cdfs 
# and a copula which we can combine to simulate more data.
# Together with the weights obtained through our previous optimization,
# we can simulate the index for the next month
library(copula)
library(sn)
#library(ks)
set.seed(999)
source('FINAL DATA IMPUTE AND CLEANUP CODE.R')#this is index3

# container for simulated data
sim_data_guass_cop = data.frame(array(NA,c(sum(my_data_3_xdr$quarters>0),14)))
names(sim_data_guass_cop) = names(my_data_3_xdr)
# populate with correct quarters
sim_data_guass_cop$quarters=my_data_3_xdr$quarters[my_data_3_xdr$quarters>0]
# populate with correct dates
sim_data_guass_cop$Date=my_data_3_xdr$Date[my_data_3_xdr$quarters>0]

for(i in 1:11){
  # convert log diffs to empirical cdf of log diffs
  # get stationary values
  x = my_data_4_log_diff[my_data_4_log_diff$quarters==i-1,c(-1,-dim(my_data_4_log_diff)[2])]
  #print(head(u))
  #print(dim(u)[1])
  # get u's which are empiracle observations of log diff distribution
  u = x
  for(j in 1:dim(u)[2]){
    u[,j] = pobs(x[,j])
  }
  u = as.matrix(u)
  # now fit copula
  ## Inverting Kendall's tau
  fit.tau <- fitCopula(normalCopula(dim=12, dispstr="un"), u, method="itau")
  # get number of log diffs we need to simulate
  n = sum(my_data_3_xdr$quarters==i)
  # simulate
  sim_log_diffs = rCopula(n, fit.tau@copula)
  #print(n)
  # change from uniform back to draws from log diff distributions
  # through using quantiles (cdf approximation)
  for(j in 1:dim(sim_log_diffs)[1]){
    for(k in 1:dim(sim_log_diffs)[2]){
      sim_log_diffs[j,k] = quantile(x[,k],sim_log_diffs[j,k])
    }
  }
  #print(sim_log_diffs[1:3])
  # index of last date in training set (need this to go back from stationary to (simulated) price)
  ind = which(max(my_data_3_xdr$Date[my_data_3_xdr$quarters==i-1])==my_data_3_xdr$Date)
  #print(ind)
  # now get back simulated price! complicated forumla
  # initialize return matrix
  sim_prices = sim_log_diffs
  for(j in 1:dim(sim_log_diffs)[2]){
    sim_prices[,j] = exp(cumsum(sim_log_diffs[,j])+log(my_data_3_xdr[ind,1+j]))
  }
  # add back to matrix!!!
  sim_data_guass_cop[sim_data_guass_cop$quarters==i,c(-1,-dim(sim_data_guass_cop)[2])] = sim_prices
}
#head(sim_prices)
# get weights
weights = read.csv('realizedweights.csv')
# add to model
cop_sim_with_weights = cbind(sim_data_guass_cop,
                             data.frame(t(weights[,-1][,sim_data_guass_cop$quarters])))
# dont let rownames be weird
rownames(cop_sim_with_weights) = 1:(dim(cop_sim_with_weights)[1])
# check
#cop_sim_with_weights[80:100,]


#ADD COP INDEX
cop_sim_with_weights['cop_pred_index'] = 1

#rowSums(cop_sim_with_weights[,2:(2+11)] * cop_sim_with_weights[,15:(15+11)])

for(i in 1:11){
  last_date_prev_quarter = max(my_data_3_xdr$Date[my_data_3_xdr$quarters==i-1])
  # ind is last date before current quarter from original data
  ind = (which(my_data_3_xdr$Date==last_date_prev_quarter))
  # ind2 is last date before current quarter from cop simulated data
  ind2 = (which(cop_sim_with_weights$Date==last_date_prev_quarter))
  if(length(ind2)==0){
    last_val = 100
  } else{
    last_val = cop_sim_with_weights$cop_pred_index[ind2]
  }
  # make cop_index
  # count how many elements in this month
  n = sum(cop_sim_with_weights$quarters==i)
  # get the subset of the data frame with only this month
  temp = cop_sim_with_weights[cop_sim_with_weights$quarters==i,]
  # index value!!! calculated with formula
  cop_sim_with_weights$cop_pred_index[cop_sim_with_weights$quarters==i] =
    last_val*rowSums(temp[,2:(2+11)] * temp[,15:(15+11)] /
                       (as.matrix(array(rep(1,n),c(n,1)))%*%
                          as.matrix(array(my_data_3_xdr[ind,2:13],c(1,12)))))
}

# plot comparing cop index to bitcoin
btc = my_data_3_xdr[my_data_3_xdr$quarters>0,]$btc
btc = (btc-mean(btc))/sqrt(var(btc))
plot(cop_sim_with_weights$Date,btc,col="blue",type='l')
points(cop_sim_with_weights$Date,(cop_sim_with_weights$cop_pred_index-mean(cop_sim_with_weights$cop_pred_index))/
         sqrt(var(cop_sim_with_weights$cop_pred_index)),col='red',type='l')

points(cop_sim_with_weights$Date,10000*(0.5-c(0,diff(cop_sim_with_weights$quarters))),col="black",type='o')
# var of the index
var(cop_sim_with_weights$cop_pred_index)


# plot comparing bitcoin price simulated by copula to actual price
plot(cop_sim_with_weights$Date,cop_sim_with_weights$btc,col='red',type='l')
points(cop_sim_with_weights$Date,my_data_3_xdr[my_data_3_xdr$quarters>0,]$btc,col="blue",type='l')
points(cop_sim_with_weights$Date,10000*(-c(0,diff(cop_sim_with_weights$quarters))),col="black",type='o')




###########################################
## 5: Garch Model
###########################################
cat("\014")  
rm(list = ls())
graphics.off()

index3 <- subset(Index3,Index3$Quarter<=10)
q = as.character(index3$Quarter)
q[1] = '17/1'
q[277] = '18/4'
q[642] = '19/8'
a =  c(1,which.max(index3$Quarter==2),which.max(index3$Quarter==3),277,which.max(index3$Quarter==5),
       which.max(index3$Quarter==6),which.max(index3$Quarter==7),642,which.max(index3$Quarter==9),which.max(index3$Quarter==10))

#index property
Data = index3$index
kpss.test(Data)
log.Data = log(Data)
kpss.test(log.Data)
sqrt.Data = sqrt(Data)
kpss.test(sqrt.Data)
boxcox(Data~1)
bc = boxcox(Data~1, lambda = seq(1.5, 2, 1/100), interp = FALSE)
lamba = 1.95
y = (Data^(lamba)-1)/lamba
kpss.test(y)
diff.Data = diff(Data)
kpss.test(diff.Data)
diff.diff.Data = diff(diff(Data))
kpss.test(diff.diff.Data)
par(mfrow=c(1,2))
plot(Data, xaxt = "n",xlab='year/quarter',ylab='Index')
axis(1, at=a, labels=q[a])
plot(diff.diff.Data, xaxt = "n",xlab='year/quarter',ylab='transform Index')
axis(1, at=a, labels=q[a])

par(mfrow=c(1,1))
boxplot(diff.diff.Data, main = "boxplot of diff(diff(Index))")
plot(density(diff.diff.Data), main = "kernel density plot of diff(diff(Index))")
tdist_par = fitdistr(diff.diff.Data, "t", start = list(m=mean(diff.diff.Data),s=sd(diff.diff.Data), df=3), lower=c(-1, 0.001,1))
tdist_ind = tdist_par$estimate[2]*rt(length(diff.diff.Data),df=tdist_par$estimate[3]) + tdist_par$estimate[1]
par(mfrow=c(1,2))
qqnorm(diff.diff.Data,datax = TRUE, main = "normal Q-Q plot")
qqplot(diff.diff.Data, tdist_ind, main = "tdist Q-Q plot")

#fitting model
GARCHRMSE <- vector()  
GARCHMAPE <-vector()
GARCHMSE <-vector()
APARCHRMSE <- vector() 
APARCHMAPE <-vector()
APARCHMSE <-vector()
for (i in 1:2)
{
  IndexRealized <- index3 # Rename the realized index dataframe
  TrainIndexi   <- subset(IndexRealized, IndexRealized$Quarter <= i) # Our current training set
  TrainIndexi   <- TrainIndexi$index # Strip out date columns
  TestIndexi    <- subset(IndexRealized, IndexRealized$Quarter == i+1) # Our current testing set
  TestIndexi    <- TestIndexi$index # Strip out date columns
  ARIMAi        <- auto.arima(TrainIndexi
                              ,max.p      = 10
                              ,max.q      = 5
                              ,trace      = TRUE
                              ,ic         = "aic") # Fit ARIMA Model
  Armavar = ARIMAi$arma
  if(Armavar[6]==1){
    TrainIndexi = diff(TrainIndexi)
    TestIndexi = diff(TestIndexi)
  } else if(Armavar[6]==2){
    TrainIndexi = diff(diff(TrainIndexi))
    TestIndexi = diff(diff(TestIndexi))
  } else if(Armavar[6]==3){
    TrainIndexi = diff(diff(diff(TrainIndexi)))
    TestIndexi = diff(diff(diff(TestIndexi)))
  }
  specgarch = ugarchspec(mean.model = list(armaOrder = c(Armavar[1],Armavar[2]), include.mean = TRUE), distribution.model = "std")
  fit = ugarchfit(specgarch, data = TrainIndexi)
  bootpred = ugarchboot(fit, method = "Partial", n.ahead = length(TestIndexi), n.bootpred = 1000)
  u = bootpred@forc@forecast
  meanboot = u$seriesFor
  
  specaparch = ugarchspec(mean.model = list(armaOrder = c(Armavar[1],Armavar[2]), include.mean = TRUE),
                          variance.model = list(model = "apARCH"), distribution.model = "std")
  fit = ugarchfit(specaparch, data = TrainIndexi)
  show(fit)
  bootpredaparch = ugarchboot(fit, method = "Partial", n.ahead = length(TestIndexi), n.bootpred = 1000)
  u = bootpredaparch@forc@forecast
  meanbootaparch = u$seriesFor
  
  RMSEi         <- rmse(meanboot, TestIndexi)
  GARCHRMSE     <- cbind(GARCHRMSE, RMSEi)
  RMSEi         <- rmse(meanbootaparch, TestIndexi)
  APARCHRMSE    <- cbind(APARCHRMSE, RMSEi)
  MAPEi         <- mape(meanboot, TestIndexi)
  GARCHMAPE     <- cbind(GARCHMAPE, MAPEi)
  MAPEi         <- mape(meanbootaparch, TestIndexi)
  APARCHMAPE     <- cbind(APARCHMAPE, MAPEi)
  MSEi          <- mse(meanboot, TestIndexi)
  GARCHMSE     <- cbind(GARCHMSE, MSEi)
  MSEi          <- mse(meanbootaparch, TestIndexi)
  APARCHMSE     <- cbind(APARCHMSE, MSEi)
}
GARCHRMSE
APARCHRMSE
GARCHMAPE
APARCHMAPE
GARCHMSE
APARCHMSE

#backtesting: put seven quarters as historical data
TrainIndex   <- subset(index3, IndexRealized$Quarter <= 5) 
TrainIndex   <- TrainIndex$index
TestIndex    <- subset(index3, IndexRealized$Quarter > 5) 
TestIndex    <- TestIndex$index 
specgarch = ugarchspec(mean.model = list(armaOrder = c(3,0), include.mean = TRUE), distribution.model = "std")
fit = ugarchfit(specgarch, data = (diff(TrainIndex)))
bootpred = ugarchboot(fit, method = "Partial", n.ahead = (length(TestIndex)-1),
                      n.bootpred = 1000)
u = bootpred@forc@forecast
meanboot = u$seriesFor
zr = t(as.data.frame(bootpred, which = "series", type = "summary"))
upperboot = zr[,2]
lowerboot = zr[,4]
fitdata = ((fit@fit$fitted.values))
indexp = diffinv(c(diff((TrainIndex)),meanboot),xi = 100)
indexpu = diffinv(c(diff((TrainIndex)),upperboot),xi = 100)
indexpl = diffinv(c(diff((TrainIndex)),lowerboot),xi = 100)
par(mfrow=c(1,1))
plot(diff(index3$index), xaxt = "n",xlab='year/quarter',ylab='Index',main='GARCH model for diff(index)')
axis(1, at=a, labels=q[a])
lines(c(fitdata,meanboot),col='red')
lines(c(fitdata,upperboot),col='blue')
lines(c(fitdata,lowerboot),col='green')
plot((index3$index), xaxt = "n",xlab='year/quarter',ylab='Index',main='GARCH model for index')
axis(1, at=a, labels=q[a])
lines(indexp,col='red')
lines(indexpu,col='blue')
lines(indexpl,col='green')
plot(fit,which=8)
Res <- residuals(fit, standardize = TRUE)
skewness(Res) #left skew
kurtosis(Res) #thick tail
shapiro.test(as.numeric(Res))

specaparch = ugarchspec(mean.model = list(armaOrder = c(3,0), include.mean = TRUE),
                        variance.model = list(model = "apARCH"), distribution.model = "std")
fit = ugarchfit(specaparch, data = (diff(TrainIndex)))
bootpredaparch = ugarchboot(fit, method = "Partial", n.ahead = (length(TestIndex)-1), n.bootpred = 1000)
u = bootpredaparch@forc@forecast
meanbootaparch = u$seriesFor
zr = t(as.data.frame(bootpredaparch, which = "series", type = "summary"))
upperbootaparch = zr[,2]
lowerbootaparch = zr[,4]
fitdata = fit@fit$fitted.values
indexpaparch = diffinv((c(diff((TrainIndex)),meanbootaparch)),xi = 100)
indexpuaparch = diffinv((c(diff((TrainIndex)),upperbootaparch)),xi = 100)
indexplaparch = diffinv((c(diff((TrainIndex)),lowerbootaparch)),xi = 100)
plot(diff(index3$index), xaxt = "n",xlab='year/quarter',ylab='Index',main='APARCH model for diff(index)')
axis(1, at=a, labels=q[a])
lines(c(fitdata,meanbootaparch),col='red')
lines(c(fitdata,upperbootaparch),col='blue')
lines(c(fitdata,lowerbootaparch),col='green')
plot((index3$index), xaxt = "n",xlab='year/quarter',ylab='Index',main='APARCH model for index')
axis(1, at=a, labels=q[a])
lines(indexpaparch,col='red')
lines(indexpuaparch,col='blue')
lines(indexplaparch,col='green')
plot(fit,which=8)
Res <- residuals(fit, standardize = TRUE)
skewness(Res) #left skew
kurtosis(Res) #thick tail
shapiro.test(as.numeric(Res))

