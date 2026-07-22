library(abind)
library(CovTools)
library(ks)
library(Riemann)
library(shapes)
library(StatPerMeCo)

load("C:/Users/andre/OneDrive - Università degli Studi di Macerata/Articolo Zhang Palma/Data_50.RData")

ytest = rcovariance[(window+1):nrow(rcovariance),]

a1 = read.csv2('C:/Users/andre/OneDrive - Università degli Studi di Macerata/Articolo Zhang Palma/1.Approx50.csv')
a2 = read.csv2('C:/Users/andre/OneDrive - Università degli Studi di Macerata/Articolo Zhang Palma/2.Cholesky50_50PC.csv')
a3 = read.csv2('C:/Users/andre/OneDrive - Università degli Studi di Macerata/Articolo Zhang Palma/4.DCC50.csv')
a3b = read.csv2('C:/Users/andre/OneDrive - Università degli Studi di Macerata/Articolo Zhang Palma/7.DCCNL50_v2.csv')
a4 = read.csv('C:/Users/andre/OneDrive - Università degli Studi di Macerata/Articolo Zhang Palma/5.SPDNetM50.csv', header = T)
a5 = read.csv('C:/Users/andre/OneDrive - Università degli Studi di Macerata/Articolo Zhang Palma/5.SPDNetF50.csv')
a6 = read.csv('C:/Users/andre/OneDrive - Università degli Studi di Macerata/Articolo Zhang Palma/5.SPDNetLE50.csv')
a7 = read.csv('C:/Users/andre/OneDrive - Università degli Studi di Macerata/Articolo Zhang Palma/5.SPDNetM50_3lag.csv')
a8 = read.csv('C:/Users/andre/OneDrive - Università degli Studi di Macerata/Articolo Zhang Palma/5.SPDNetF50_3lag.csv')
a9 = read.csv('C:/Users/andre/OneDrive - Università degli Studi di Macerata/Articolo Zhang Palma/5.SPDNetLE50_3lag.csv')
a10 = read.csv('C:/Users/andre/OneDrive - Università degli Studi di Macerata/Articolo Zhang Palma/5.SPDNetM50_5lag.csv')
a11 = read.csv('C:/Users/andre/OneDrive - Università degli Studi di Macerata/Articolo Zhang Palma/5.SPDNetF50_5lag.csv')
a12 = read.csv('C:/Users/andre/OneDrive - Università degli Studi di Macerata/Articolo Zhang Palma/5.SPDNetLE50_5lag.csv')
a13 = read.csv('C:/Users/andre/OneDrive - Università degli Studi di Macerata/Articolo Zhang Palma/5.SPDNetM50_10lag.csv')
a14 = read.csv('C:/Users/andre/OneDrive - Università degli Studi di Macerata/Articolo Zhang Palma/5.SPDNetF50_10lag.csv')
a15 = read.csv('C:/Users/andre/OneDrive - Università degli Studi di Macerata/Articolo Zhang Palma/5.SPDNetLE50_10lag.csv')
a16 = read.csv('C:/Users/andre/OneDrive - Università degli Studi di Macerata/Articolo Zhang Palma/5.SPDNetMHAR.csv')
a17 = read.csv('C:/Users/andre/OneDrive - Università degli Studi di Macerata/Articolo Zhang Palma/5.SPDNetFHAR.csv')
a18 = read.csv('C:/Users/andre/OneDrive - Università degli Studi di Macerata/Articolo Zhang Palma/5.SPDNetLEHAR.csv')
a19 = read.csv('C:/Users/andre/OneDrive - Università degli Studi di Macerata/Articolo Zhang Palma/5.SPDNetLEHARProc_v2.csv')
a20 = read.csv('C:/Users/andre/OneDrive - Università degli Studi di Macerata/Articolo Zhang Palma/5.SPDNetMHARProc_v2.csv')
a21 = read.csv('C:/Users/andre/OneDrive - Università degli Studi di Macerata/Articolo Zhang Palma/6.NNetMHARle.csv')
a22 = read.csv('C:/Users/andre/OneDrive - Università degli Studi di Macerata/Articolo Zhang Palma/6.NNetMHARP_v2.csv')
a23 = read.csv('C:/Users/andre/OneDrive - Università degli Studi di Macerata/Articolo Zhang Palma/7.GEONETle.csv')
a24 = read.csv('C:/Users/andre/OneDrive - Università degli Studi di Macerata/Articolo Zhang Palma/7.GEONETP_v2.csv')
a25 = as.data.frame(rcovariance[(window):(nrow(rcovariance)-1),])

predlist = list(a1, a2, a3, a3b, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13, a14, a15, a16, a17, a18, a19, a20, a21, a22, a23, a24, a25)

npred = nrow(a5)
nassets = 50
####Computing losses####
euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))

approximant <- function(A){
  B <- (A+t(A))/2
  C <- (A-t(A))/2
  BB <- tensr::polar(B)
  H <- BB$Z
  speB <- eigen(B)
  L <- speB$values
  Z <- speB$vectors
  L[L<0] <- 0
  Xf <- Z %*% diag(L) %*% t(Z)
  return(Xf)
}

LERM <- function(A, B){
  Aspd = as.matrix(Matrix::nearPD(A)$mat)
  Bspd = as.matrix(Matrix::nearPD(B)$mat)
  lerm = norm(expm::logm(Aspd) - expm::logm(Bspd))
  return(lerm)
}

LERM2 <- function(A, B){
  Aspd = as.matrix(Matrix::nearPD(A)$mat)
  Bspd = as.matrix(Matrix::nearPD(B)$mat)
  lerm = shapes::distcov(Aspd, Bspd, method = 'LogEuclidean')
  return(lerm)
}

Stein <- function(Sigma, Sigmahat){
  Sigmaspd = as.matrix(Matrix::nearPD(Sigma)$mat)
  Sigmahatspd = as.matrix(Matrix::nearPD(Sigmahat)$mat)
  stein = matrix.trace(Sigmahatspd %*% MASS::ginv(Sigmaspd)) - log(det(Sigmahatspd %*% MASS::ginv(Sigmaspd))) - nrow(Sigmaspd)
  return(stein)
}

mcov = list()
m_list <- vector("list", length(predlist))
names(m_list) <- paste0("m", seq_along(predlist))
lossF = matrix(NA, ncol = length(predlist), nrow = npred)
lossE = matrix(NA, ncol = length(predlist), nrow = npred)
lossP = matrix(NA, ncol = length(predlist), nrow = npred)
lossLE = matrix(NA, ncol = length(predlist), nrow = npred)
lossSt = matrix(NA, ncol = length(predlist), nrow = npred)

for (i in 1:npred){
  mcov[[i]] =  invvech(as.matrix(ytest[i,]))
  for (k in seq_along(predlist)){
    m_list[[k]][[i]] <- invvech(as.matrix(predlist[[k]][i, ]))
    lossF[i,k] = Frobenius(mcov[[i]], m_list[[k]][[i]])
    lossE[i,k] = euc.dist(t(as.matrix(ytest[i,])),as.matrix(predlist[[k]][i,]))
    lossP[i,k] = procOPA(mcov[[i]], m_list[[k]][[i]], scale = F)$rmsd
    lossLE[i,k] = LERM2(mcov[[i]], m_list[[k]][[i]])
    lossSt[i,k] = Stein(mcov[[i]], m_list[[k]][[i]])
  }
  print(i)
}
###Frobenius####
colnames(lossF) <- c('Approx', 'Cholesky', 'GOGARCH', 'DCCNL', 'SPDNetM','SPDNetF','SPDNetLE','SPDNetM3l','SPDNetF3l',
                     'SPDNetLE3l','SPDNetM5l','SPDNetF5l','SPDNetLE5l','SPDNetM10l','SPDNetF10l',
                     'SPDNetLE10l','GeoHARMLE','GeoHARF','GeoHARLELE','GeoHARLEProc', 'GeoHARMProc',
                     'MLPHARLE', 'MLPHARP', 'GEOREGHARLE', 'GEOREGHARP','RW')
colMeans(lossF[,c('SPDNetLE','SPDNetLE3l','SPDNetLE5l','SPDNetLE10l',
               'SPDNetM','SPDNetM3l','SPDNetM5l','SPDNetM10l',
               'GeoHARLELE','GeoHARLEProc',
               'GeoHARMLE', 'GeoHARMProc',
               'MLPHARLE', 'MLPHARP', 
               'GEOREGHARLE', 'GEOREGHARP',
               'Cholesky','GOGARCH', 'DCCNL', 'RW')])
###Euclidean
colnames(lossE) <- c('Approx', 'Cholesky', 'GOGARCH', 'DCCNL', 'SPDNetM','SPDNetF','SPDNetLE','SPDNetM3l','SPDNetF3l',
                     'SPDNetLE3l','SPDNetM5l','SPDNetF5l','SPDNetLE5l','SPDNetM10l','SPDNetF10l',
                     'SPDNetLE10l','GeoHARMLE','GeoHARF','GeoHARLELE','GeoHARLEProc', 'GeoHARMProc',
                     'MLPHARLE', 'MLPHARP', 'GEOREGHARLE', 'GEOREGHARP','RW')
colMeans(lossE[,c('SPDNetLE','SPDNetLE3l','SPDNetLE5l','SPDNetLE10l',
                  'SPDNetM','SPDNetM3l','SPDNetM5l','SPDNetM10l',
                  'GeoHARLELE','GeoHARLEProc',
                  'GeoHARMLE', 'GeoHARMProc',
                  'MLPHARLE', 'MLPHARP', 
                  'GEOREGHARLE', 'GEOREGHARP',
                  'Cholesky','GOGARCH', 'DCCNL', 'RW')])

###Stein
colnames(lossSt) <- c('Approx', 'Cholesky', 'GOGARCH', 'DCCNL', 'SPDNetM','SPDNetF','SPDNetLE','SPDNetM3l','SPDNetF3l',
                      'SPDNetLE3l','SPDNetM5l','SPDNetF5l','SPDNetLE5l','SPDNetM10l','SPDNetF10l',
                      'SPDNetLE10l','GeoHARMLE','GeoHARF','GeoHARLELE','GeoHARLEProc', 'GeoHARMProc',
                      'MLPHARLE', 'MLPHARP', 'GEOREGHARLE', 'GEOREGHARP','RW')
colMeans(lossSt[,c('SPDNetLE','SPDNetLE3l','SPDNetLE5l','SPDNetLE10l',
                   'SPDNetM','SPDNetM3l','SPDNetM5l','SPDNetM10l',
                   'GeoHARLELE','GeoHARLEProc',
                   'GeoHARMLE', 'GeoHARMProc',
                   'MLPHARLE', 'MLPHARP', 
                   'GEOREGHARLE', 'GEOREGHARP',
                   'Cholesky','GOGARCH', 'DCCNL', 'RW')])

###Procrustes
colnames(lossP) <- c('Approx', 'Cholesky', 'GOGARCH', 'DCCNL', 'SPDNetM','SPDNetF','SPDNetLE','SPDNetM3l','SPDNetF3l',
                     'SPDNetLE3l','SPDNetM5l','SPDNetF5l','SPDNetLE5l','SPDNetM10l','SPDNetF10l',
                     'SPDNetLE10l','GeoHARMLE','GeoHARF','GeoHARLELE','GeoHARLEProc', 'GeoHARMProc',
                     'MLPHARLE', 'MLPHARP', 'GEOREGHARLE', 'GEOREGHARP','RW')
colMeans(lossP[,c('SPDNetLE','SPDNetLE3l','SPDNetLE5l','SPDNetLE10l',
                  'SPDNetM','SPDNetM3l','SPDNetM5l','SPDNetM10l',
                  'GeoHARLELE','GeoHARLEProc',
                  'GeoHARMLE', 'GeoHARMProc',
                  'MLPHARLE', 'MLPHARP', 
                  'GEOREGHARLE', 'GEOREGHARP',
                  'Cholesky','GOGARCH', 'DCCNL', 'RW')])

###LogE
colnames(lossLE) <- c('Approx', 'Cholesky', 'GOGARCH', 'DCCNL', 'SPDNetM','SPDNetF','SPDNetLE','SPDNetM3l','SPDNetF3l',
                      'SPDNetLE3l','SPDNetM5l','SPDNetF5l','SPDNetLE5l','SPDNetM10l','SPDNetF10l',
                      'SPDNetLE10l','GeoHARMLE','GeoHARF','GeoHARLELE','GeoHARLEProc', 'GeoHARMProc',
                      'MLPHARLE', 'MLPHARP', 'GEOREGHARLE', 'GEOREGHARP','RW')
colMeans(lossLE[,c('SPDNetLE','SPDNetLE3l','SPDNetLE5l','SPDNetLE10l',
                   'SPDNetM','SPDNetM3l','SPDNetM5l','SPDNetM10l',
                   'GeoHARLELE','GeoHARLEProc',
                   'GeoHARMLE', 'GeoHARMProc',
                   'MLPHARLE', 'MLPHARP', 
                   'GEOREGHARLE', 'GEOREGHARP',
                   'Cholesky','GOGARCH', 'DCCNL', 'RW')])



####MCS####
source('C:/Users/andre/OneDrive - Università Politecnica delle Marche/Articolo VLSTAR/Dati New/Predictions_Final/MCS1.R')
set.seed(123)
lossf = lossF[,c('SPDNetLE','SPDNetLE3l','SPDNetLE5l','SPDNetLE10l',
                 'SPDNetM','SPDNetM3l','SPDNetM5l','SPDNetM10l',
                 'GeoHARLELE','GeoHARLEProc',
                 'GeoHARMLE', 'GeoHARMProc',
                 'MLPHARLE', 'MLPHARP', 
                 'GEOREGHARLE', 'GEOREGHARP',
                 'Cholesky','GOGARCH', 'DCCNL', 'RW')]
MCSF1 = MCS1(lossf, alpha = 0.1,B=10000,statistic='TR')

losse = lossE[,c('SPDNetLE','SPDNetLE3l','SPDNetLE5l','SPDNetLE10l',
                 'SPDNetM','SPDNetM3l','SPDNetM5l','SPDNetM10l',
                 'GeoHARLELE','GeoHARLEProc',
                 'GeoHARMLE', 'GeoHARMProc',
                 'MLPHARLE', 'MLPHARP', 
                 'GEOREGHARLE', 'GEOREGHARP',
                 'Cholesky','GOGARCH', 'DCCNL', 'RW')]
MCSE1 = MCS1(losse, alpha = 0.1,B=10000,statistic='TR')

lossst = lossSt[,c('SPDNetLE','SPDNetLE3l','SPDNetLE5l','SPDNetLE10l',
                   'SPDNetM','SPDNetM3l','SPDNetM5l','SPDNetM10l',
                   'GeoHARLELE','GeoHARLEProc',
                   'GeoHARMLE', 'GeoHARMProc',
                   'MLPHARLE', 'MLPHARP', 
                   'GEOREGHARLE', 'GEOREGHARP',
                   'Cholesky','GOGARCH', 'DCCNL', 'RW')]
MCSST1 = MCS1(lossst, alpha = 0.1,B=10000,statistic='TR')

lossproc = lossP[,c('SPDNetLE','SPDNetLE3l','SPDNetLE5l','SPDNetLE10l',
                    'SPDNetM','SPDNetM3l','SPDNetM5l','SPDNetM10l',
                    'GeoHARLELE','GeoHARLEProc',
                    'GeoHARMLE', 'GeoHARMProc',
                    'MLPHARLE', 'MLPHARP', 
                    'GEOREGHARLE', 'GEOREGHARP',
                    'Cholesky','GOGARCH', 'DCCNL', 'RW')]
MCSP1 = MCS1(lossproc, alpha = 0.1,B=10000,statistic='TR')

lossle = lossLE[,c('SPDNetLE','SPDNetLE3l','SPDNetLE5l','SPDNetLE10l',
                   'SPDNetM','SPDNetM3l','SPDNetM5l','SPDNetM10l',
                   'GeoHARLELE','GeoHARLEProc',
                   'GeoHARMLE', 'GeoHARMProc',
                   'MLPHARLE', 'MLPHARP', 
                   'GEOREGHARLE', 'GEOREGHARP',
                   'Cholesky','GOGARCH', 'DCCNL', 'RW')]
MCSLE1 = MCS1(lossle, alpha = 0.1,B=10000,statistic='TR')



####Robustness###
library(quantmod)
library(xlsx)
library(highfrequency)
library(xts)
library(ggplot2)
start_date <- as.Date("2007-06-26")
end_date <- as.Date("2021-07-02")

data2017 = read.csv2('SPXUSD2017.csv', header = F)
data2018 = read.csv2('SPXUSD2018.csv', header = F)
data2019 = read.csv2('SPXUSD2019.csv', header = F)
data2020 = read.csv2('SPXUSD2020.csv', header = F)
data2021 = read.csv2('SPXUSD2021.csv', header = F)

data = rbind(data2017, data2018, data2019, data2020, data2021)
colnames(data) = c('date','open', 'high', 'low', 'close', 'adj')
data$date = as.POSIXct(data$date, format = "%Y-%m-%d %H:%M")

df = data[,c('date', 'close')]
df.xts = xts(df$close, order.by = df$date)
rv <- rRVar(df.xts, makeReturns = T)

Date1 <- tail(row.names(dairet), npred)

SP500test = rv[as.Date(index(rv)) %in% Date1,]

thresh = quantile(SP500test, 0.90)

SP500test = as.data.frame(SP500test)
SP500test$regime[SP500test$V1 < thresh] = 'Low'
SP500test$regime[SP500test$V1 >= thresh] = 'High'
SP500test$regime[SP500test$V1 < thresh] = 'Low'

plot(index(SP500test), SP500test$V1, type = "l", col = "black", ylim = range(SP500test$V1), ylab = "y2", xlab = "")
for (j in 1:nrow(SP500test)) {
  if (SP500test$regime[j] == 'High') {
    points(Date1[j],  SP500test$V1[j], pch = 2, col = "blue")
  } else if (SP500test$regime[j] == 'Low') {
    points(Date1[j],  SP500test$V1[j], pch = 15, col = "red")
  }
}

####Low regime###
set.seed(123)
#Frobenius
source('C:/Users/andre/OneDrive - Università Politecnica delle Marche/Articolo VLSTAR/Dati New/Predictions_Final/MCS1.R')
lossf_low = lossf[SP500test$regime == 'Low',]
round(colMeans(lossf_low),3)
MCSFlow = MCS1(lossf_low, alpha = 0.1,B=10000,statistic='TR')

#Euclidean
losse_low = losse[SP500test$regime == 'Low',]
round(colMeans(losse_low),3)
MCSElow = MCS1(losse_low, alpha = 0.1,B=10000,statistic='TR')

#Stein
lossst_low = lossst[SP500test$regime == 'Low',]
round(colMeans(lossst_low),3)
MCSSTlow = MCS1(lossst_low, alpha = 0.1,B=10000,statistic='TR')

#Procrustes
lossp_low = lossproc[SP500test$regime == 'Low',]
round(colMeans(lossp_low),3)
MCSPlow = MCS1(lossp_low, alpha = 0.1,B=10000,statistic='TR')

#Log-Euclidean
lossle_low = lossle[SP500test$regime == 'Low',]
round(colMeans(lossle_low),3)
MCSLElow = MCS1(lossle_low, alpha = 0.1,B=10000,statistic='TR')


####High regime####
#Frobenius
source('C:/Users/andre/OneDrive - Università Politecnica delle Marche/Articolo VLSTAR/Dati New/Predictions_Final/MCS1.R')
set.seed(123)
lossf_high = lossf[SP500test$regime == 'High',]
round(colMeans(lossf_high),3)
MCSFhigh = MCS1(lossf_high, alpha = 0.1,B=10000,statistic='TR')

#Euclidean
losse_high = losse[SP500test$regime == 'High',]
round(colMeans(losse_high),3)
MCSEhigh = MCS1(losse_high, alpha = 0.1,B=10000,statistic='TR')

#Stein
lossst_high = lossst[SP500test$regime == 'High',]
round(colMeans(lossst_high),3)
MCSSThigh = MCS1(lossst_high, alpha = 0.1,B=10000,statistic='TR')

#Procrustes
lossp_high = lossproc[SP500test$regime == 'High',]
round(colMeans(lossp_high),3)
MCSPhigh = MCS1(lossp_high, alpha = 0.1,B=10000,statistic='TR')

#Log-Euclidean
lossle_high = lossle[SP500test$regime == 'High',]
round(colMeans(lossle_high),3)
MCSLEhigh = MCS1(lossle_high, alpha = 0.1,B=10000,statistic='TR')



####Portfolio optimization####
library(data.table)
library(dplyr)
library(ks)
library(lessR)
library(lubridate)
library(Matrix)
library(matrixcalc)
library(pracma)
library(quantmod)
library(reshape2)
library(rmgarch)
library(rockchalk)
library(rugarch)
library(shapes)
library(starvars)
library(StatPerMeCo)
library(tensr)
library(tidyr)
library(vars)
require(nloptr)
library(quadprog)

dairet_cc = read.csv('C:/Users/andre/OneDrive - Università degli Studi di Macerata/Articolo Zhang Palma/Data Chao/SP100_return.csv')
row.names(dairet_cc) = dairet_cc[,1]
dairet_cc = dairet_cc[,-1]
dairet_cc = dairet_cc[,colnames(dairet_cc) %in% stock]

ref_dates = as.Date(row.names(dairet_cc))

dairet_oc = data.frame(matrix(NA_real_, nrow = length(ref_dates), ncol = length(stock)))
colnames(dairet_oc) = stock
rownames(dairet_oc) = as.character(ref_dates)

missing_report = list()

for (k in 1:length(stock)) {
  tmp = getSymbols(stock[k], from = '2007-06-27', to = "2021-07-02",
                   warnings = FALSE, auto.assign = FALSE)
  
  # open-to-close returns in %
  r_oc = 100 * (log(as.numeric(tmp[, 4])) - log(as.numeric(tmp[, 1])))
  names(r_oc) = as.character(index(tmp))
  
  common_dates = intersect(as.character(ref_dates), names(r_oc))
  missing_report[[stock[k]]] = setdiff(as.character(ref_dates), names(r_oc))
  
  dairet_oc[common_dates, k] = r_oc[common_dates]
}

dairet = dairet_oc
nassets = ncol(dairet)


GMV <- function(Sigma){
  D_matrix <- 2 * as.matrix(nearPD(Sigma)$mat)
  n <- nrow(Sigma)
  d_vector <- rep(0, n)
  A_matrix <- cbind(rep(1, n), diag(n))
  b_vector <- c(1, rep(0, n))
  # use solve.QP to minimize portfolio variance
  quad_prog <- solve.QP(Dmat = D_matrix, dvec = d_vector, Amat = A_matrix, bvec = b_vector, meq = 1)
  return(quad_prog$solution)
}
mcov <- list()

dairet <- as.matrix(dairet[(window+1):nrow(dairet),])
dairetp = dairet/100
pesina <- rep((1/nassets),nassets)
por_ret = matrix(0, ncol = ncol(lossF), nrow = npred)
por_retlong = matrix(0, ncol = ncol(lossF), nrow = npred)
naive <- NULL
turnover = matrix(0, ncol = ncol(lossF), nrow = npred)
turnoverlong = matrix(0, ncol = ncol(lossF), nrow = npred)
turnoverna = rep(0, npred)
turnoverlongna = rep(0, npred)
iota = as.matrix(rep(1, nassets))

failed_calls = list()  

for (i in 1:ntest){
  for (k in seq_along(predlist)){
    tryCatch({
      Sigma_hat = Matrix::nearPD(m_list[[k]][[i]])$mat  
      Sigma_inv = solve(Sigma_hat)
      
      pesilong = GMV(m_list[[k]][[i]])
      pesi = as.matrix(Sigma_inv %*% iota) / as.numeric(t(iota) %*% Sigma_inv %*% iota)
      
      por_ret[i,k] = t(dairet[i,]) %*% pesi
      por_retlong[i,k] = t(dairet[i,]) %*% pesilong
      
      if (i > 1){
        for (j in 1:nassets){
          turnover[i,k] = turnover[i,k] + abs(pesi[j] - prev_pesi[j] * (1 + dairetp[i-1, j]) /
                                                as.numeric(1 + t(prev_pesi) %*% dairetp[i-1,]))
          turnoverlong[i,k] = turnoverlong[i,k] + abs(pesilong[j] - prev_pesilong[j] * (1 + dairetp[i-1, j]) /
                                                        as.numeric(1 + t(prev_pesilong) %*% dairetp[i-1,]))
        }
      }
      prev_pesi <<- pesi
      prev_pesilong <<- pesilong
    }, error = function(e){
      failed_calls[[length(failed_calls)+1]] <<- list(i = i, k = k, error = conditionMessage(e))
    })
  }
  
  naive[i] <- t(dairet[i,]) %*% pesina
  if (i > 1){
    for (j in 1:nassets){
      turnoverna[i] = turnoverna[i] + abs(pesina[j] - pesina[j] * (1 + dairetp[i-1, j]) /
                                            as.numeric(1 + t(pesina) %*% dairetp[i-1,]))
      turnoverlongna[i] = turnoverlongna[i] + abs(pesina[j] - pesina[j] * (1 + dairetp[i-1, j]) /
                                                    as.numeric(1 + t(pesina) %*% dairetp[i-1,]))
    }
  }
  print(i)
}

if (length(failed_calls) > 0) {
  cat(sprintf("%d failed calls (i,k) during portfolio optimization:\n", length(failed_calls)))
  print(do.call(rbind, lapply(failed_calls, as.data.frame)))
}

por_return = as.data.frame(cbind(naive, por_ret))
colnames(por_return) <- c('Naive', 'Approx', 'Cholesky', 'DCC', 'DCCNL', 'SPDNetM','SPDNetF','SPDNetLE','SPDNetM3l','SPDNetF3l',
                          'SPDNetLE3l','SPDNetM5l','SPDNetF5l','SPDNetLE5l','SPDNetM10l','SPDNetF10l',
                          'SPDNetLE10l','GeoHARM','GeoHARF','GeoHARLE','GeoHARLEProc', 'GeoHARMProc',
                          'MLPLE', 'MLPP', 'GEONETLE', 'GEONETP','RW')
porreturn = por_return[,c('SPDNetLE','SPDNetLE3l','SPDNetLE5l','SPDNetLE10l',
                          'SPDNetM','SPDNetM3l','SPDNetM5l','SPDNetM10l',
                          'GeoHARLE','GeoHARLEProc',
                          'GeoHARM', 'GeoHARMProc',
                          'MLPLE', 'MLPP', 'GEONETLE', 'GEONETP',
                          'Cholesky','DCC', 'DCCNL', 'RW', 'Naive')]
porreturnlong = as.data.frame(cbind(naive, por_retlong))
colnames(porreturnlong) <- c('Naive', 'Approx', 'Cholesky', 'DCC', 'DCCNL', 'SPDNetM','SPDNetF','SPDNetLE','SPDNetM3l','SPDNetF3l',
                             'SPDNetLE3l','SPDNetM5l','SPDNetF5l','SPDNetLE5l','SPDNetM10l','SPDNetF10l',
                             'SPDNetLE10l','GeoHARM','GeoHARF','GeoHARLE','GeoHARLEProc', 'GeoHARMProc',
                             'MLPLE', 'MLPP', 'GEONETLE', 'GEONETP','RW')
porreturnlong = porreturnlong[,c('SPDNetLE','SPDNetLE3l','SPDNetLE5l','SPDNetLE10l',
                                 'SPDNetM','SPDNetM3l','SPDNetM5l','SPDNetM10l',
                                 'GeoHARLE','GeoHARLEProc',
                                 'GeoHARM', 'GeoHARMProc',
                                 'MLPLE', 'MLPP', 'GEONETLE', 'GEONETP',
                                 'Cholesky','DCC', 'DCCNL', 'RW', 'Naive')]
turnover = as.data.frame(cbind(turnoverna, turnover))
colnames(turnover) <- c('Naive', 'Approx', 'Cholesky', 'DCC', 'DCCNL', 'SPDNetM','SPDNetF','SPDNetLE','SPDNetM3l','SPDNetF3l',
                        'SPDNetLE3l','SPDNetM5l','SPDNetF5l','SPDNetLE5l','SPDNetM10l','SPDNetF10l',
                        'SPDNetLE10l','GeoHARM','GeoHARF','GeoHARLE','GeoHARLEProc', 'GeoHARMProc',
                        'MLPLE', 'MLPP', 'GEONETLE', 'GEONETP','RW')
turnover = turnover[,c('SPDNetLE','SPDNetLE3l','SPDNetLE5l','SPDNetLE10l',
                       'SPDNetM','SPDNetM3l','SPDNetM5l','SPDNetM10l',
                       'GeoHARLE','GeoHARLEProc',
                       'GeoHARM', 'GeoHARMProc',
                       'MLPLE', 'MLPP', 'GEONETLE', 'GEONETP',
                       'Cholesky','DCC', 'DCCNL', 'RW', 'Naive')]
turnoverlong = as.data.frame(cbind(turnoverlongna, turnoverlong))
colnames(turnoverlong) <- c('Naive', 'Approx', 'Cholesky', 'DCC', 'DCCNL', 'SPDNetM','SPDNetF','SPDNetLE','SPDNetM3l','SPDNetF3l',
                            'SPDNetLE3l','SPDNetM5l','SPDNetF5l','SPDNetLE5l','SPDNetM10l','SPDNetF10l',
                            'SPDNetLE10l','GeoHARM','GeoHARF','GeoHARLE','GeoHARLEProc', 'GeoHARMProc',
                            'MLPLE', 'MLPP', 'GEONETLE', 'GEONETP','RW')
turnoverlong = turnoverlong[,c('SPDNetLE','SPDNetLE3l','SPDNetLE5l','SPDNetLE10l',
                               'SPDNetM','SPDNetM3l','SPDNetM5l','SPDNetM10l',
                               'GeoHARLE','GeoHARLEProc',
                               'GeoHARM', 'GeoHARMProc',
                               'MLPLE', 'MLPP', 'GEONETLE', 'GEONETP',
                               'Cholesky','DCC', 'DCCNL', 'RW', 'Naive')]
round(colMeans(turnover),3)
round(colMeans(turnoverlong),3)

annsd = function(x){
  return(sd(x)*sqrt(252))
}
round(apply(porreturn, 2, annsd),3)
round(apply(porreturnlong, 2, annsd),3)


####Return loss####

lossreturns = function(returns, turnover, c1){
  nmodel = ncol(returns)
  ntest = nrow(returns)
  rnet = (1 + returns) * (1 - c1 * turnover) - 1
  colnames(rnet) = colnames(returns)
  rew = rnet$Naive
  muew = mean(rew)*252
  sigmaew = sd(rew)*sqrt(252)
  returnloss = rep(NA, nmodel)
  for(i in 1:(nmodel)){
    mui = mean(rnet[,i])*252
    sigmai = sd(rnet[,i])*sqrt(252)
    returnloss[i] = (muew/sigmaew)*sigmai - mui
  }
  names(returnloss) = colnames(returns)
  return(returnloss)
}
c1 = 0.005
###GMV
round(lossreturns(porreturn, turnover, c1), 3)

###GMV+
round(lossreturns(porreturnlong, turnoverlong, c1), 3)
