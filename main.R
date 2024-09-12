library('forecast')
library('tidyverse')
library(stats)
library('latex2exp')
library('tseries')
library("rugarch")
library(rmgarch)
library(qqplotr)
#1 
df <- read_csv("denmark_covid_logrets.csv")

model <- arima(df$x, order=c(6,0,5), include.mean = FALSE)

model.resid <- resid(model)
residuals.2 <- model.resid^2

#2

acf(model.resid, lag=20, main="ACF Plot of Residuals")
acf(model.resid, lag = 20, type=c('partial'), main="PACF Plot of Residuals")
acf(residuals.2, lag=20, main="ACF Plot of  Squared Residuals")

Box.test(model.resid, lag = 20, type=c('Ljung-Box'))


#3 

garch.models <- NULL
GFSpec=ugarchspec(variance.model = list(model="sGARCH",garchOrder=c(1,1)),
                  mean.model = list(armaOrder=c(0,0),include.mean=TRUE))
GF=ugarchfit(spec=GFSpec, data = model.resid,
             solver="gosolnp",solver.control = list(tol=1e-7))
GF@fit$coef
-2*-103.9613 +4*log(length(model.resid))


GFSpec=ugarchspec(variance.model = list(model="sGARCH",garchOrder=c(6,8)),
                  mean.model = list(armaOrder=c(0,0),include.mean=TRUE))
GF=ugarchfit(spec=GFSpec, data = model.resid,
             solver="gosolnp",solver.control = list(tol=1e-7))
infocriteria(GF)[2]*nrow(df)

Garch = function(df){
  tryCatch(GF=ugarchfit(spec=GFSpec, data = model.resid,
                        solver="gosolnp",solver.control = list(tol=1e-7)),
           error = function(e) { print("inf or 0 eigen");
             NULL
           } )
}

m <- matrix(NULL*10,10,10)


for (i in 1:10){
  for (j in 1:10){
    
    GFSpec=ugarchspec(variance.model = list(model="sGARCH",garchOrder=c(i,j)),
                      mean.model = list(armaOrder=c(0,0),include.mean=TRUE))
    
    GF <- Garch(model.resid)
    
    if(is.null(GF)){
      BIC <- 0
    } else {
      BIC <- -2*GF@fit$LLH + (i+j+2*log(nrow(df)))
    }
    
    m[i,j] <- BIC
    
  }
}

GFSpec=ugarchspec(variance.model = list(model="sGARCH",garchOrder=c(1,1)),
                  mean.model = list(armaOrder=c(0,0),include.mean=TRUE))

GF=ugarchfit(spec=GFSpec, data = df,
             solver="gosolnp",solver.control = list(tol=1e-7))

h <- Garch(df)

m <- matrix(NULL*10,10,10)


for (i in 1:10){
  for (j in 1:10){
    
    GFSpec=ugarchspec(variance.model = list(model="sGARCH",garchOrder=c(i,j)),
                      mean.model = list(armaOrder=c(0,0),include.mean=TRUE))
    
    GF=ugarchfit(spec=GFSpec, data = model.resid,
                 solver="gosolnp",solver.control = list(tol=1e-7))
    
    
    BIC <- -2*GF@fit$LLH + (i+j+2*log(nrow(df)))
    m[i,j] <- BIC
    
  }
}

m[is.na(m)] <- 5000

c <-  which(m == min(m), arr.ind = TRUE)
m[c] #best is therefore garch(2,3)

GFSpec=ugarchspec(variance.model = list(model="sGARCH",garchOrder=c(5,5)),
                  mean.model = list(armaOrder=c(0,0),include.mean=TRUE))

GF=ugarchfit(spec=GFSpec, data = model.resid,
             solver="gosolnp",solver.control = list(tol=1e-7))





# 4

#including the mean is the offset term that should be included
GFSpec=ugarchspec(variance.model = list(model="sGARCH",garchOrder=c(2,3)),
                  mean.model = list(armaOrder=c(0,0),include.mean=TRUE), distribution.model = "norm")

GF=ugarchfit(spec=GFSpec, data = model.resid,
             solver="gosolnp",solver.control = list(tol=1e-7))
resids <- GF@fit$residuals
resids2 <- (GF@fit$residuals)^2


data.frame(resids) %>% ggplot(aes(resids)) +
  geom_density()

kurtosis(resids)


shapiro.test(resids)

Box.test(resids2, type = c('Ljung'))
acf(resids, lag = 20, type = c('correlation'), main="Autocorrelation Plot of GARCH(2,3) Residuals")
acf(resids2, lag = 20, type = c('correlation'), main="Autocorrelation Plot of GARCH(2,3) Squared Residuals")

pacf(resids, lag = 20)
pacf(resids2, lag = 20)

resids.table <- tibble(resids)

library('ggplot2')

data.frame(resids) %>% ggplot(mapping = aes(sample = resids)) + 
  stat_qq_point(size = 2, colour="Royalblue")  + 
  stat_qq_line(colour="grey") + 
  xlab("Theoretical Quantiles") + 
  ylab("Sample Quantiles") 


# 6

GFSpec=ugarchspec(variance.model = list(model="sGARCH",garchOrder=c(2,3)),
                  mean.model = list(armaOrder=c(0,0),include.mean=TRUE), distribution.model = "std")

GF=ugarchfit(spec=GFSpec, data = model.resid,
             solver="gosolnp",solver.control = list(tol=1e-7))
resids <- GF@fit$residuals

resids2 <- resid(GF@fit)^2
acf(resids, lag = 20, type = c('correlation'))
acf(resids2, lag = 20, type = c('correlation'))



pacf(resids, lag = 20)

resids.table <- tibble(resids)

library('ggplot2')

tibble(resids) %>% ggplot(mapping = aes(sample = resids)) + 
  stat_qq_point( colour="Royalblue", distribution = "t.scaled")  + 
  stat_qq_line(colour="grey", distribution = "t.scaled") + 
  xlab("Theoretical Quantiles") + 
  ylab("Sample Quantiles") 


params <- as.list(MASS::fitdistr(resids, "t")$estimate)
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
#> Warning: NaNs produced
tibble(resids) %>% ggplot(aes(sample = resids$resids)) +
  stat_qq(distribution = qt, dparams = params["df"]) +
  stat_qq_line(distribution = qt, dparams = params["df"])


####### Part 2 ########

#assumptions
psi <- 0.9827
sigma2 <- 4.718

#1

testdata <- read_csv("test_data.csv")
St <- colnames(testdata)
St <- as.numeric(St)
St <- append(St, testdata$`24.2876668779738`)
testdata <- as.data.frame(St)
view(testdata)
plot(S$St)

testdata %>% ggplot(aes(x = 1:nrow(testdata),y = St)) + 
  geom_line() + 
  labs(x = "Time", y = "Values") +
  ggtitle('Test Data')

lag.s <- lag(testdata$St, 1)
testdata['S_(t-1)'] <- lag.s

Z <- rnorm(143-2, mean=0, sd=sqrt(sigma2))

S <- NULL
Si <- df

#S_T <- function(S_T_1){
#  return (psi*S_T_1 + rnorm(1,mean=0, sd=Z.sigma2))
#}


#average.preds <- 0
#for (n in seq(2,142)){
#  S <- NULL
#  for (j in seq(1,1000)){
#    S <- append(S, testdata[n,'S_(t-1)'])
#  }
#  average.preds <- append(average.preds, mean(S))
#}

#S <- NULL
#for (n in seq(2,142)){
#  for(i in seq(1,1000)){
#    Z <- rnorm(1,mean=0, sd=sigma2)
#    S <- append(S, psi*testdata[n, 'S_(t-1)'] + Z)
#  }
#}

S.sim <- testdata[2:142,]
for (n in seq(1,1000)){
  Z <- rnorm(141, mean = 0, sd = sqrt(sigma2))
  S.sim[n] <- psi* S.sim$`S_(t-1)` + Z
}


forecastAverage <- function(i,df,predsstart){
  lst <- NULL
  #df <- as.data.frame(df)
  for (j in seq(predsstart,ncol(df))){
    lst <- append(lst, df[i,predsstart])
  }
  return(mean(lst))
}

averages <- NULL

for(i in seq(1,141)){
  averages <- append(averages, forecastAverage(i,S.sim, predsstart = 3))
}

averages <- data.frame(averages)

averages['S_(t-1)'] <- S.sim$`S_(t-1)`

averages %>%  ggplot(aes(x = 1:141, y = averages)) + 
  geom_line(aes(colour="Forecasted"), size = 1, alpha = 0.6) + 
  geom_line(aes(x = 1:141, y=`S_(t-1)`, colour = "True Values"), alpha=0.4 , size = 1) +
  labs(x = "Time", y = TeX("$S_i$")) +
  ggtitle(TeX("$S_i$ vs $\\hat{S}_{i}$"))

getAllPredictions <- function(df,i,predsstart){
  lst <- NULL
  for (j in seq(predsstart, ncol(df))){
    lst <- append(lst, df[i,j])
  }
  return(lst)
}

table.of.predictions <- getAllPredictions(S.sim,1,predsstart = 3)
table.of.predictions <- data.frame(table.of.predictions)
for (i in seq(2,141)){
  table.of.predictions[glue("{i}")] <- getAllPredictions(S.sim, i, predsstart = 3)
}


getQuantiles <- function(lst.of.quantiles,df,i){
  lst <- NULL
  for (quant in lst.of.quantiles){
    # this is wrong, getting multiple quantiles
    lst <- append(lst,quantile(df[,i], quant))
  }
  return(lst)
}

lst <- NULL
for (i in seq(1,ncol(table.of.predictions))){
  lst <- append(lst, getQuantiles(c(0.025, 0.975), table.of.predictions, i))
}

quantile.25 <- NULL
for (i in seq(1,282,2)){
  quantile.25 <- append(quantile.25,lst[i])
}

quantile.75 <- NULL
for (i in seq(2,282,2)){
  quantile.75 <- append(quantile.75,lst[i])
}

averages %>% ggplot(aes(x = 1:141, y = averages)) + 
  geom_line(aes(colour="Forecasted")) + 
  geom_line(aes(x = 1:141, y =`S_(t-1)`, colour = "True Values"), alpha = 0.7) +
  geom_line(aes(x= 1:141, y = quantile.25 , colour = "2.5% Quantile"), 
            size = 0.5,
            linetype = "dashed") +
  geom_line(aes(x = 1:141, y = quantile.75, colour = "97.5% Quantile"),
            size = 0.5,
            linetype = "dashed") + 
 # geom_point(aes(x = 1:141, y = quantile.25, colour = "2.5% Quantile")) +
 # geom_point(aes(x = 1:141, y = quantile.75, colour = "97.5% Quantile")) +
  labs(x = "Time", y = TeX("$S_i$")) +
  ggtitle(TeX("$S_i$ vs $\\hat{S}_{i}$ 1-Step-Ahead Forecast")) 

for (i in seq(1,nrow(table.ofpredictions))){
  lst <- append(lst, getQuantiles(0.025, table.of.predictions, i))
}

a <- getQuantiles(c(0.025, 0.975), table.of.predictions, 1)
getQuantiles(c(0.025,0.975), table.of.predictions, 2)


# 2. a

MSE <- function(h,forecast,truevalues){
  s <- NULL
  for(i in seq(1,length(truevalues)-h)){
    s <- append(s, (truevalues[i] - forecast[i])^2)
  }
  s <- sum(s)
  return((1/length(truevalues) - h)*s)
}

#s <- NULL
#for (i in seq(1:length(testdata))){
#  s <- append(s, (testdata$St[] - ))
#}

(1/(nrow(averages)-1))*sum((averages$`S_(t-1)` - averages$averages)^2)

MSE(1,forecast = averages, truevalues = testdata$St[2:length(testdata)])

getDifferential <- function(vector1, i){
  difference <- vector1[i] - vector1[i-1]
  if(difference<0){
    return('-')
  } else if (difference == 0){
    return('0')
  } else {
    return('+')
  }
}



M11 <- function(predicted,true){
  s <- 0
  for(i in seq(2,length(predicted))){
    if((getDifferential(predicted, i) == "+") & (getDifferential(true,i) == "+")){
      s <- s + 1
    } 
  }
  return(s)
}

m11 <- M11(averages$averages, averages$`S_(t-1)`)

M21 <- function(predicted, true){
  s <- 0
  for(i in seq(2,length(predicted))){
    if((getDifferential(predicted,i)=="+") & (getDifferential(true,i) == "-")){
      s <- s + 1
    }
  }
  return(s)
}

m21 <- M21(averages$averages, averages$`S_(t-1)`)

M12 <- function(predicted, true){
  s <- 0
  for(i in seq(2,length(predicted))){
    if((getDifferential(predicted,i)=="-") & (getDifferential(true,i) == "+")){
      s <- s + 1
    }
  }
  return(s)
}

m12 <- M12(averages$averages, averages$`S_(t-1)`)

M22 <- function(predicted, true){
  s <- 0
  for(i in seq(2,length(predicted))){
    if((getDifferential(predicted,i)=="-") & (getDifferential(true,i) == "-")){
      s <- s + 1
    }
  }
  return(s)
}
m22 <- M22(averages$averages, averages$`S_(t-1)`)

m20 <- m21 + m22
m10 <- m11 + m12
m01 <- m11 + m21
m02 <- m12 + m22

contigency.table <- data.frame(cbind(c(m11,m21),c(m12,m22)))
contigency.table <- contigency.table %>% rename('Up'='X1', 'Down'='X2')
contigency.table

createDummyvariables <- function(vector){
  lst <- NULL
  for (i in seq(2,length(vector))){
    if(getDifferential(vector,i)=="+"){
      lst <- append(lst,1)
    } else {
      lst <- append(lst,0)
    }
  }
  return(lst)
}

predictions <- createDummyvariables(averages$averages)
true <- createDummyvariables(averages$`S_(t-1)`)

chisq.test(predictions, true)


# 4


model <- arima(testdata$St, order = c(1,0,0), include.mean = FALSE)
pred <- predict(model, n.ahead = 10)


data.frame(pred) %>% ggplot(aes(x = 1:141, y = testdata$St)) + 
  geom_line(aes(colour="True Values")) + 
  geom_line(aes(x = 1:141, y =`S_(t-1)`, colour = "Forecasted")) +
  geom_line(aes(x= 1:141, y = quantile.25 , colour = "2.5% Quantile"), 
            size = 0.5,
            linetype = "dashed") +
  geom_line(aes(x = 1:141, y = quantile.75, colour = "97.5% Quantile"),
            size = 0.5,
            linetype = "dashed") + 
  # geom_point(aes(x = 1:141, y = quantile.25, colour = "2.5% Quantile")) +
  # geom_point(aes(x = 1:141, y = quantile.75, colour = "97.5% Quantile")) +
  labs(x = "Time", y = TeX("$S_i$")) +
  ggtitle(TeX("$S_i$ vs $\\hat{S}_{i} 1-Step-Ahead Forecast$")) 


#####

multstepForecast <- function(h,timeseries,i,psi,deviation){
  z <- rnorm(1,mean=0,sd = deviation)
  s <- timeseries[i]
  for(k in seq(1,h)){
    z <- rnorm(1, mean=0, sd = deviation)
    s <- psi*s + z
  }
  return(s)
}

for(i in seq(10:length(testdata))){
  for (n in seq(1,1000)){
    
  }
}

multistepahead.preds <- NULL
for (i in seq(1:length(testdata))){
  for (n in seq(2,1001)){
    j <- n+1
    multistepahead.preds <- append(multistepahead.preds,
                                  multstepForecast(10,
                                                    timeseries = testdata$St,
                                                    psi = psi, 
                                                    i = i,
                                                    deviation = sigma2 ))
  }
}




multistep <- function(h,series,start, psi, sigma2){
  s <- (psi)*series[start-9] + rnorm(1,mean=0, sd = sqrt(sigma2))
  for (j in seq(0,h)){
    s <- (psi)*s + rnorm(1,mean = 0, sd = sqrt(sigma2))
  }
  return(s)
}


series <- testdata$St

multistep(10,series=series,start=10,psi = psi, sigma2 = sigma2)

df.series <- testdata
series <- testdata$St


nstepaheadpreds <- data.frame(rnorm(134,0,0))
for (n in seq(1,1000)){
  preds <- NULL
  for (j in seq(10,length(series))){
    preds <- append(preds,multistep(10, series=series,start=j,psi = psi, sigma2 = sigma2 ))
  }
  nstepaheadpreds[glue("{n}")] <- preds
}

nstepaheadpreds <- nstepaheadpreds[, 2:134]

forecastAverage <- function(i,df){
  lst <- NULL
  #df <- as.data.frame(df)
  for (j in seq(1,ncol(df))){
    lst <- append(lst, df[i,j])
  }
  return(mean(lst))
}

averages.nstepahead <- NULL



for(i in seq(1,134)){
  averages.nstepahead <- append(averages.nstepahead, forecastAverage(i,nstepaheadpreds))
}

averages.nstepahead <- data.frame(averages.nstepahead)

df.series %>%  ggplot(aes(x = 1:143, y = St)) + 
  geom_line() + 
  geom_line(aes(x=11:143, y = averages.nstepahead$averages.nstepahead))

### FINISH THIS OFF AND WRITE ALL IN REPORT
### REMEMBER TO CHANGE THE PREDICTIONS START COLUMN

vec <- rep(NA,9)
vec <- append(vec,averages.nstepahead)
averages.nstepahead <- data.frame(vec)


averages.nstepahead['St'] <- testdata$St


df.series %>%  ggplot(aes(x = 1:143, y = St)) + 
  geom_line() + 
  geom_line(aes(x=1:143, y = averages.nstepahead$vec))


table.multi <- getAllPredictions(nstepaheadpreds,1,predsstart = 1)
table.multi <- data.frame(table.multi)
for (i in seq(2,134)){
  table.multi[glue("{i}")] <- getAllPredictions(nstepaheadpreds, i, predsstart = 1)
}

lst <- NULL
for (i in seq(1,ncol(table.multi))){
  lst <- append(lst, getQuantiles(c(0.025, 0.975), table.multi, i))
}

quantile.25 <- NULL
for (i in seq(1,268,2)){
  quantile.25 <- append(quantile.25,lst[i])
}


temp <- rep(NA,9)
temp <- append(temp, quantile.25[1:134])
quantile.25 <- temp 


quantile.75 <- NULL
for (i in seq(2,268,2)){
  quantile.75 <- append(quantile.75,lst[i])
}

temp <- rep(NA,9)
temp <- append(temp, quantile.75[1:134])
quantile.75 <- temp 


df.series %>% ggplot(aes(x = 1:143, y = averages.nstepahead$vec)) + 
  geom_line(aes(colour="Forecasted")) + 
  geom_line(aes(x = 1:143, y = St, colour = "True Values"), alpha = 0.7) +
  geom_line(aes(x= 1:143, y = quantile.25 , colour = "2.5% Quantile"), 
            size = 0.5,
            linetype = "dashed") +
  geom_line(aes(x = 1:143, y = quantile.75, colour = "97.5% Quantile"),
            size = 0.5,
            linetype = "dashed") + 
  # geom_point(aes(x = 1:141, y = quantile.25, colour = "2.5% Quantile")) +
  # geom_point(aes(x = 1:141, y = quantile.75, colour = "97.5% Quantile")) +
  labs(x = "Time", y = TeX("$S_i$")) +
  ggtitle(TeX("$S_i$ vs $\\hat{S}_{i}$ 10-Step-Ahead Forecast")) 



(1/(133-1))*sum((averages.nstepahead$St[10:143] - averages.nstepahead$vec[10:143])^2)



m11 <- M11(averages.nstepahead$vec[10:nrow(averages.nstepahead)], 
           averages.nstepahead$St[10:nrow(averages.nstepahead)])

m21 <- M21(averages.nstepahead$vec[10:nrow(averages.nstepahead)], 
            averages.nstepahead$St[10:nrow(averages.nstepahead)])

m12 <- M12(averages.nstepahead$vec[10:nrow(averages.nstepahead)], 
           averages.nstepahead$St[10:nrow(averages.nstepahead)])

m22 <- M22(averages.nstepahead$vec[10:nrow(averages.nstepahead)], 
           averages.nstepahead$St[10:nrow(averages.nstepahead)])

m20 <- m21 + m22
m10 <- m11 + m12
m01 <- m11 + m21
m02 <- m12 + m22

contigency.table <- data.frame(cbind(c(m11,m21),c(m12,m22)))
contigency.table <- contigency.table %>% rename('Up'='X1', 'Down'='X2')
contigency.table

createDummyvariables <- function(vector){
  lst <- NULL
  for (i in seq(2,length(vector))){
    if(getDifferential(vector,i)=="+"){
      lst <- append(lst,1)
    } else {
      lst <- append(lst,0)
    }
  }
  return(lst)
}

predictions <- createDummyvariables(averages.nstepahead$vec[10:nrow(averages.nstepahead)])
true <- createDummyvariables(averages.nstepahead$vec[10:nrow(averages.nstepahead)])

chisq.test(predictions, true)


# 4


set.seed(12345678)
library(TSA)
garch01.sim = garch.sim(alpha=c(.01,.9), n = 500)
plot(garch01.sim, type="l", ylab = expression(r[t]), xlab="t")
mean(garch01.sim)
var(garch01.sim)

