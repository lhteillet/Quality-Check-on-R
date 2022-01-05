#' ---
#' title: "Report TD4"
#' author: "Louis TEILLET"
#' ---
knitr::opts_chunk$set(message = FALSE)
#' # Exercise 1 : TS and DS process
#' 
#' ### 1. Computing 2 random walk with a drift and store them within a matrix
library(tseries)
N = 1000
u1 = rnorm(N) # Simulation of 1000 Gaussian observations
y1=c()
y2=c()
y1[1]=u1[1] # First random walk
y2[1]=u1[1] # Second random walk

for (k in 2:N)
{
 y1[k] = y1[k-1] + 0.5 +u1[k] #With a drift=0.5
 y2[k]= y2[k-1]+1+u1[k] #With a drift=1
 
}
#' ### 2. Plot them all in order to verify(visually) the existence of a trend
par(mfrow=c(2,1)) 
plot(y1,type="l",col="blue", main="Random Walk with drift = 2") #Plotting the first one
plot(y2,type="l",col="red",main ="Random Xalk with drift=1") #Plotting the second one

#' We observe a trend but it is hard to say if it's a deterministic trend or a stochastic one
#' 
#' ### 3. Select one of the random walk you have generated and compute a linear regression where the explanatory variable is a trend t
time = seq(1:N)
lr = lm(y1~time)
summary(lr)
plot(y1,type="l",main="y1(t) and linear regression",xlab="t",ylab="y1(t)")
lines(time,lr$fitted.values,col="red")

#' ### 4. Is it the right way to neutralize the source of non stationary ? 
#' 
#' We try to neutralize the source of non stationary upon making the difference between y(t) and lr(t)
#' 
y1.minus.lr = y1-lr$fitted.values
plot(y1.minus.lr,type="l") #There is still a trend
acf(y1.minus.lr) #The acf show a really high persistence, the process is not stationary

#' ### 5. Differentiate the selected variable and check the autocorrelation function
#' 
#' 
y1.diff = diff(y1)
plot(y1.diff,type="l")
acf(y1.diff) # The acf doesn't show persistence anymore 

#' # Exercise 2 : ADF, PP and KPSS unit root tests
#' ### Dickey Fuller

library("urca") #Loading "urca" package
#' The Dickey Fuller tests the type of non stationary of a process.
#' 
#' The first type of non stationary is a model without drift and temporal trend:
#' 
#' $y_{t} = \rho_{1}.y_{t-1} + \epsilon_{t}$
#' 
#' The second one is a model with drift but without temporal trend:
#' 
#' $y_{t} = \rho_{1}.y_{t-1} + c +  \epsilon_{t}$
#' 
#' The last one is a model with drift and temporal trend:
#' 
#' $y_{t} = \rho_{1}.y_{t-1} + c + a.t +  \epsilon_{t}$
#' 
#' 
#' The null hypothesis is $H_{0} : \rho = 0$
#' 
#' We start to test the third model, if the non-nullity of the temporal trend coefficient is not significant 
#' we go trough the second one, else we check if $\rho_{1}$ is significantly different from 0.  
#' 
#' We do the same thing for the model 2, and then for the model 1. 
#' 

data = read.csv("usgdp.csv",sep=";")
usgdp.ts = ts(data$usgdp,frequency=4,start=c(1950,1))
plot(usgdp.ts)

summary(ur.df(usgdp.ts,type="drift",selectlags="AIC"))
#' We can't reject Ho : The presence of a unit root, hence the serie is a random walk
par(mfrow=c(2,1))
acf(usgdp.ts)
pacf(usgdp.ts) 
#' The dickey fuller result is confirmed by the autocorrelograms 

usgdp.diff= diff(usgdp.ts) # We differentiate the serie one time
acf(usgdp.diff)
pacf(usgdp.diff)
#' Autocorrelograms seems to indicate the stationnarity of the time serie
summary(ur.df(usgdp.diff,type="drift",selectlags="AIC"))
#' We now can reject Ho : there is a unit root and accept H1 : the serie is stationnary
#' 
#' Conclusion : usGDP is I(1)
#' 
#' 

#' ### Phillips-Perron
#' 
#' The Phillips-Perron (PP) unit root tests differ from the ADF tests mainly in how they deal with serial correlation and heteroskedasticity in the errors. 
#' In particular, where the ADF tests use a parametric autoregression to approximate the ARMA structure of the errors in the test regression, the PP tests ignore any serial correlation in the test regression.
library(aTSA)
pp.test(data$usgdp)
pp.test(usgdp.diff)
#' p value <=0.01 so we reject Ho : the serie has a unit root.
#' 
#' Conclusion: The Phillips-Perron unit root tests return the same result than Dickey Fuller.
#' **Usgdp is I(1) for this test**
#' 

#' ### KPSS 
#' 
#' Here we just have to be care of the fact that the null hypothesis is the stationnarity of the serie
kpss.test(data$usgdp)
kpss.test(usgdp.diff)
kpss.test(diff(data$usgdp,2))
kpss.test(diff(data$usgdp,3))
kpss.test(diff(data$usgdp,4))
kpss.test(diff(data$usgdp,5)) #p value >0.05 we can't reject the null hypothesis 

#' Kpss returns an order of integration of 5.
#' **Conclusion : With KPSS test, usgdp is I(5)
#' 






#' # Exercise 3 : Estimating ARIMA(p,d,q)
#' 
library(tidyquant)
jnj = tq_get("JNJ",get="stock.prices",from="1997-01-01") %>% tq_transmute(mutate_fun = to.period,period="months")

#' ### 1. Determine the degree of integration of the Johnson & Johnson stock price.
#' 
#' 
plot(jnj$close,type="l",main="Close price P(t)",ylab="Close price",xlab="time")
df.test=ur.df(jnj$close,type="trend",selectlags="AIC")
summary(df.test)
#' All of our test value are within the "fail to reject" zone, 
#' so jnj is a time serie which has a unit root (it's a random walk), but neither trend nor constant drift
#'  
#'  We try to differentiate the time serie, one time

jnj.diff = diff(jnj$close)
df.test.diff=ur.df(jnj.diff,type="trend",selectlags="AIC")
summary(df.test.diff)
#' The three test-statistic are within the "reject" zone, so the serie is staionnary but with less than 1% of error, 
#' The p-value of the linear trend is 0.27 so we fail to reject the null hypothesis of the coeff, the model has not linear trend
#' same thing for the constant, we can't reject the null hypothesis. 
#' 
#' **The time serie is I(1)**
#' 

#' ### 2. Determine the order of the required ARIMA model, i.e the values of p,d,q.
#' 

par(mfrow=c(1,1))
acf(jnj.diff) 
pacf(jnj.diff) 

#' We observe 2 pics out of the significance limit, for the lag 2 and 6 
#' 
#'  **p = 2 or p =6, we will choose 2 to respect the parsimony of the model**
#' 
#' We observe 3 pics out of the significance limit, for the lag 2,6 and 7 
#' 
#'  **q = 2 or q = 9 we will choose 2 to respect the parsimony of the model**

#determiner (p,q) avec les criteres d'informations
library(qpcR)
model = arima(jnj.diff,order=c(1,0,0))
min_AIC = c(AIC(model),1,0)
min_BIC=c(BIC(model),1,0)
for(p in 1:8)
{
  for(q in 1:8)
  {
    model = arima(jnj.diff,order=c(p,0,q))
    ci1 = AIC(model)
    ci2 = BIC(model)
    if (ci1<min_AIC[1])
    {
      min_AIC=c(ci1,p,q)
    }
    if (ci2<min_BIC[1])
    {
      min_BIC=c(ci2,p,q)
    }
    
    
  }
}
cat("Akaike criteria : p =",min_AIC[2], "and q =", min_AIC[3])

cat("Baysian criteria : p =",min_BIC[2], "and q =", min_BIC[3])


#' We choose to keep p=1 and q=1 to respect the parsimony of the model
#' The final model is an ARIMA(1,1,1)
#' 

#' ### 3. Estimating our ARIMA(1,1,1)
#' 
#' 
model = arima(jnj$close,order=c(2,1,2))
model
plot(jnj$close,type="l",ylab="prices",xlab="time",col="blue")
points(jnj$close+model$residuals,type="l",col="red",lty=2)
legend("topleft",legend=c("real value","fitted values"),col=c("blue","red"),lty=1)
plot(model$residuals)


#' ### 4. Check if the residuals are gaussian
#' 

jarque.bera.test(model$residuals) #p-value < $2^{-16}$ we reject Ho :residuals have a normal kurtosis

#' Maybe because of the structural break of the time serie
#' 



#' ### 5. Forecast
#' 
jnj.forecasting = predict(model,n.ahead=3)
jnj.forecasting = c(jnj[300,]$close,jnj.forecasting$pred)
plot(jnj$close,type="l",col="red")
points(c(300,301,302,303),jnj.forecasting,col="blue",type="l")



#' # Exercise 4 : Unit root test for another one
#' 
#' 
#' ### 2. Generating 3 differents random walk
#' 
N = 1000
u1 = rnorm(N) # Simulation of 1000 Gaussian observations
u2 = rnorm(N) # Simulation of 1000 Gaussian observations
u3 = rnorm(N) # Simulation of 1000 Gaussian observationsy1=c()
y1=c()
y2=c()
y3=c()
y1[1]=u1[1] # First random walk without break
y2[1]=u2[1] # with level break
y3[1]=u3[1] # with trend and level break

for (k in 2:(N/2))
{
  y1[k] = y1[k-1] + 0.5 +u1[k] 
  y2[k]= y2[k-1]+0.3+u2[k] 
  y3[k]=y3[k-1]+0.2+0.3*k + u3[k]
  
}
for (k in (N/2):N)
{
  y1[k] = y1[k-1] + 0.5 +u1[k] 
  y2[k]= y2[k-1]+1+u2[k] 
  y3[k]=y3[k-1]-0.1+0.6*k + u3[k]
  
}
par(mfrow=c(1,1))
plot(y1,type="l")
plot(y2,type="l")
plot(y3,type="l")


#' ### 3. Compute the Zivot and Andrews test

summary(ur.za(y1))
summary(ur.za(y2,model="intercept"))
summary(ur.za(y3,model="both"))

#' ### 4. Is it relevant to use such test for the US GDP? Compute the Zivot and Andrews unit root test using the US GDP
#' 


#' # Exercise 4.2 : Modeling

#' ### Order of integration
plot(usgdp.ts)

summary(ur.df(usgdp.ts,type="drift",selectlags="AIC"))
#' We can't reject Ho : The presence of a unit root, hence the serie is a random walk
par(mfrow=c(2,1))
acf(usgdp.ts) # q=0
pacf(usgdp.ts) # p  =1
#' The dickey fuller result is confirmed by the autocorrelograms 

usgdp.diff= diff(usgdp.ts) # We differentiate the serie one time
acf(usgdp.diff)
pacf(usgdp.diff)
#' Autocorrelograms seems to indicate the stationnarity of the time serie
summary(ur.df(usgdp.diff,type="drift",selectlags="AIC"))
#' We now can reject Ho : there is a unit root and accept H1 : the serie is stationnary
#' 
#' Conclusion : usGDP is I(1)
#' 
#' ### ARMA(p,q) analyse

model = arima(usgdp.ts,order=c(1,1,0),method="ML")
min_AIC = c(AIC(model),1,0)
min_BIC=c(BIC(model),1,0)
for(p in 1:4)
{
  for(q in 1:4)
  {
    model = arima(usgdp.ts,order=c(p,1,q),method="ML")
    ci1 = AIC(model)
    ci2 = BIC(model)
    if (ci1<min_AIC[1])
    {
      min_AIC=c(ci1,p,q)
    }
    if (ci2<min_BIC[1])
    {
      min_BIC=c(ci2,p,q)
    }
    
    
  }
}
cat("Akaike criteria : p =",min_AIC[2], "and q =", min_AIC[3])

cat("Baysian criteria : p =",min_BIC[2], "and q =", min_BIC[3])

#' **Conclusion:** we choose an ARIMA(1,1,2) to modelize GDP
#' 
#' 
#' 

#' ### Quality Check
#' 

par(mfrow=c(1,1))
gdp_model = arima(usgdp.ts,order=c(1,1,2),method="ML")
plot(usgdp.ts,ylab="us gdp",col="blue")
points(usgdp.ts+gdp_model$residuals,col="red",type="l")
legend("topleft",legend=c("real value","fitted values"),col=c("blue","red"),lty=1)

plot(gdp_model$residuals)
#' our model seems to be not really robust when an economic crash appear (problem of structural break, but I have not really understand this concept and the test of the presence of these breaks.)

jarque.bera.test(model$residuals)
jarque.bera.test(model$residuals[1:100])

#' ### Conclusion
#' 
#' The final test of jarque bera is the proof that our model is good in a period of non-break because residuals have
#' a normal kurtosis but the same model is not robust when crashes appear. Actually, is period of crash our model make high
#' error.
#' 
#' Maybe we should test another model with more differentiation, or another type of model which treat the influence of break on our model.
#' 



