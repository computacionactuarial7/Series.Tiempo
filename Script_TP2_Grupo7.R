######################
##Trabajo Practico 2##
## Series de tiempo ##
##    Grupo 7       ##
######################

library (FinTS)
library(vars)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(lubridate)
library(tseries)
library(astsa)
library(moments)
library(forecast)
library(foreign)
library(quantmod)
library (fGarch)
library(dynlm)
library(rmgarch)
library(rugarch)

#EJEMPLO 1. RETORNOS VISA

#Se buscan datos de las acciones de Visa en Yahoo finance
startDate = as.Date ("2010-01-04")
endDate = as.Date ("2020-07-23")

getSymbols("V", from = startDate, to= endDate)
head(V)
plot(V$V.Adjusted, ylab = "Precios",main = "Precios de cierre ajustados Visa", col = "red")



# Se obtienen los retornos
rV <- dailyReturn(V)
head(rV)

#Se utilizan los datos de los retornos para crear una serie de tiempo
Visa.ts <- ts(rV,
              frequency = 252,      
             start= c(2010,01,04),
             end= c(2020, 07, 22))

class(Visa.ts)

str(Visa.ts)


#Gráfico de la serie de tiempo
plot(Visa.ts, ylab="Rendimientos diarios", main= "Serie de Tiempo de los Rendimientos de Visa")


#Verificar si es estacionario. P VALUE < 0.05
adf.test(Visa.ts)

#Funcion de autocorrelacion y autocorrelacion parcial 
#para ver cuantas medias moviles y cuantos regresivos vamos a utilizar en el modelo
par(mfrow=c(1,2))
acf(Visa.ts)
pacf(Visa.ts)

#Tenemos que hacer que el resago coincida con las frecuencia
acf(ts(Visa.ts, frequency=1))    # numero de medias moviles 
pacf(ts(Visa.ts, frequency=1))   # numero de autorregresivos

modelo1 <- arima(Visa.ts, order = c(1,0,1))  
coeftest(modelo1)
modelo1
AIC(modelo1)

#Diagnostico para ver si el modelo es bueno
tsdiag (modelo1)

#Calculamos el error
par(mfrow=c(1,1))
error <- residuals (modelo1)
plot(error)

hist(error)
qqnorm(error)
qqline(error, col = "yellow")
acf(error)
pacf(error)
kurtosis(error)

#Validar si el error es ruido blanco, PVALUE >0.05
Box.test(error, type = "Ljung-Box")

#Por lo tanto evaluamos si hay efectos GARCH

#Calculo los residuales al cuadrado
residuos <- residuals(modelo1)^2
plot(residuos,main='Residuales al cuadrado')

#Muestra que los residuos al cuadrado no son constantes
#Hay varianza heterocedastica. Modelo GARCH

#Hay que agregar la varianza al proceso ARIMA porque se ve que la varianza presenta heterocedasticidad
V.arch<- dynlm (residuos ~ L (residuos), data = V)
summary(V.arch)
V.arch

#H0: No hay efectos GARCH >0.05
#H1: Hay efectos GARCH <0.05
#Los errores al cuadrado son significativos.

Box.test(residuals(modelo1)^2, type = "Ljung-Box")  #No hay Ruido Blanco. Los errores estan correlacionados 

#Tambien podemos probarlo con ARCHTEST
prueba <- ArchTest(Visa.ts, lags=1, demean= TRUE)
prueba

prueba2 <-ArchTest(Visa.ts, lags=2, demean = T)
prueba2 

#Si hay efectos GARCH
par(mfrow=c(1,2))
acf(residuos)
pacf(residuos)

ug_spec <- ugarchspec(mean.model = list(armaOrder = c(1,1)))

ugfit <- ugarchfit(spec = ug_spec, data = rV)

ug_var <- ugfit@fit$var
ug_res <- (ugfit@fit$residuals)^2
plot(ug_res, type = "l")
lines(ug_var, col = "green")

#Pronóstico del modelo ARIMA(1,0,1) + GARCH(1,1)

ugforce <- ugarchforecast(ugfit, n.ahead = 10)

plot(ugforce)


##############################################################

#EJEMPLO 2. PRECIO ACCIONES MERCADO LIBRE

#Buscamos datos de las acciones de Mercado Libre en Yahoo finance
startDate = as.Date ("2010-01-04")
endDate = as.Date ("2020-07-23")

getSymbols("MELI", from = startDate, to= endDate)
head(MELI)

#Se arma una serie de tiempo con los precios de las acciones
Meli.ts <- ts(MELI$MELI.Adjusted,
              frequency = 252,      
              start= c(2010,01,04),
              end= c(2020, 07, 22))

class(Meli.ts)

str(Meli.ts)

#Grafico de la serie de tiempo
plot(Meli.ts, ylab = "Precios",main = "Precios de cierre ajustados Mercado Libre", col = "red")

#Evaluo si es modelo aditivo o multiplicativo
Modelo.aditivoML <- decompose(Meli.ts)
plot(Modelo.aditivoML, col= "blue")

Modelo.multiplicativoML <- decompose (Meli.ts, type="mult")
plot(Modelo.multiplicativoML, col="purple")

#Busco la tendencia
TendenciaML <- Modelo.multiplicativoML$trend
print (TendenciaML)   

#Busco estacionalidad
EstacionalidadML <- Modelo.multiplicativoML$seasonal
print(EstacionalidadML)

#La tendencia y la estacionalidad se deben agregar al modelo de regresion 
ts.plot (cbind(TendenciaML, TendenciaML*EstacionalidadML), lty=1:2)

#Chequeamos la transformacion logaritmica para convertirlo 
#en un modelo estacionario
Meli.ts2 <- log(Meli.ts)
plot(Meli.ts2, col="red")
plot(decompose(Meli.ts2), col="red")

#Hago trasformacion con diferencial para convertirlo 
#en un modelo estacionario
Meli.ts3 <- diff(Meli.ts)
plot(Meli.ts3, col="orange")
plot(decompose(Meli.ts3), col="orange")

#Test de Dickey Fuller para ver si hay raices unitarias. 
#PVALUE <0.05 ES ESTACIONARIA
adf.test(Meli.ts3)

# Es estacionario con una diferencia

#Hacemos ARIMA con una diferencia de los precios de MELI

#Funcion de autocorrelacion y autocorrelacion parcial 
#para ver cuantas medias moviles y cuantos regresivos vamos a utilizar en el modelo
par(mfrow=c(1,2))
acf(Meli.ts3)
pacf(Meli.ts3)

#Tenemos que hacer que el resago coincida con las frecuencia
acf(ts(Meli.ts3, frequency=1))    # numero de medias moviles 
pacf(ts(Meli.ts3, frequency=1))   # numero de autorregresivos 

#Se hace sobre el modelo de serie de tiempo original poniendo 1 por diferencia
modeloML <- arima(Meli.ts, order = c(1,1,1))   
modeloML
coeftest(modeloML)
AIC(modeloML)

#Hago diagnostico para ver si el modelo es bueno
tsdiag (modeloML)

#Prueba de Ruido Blanco
#si >0.05 el modelo ajusta bien hay ruido blanco

Box.test(residuals(modeloML), type = "Ljung-Box")  

#Hay Ruido Blanco. Error media=0

#Calculamos el error

errorML <- residuals (modeloML)
par(mfrow=c(1,1))
plot(errorML)
summary(errorML)

hist(errorML)
qqnorm(errorML)
qqline(errorML, col = "pink")
acf(errorML)
pacf(errorML)
kurtosis(errorML)

#MODELO GARCH

#Calculamos los residuales al cuadrado
residuosML <- resid(modeloML)^2
plot(residuosML,main='Residuales')

qqnorm(residuosML)
qqline(residuosML, col = "green")
acf(residuosML)
pacf(residuosML)

#Muestra que los residuos al cuadrado no son constantes
#Hay varianza heterocedastica. Modelo GARCH

#Hay que agregar la varianza al proceso ARIMA porque se ve que la varianza presenta heterocedasticidad
modelo2ML<- dynlm (residuosML ~ L (residuosML), data = MELI)
summary(modelo2ML)

Box.test(residuosML, type = "Ljung-Box") #No hay ruido blanco

#Es significativo
#H0: No hay efectos GARCH >0.05
#H1: Hay efectos GARCH <0.05
#Los errores al cuadrado son significativos.

#Tambien podemos probarlo con ARCHTEST
pruebaML <- ArchTest(Meli.ts3, lags=1, demean= TRUE)
pruebaML

# Probamos el ARCHTEST con 2 rezagos
prueba2ML <- ArchTest(Meli.ts3, lags=2, demean= TRUE)
prueba2ML
# Corroboramos que hay efectos GARCH

#Buscamos el modelo que se ajuste

plot.ts(residuals(modeloML),main = 'Residuos, Modelo en Varianza, Garch(1,1)')
qqnorm(residuals(modelo2ML))
qqline(residuals(modelo2ML))
acf(residuals(modelo2ML)^2,na.action = na.omit)
pacf(residuals(modelo2ML)^2,na.action = na.omit)

#Entonces modelo GARCH(1,1)

ModeloGarch11 <-   garch(na.omit(residuals(modeloML)),order = c(1,1),na.action=na.omit)
coeftest(ModeloGarch11)

plot.ts(residuals(ModeloGarch11),main = 'Residuos, Modelo en Varianza, Garch(1,1)')
qqnorm(residuals(ModeloGarch11))
qqline(residuals(ModeloGarch11), col = "cyan")
acf(residuals(ModeloGarch11)^2,na.action = na.omit)
pacf(residuals(ModeloGarch11)^2,na.action = na.omit)

fit <- fitted.values(ModeloGarch11)

plot(MELI$MELI.Adjusted,
     main = 'Mercado Libre: Precios de cierre ajustados
Modelo ARIMA(1,1,1)-GRACH(1,1)',
     type = 'l')

