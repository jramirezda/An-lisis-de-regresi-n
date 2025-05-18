## Parcial 2 tercera parte
## Jhon Ramirez Daza
## Miguel Angel Martin
install.packages("lmtest")
library(car)
library(lmtest)

## Carga de los datos
datos <- read.table("ExpImp.txt", header = TRUE, sep = "\t")
## aplicamos el attach para trabajar más cómodamente.
attach(datos)
## un breve resumen de los datos donde usaremos la variable Invext como variable dependiente y como 
## variables regresoras usaremos a mes, Export, Import, IPC
str(datos)

## Empezamos definiendo variables de interés
n <- length(Invext)
datos=as.matrix(datos)
colnames(datos) <- c("x1","y","x2","x3","x4")
X=as.matrix(cbind(1,datos[,1],datos[,3],datos[,4]))
H=X %*% solve( t(X) %*% X) %*% t(X)
diagH=diag(H)

## Creación de la matriz de datos
X=as.matrix(cbind(1,datos[,1],datos[,3],datos[,4],datos[,5]))

# 

## Para definir en qué orden incluir las variables en el modelo usamos el valor absoluto del t value.
summary(lm(Invext ~ mes + Import + Export + Ipc, data = datos))

## Empezamos desde el escenario base de Invext explicada solo por su promedio.
lm.0 <- lm(Invext ~ 1, data = datos)

summary(lm.0) ## resumen del modelo base
## vamos a usar este modelo teniendo en cuenta que no es un buen modelo pero lo haremos
## para comparar los modelos posteriores.
## como medida de comparación usaremos el AIC para la selección de las variables bajo este criterio
## expliquen más del modelo
AIC0 <- AIC(lm.0); AIC0
BIC0 <- BIC(lm.0); BIC0

-2 * sum(log(dnorm(Invext, mean = mean(Invext), sd = sd(Invext)))) + 4 ## Esta es la forma manual de calcular el AIC 

## también podemos ver el R cuadrado ajustado de forma equivalente al AIC ya que cuando el AIC disminuye el 
## R cuadrado ajustado aumentará de forma proporcional.

## por otro lado también usaremos la varianza estimada del modelo como criterio para usar una variable o no.
sigM0 <- summary(lm.0)$sigma; sigM0; sqrt(var(Invext)); sd(Invext)

## Realizaremos este mismo proceso con las siguientes variables según su nivel de significancia.
lm.1 <- lm(Invext ~ mes, data = datos)
summary(lm.1)
## como vemos en este caso todos nuestros indicadores mejoran por lo que decimos que la variable mes sí debe estar en el modelo.
AIC1 <- AIC(lm.1); AIC1
BIC1 <- BIC(lm.1); BIC1
R2Aj1 <- summary(lm.1)$adj.r.squared; R2Aj1
sigM1 <- summary(lm.1)$sigma; sigM1

## Segundo paso
lm.2 <- lm(Invext ~ mes + Import, data = datos)
summary(lm.2)
## como vemos también mejoran al agregar la variable Import por lo que debe estar en el modelo.
AIC2 <- AIC(lm.2); AIC2
BIC2 <- BIC(lm.2); BIC2
R2Aj2 <- summary(lm.2)$adj.r.squared; R2Aj2
sigM2 <- summary(lm.2)$sigma; sigM2


## Tercer paso 
lm.3 <- lm(Invext ~ mes + Import + Export, data = datos)
summary(lm.3)
## como vemos aún mejoran nuestros indicadores por lo que agregamos la variable Export
AIC3 <- AIC(lm.3); AIC3
BIC3 <- BIC(lm.3); BIC3
R2Aj3 <- summary(lm.3)$adj.r.squared; R2Aj3
sigM3 <- summary(lm.3)$sigma; sigM3

## Cuarto paso 
lm.4 <- lm(Invext ~ mes + Import + Export + Ipc, data = datos)
summary(lm.4)

## es claro que en este caso no mejoraron nuestros indicadores por lo que no se va a incluir la variable IPC.
AIC4 <- AIC(lm.4); AIC4
BIC4 <- BIC(lm.4); BIC4
R2Aj4 <- summary(lm.4)$adj.r.squared; R2Aj4
sigM4 <- summary(lm.4)$sigma; sigM4

## Por lo anterior, el modelo definitivo incluye a las variables mes, Import y Export.
lm.def <- lm(Invext ~ mes + Import + Export, data = datos)

## Ahora haremos la revisión de los supuestos.
## Primero revisaremos el estado de la multicolinealidad mediante el VIF.

lm.3 <- lm(Invext ~ mes + Import + Export, data = ExpImp)
lm.3a <- lm(mes ~ Import + Export, data = ExpImp)
RC3a <- summary(lm.3a)$r.squared; RC3a
vif_mes = 1 / (1 - RC3a)

lm.3b <- lm(Import ~ mes + Export, data = ExpImp)
RC3b <- summary(lm.3b)$r.squared; RC3b
vif_Import = 1 / (1 - RC3b)

lm.3c <- lm(Export ~ mes + Import, data = ExpImp)
RC3c <- summary(lm.3c)$r.squared; RC3c
vif_Export = 1 / (1 - RC3c)

## Como vemos, ninguno toma valores mayores a 2, lo cual no señala problemas de multicolinealidad.
vif = c(vif_mes, vif_Import, vif_Export); vif
vif(lm(Invext ~ mes + Import + Export, data = ExpImp))
## para verificar nuestros resultados usamos la función vif de la librería car y vemos que son
## iguales.

## En este caso, no necesitamos hacer correcciones por multicolinealidad, por lo que no haremos la 
## regresión Ridge.

## para verificar autocorrelación vamos a usar la función durbinWatsonTest de la librería car.
durbinWatsonTest(lm.def)
## No existe evidencia estadística para rechazar la hipótesis de que no existe autocorrelación de primer orden.

## Ahora veamos el supuesto de normalidad en los residuales empezando el análisis de residuales.
## Usando el test de Kolmogorov-Smirnov.
ks_test_result <- ks.test(lm.def$residuals, "pnorm", mean = mean(lm.def$residuals), sd = sd(lm.def$residuals))
ks_test_result
## Aplicando el test vemos que no existe evidencia estadística para rechazar la hipótesis de normalidad en los
## residuales. También veremos la gráfica de los Invext estimados vs los residuales del modelo.
plot(lm.def$fitted.values, lm.def$residuals)
## Tampoco se evidencian problemas con la normalidad ni patrones aparentes que nos indicaran 
## deficiencias en el modelo.

## vamos a repetir esta sección de análisis de residuales con los residuales estudentizados.
r = lm.def$residuals / (summary(lm.def)$sigma * (1 - diagH)); r
## estos son los errores estudentizados.

## Ahora veamos el supuesto de normalidad en los residuales estudentizados empezando el análisis de residuales.
## Usando el test de Kolmogorov-Smirnov.
ks_test_result <- ks.test(r, "pnorm", mean = mean(r), sd = sd(r))
ks_test_result
## Aplicando el test vemos que no existe evidencia estadística para rechazar la hipótesis de normalidad en los
## residuales estudentizados. También veremos la gráfica de los Invext estimados vs los residuales estudentizados
## del modelo.
plot(lm.def$fitted.values, r)
## Tampoco se evidencian problemas con la normalidad ni patrones aparentes que nos indicaran 
## deficiencias en el modelo.

## Ahora para verificar el supuesto de homocedasticidad usaremos la prueba Goldfeld-Quandt 
goldfeld_quandt_test <- gqtest(lm.def); goldfeld_quandt_test
## en consecuencia no existe evidencia estadística para rechazar el supuesto de homocedasticidad en los residuales.

## Llegados a este punto encontramos luego de las numerosas pruebas que el modelo es válido 
## ya que cumple con todos los supuestos necesarios para aplicar un modelo de regresión lineal múltiple 
## sin la necesidad de hacer ninguna corrección, por lo que es un buen modelo que explica el 93% de la variabilidad 
## total del modelo basándonos en el R cuadrado ajustado. Cabe recalcar que este modelo tiene en cuenta solo las
## variables necesarias para explicar el modelo y fueron agregadas al modelo según su significancia.
