rm(list = ls())
##############################################################
######## Ajuste del Modelo de Regresión Multiple ##############
##############################################################

ExpImp<-read.table("ExpImp.txt", header=T)
attach(ExpImp)

#### Generación de vectores
y=ExpImp$Invext
mes=ExpImp$mes;Import= ExpImp$Import;Export=ExpImp$Export;Ipc= ExpImp$Ipc

########## ESTINACIÓN DEL MODELO DE REGRESIÓN USANDO LA INSTRUCCIÓN "lm" ########
lm.4<-lm(Invext ~ mes + Import + Export + Ipc, data=ExpImp)
summary(lm.4)


######## Construcción del Modelo sin el Uso Directo de la Instrucción "lm"  ############

######## Matrix de información.
# a) Obtención mediante el uso de la instrucción "lm"
x<-model.matrix(lm.4)
hat(x)            #Diagonal de H  #Salida de la instrucción "lm"
hatvalues(lm.4)   #Diagonal de H  #Salida de la instrucción "lm"

# b) Usando instrucciones básicas matriciales del R
x=cbind(1, mes, Import, Export, Ipc)
x=as.matrix(x);x

####### Matriz Hat y Matriz de Proyección
n=length(y)
XtX=t(x)%*%x;eigen(XtX)
IXtX=solve(t(x)%*%x)
XtY=t(x)%*%y
H=(x)%*%(IXtX)%*%t(x)
#H%*%H-H; H- t(H)
I=diag(n)
M=I-H
round(M - M%*%M,10); 
round(M - t(M),10)
sum(diag(M)) # n-p-1=24-4-1=19
hat(x) # Diagonal de H

####### Estimación de beta
Be=IXtX%*%XtY;Be
lm.4$coefficients

####### Estimaciones de Y (Y estimado)
YE=H%*%y
cbind(YE , x%*%Be , as.vector(lm.4$fitted.values))

plot(y,main="Modelo Ajustado",col="red")
lines(YE,col="blue")

####### Errores del modelo
e=y-YE
mean(e)
cbind(e , M%*%y , lm.4$residuals )

####### Suma de cuadrados del error, grados de libertad y varianza del error
SCE=as.numeric(t(e)%*%e);SCE
var(e)*(n-1)
gle=n-length(Be)
Sigma2=SCE/gle;Sigma2
Sigma=sqrt(Sigma2);Sigma
summary(lm.4)$sigma                 # Salida de la instrucción "lm"

### matriz de varianzas y covarianzas de beta  y varianza de Y estimado
SIGMA = IXtX * Sigma2;SIGMA
vcov(lm.4)                           # Salida de la instrucción "lm"
seBeta=sqrt(diag(SIGMA));seBeta      # errores estándar
summary(lm.4)$coefficients[,2]
summary(lm.4)				 # Salida de la instrucción "lm"
V=H*Sigma2

####### Sumas de cuadrados total, grados de libertad total y varianza de Y
SCT = as.numeric(t(y)%*%y - n*(mean(y))^2);SCT   # Syy
var(y)*(n-1)
glt=n-1
S2y=SCT/glt;S2y;var(y)

####### Sumas de cuadrados de la regresión, "gl" y cuadrado medio de la regresión
SCR = t(Be) %*% XtY - n*(mean(y))^2 ;
SCR =as.numeric(SCR);SCR
sum((YE-mean(y))^2)                      # SCR = suma ( y estimado - M[y] )^2
glr=length(Be)-1
CMR=SCR/glr;CMR

####### Estadística F, y prueba global
Fc = CMR/Sigma2 ; Fc  		    # CMR/CME
1-pf(Fc,glr, gle) 		    # valor-p
summary(lm.4)$fstatistic
summary(lm.4)                     #Salida de la instrucción "lm"

####### Pruebas individuales
tc=Be/seBeta;tc                   #estadística
cbind(tc, summary(lm.4)$coefficients[,3])
2*(1-pt(abs(tc), gle))            #valores-p
cbind(2*(1-pt(abs(tc), gle)) , summary(lm.4)$coefficients[,4])

####### R cuadrado
R2 = SCR/SCT ; R2
summary(lm.4)$r.squared         #Salida de la instrucción "lm"

####### R cuadrado ajustado
R2Aj = 1-Sigma2/S2y ; R2Aj
summary(lm.4)$adj.r.squared     #Salida de la instrucción "lm"


####### Inferencia simultánea de parámetros
#Para beta1 a beta4
alpha=0.05;
qt(1-alpha/2, gle) 			# percentil para IC individuales
alphaT=alpha/(2*4)
DBonferroni=qt(1-alphaT, gle);
DBonferroni        			# percentil para IC simultáneos
LIs=Be[2:5]-DBonferroni*seBeta[2:5]
LSs=Be[2:5]+DBonferroni*seBeta[2:5]
cbind(LIs, Be[2:5], LSs)
#Si no es simultánea
LI=Be[2:5]-qt(1-alpha/2, gle)*seBeta[2:5]
LS=Be[2:5]+qt(1-alpha/2, gle)*seBeta[2:5]
cbind(LI,LS)
confint(lm.4,level=0.95)[2:5,]   #Salida de la instrucción "lm"


####### Inferencia simultánea de la respuesta media
#Para las observaciones 1 y 2
alphaT=alpha/(2*2)
DBonferroni=qt(1-alphaT, gle)
x1 = x[1,] ; x2 = x[2,]

h11 = t(x1) %*% IXtX %*% x1 ; h11 ; H[1,1] ; hat(x)[1]
h22 = t(x2) %*% IXtX %*% x2 ; h22 ; H[2,2] ; hat(x)[2]

VarY1 = Sigma2 * h11 ; VarY1; V[1,1]
VarY2 = Sigma2 * h22 ; VarY2; V[2,2]

eeY1 = sqrt(VarY1)
eeY2 = sqrt(VarY2)

Ye=YE[1:2];Ye
LiY=Ye-DBonferroni*c(eeY1,eeY2)
LuY=Ye+DBonferroni*c(eeY1,eeY2)
cbind(LiY,Ye,LuY)
#y[1:2]

#Si no es simultánea
LiY=Ye-qt(1-alpha/2, gle)*c(eeY1,eeY2)
LuY=Ye+qt(1-alpha/2, gle)*c(eeY1,eeY2)
cbind(LiY,Ye,LuY)


LLl=YE-qt(1-alpha/2, gle)* Sigma * sqrt(hat(x))
LLu=YE+qt(1-alpha/2, gle)* Sigma * sqrt(hat(x))
plot(y,main="Modelo Ajustado",col="red")
lines(YE,col="blue")
lines(LLl)
lines(LLu)

LLlb=YE-qt(1-alpha/2, gle)* Sigma * sqrt(1+hat(x))
LLub=YE+qt(1-alpha/2, gle)* Sigma * sqrt(1+hat(x))
plot(y,main="Modelo Ajustado",col="red")
lines(YE,col="blue")
lines(LLl)
lines(LLu)
lines(LLlb, col="red")
lines(LLub, col="red")

