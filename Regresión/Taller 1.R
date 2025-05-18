##Jhon Ramirez Daza
##Miguel Angel Matrin

## Carga de los datos
datos=read.table("ExpImp.txt", header = TRUE, sep = "\t")

#manipulación de los datos
datos=as.matrix(datos)
datos
colnames(datos) <- c("x1","y","x2","x3","x4")

## Creación de la matriz de datos
X=as.matrix(cbind(1,datos[,1],datos[,3],datos[,4],datos[,5]))
X
# creacion de el vector respuesta  
y=datos[,2]
y

n=dim(datos)[1]
n
p=dim(X)[2]
p
k=p-1

# creacion de la matriz hat
H=X %*% solve( t(X) %*% X) %*% t(X)
H
# Creadion del vector de Betas estimados 
BetasE=solve( t(X) %*% X) %*% t(X) %*% y
BetasE

# Creacion del vector y estimado
yE= as.matrix(H %*% y)
yE

## calculo de los errores del modelo
errores=y-yE
errores
sum(errores)

#calculo de la suma de los errores al cuadrado
SSe=sum(errores^2)
SSe
gle=n-k-1;gle
glr=k;glr
# calculo de la varianza estimada
MSe=SSe/gle;MSe
#Matriz de varianzas y covarianzas de BetaE
covBetaE=MSe*solve(t(X) %*% X);covBetaE

#varianza de y estimado 


##Sumas de cuadrados total, grados de libertad total y varianza de Y
##Sumas de cuadrados de la regresión, "gl" y cuadrado medio de la regresión

SSr <- t(BetasE) %*% t(X) %*% y - ((sum(y))^2 / n)
MSr <- SSr / glr;MSr

SSt=SSr+SSe
MSt=SSt/glr+gle;MSt

#estadistica F calculada y p valor
fcal <- MSr / MSe;fcal
pval <- 1 - pf(fcal, k, n - k - 1);pval


# Pruebas individuales


#R2 y  R2 ajustado
R_2=1-SSe/SSt;R_2
R_2A=1-MSr/MSt;R_2A

###########################################################3
#sumas de cuadrados 
############################################################
## Carga de los datos
ExpImp=read.table("ExpImp.txt", header = TRUE, sep = "\t")
attach(ExpImp)
n=length(Invext)
cor(Invext,Import)
cor(Invext,Export)
cor(Import,Export)

summary(lm(Invext~Import+Export,data = ExpImp))

l1=lm(Invext~Export,data = ExpImp)#export:x2 invext:y
l2=lm(Import~Export,data = ExpImp)#x1:import

R_InvEImp.Exp=cor(l1$residuals,l2$residuals);R_InvEImp.Exp
cor(Invext,Import)/cor(l1$residuals,l2$residuals)

(-2.776)^2/((-2.776)^2+21);R_InvEImp.Exp^2

##############Selección de  un modelo 


summary(lm(Invext~mes+Import+Export+Ipc,data = ExpImp))

#paso 0
lm.0=lm(Invext~1,data=ExpImp)
summary(lm.0)
mean(Invext)
AIC0=AIC(lm.0);AIC0
-2*sum(log(dnorm(Invext,mean = mean(Invext),sd=sd(Invext))))+4

R2Aj0=summary(lm.0)$adj.r.squared;R2Aj0
sigM0=summary(lm.0)$sigma;sigM0;sqrt(var(Invext));sd(Invext)

gle0=lm.0$df.residual;gle0
SCT=var(Invext)*(n-1);SCT
SCE0=sigM0^2*gle0;SCE0
SCR0=sum((lm.0$fitted.values-mean(Invext))^2);SCR0
SCT-SCE0


####paso 1
##la variable mas importante es mes
lm.1=lm(Invext~mes,data=ExpImp)
summary(lm.1)
AIC1=AIC(lm.1);AIC1
-2*sum(log(dnorm(Invext,mean =15310.50+247.90*mes ,sd=518.2)))+2*(3)



R2Aj1=summary(lm.1)$adj.r.squared;R2Aj1
sigM1=summary(lm.1)$sigma;sigM1

gle1=lm.1$df.residual;gle1
SCT=var(Invext)*(n-1);SCT
SCE1=sigM1^2*gle1;SCE1;sum(lm.1$residuals^2)
SCR1=sum((lm.1$fitted.values-mean(Invext))^2);SCR1
SCT-SCE1

summary(lm.1)$fstatistic

SCRs1=SCR1

#######################
###paso 2
##la variable mas importante sin  mes
summary(lm(Invext~Import+Export+Ipc,data = ExpImp))
lm.2=lm(Invext~mes+Import,data=ExpImp)
summary(lm.2)
AIC2=AIC(lm.2);AIC2
-2*sum(log(dnorm(Invext,mean =16504.7465+229.8082*mes-0.6020*Import ,sd=468.2)))+2*(4)



R2Aj2=summary(lm.2)$adj.r.squared;R2Aj2
sigM2=summary(lm.2)$sigma;sigM2

gle2=lm.2$df.residual;gle2
SCT=var(Invext)*(n-1);SCT
SCE2=sigM2^2*gle2;SCE2;sum(lm.2$residuals^2)
SCR2=sum((lm.2$fitted.values-mean(Invext))^2);SCR2
SCT-SCE2

summary(lm.2)$fstatistic

SCRs2=SCR2-SCR1;SCRs2
(SCRs2/1)/(sigM2^2)
(1-pt(sqrt((SCRs2/1)/(sigM2^2)),gle2))*2
1-pf((SCRs2/1)/(sigM2^2),1,gle2)

#######################
###paso 3
##la variable mas importante sin  mes
summary(lm(Invext~Export+Ipc,data = ExpImp))
lm.3=lm(Invext~mes+Import+Export,data=ExpImp)
summary(lm.3)
AIC3=AIC(lm.3);AIC3


R2Aj3=summary(lm.3)$adj.r.squared;R2Aj3
sigM3=summary(lm.3)$sigma;sigM3

gle3=lm.3$df.residual;gle3
SCT=var(Invext)*(n-1);SCT
SCE3=sigM3^2*gle3;SCE3;sum(lm.3$residuals^2)
SCR3=sum((lm.3$fitted.values-mean(Invext))^2);SCR3
SCT-SCE3

summary(lm.3)$fstatistic

SCRs3=SCR3-SCR2;SCRs3
(SCRs3/1)/(sigM3^2)
(1-pt(sqrt((SCRs3/1)/(sigM3^2)),gle3))*2
1-pf((SCRs3/1)/(sigM3^2),1,gle3)


#######################
###paso 4
##la variable mas importante sin  mes
summary(lm(Invext~Ipc,data = ExpImp))
lm.4=lm(Invext~mes+Import+Export+Ipc,data=ExpImp)
summary(lm.4)
AIC4=AIC(lm.4);AIC4


R2Aj4=summary(lm.4)$adj.r.squared;R2Aj4
sigM4=summary(lm.4)$sigma;sigM4

gle4=lm.4$df.residual;gle4
SCT=var(Invext)*(n-1);SCT
SCE4=sigM4^2*gle4;SCE4;sum(lm.4$residuals^2)
SCR4=sum((lm.4$fitted.values-mean(Invext))^2);SCR4
SCT-SCE4

summary(lm.3)$fstatistic

SCRs4=SCR4-SCR3;SCRs4
(SCRs4/1)/(sigM4^2)
(1-pt(sqrt((SCRs4/1)/(sigM4^2)),gle4))*2
1-pf((SCRs4/1)/(sigM4^2),1,gle4)

SCRs1+SCRs2+SCRs3+SCRs4;SCR4






###################
#######################
###paso 4
##la variable mas importante sin  mes
lm.5=lm(Invext~mes+Import+Ipc,data=ExpImp)
summary(lm.5)
AIC5=AIC(lm.5);AIC5


R2Aj4=summary(lm.4)$adj.r.squared;R2Aj4
sigM4=summary(lm.4)$sigma;sigM4

gle4=lm.4$df.residual;gle4
SCT=var(Invext)*(n-1);SCT
SCE4=sigM4^2*gle4;SCE4;sum(lm.4$residuals^2)
SCR4=sum((lm.4$fitted.values-mean(Invext))^2);SCR4
SCT-SCE4

summary(lm.3)$fstatistic

SCRs4=SCR4-SCR3;SCRs4
(SCRs4/1)/(sigM4^2)
(1-pt(sqrt((SCRs4/1)/(sigM4^2)),gle4))*2
1-pf((SCRs4/1)/(sigM4^2),1,gle4)

SCRs1+SCRs2+SCRs3+SCRs4;SCR4
############################
###parentesis de maxima logverosimilitud
fynorm=function(v){
  sig2=v[1];b0=v[2];b1=v[3]
  ind=b0+b1*mes
  l=log((2*pi*sig2)^(n/2))+0.5/sig2*sum((Invext-ind)^ 2)
}
b=c(1,1,1)
nlminb(b,fynorm)

summary(lm(Invext~mes))










#####################################
#####multicolinealidad
install.packages("car")
library(car)

vif(lm(Invext~mes+Import+Export,Data=ExpImp))




install.packages("lmridge")
library(lmridge)

lmridge(Invext~mes+Import+Export+Ipc, data=ExpImp,K=0.1)


#####descripcion  de la correlacion 
pairs(ExpImp,gap=0.5)
x=cbind(mes,Import,Export,Ipc);x=as.matrix(x)
cor(x)
x<-model.matrix(lm(Invext~mes+Import+Export+Ipc-1))

#scale(x,center=TRUE,scale=TRUE)
Cmes=scale(mes,F,T);CImport=scale(Import,F,T);CExport=scale(Export,F,T);cIpc=scale(Ipc,F,T)
x=cbind(Cmes,CImport,CExport);x=as.matrix(x)
sqrt(max(eigen(t(x)%*%x)$values)/eigen(t(x)%*%x)$values)

#convirtiendo  el vector unitario.
Cmes = mes/sqrt(sum(mes^2))
CImport = Import/sqrt(sum(Import^2))
CExport = Export/sqrt(sum(Export^2))
cIpc = Ipc/sqrt(sum(Ipc^2))
x=cbind(Cmes,CImport,CExport);x=as.matrix(x)
sqrt(max(eigen(t(x)%*%x)$values)/eigen(t(x)%*%x)$values)

##factor inflacion de la vaianza
library(car)
vif(lm(Invext~mes+Import+Export,data = ExpImp))

library(MASS)

library(lmridge)
betRidg3=lmridge(Invext~mes+Import+Export,data = ExpImp,scaling = "sc",k=0.1);betRidg3  
lm(Invext~mes+Import+Export,data = ExpImp)$coefficients
plot.lmridge(lmridge(Invext~mes+Import+Export,data = ExpImp,k=seq(0,1,0.1)))


Cmes=(mes-mean(mes))/sqrt(var(mes)*(n-1))
CImport=(Import-mean(Import))/sqrt(var(Import)*(n-1))
CExport=(Export-mean(Export))/sqrt(var(Export)*(n-1))


Cmes=(mes-mean(mes))/sqrt(sum(mes^2)/(n-1))
lm(Invext~Cmes)$coefficients




lm.3=lm(Invext~mes+Import+Export,data = ExpImp)
lm.3a=lm(mes~Import+Export,data = ExpImp)
RC3a<-summary(lm.3a)$r.squared;RC3a
vifmes=1/(1-RC3a)


lm.3b=lm(Import~mes+Export,data = ExpImp)
RC3b<-summary(lm.3b)$r.squared;RC3b
vifImport=1/(1-RC3b)

lm.3c=lm(Export~mes+Import,data = ExpImp)
RC3c<-summary(lm.3c)$r.squared;RC3c
vifExport=1/(1-RC3c)
vif=c(vifmes,vifImport,vifExport);vif




## graficas de dispercion
pairs(ExpImp,gap=0.5)

