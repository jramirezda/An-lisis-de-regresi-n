#Clase regresión.Jhon Ramirez 
#parametros bivariada
Mux=5
Muy=0.5
Sigmax=3
Sigmay=0.3
Sigmaxy=-0.72
cvy=Sigmay/Muy;cvy

## parametros regresión
Rhoxy=Sigmaxy/(Sigmax*Sigmay);Rhoxy
Beta1=Sigmaxy/(Sigmax^2);Beta1
Beta0=Muy-Beta1*Mux;Beta0
Sigma2=Sigmay^2*(1-Rhoxy^2);Sigma2
Sigma=sqrt(Sigma2);Sigma

## estimmación de ybarra
varybarra=Sigmay^2/n; varybarra
coefvarybarra=sqrt(varybarra)/Muy;coefvarybarra

##Simulaciín de  extraacción de muestras marginales de Y
m=1000
Ybarra=NULL;varYbarra=NULL
varmuestraly=NULL#S^2_y
muestras=NULL
coefvarybarrae=NULL
varybarrae=NULL

for (j in 1:m) {
  ys=rnorm(n,mean=Muy,sd=Sigmay)
  muestras=rbind(muestras,ys[])#concadena horizontal
  Ybarra[j]=mean(ys)
  varmuestraly=var(ys)
  varybarrae[j]=var(ys)/n
  coefvarybarrae=sqrt(varybarrae)/Ybarra[j]
}

apply(muestras,2,mean)
apply(muestras,2,var)
cbind(muestras,Ybarra)
cbind(Muy,Ybarra,Muy/Ybarra)
cbind(varybarra,varybarrae,varybarra/varybarrae)
cbind(coefvarybarra,coefvarybarrae,coefvarybarra/coefvarybarrae)

par(mfrow=c(3,3))
apply(muestras, 2,hist)
par(mfrow=c(1,1))


#Esenario  de muestreo
##fijar las x
X=1:9
n=length(X);n
xbar=mean(X);Sxx=sum((X-xbar)^2)
VarBeta0e=Sigma2*(1/n+xbar^2/Sxx);VarBeta0e
VarBeta1e=Sigma2/Sxx;VarBeta1e
covBeta0eBeta1e=-Sigma2*xbar/Sxx;covBeta0eBeta1e

#calculamos el coeficiente de varianción

coefvarBeta1e=sqrt(VarBeta1e)/Beta1;coefvarBeta1e
coefVarBeta0e=sqrt(VarBeta0e)/Beta0;coefVarBeta0e

## con muchas muestras vamos a ver como se comporta
m=1000
Beta0e=NULL
Beta1e=NULL
VareBeta0e=NULL
VareBeta1e=NULL
Sigma2e=NULL
Ybarra=NULL
vareY=NULL
varYbarrae=NULL
coefvarybarra=NULL
muestras=NULL
varmuestraly=NULL

for (j in 1:m){
  E=rnorm(n,mean = 0,sd=Sigma)
  ys=Beta0+Beta1*X+E
  muestras=rbind(muestras,ys[])
  reg=summary(lm(ys ~X))
  Sigma2e[j]=(reg$sigma)^2
  Beta0e[j]=reg$coefficients[1]
  Beta1e[j]=reg$coefficients[2]
  Ybarra[j]=mean(ys)
  varmuestraly[j]=var(ys)
  varYbarrae[j]=var(ys)/n
}

apply(muestras, 2, var)
Beta0+Beta1*X

apply(muestras,2,mean)

par(mfrow=c(3,3))
apply(muestras, 2,hist)
par(mfrow=c(1,1))

cbind(muestras)
cbind(Ybarra,vareY,varYbarrae,sqrt(varYbarrae)/Ybarra)

mean(Ybarra)
var(Ybarra)

mean(Beta0e)
Beta0

var_Beta0=Sigma2*(1/n+(xbar^2)/Sxx)
var(Beta0e)


## con muchas muestras vamos a ver como se comporta
m=1000
Beta0e=NULL
Beta1e=NULL
VareBeta0e=NULL
VareBeta1e=NULL
Sigma2e=NULL
Ybarra=NULL
vareY=NULL
varYbarrae=NULL
coefvarybarra=NULL
muestras=NULL
varmuestraly=NULL

for (j in 1:m){
  E=rnorm(n,mean = 0,sd=Sigma)
  ys=Beta0+Beta1*X+E
  muestras=rbind(muestras,ys[])
  reg=summary(lm(ys ~X))
  Sigma2e[j]=(reg$sigma)^2
  Beta0e[j]=reg$coefficients[1]
  Beta1e[j]=reg$coefficients[2]
  Ybarra[j]=mean(ys)
  varmuestraly[j]=var(ys)
  varYbarrae[j]=var(ys)/n
  VareBeta0e[j]=(reg$sigma)^2*(1/n+xbar^2/Sxx)
  VareBeta1e[j]=(reg$sigma)^2/Sxx
}

## Calculo y comparacion entre  los parametros de  Beta0

c(mean(Beta0e) ,Beta0)

var_Beta0=Sigma2*(1/n+(xbar^2)/Sxx)
var(Beta0e)
cbind(VareBeta0e,VarBeta0e)
cbind(Beta0,Beta0e)
cv_Beta0e=sqrt(var_Beta0)/Beta0;cv_Beta0e

cbind(Beta0,Beta0e,Beta0/Beta0e)
cbind(cv_Beta0e, sqrt(VareBeta0e)/Beta0e)

## Calculo con Beta1
c(mean(Beta0e) ,Beta0)

var_Beta1=Sigma2/Sxx;var_Beta1
var(Beta1e)
cbind(VareBeta1e,VarBeta1e)
cbind(Beta1,Beta1e)
cv_Beta1e=sqrt(var_Beta1)/Beta1;cv_Beta1e

cbind(Beta1,Beta1e,Beta1/Beta1e)
cbind(cv_Beta1e, min(abs(sqrt(VareBeta1e)/Beta1e)))

## Estadistico  de prueba de beta1=0
cbind(Beta1e/sqrt(VareBeta1e))


##prueba de hipótesis para Beta1 
d=rep(0,m)
tc=abs(Beta1e+0.08)/sqrt(VarBeta1e)      
d[tc>qt(0.95,7)]=1
pval=2*(1-pt(tc,7));pval
cbind(Beta1,Beta1e,cv_Beta1e,tc,d,pval)
sum(d)
mean(pval)
?pt

mean(pval[d==1])



## prueba de la significancia de  beta1
d=rep(0,m)
tc=abs(Beta1e)/sqrt(VarBeta1e)      
d[tc>qt(0.95,7)]=1
pval=2*(1-pt(tc,7));pval
cbind(Beta1,Beta1e,cv_Beta1e,tc,d,pval)
sum(d)
mean(pval)
?pt

mean(pval[d==1])

install.packages("sqldf")
library(sqldf)


x = c(4, 3, 6, 5, 12, 22, 2, 2, 10)
y = c(3, 2, 4, 4, 5, 6, 2, 1, 5)

xy = data.frame(x, y)

xy_ordered = sqldf("SELECT * 
                    FROM xy 
                    ORDER BY x ASC")


xy_ordered

lx=log(x)
reg.1=lm(y~lx)
summary(reg.1)

plot(x,y)
?qt


beta0e=0.2439
beta1e=1.9672

muestyx0=beta0e+beta1e*xy_ordered$x
tq=qt(0.025,7)
sigma2e=