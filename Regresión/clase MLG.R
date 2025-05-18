library(survival)

pulmon<-read.table("DatosPulmon.txt", header=TRUE)
attach(pulmon)
# Si W es N(0,1), entonces Y=ln T = muY +sigmaY * W, es N(muY,sigmaY)
Y=log(Tiempo)
n=length(Tiempo)
tt=1:1500; yy = log(tt)


##################################################
##########  Ajuste de una Densidad ###############
##################################################

# Gráfica del Tiempo con base en una densidad kernell. No paramétrica
hist(Tiempo,prob=TRUE)
lines(density(Tiempo,bw=100),col=2,lwd=3)

# Gráfica del Tiempo con base en una normal y una lognormal
# f LogNormal. f(t) = 1/( t*Sig* sqrt(2*pi) ) * exp( -0.5*( (ln t - Mu)/Sig )^2 )
hist(Tiempo,prob=TRUE,ylim=c(0,0.003))
flognorm = ( 1/sqrt(2*pi*sd(Y)) ) * exp (-0.5*(yy-mean(Y))^2/var(Y) )/tt
fnormT = (1/sqrt(2*pi*var(Tiempo))) * exp (-0.5*(tt-mean(Tiempo))^2/var(Tiempo) )
lines(tt,flognorm,col=3,lwd=3)
lines(tt,fnormT,col=2,lwd=3)

# Gráfica de Y = ln T, asumiento Y normal
fnorm = (1/sqrt(2*pi*sd(Y))) * exp (-0.5*(yy-mean(Y))^2/var(Y) )
hist(Y,prob=TRUE,xlim=c(0,9))
lines(yy,fnorm,col=2, lwd=3)


#### Ajustar Y a una logística
# Si W es logis(0,1), entonces Y=ln T = muY +sigmaY * W, es logis(muY,sigmaY)
# f(y)=exp( (y - mu)/sigma ) / (sigma*( 1+exp( (y - mu)/sigma ) )^2  )
# Estimar parámetros distribución logística de Y = ln T
fLog=function(b){
  muY=b[1];sigY=b[2];est=(Y-muY)/sigY
  flog= -sum( est - log( (1+exp(est))^2*sigY ) )
}
optim(c(1,2.8) , fLog)
muL=5.51;sigL=0.467;est=(yy-muL)/sigL
flog = exp(est)/( (1+exp(est))^2*sigL )
# Gráfica de una logística para Y comparada frente a una normal
hist(Y,prob=TRUE)
lines(yy,flog,col=2,lwd=4)
lines(yy,fnorm,col=3, lwd=3) # Comparar frente a normal


# Estimar parámetros distribución log-logística de T = exp (Y)
# f(t)=exp( (ln(t) - mu)/sigma ) / (t*sigma*( 1+exp( (ln(t) - mu)/sigma ) )^2  )
fLoglog=function(b){
  muY=b[1];sigY=b[2];est=(Y-muY)/sigY
  flog= -sum( est - log( (1+exp(est))^2*sigY ) )
}
optim(c(10,10) , fLoglog)
muLogl=5.5;sigLogl=0.46;est=(yy-muLogl)/sigLogl
floglog = exp(est)/( (1+exp(est))^2*sigL *tt )
hist(Tiempo,prob=TRUE,ylim=c(0,0.003))
lines(tt,floglog,col=3,lwd=4)

# Gráfica de LogLogística y LogNormal
hist(Tiempo,prob=TRUE,ylim=c(0,0.003))
lines(tt,floglog,col=2,lwd=4)
lines(tt,flognorm,col=3,lwd=4)


##################################################
################  Regresión ######################
##################################################


##### Regresión LogLogística. 
# f(t)=exp( (ln(t) - mu)/sigma ) / (t*sigma*( 1+exp( (ln(t) - mu)/sigma ) )^2  )
# mu se puete tomar para el modelo de regresión, como mu=bet0 + bet1*X1+...+betkXk.
# Modelo con una sola variable independiente
fLogL<-function(b){
   bet0=b[1];bet1=b[2];sig=b[3]
   ind = bet0 + bet1*Sexo
   c = ( Y - bet0 - bet1*Sexo)/sig 
   -sum( c - log(sig) - 2*(  log(1+exp(c) ) ) )
}
RlogL=optim(c(1,1,1),fLogL,hessian=T);RlogL
RlogL$par
sqrt(diag(solve(RlogL$hessian)))
RlogL$par / sqrt(diag(solve(RlogL$hessian)))
summary(survreg(Surv(Tiempo,rep(1,n))~Sexo,dist="loglogistic"))
# Los hombres tienen un tiempo de muerte un 25,3% menor que el de las mujeres

# Grágica de la densidad
# Función para la densidad loglogística por sexo y total
fLogL=function(Sexo){
	ind = 5.661 -0.2535*Sexo; ex1 = exp( (yy - ind)/0.462 )
	fLogL1 = ex1 / (tt*0.462 * ( ( 1 + ex1 )^2 ) )
}
hist(Tiempo,prob=T,ylim=c(0,0.003))
lines(tt,fLogL(1),col=2,lwd=3,lty=1)
lines(tt,fLogL(0),col=2,lwd=3,lty=4)
lines(tt,floglog,col=3,lwd=4)


# Regresión MCO para Log(Tiempo) y regresión LOGNORMAL para T
summary(lm(Y~Sexo))
summary(survreg(Surv(Tiempo,rep(1,n))~Sexo,dist="lognorm"))

# Grágica de la densidad LogNormal por sexo y total.
# Función para la regresión lognormal
fLogNorm=function(Sexo){
	ind = 5.6092 -0.3082*Sexo;
	fLogN = 1/( tt*0.892* sqrt(2*pi) ) * exp( -0.5*( (yy - ind)/0.892 )^2 )
}
hist(Tiempo,prob=T,ylim=c(0,0.0035))
lines(tt,fLogNorm(1),col=2,lwd=3,lty=1)
lines(tt,fLogNorm(0),col=2,lwd=3,lty=4)
lines(tt,flognorm,col=3,lwd=4)


# Gráfica de la densidad de T para hombres. LogLogística y LogNormal
hist(Tiempo,prob=T,ylim=c(0,0.0035))
lines(tt,fLogL(1),col=2,lwd=3,lty=1)
lines(tt,fLogNorm(1),col=3,lwd=3,lty=1)


# Esperanza condicional estimada 
#mean(Tiempo)
by(Tiempo,Sexo,mean)
# LogLogística
exp(5.6612245 - 0.2534744*0)
exp(5.6612245 - 0.2534744*1)
# Lognormal
exp(5.6092 - 0.3082*0)
exp(5.6092 - 0.3082*1)



####### Regresión GAMMA. f(y)=1/Gamma(alpha) *(alpha*y/mu)âlpha*exp(-alpha*y/mu)*1/y
# Función de enlace "inv"
fGamInv<-function(b){
   bet0=b[1];bet1=b[2];alpha=b[3]
   Mu = bet0 + bet1*Sexo
   -sum( -log(gamma(alpha)) + alpha*( log( alpha*Tiempo/ (1/Mu)  ) ) - alpha*Tiempo / (1/Mu) )
}
RGamInv=optim(c(0.005,0.005,0.5),fGamInv,hessian=T);RGamInv$par
sqrt(diag(solve(RGamInv$hessian)))
summary(glm(Tiempo~Sexo,family=Gamma(link="inverse"), data=pulmon))$coefficients


# Gráfica de la densidad Gamma de T.
hist(Tiempo,prob=T,ylim=c(0,0.0035))
MuInv = 0.00295 + 0.000580 * (0)
MuInv1 = 0.00295 + 0.000580 * (1)
fGamIn = 1/gamma(1.847) / tt * (1.847*tt / (1/MuInv) )^1.84 * exp(-1.847*tt / (1/MuInv) )
fGamIn1 = 1/gamma(1.847) / tt * (1.847*tt / (1/MuInv1) )^1.84 * exp(-1.847*tt / (1/MuInv1) )
lines(tt,fGamIn,col=1,lwd=4)
lines(tt,fGamIn1,col=1,lwd=4,lty=2)

# Gráfica de la densidad de T para hombres. LogLogística, LogNormal y Gamma
hist(Tiempo,prob=T,ylim=c(0,0.0035))
lines(tt,fLogL(1),col=2,lwd=3,lty=1)
lines(tt,fLogNorm(1),col=3,lwd=3,lty=1)
lines(tt,fGamIn,col=1,lwd=4)


# Esperanza condicional estimada 
#mean(Tiempo)
by(Tiempo,Sexo,mean)
# LogLogística
exp(5.6613805 - 0.2535996*0)
exp(5.6613805 - 0.2535996*1)
# Lognormal
exp(5.609 - 0.3084*0)
exp(5.609 - 0.3084*1)

# Gamma. ling ="inv"
(0.0029502158 + 0.0005800026*0)^(-1)
(0.0029502158 + 0.0005800026*1)^(-1)





# Regresión Poisson

Muertes=c(32,104,206,186,102,2,12,28,28,31)
Pob=c(52407,43248,28612,12663,5317,18790,10673,5710,2585,1462)
Fuma=c(1,1,1,1,1,0,0,0,0,0)
Edad=c(40,50,60,70,80,40,50,60,70,80)
EdadC=Edad^2
Datos=data.frame(Edad,Muertes,Pob,Fuma)
Datos