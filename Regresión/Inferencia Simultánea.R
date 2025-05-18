# Programa de simulación sobre inferencia simultánea de los parámetros del modelo MCO

# Tamaño de muestra
n=100

# Generación de Variables Independientes
x1=rnorm(n,mean=0, sd=1)
x2=rnorm(n,mean=15, sd=2)
x3=rnorm(n,mean=100, sd=10)

# Fijar Parámetros
beta0=5;beta1=10;beta2=4;beta3=2;Sig=5;

#Número réplicas
m=20000
bet0=NULL;bet1=NULL;bet2=NULL;bet3=NULL;
seb0=NULL;seb1=NULL;seb2=NULL;seb3=NULL;

for (j in 1:m){
   
   e=rnorm(n,mean=0, sd=Sig)
   y=beta0+beta1*x1+beta2*x2+beta3*x3+e
   lm1=lm(y~x1+x2+x3)

   bet0[j]=lm1$coefficients[1] # Extraer beta0 estimado
   bet1[j]=lm1$coefficients[2] # Extraer beta1 estimado
   bet2[j]=lm1$coefficients[3] # Extraer beta2 estimado
   bet3[j]=lm1$coefficients[4] # Extraer beta3 estimado

   seb0=summary(lm1)$coefficients[5];  # extraer error estándar de beta0
   seb1=summary(lm1)$coefficients[6];  # extraer error estándar de beta1
   seb2=summary(lm1)$coefficients[7];  # extraer error estándar de beta2
   seb3=summary(lm1)$coefficients[8];  # extraer error estándar de beta3
      
}


# Inferencia simultánea para beta1 y beta2

# Nivel de significancia y percentil
alpha=0.05
q=qt(1-alpha/4,n)

# Límetes para beta1
LI1=bet1-q*seb1
LS1=bet1+q*seb1

# Límetes para beta2
LI2=bet2-q*seb2
LS2=bet2+q*seb2

# Nivel de significancia empírico para beta1
alp1=rep(1,m)
alp1[LI1>beta1 | LS1<beta1]<-0
cbind(beta1,LI1,LS1,alp1)
1-sum(alp1)/m

# Nivel de significancia empírico para beta2
alp2=rep(1,m)
alp2[LI2>beta2 | LS2<beta2]<-0
cbind(beta2,LI2,LS2,alp2)
1-sum(alp2)/m

# Nivel de significancia empírico para beta1 y beta2
cbind(beta1,LI1,LS1, beta2,LI2,LS2,alp1, alp2, alp1*alp2)
1-sum(alp1*alp2)/m







