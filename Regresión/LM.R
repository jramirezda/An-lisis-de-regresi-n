Lm <- function(X, y) {
  n <- length(y)
  p <- dim(X)[2]
  k <- p - 1
  I <- diag(n)
  XTX <- t(X) %*% X  # Calcular la matriz t(X) %*% X una vez
  
  if (det(XTX) == 0) {
    stop("La matriz XTX es singular. No se puede calcular la inversa.")
  }
  
  C <- solve(XTX)
  BetaE <- C %*% t(X) %*% y
  SSer <- t(y) %*% y - t(BetaE) %*% t(X) %*% y
  sigma2E <- SSer / (n - p)
  
  ## Tabla ANOVA
  SSr <- t(BetaE) %*% t(X) %*% y - ((sum(y))^2 / n)
  MSr <- SSr / k
  SSer <- t(y) %*% y - t(BetaE) %*% t(X) %*% y
  sigma2E <- SSer / (n - p)
  
  fcal <- MSr / sigma2E
  pval <- 1 - pf(fcal, k, n - k - 1)
  # Formatear el valor p en notación científica con 3 decimales
  pval<- sprintf("%.3e", pval)
  
  resultados <- data.frame(
    Fuente_de_variación = c("Suma de cuadrados", "Grados de libertad", "Cuadrado medio", "F_0", "Valor p"),
    regresion = c( SSr, k, MSr, fcal, pval),
    errores = c( SSer, n - k - 1, sigma2E, "", ""),
    total = c( SSer + SSr, n - 1, "", "", "")
  )
  resultados <- t(resultados)
  return(resultados)
}

# Crear un conjunto de datos de ejemplo
set.seed(123)
X <- matrix(rnorm(100), ncol = 2)
y <- X[, 1] + 2 * X[, 2] + rnorm(50)

# Llamar a la función Lm
resultados_anova <- Lm(X, y)

# Imprimir la tabla ANOVA
print(resultados_anova)


# Crear un conjunto de datos de ejemplo
datos <- data.frame(
  Publicidad = c(100, 200, 300, 400, 500),
  Precio = c(10, 15, 20, 25, 30),
  Ventas = c(500, 600, 750, 900, 1100)
)

# Llamar a la función Lm
resultados_anova2 <- Lm(as.matrix(datos[, c("Publicidad", "Precio")]), as.numeric(datos$Ventas))

# Imprimir la tabla ANOVA
print(resultados_anova2)


Coef_regresion<- function(X, y) {
  n <- length(y)
  p <- dim(X)[2]
  k <- p - 1
  I <- diag(n)
  XTX <- t(X) %*% X  # Calcular la matriz t(X) %*% X una vez
  
  if (det(XTX) == 0) {
    stop("La matriz XTX es singular. No se puede calcular la inversa.")
  }
  
  C <- solve(XTX)
  BetaE <- C %*% t(X) %*% y
  
  return(BetaE)
}


# Crear un conjunto de datos de ejemplo
set.seed(123)
X <- matrix(rnorm(100), ncol = 2)
y <- X[, 1] + 2 * X[, 2] + rnorm(50)

resultados_betas <- Coef_regresion(X, y)

# Imprimir la tabla betas
print(resultados_betas)



# Crear un conjunto de datos de ejemplo
datos <- data.frame(
  Publicidad = c(100, 200, 300, 400, 500),
  Precio = c(10, 15, 20, 25, 30),
  Ventas = c(500, 600, 750, 900, 1100)
)

# Llamar a la función Lm
resultados_betas2 <- Coef_regresion(as.matrix(datos[, c("Publicidad", "Precio")]), as.numeric(datos$Ventas))

# Imprimir la tabla ANOVA
print(resultados_betas2)
