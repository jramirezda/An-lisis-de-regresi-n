## Jhon Alejandro Ramírez Daza C.C. 1000118906
# Definir los valores
x <- c(1:9)
sigma <- 0.18
B1h0 <- -0.08
alpha <- 0.05
effect_size <- NULL  # Tamaño del efecto, lo calcularemos
sample_size <- NULL  # Tamaño de la muestra, lo calcularemos
num_predictors <- 2  # Número de predictores (incluyendo B1)

# Función para calcular la potencia
calculate_power <- function(effect_size, sample_size, alpha, num_predictors) {
  # Simular datos bajo la hipótesis alternativa
  set.seed(123)  # Establecer semilla para reproducibilidad
  x <- c(1:9)
  error <- rnorm(9, mean = 0, sd = sigma)  # Error aleatorio
  y <- B1h0 + effect_size * x + error  # Modelo de regresión lineal con efecto en B1
  
  # Ajustar el modelo de regresión lineal
  lm_model <- lm(y ~ x)
  
  # Realizar la prueba de hipótesis sobre B1
  summary_lm <- summary(lm_model)
  t_stat <- summary_lm$coefficients["x", "t value"]
  
  # Calcular el valor crítico
  critical_value <- qt(1 - alpha/2, df = sample_size - num_predictors)
  
  # Calcular la potencia
  power <- 1 - pt(critical_value - abs(t_stat), df = sample_size - num_predictors)
  
  return(power)
}

# Establecer un tamaño de muestra inicial (puedes ajustarlo según tus necesidades)
sample_size <- 100

# Calcular la potencia para diferentes tamaños del efecto
effect_sizes <- seq(0.01, 0.5, by = 0.01)
power_values <- numeric(length(effect_sizes))

for (i in 1:length(effect_sizes)) {
  power_values[i] <- calculate_power(effect_sizes[i], sample_size, alpha, num_predictors)
}

# Encontrar el tamaño de muestra necesario para una potencia específica (por ejemplo, 0.8)
target_power <- 0.8
required_sample_size <- NULL

for (i in 1:length(effect_sizes)) {
  if (power_values[i] >= target_power) {
    required_sample_size <- sample_size
    break
  }
  sample_size <- sample_size + 1
}

# Imprimir el resultado
cat("Tamaño de muestra necesario para alcanzar una potencia de", target_power, "es:", required_sample_size, "\n")

# Graficar la relación entre el tamaño del efecto y la potencia
plot(effect_sizes, power_values, type = "l", xlab = "Tamaño del Efecto", ylab = "Potencia")
abline(h = target_power, col = "red", lty = 2)
