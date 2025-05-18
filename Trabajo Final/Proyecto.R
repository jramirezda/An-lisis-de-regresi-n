# Cargar la librería readxl
library(readxl)

# Ruta al archivo de Excel
ruta <- "Parcial 3 - II 2023.xlsx"

# Leer la primera hoja del archivo Excel
data <- read_excel(ruta, sheet = 1)
data <- data[1:50, ]
# Mostrar los primeros registros para verificar que se ha cargado correctamente
head(data)

# Corregir la fórmula del modelo lineal
modelo <- lm(Demanda ~ pedidos + PUP + `Venta Neta Final` + `% Dsc. Catalogo`, data = data)

resumen <- summary(modelo)

# Obtener los valores p
p_values <- resumen$coefficients[, 4]
significancia_10 <- p_values <= 0.1
summary(resumen$coefficients[significancia_10, ])

# Obtener los residuos del modelo
residuos <- residuals(modelo)
plot(residuos, type = "l")

# Función para realizar la prueba de Park
prueba_park <- function(variable) {
  # Eliminar valores problemáticos antes de realizar la prueba de Park
  data_clean <- data[!is.na(data$`% Dsc. Revista`), ]
  
  # Regresar los residuos al cuadrado contra la variable
  plot(variable, residuos^2, main = paste("Gráfico de errores al cuadrado vs", deparse(substitute(variable))))
  
  # Realizar la prueba de Park
  park_test <- lm(log(residuos^2) ~ log(variable), data = data_clean)
  summary(park_test)
}

# Realizar la prueba de Park para cada variable independiente
prueba_park(data$pedidos)
prueba_park(data$PUP)
prueba_park(data$`Venta Neta Final`)
prueba_park(data$`% Dsc. Catalogo`)



# Definir la variable de ponderación
ponderacion <- data$`% Dsc. Catalogo`

# Realizar la regresión de mínimos cuadrados ponderados
modelo_ponderado <- lm(Demanda ~ pedidos + PUP + `Venta Neta Final` + `% Dsc. Catalogo`, data = data, weights = ponderacion)
modelo_ponderado <- glm(Demanda ~ pedidos + PUP + `Venta Neta Final` + `% Dsc. Catalogo`, data = data, weights = ponderacion , family = gaussian)

resumen_ponderado <- summary(modelo_ponderado)

# Obtener los valores p
p_values_ponderado <- resumen_ponderado$coefficients[, 4]
significancia_10_ponderado <- p_values_ponderado <= 0.1
summary(resumen_ponderado$coefficients[significancia_10_ponderado, ])

# Obtener los residuos del modelo ponderado
residuos_ponderados <- residuals(modelo_ponderado)
plot(ponderacion, residuos_ponderados^2, main = "Gráfico de errores al cuadrado vs Porcentaje Descuento Catalogo")

# Realizar la prueba de Park solo para la variable "Porcentaje Descuento Catalogo"
park_test_ponderado <- lm(log(residuos_ponderados^2) ~ log(data$`% Dsc. Catalogo`))
summary(park_test_ponderado)

#############
########################################
########################################

# Crear el conjunto de datos
datos <- data.frame(
  poblacion = c(500, 1200, 100, 400, 500, 300),
  reclamaciones = c(42, 37, 1, 101, 73, 14),
  tamano_coche = c("Pequeño", "Mediano", "Grande", "Pequeño", "Mediano", "Grande"),
  grupo_edad = c(1, 1, 1, 2, 2, 2)
)
mean(datos$reclamaciones);var(datos$reclamaciones)
# Ajustar el modelo Poisson sin offset con poblacion
modelo_poisson <- glm(reclamaciones ~  as.factor(tamano_coche) + as.factor(grupo_edad)+poblacion, data = datos, family = "poisson")
summary(modelo_poisson)

#sin poblacion 
modelo_poisson <- glm(reclamaciones ~  as.factor(tamano_coche) + as.factor(grupo_edad), data = datos, family = "poisson")
summary(modelo_poisson)


# Ajustar el modelo Poisson con offset
modelo_poisson_offset <- glm(reclamaciones ~ as.factor(tamano_coche) + as.factor(grupo_edad)+ offset(log(poblacion)), data = datos, family = "poisson")
summary(modelo_poisson_offset)

nueva_data <- data.frame(
  tamano_coche = "Pequeño",
  grupo_edad = 1,
  poblacion = 1000
)

# Sin offset
prob_sin_offset <- predict(modelo_poisson, newdata = nueva_data, type = "response")

# Con offset
prob_con_offset <- predict(modelo_poisson_offset, newdata = nueva_data, type = "response")



ppois(50,prob_sin_offset)
ppois(50,prob_con_offset)
lam=prob_con_offset*log(1000)
ppois(50/log(1000),lam)




# Ruta al archivo de Excel
ruta <- "C:/Users/SANTIAGO/Desktop/Parcial 3 - II 2023.xlsx"

# Leer la primera hoja del archivo Excel
P3 <- read_excel(ruta, sheet = 3)
#modificacion de la tabla
P3$TipoRiesgo1 <- ifelse(P3$TipoRiesgo == "A" | P3$TipoRiesgo == "B", 1, 0)
P3$TipoRiesgo1<-as.factor(P3$TipoRiesgo1)
P3$Comp<- as.factor(P3$Comp)
P3$Competencia<-as.factor(P3$Competencia)
#eliminar TipoRiesgo
P3 <- P3[, -which(names(P3) == "TipoRiesgo")]
P3 <- P3[, -which(names(P3) == "Comp")]
#lectura de la tabla
str(P3)
coef3<-glm(TipoRiesgo1~.,data = P3,family = binomial)$coef
summary(coef3)
#predicion
wa<-c(1,10000,0.2,0.05,1,0,15)

pi<-t(wa)%*%coef3

1/(1+exp(-pi))