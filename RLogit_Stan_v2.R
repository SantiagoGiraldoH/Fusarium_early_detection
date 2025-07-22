# Cargar librerías
library(rstan)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(tidyr)
library(bayesplot)
library(coda)
library(gridExtra)


# Cargar datos
data_raw <- read.csv("datos_RLB.csv")
data_raw <- data_raw[, c("etiqueta", setdiff(names(data_raw), "etiqueta"))]

# Extraer componentes
etiquetas <- data_raw[,1]          # Fusarium (0/1)
dpi_values <- data_raw[,2]         # Días post-inoculación
spectral_data <- data_raw[,3:ncol(data_raw)]  # Datos espectrales

dia_actual <- 10 #modificar para cada día
indices_dia <- which(dpi_values == dia_actual)

# Datos del día específico
y <- etiquetas[indices_dia]
X <- as.matrix(spectral_data[indices_dia, ])
N <- length(y)
P <- ncol(X)

# Verificar
cat("Día:", dia_actual, "\n")
cat("Observaciones:", N, "\n") 
cat("Longitudes de onda:", P, "\n")
cat("Casos Fusarium:", sum(y), "\n")
cat("Casos No-Fusarium:", sum(1-y), "\n")

wavelengths <- seq(350, 2500, length.out = P)

#================== === CONFIGURACIÓN DE PRIORS PARA ETAPA TEMPRANA DE FUSARIUM ===

setup_early_fusarium_priors <- function(wavelengths, P) {
  # Inicializar con priors restrictivos (sesgado hacia exclusión)
  alpha_pi <- rep(1, P)    # Base neutral  
  beta_pi <- rep(3, P)     # Sesgado hacia exclusión (E[π] = 0.25)
  
  # === BANDAS CRÍTICAS BASADAS EN LITERATURA CIENTÍFICA ===
  
  # CLOROFILA (680-733 nm) - CRÍTICA para detección temprana
  #    Justificación: Degradación temprana de clorofila antes de síntomas visibles
  chlorophyll <- which(wavelengths >= 680 & wavelengths <= 733)
  alpha_pi[chlorophyll] <- 4    # E[π] = 0.8 (muy alta probabilidad)
  beta_pi[chlorophyll] <- 1
  
  # ABSORCIÓN DE AGUA (1400-1450, 1900-1950 nm) - MUY IMPORTANTE
  #    Justificación: Estrés hídrico temprano por infección vascular
  water_1400 <- which(wavelengths >= 1400 & wavelengths <= 1450)
  water_1900 <- which(wavelengths >= 1900 & wavelengths <= 1950)
  water_bands <- c(water_1400, water_1900)
  alpha_pi[water_bands] <- 3    # E[π] = 0.75
  beta_pi[water_bands] <- 1
  
  # CAROTENOIDES (500-533 nm) - Respuesta de estrés temprana
  #    Justificación: Aumento de carotenoides como respuesta a infección FOC
  carotenoids <- which(wavelengths >= 500 & wavelengths <= 533)
  alpha_pi[carotenoids] <- 2    # E[π] = 0.67
  beta_pi[carotenoids] <- 1
  
  # COMPUESTOS FENÓLICOS (300-400 nm) - Defensa temprana
  #    Justificación: Respuesta de defensa molecular temprana
  phenolic <- which(wavelengths >= 300 & wavelengths <= 400)
  alpha_pi[phenolic] <- 2       # E[π] = 0.67
  beta_pi[phenolic] <- 1
  
  # AGUA SWIR (2200-2300 nm) - Complementaria para estrés hídrico
  #    Justificación: Cambios en contenido hídrico celular
  swir_water <- which(wavelengths >= 2200 & wavelengths <= 2300)
  alpha_pi[swir_water] <- 2     # E[π] = 0.67
  beta_pi[swir_water] <- 1
  
  # RED-EDGE (700-750 nm) - Transición clorofila-NIR
  #    Justificación: Cambios en estructura foliar temprana
  red_edge <- which(wavelengths >= 700 & wavelengths <= 750)
  alpha_pi[red_edge] <- 2       # E[π] = 0.67
  beta_pi[red_edge] <- 1
  
  # Calcular expectativas para reporte
  expected_inclusion <- alpha_pi / (alpha_pi + beta_pi)
  
  cat("=== CONFIGURACIÓN DE PRIORS PARA ETAPA TEMPRANA ===\n")
  cat("Bandas con alta expectativa (>0.7):", sum(expected_inclusion > 0.7), "\n")
  cat("Bandas moderadas (0.5-0.7):", sum(expected_inclusion >= 0.5 & expected_inclusion <= 0.7), "\n")
  cat("Bandas con baja expectativa (<0.5):", sum(expected_inclusion < 0.5), "\n")
  cat("Expectativa media de inclusión:", round(mean(expected_inclusion), 3), "\n")
  
  # Detalles por región
  cat("\nDETALLE POR REGIÓN:\n")
  cat("Clorofila (680-733 nm):", length(chlorophyll), "bandas, E[π] =", 
      round(mean(expected_inclusion[chlorophyll]), 2), "\n")
  cat("Agua (1400-1450, 1900-1950 nm):", length(water_bands), "bandas, E[π] =", 
      round(mean(expected_inclusion[water_bands]), 2), "\n")
  cat("Carotenoides (500-533 nm):", length(carotenoids), "bandas, E[π] =", 
      round(mean(expected_inclusion[carotenoids]), 2), "\n")
  cat("Fenólicos (300-400 nm):", length(phenolic), "bandas, E[π] =", 
      round(mean(expected_inclusion[phenolic]), 2), "\n")
  
  return(list(
    alpha_pi = alpha_pi, 
    beta_pi = beta_pi,
    expected_inclusion = expected_inclusion,
    regions = list(
      chlorophyll = chlorophyll,
      water = water_bands,
      carotenoids = carotenoids,
      phenolic = phenolic,
      red_edge = red_edge,
      swir_water = swir_water
    )
  ))
}

prior_config <- setup_early_fusarium_priors(wavelengths, P)

##=======  CONFIGURACIÓN DE PRIORS PARA ETAPA MEDIA DE FUSARIUM ===

setup_middle_stage_priors <- function(wavelengths, P) {
  # Base moderadamente restrictiva - cambios fisiológicos se intensifican
  alpha_pi <- rep(1, P)    
  beta_pi <- rep(2, P)     # E[π] = 0.33 (menos restrictivo que etapa temprana)
  
  # === CAMBIOS FISIOLÓGICOS PROGRESIVOS ===
  
  #  CLOROFILA (680-733 nm) - CRÍTICA para amarillamiento progresivo
  # Justificación: Un síntoma inicial común es la aparición de una raya amarilla pálida débil en la base del pecíolo de la hoja más antigua. Esto es seguido por la clorosis foliar que progresa de las hojas inferiores a las superiores
  chlorophyll <- which(wavelengths >= 680 & wavelengths <= 733)
  alpha_pi[chlorophyll] <- 5    # E[π] = 0.83 (mayor que etapa temprana)
  beta_pi[chlorophyll] <- 1
  
  #  AGUA (1400-1450, 1900-1950, 2200-2300 nm) - Estrés hídrico severo
  # Justificación: Después de que el patógeno haya invadido sistémicamente los elementos del vaso xilemático con invasión apreciable del rizoma, se desarrolla una severa escasez de agua debido al taponamiento vascular
  water_1400 <- which(wavelengths >= 1400 & wavelengths <= 1450)
  water_1900 <- which(wavelengths >= 1900 & wavelengths <= 1950) 
  water_2200 <- which(wavelengths >= 2200 & wavelengths <= 2300)
  water_bands <- c(water_1400, water_1900, water_2200)
  alpha_pi[water_bands] <- 4    # E[π] = 0.8
  beta_pi[water_bands] <- 1
  
  #  ÁCIDO FUSÁRICO - Efectos en fotosíntesis (560-680 nm)
  # Justificación: El ácido fusárico (FA), una fitotoxina que es producida por F. oxysporum f. sp. cubense... ha sido señalado como la causa del síntoma de clorosis foliar. El daño del cloroplasto, la eficiencia fotoquímica reducida del fotosistema II (FV/Fmax), y una reducción concomitante en la fotosíntesis (asimilación neta de CO2) estuvieron asociados con la progresión del amarillamiento foliar
  photosynthesis <- which(wavelengths >= 560 & wavelengths <= 680)
  alpha_pi[photosynthesis] <- 4    # E[π] = 0.8
  beta_pi[photosynthesis] <- 1
  
  #  CAROTENOIDES (500-533 nm) - Respuesta de estrés intensificada
  # Justificación: En la etapa media, los carotenoides aumentan como respuesta al estrés
  carotenoids <- which(wavelengths >= 500 & wavelengths <= 533)
  alpha_pi[carotenoids] <- 3    # E[π] = 0.75
  beta_pi[carotenoids] <- 1
  
  #  ESTRUCTURA VASCULAR - NIR (800-1100 nm)
  # Justificación: En el día 12, después de que las plántulas de banano fueron inoculadas con Foc 4, los bordes de las hojas de banano se volvieron amarillos
  structural <- which(wavelengths >= 800 & wavelengths <= 1100)
  alpha_pi[structural] <- 3    # E[π] = 0.75
  beta_pi[structural] <- 1
  
  #  RED-EDGE (700-750 nm) - Transición crítica
  # Justificación: Cambios en la estructura foliar durante la transición
  red_edge <- which(wavelengths >= 700 & wavelengths <= 750)
  alpha_pi[red_edge] <- 3      # E[π] = 0.75
  beta_pi[red_edge] <- 1
  
  #  METABOLITOS SECUNDARIOS (300-400 nm) - Respuesta de defensa
  # Justificación: Intensificación de la respuesta de defensa de la planta
  phenolic <- which(wavelengths >= 300 & wavelengths <= 400)
  alpha_pi[phenolic] <- 2      # E[π] = 0.67
  beta_pi[phenolic] <- 1
  
  # Calcular expectativas para reporte
  expected_inclusion <- alpha_pi / (alpha_pi + beta_pi)
  
  cat("=== CONFIGURACIÓN DE PRIORS PARA ETAPA MEDIA (5-10 DPI) ===\n")
  cat("Bandas con alta expectativa (>0.7):", sum(expected_inclusion > 0.7), "\n")
  cat("Bandas moderadas (0.5-0.7):", sum(expected_inclusion >= 0.5 & expected_inclusion <= 0.7), "\n")
  cat("Bandas con baja expectativa (<0.5):", sum(expected_inclusion < 0.5), "\n")
  cat("Expectativa media de inclusión:", round(mean(expected_inclusion), 3), "\n")
  
  # Detalles por región
  cat("\nDETALLE POR REGIÓN:\n")
  cat("Clorofila (680-733 nm):", length(chlorophyll), "bandas, E[π] =", 
      round(mean(expected_inclusion[chlorophyll]), 2), "\n")
  cat("Fotosíntesis (560-680 nm):", length(photosynthesis), "bandas, E[π] =", 
      round(mean(expected_inclusion[photosynthesis]), 2), "\n")
  cat("Agua (multi-región):", length(water_bands), "bandas, E[π] =", 
      round(mean(expected_inclusion[water_bands]), 2), "\n")
  cat("Estructura vascular (800-1100 nm):", length(structural), "bandas, E[π] =", 
      round(mean(expected_inclusion[structural]), 2), "\n")
  
  return(list(
    alpha_pi = alpha_pi, 
    beta_pi = beta_pi,
    expected_inclusion = expected_inclusion,
    regions = list(
      chlorophyll = chlorophyll,
      photosynthesis = photosynthesis,
      water = water_bands,
      carotenoids = carotenoids,
      structural = structural,
      red_edge = red_edge,
      phenolic = phenolic
    )
  ))
}

prior_config <- setup_middle_stage_priors(wavelengths, P)


#============== CONFIGURACIÓN DE PRIORS PARA ETAPA TARDÍA DE FUSARIUM  ========

setup_late_stage_priors <- function(wavelengths, P) {
  # Base menos restrictiva - cambios visibles y severos
  alpha_pi <- rep(1, P)    
  beta_pi <- rep(1.5, P)   # E[π] = 0.40 (menos restrictivo que etapas anteriores)
  
  # === SÍNTOMAS VISIBLES Y DAÑO SEVERO ===
  
  #  CLOROSIS SEVERA Y MUERTE FOLIAR (680-733 nm, 560-680 nm) - EXTREMA
  # Justificación: "Los primeros síntomas externos del marchitamiento por Fusarium incluyen la clorosis y muerte de las hojas más viejas, que típicamente se doblan y colapsan contra el pseudotallo"
  chlorophyll_death <- which(wavelengths >= 560 & wavelengths <= 733)
  alpha_pi[chlorophyll_death] <- 6    # E[π] = 0.86 (máxima prioridad)
  beta_pi[chlorophyll_death] <- 1
  
  #  COLAPSO DEL SISTEMA HÍDRICO (1400-1450, 1900-1950, 2200-2300 nm) - CRÍTICO
  # Justificación: "A medida que la enfermedad progresa, las hojas más jóvenes se ven afectadas, se vuelven amarillas y se arrugan y todo el dosel comienza a consistir en hojas muertas o moribundas"
  water_collapse <- which((wavelengths >= 1400 & wavelengths <= 1450) |
                            (wavelengths >= 1900 & wavelengths <= 1950) |
                            (wavelengths >= 2200 & wavelengths <= 2300))
  alpha_pi[water_collapse] <- 5       # E[π] = 0.83
  beta_pi[water_collapse] <- 1
  
  # NECROSIS MARGINAL Y SENESCENCIA (750-850 nm) - NUEVA REGIÓN
  # Justificación: "F. oxysporum generalmente produce síntomas como marchitamiento, clorosis, necrosis, caída prematura de hojas, pardeo del sistema vascular, enanismo y damping-off"
  necrosis <- which(wavelengths >= 750 & wavelengths <= 850)
  alpha_pi[necrosis] <- 5             # E[π] = 0.83
  beta_pi[necrosis] <- 1
  
  #  SISTEMA VASCULAR COMPROMETIDO (800-1200 nm) - EXTENSO
  # Justificación: "internamente incluyen una decoloración marrón a rojo ladrillo del sistema vascular que va desde puntos discretos hasta grandes áreas confluentes"
  vascular_damage <- which(wavelengths >= 800 & wavelengths <= 1200)
  alpha_pi[vascular_damage] <- 4      # E[π] = 0.8
  beta_pi[vascular_damage] <- 1
  
  #  FOTOSÍNTESIS COLAPSO TOTAL (400-700 nm) - AMPLIO ESPECTRO
  # Justificación: "Daño del cloroplasto, eficiencia fotoquímica reducida del fotosistema II (FV/Fmax), y una reducción concomitante en la fotosíntesis"
  photosynthesis_collapse <- which(wavelengths >= 400 & wavelengths <= 700)
  alpha_pi[photosynthesis_collapse] <- 4  # E[π] = 0.8
  beta_pi[photosynthesis_collapse] <- 1
  
  #  ESTRÉS OXIDATIVO SEVERO (500-533 nm) - CAROTENOIDES
  # Justificación: Respuesta de estrés oxidativo máxima
  oxidative_stress <- which(wavelengths >= 500 & wavelengths <= 533)
  alpha_pi[oxidative_stress] <- 4     # E[π] = 0.8
  beta_pi[oxidative_stress] <- 1
  
  #  METABOLITOS TÓXICOS (300-450 nm) - EXTENSO
  # Justificación: "El ácido fusárico (FA), una fitotoxina que es producida por F. oxysporum f. sp. cubense... ha sido señalado como la causa del síntoma de clorosis foliar"
  toxin_response <- which(wavelengths >= 300 & wavelengths <= 450)
  alpha_pi[toxin_response] <- 3       # E[π] = 0.75
  beta_pi[toxin_response] <- 1
  
  #  WILTING Y COLAPSO ESTRUCTURAL (1100-1400 nm) - NUEVA REGIÓN
  # Justificación: "Las hojas marchitas también pueden romperse en el pecíolo y colgar del pseudotallo"
  structural_collapse <- which(wavelengths >= 1100 & wavelengths <= 1400)
  alpha_pi[structural_collapse] <- 3  # E[π] = 0.75
  beta_pi[structural_collapse] <- 1
  
  #  SWIR AMPLIO - Cambios en composición (2300-2500 nm)
  # Justificación: Cambios en lignina, celulosa y otros componentes estructurales
  swir_composition <- which(wavelengths >= 2300 & wavelengths <= 2500)
  alpha_pi[swir_composition] <- 2     # E[π] = 0.67
  beta_pi[swir_composition] <- 1
  
  # Calcular expectativas para reporte
  expected_inclusion <- alpha_pi / (alpha_pi + beta_pi)
  
  cat("=== CONFIGURACIÓN DE PRIORS PARA ETAPA TARDÍA (11-15 DPI) ===\n")
  cat("Bandas con alta expectativa (>0.7):", sum(expected_inclusion > 0.7), "\n")
  cat("Bandas moderadas (0.5-0.7):", sum(expected_inclusion >= 0.5 & expected_inclusion <= 0.7), "\n")
  cat("Bandas con baja expectativa (<0.5):", sum(expected_inclusion < 0.5), "\n")
  cat("Expectativa media de inclusión:", round(mean(expected_inclusion), 3), "\n")
  
  # Detalles por región
  cat("\nDETALLE POR REGIÓN:\n")
  cat("Clorosis/muerte (560-733 nm):", length(chlorophyll_death), "bandas, E[π] =", 
      round(mean(expected_inclusion[chlorophyll_death]), 2), "\n")
  cat("Fotosíntesis colapso (400-700 nm):", length(photosynthesis_collapse), "bandas, E[π] =", 
      round(mean(expected_inclusion[photosynthesis_collapse]), 2), "\n")
  cat("Sistema vascular (800-1200 nm):", length(vascular_damage), "bandas, E[π] =", 
      round(mean(expected_inclusion[vascular_damage]), 2), "\n")
  cat("Agua colapso (multi-región):", length(water_collapse), "bandas, E[π] =", 
      round(mean(expected_inclusion[water_collapse]), 2), "\n")
  cat("Necrosis (750-850 nm):", length(necrosis), "bandas, E[π] =", 
      round(mean(expected_inclusion[necrosis]), 2), "\n")
  
  return(list(
    alpha_pi = alpha_pi, 
    beta_pi = beta_pi,
    expected_inclusion = expected_inclusion,
    regions = list(
      chlorophyll_death = chlorophyll_death,
      photosynthesis_collapse = photosynthesis_collapse,
      vascular_damage = vascular_damage,
      water_collapse = water_collapse,
      necrosis = necrosis,
      oxidative_stress = oxidative_stress,
      toxin_response = toxin_response,
      structural_collapse = structural_collapse,
      swir_composition = swir_composition
    )
  ))
}

prior_config <- setup_late_stage_priors(wavelengths, P)

#----------------------------------------------------------------------------

# Estandarizar datos espectrales 
X_scaled <- scale(X)  # Media 0, SD 1 por columna

# Configurar hiperparámetros
sigma_beta_value <- 0.1    # Regularización moderada
sigma_alpha_value <- 2.5   # Prior débil para intercepto

# Lista para Stan
stan_data <- list(
  N = N,
  P = P,
  X = X_scaled,
  y = y,
  alpha_pi = prior_config$alpha_pi,
  beta_pi = prior_config$beta_pi,
  sigma_beta = sigma_beta_value,
  sigma_alpha = sigma_alpha_value
)

# Verificar datos
str(stan_data)

# Compilar modelo
model <- stan_model("ModeloRLB_v2.stan")

# Ejecutar (debería tomar 30-60 minutos vs. 6 horas)
fit <- sampling(
  model,
  data = stan_data,
  chains = 4,
  iter = 2000,
  warmup = 1000,
  cores = 4,
  control = list(adapt_delta = 0.95, max_treedepth = 12)
)

# Verificar convergencia
print(fit, pars = c("B_0", "lp__"))

#------------------------------------------------------------------------------ día 1 
fit_d1 = fit

B_0_samplesd1 <- rstan::extract(fit, "B_0")[[1]]
B_j_samplesd1 <- rstan::extract(fit, "B_j")[[1]]
pi_j_samplesd1 <- rstan::extract(fit, "pi_j")[[1]]

# === CALCULAR MÉTRICAS EN R  ===

# 1. Estadísticas posteriores
B_0_meand1 <- mean(B_0_samplesd1)
B_j_meand1 <- colMeans(B_j_samplesd1)
inclusion_probd1 <- colMeans(pi_j_samplesd1)

# 2. Coeficientes efectivos
effective_betad1 <- B_j_meand1 * inclusion_probd1

# 3. Métricas de importancia
importanced1 <- abs(effective_betad1)
selected_bandsd1 <- inclusion_probd1 > 0.5
num_selectedd1 <- sum(selected_bandsd1)

# 4. Predicciones
logit_predd1 <- B_0_meand1 + X_scaled %*% effective_betad1  # Usar X_scaled
y_pred_probd1 <- plogis(logit_predd1)  # inv_logit en R
y_predd1 <- as.numeric(y_pred_probd1 > 0.5)

# 5. Métricas de rendimiento
accuracyd1 <- mean(y_predd1 == y)
sensitivityd1 <- sum(y_predd1 == 1 & y == 1) / sum(y == 1)
specificityd1 <- sum(y_predd1 == 0 & y == 0) / sum(y == 0)

# === REPORTE (ACTUALIZADO) ===
cat("=== MÉTRICAS DEL MODELO - DÍA", dia_actual, "===\n")
cat("Precisión:", round(accuracyd1 * 100, 1), "%\n")
cat("Sensibilidad:", round(sensitivityd1 * 100, 1), "%\n")
cat("Especificidad:", round(specificityd1 * 100, 1), "%\n")
cat("Longitudes importantes (>50%):", num_selectedd1, "\n")
cat("Rango:", min(wavelengths[selected_bandsd1]), "-", max(wavelengths[selected_bandsd1]), "nm\n")


# Top 10 bandas más importantes (ACTUALIZADO)
top_bandsd1 <- order(importanced1, decreasing = TRUE)[1:10]
cat("\nTop 10 bandas:\n")
for(i in 1:10) {
  idx <- top_bandsd1[i]
  cat(sprintf("%.0f nm: π=%.3f, |β_eff|=%.4f\n", 
              wavelengths[idx], inclusion_prob[idx], importance[idx]))
}

# === ANÁLISIS ADICIONAL ===
sparsity_ratiod1 <- mean(inclusion_probd1)
cat("\nSparsity ratio:", round(sparsity_ratiod1, 3), "\n")

# Incertidumbre en selección
inclusion_sdd1 <- apply(pi_j_samplesd1, 2, sd)
high_uncertaintyd1 <- sum(inclusion_sdd1 > 0.3)
cat("Bandas con alta incertidumbre (SD > 0.3):", high_uncertaintyd1, "\n")

# Contribución por región espectral
visible_contribd1 <- sum(importanced1[wavelengths <= 700])
nir_contribd1 <- sum(importanced1[wavelengths > 700 & wavelengths <= 1300])
swir_contribd1 <- sum(importanced1[wavelengths > 1300])
total_contribd1 <- sum(importanced1)

cat("\nContribución por región:\n")
cat("VIS (≤700 nm):", round(visible_contribd1/total_contribd1*100, 1), "%\n")
cat("NIR (700-1300 nm):", round(nir_contribd1/total_contribd1*100, 1), "%\n")
cat("SWIR (>1300 nm):", round(swir_contribd1/total_contribd1*100, 1), "%\n")

# Consistencia con literatura
if(exists("prior_config")) {
  prior_expectationd1 <- prior_config$expected_inclusion
  literature_agreementd1 <- cor(importanced1, prior_expectation)
  cat("Correlación con priors:", round(literature_agreementd1, 3), "\n")
}




#------------------------------------------------------------------------------ día 2

fit_d2 = fit

B_0_samplesd2 <- rstan::extract(fit, "B_0")[[1]]
B_j_samplesd2 <- rstan::extract(fit, "B_j")[[1]]
pi_j_samplesd2 <- rstan::extract(fit, "pi_j")[[1]]

# === CALCULAR MÉTRICAS EN R ===

# 1. Estadísticas posteriores
B_0_meand2 <- mean(B_0_samplesd2)
B_j_meand2 <- colMeans(B_j_samplesd2)
inclusion_probd2 <- colMeans(pi_j_samplesd2)

# 2. Coeficientes efectivos
effective_betad2 <- B_j_meand2 * inclusion_probd2

# 3. Métricas de importancia
importanced2 <- abs(effective_betad2)
selected_bandsd2 <- inclusion_probd2 > 0.5
num_selectedd2 <- sum(selected_bandsd2)

# 4. Predicciones
logit_predd2 <- B_0_meand2 + X_scaled %*% effective_betad2  # Usar X_scaled
y_pred_probd2 <- plogis(logit_predd2)  # inv_logit en R
y_predd2 <- as.numeric(y_pred_probd2 > 0.5)

# 5. Métricas de rendimiento
accuracyd2 <- mean(y_predd2 == y)
sensitivityd2 <- sum(y_predd2 == 1 & y == 1) / sum(y == 1)
specificityd2 <- sum(y_predd2 == 0 & y == 0) / sum(y == 0)

# === REPORTE (ACTUALIZADO) ===
cat("=== MÉTRICAS DEL MODELO - DÍA", dia_actual, "===\n")
cat("Precisión:", round(accuracyd2 * 100, 1), "%\n")
cat("Sensibilidad:", round(sensitivityd2 * 100, 1), "%\n")
cat("Especificidad:", round(specificityd2 * 100, 1), "%\n")
cat("Longitudes importantes (>50%):", num_selectedd2, "\n")
cat("Rango:", min(wavelengths[selected_bandsd2]), "-", max(wavelengths[selected_bandsd2]), "nm\n")


# Top 10 bandas más importantes (ACTUALIZADO)
top_bandsd2 <- order(importanced2, decreasing = TRUE)[1:10]
cat("\nTop 10 bandas:\n")
for(i in 1:10) {
  idx <- top_bandsd2[i]
  cat(sprintf("%.0f nm: π=%.3f, |β_eff|=%.4f\n", 
              wavelengths[idx], inclusion_probd2[idx], importanced2[idx]))
}

# === ANÁLISIS ADICIONAL ===
sparsity_ratiod2 <- mean(inclusion_probd2)
cat("\nSparsity ratio:", round(sparsity_ratiod2, 3), "\n")

# Incertidumbre en selección
inclusion_sdd2 <- apply(pi_j_samplesd2, 2, sd)
high_uncertaintyd2 <- sum(inclusion_sdd2 > 0.3)
cat("Bandas con alta incertidumbre (SD > 0.3):", high_uncertaintyd2, "\n")

# Contribución por región espectral
visible_contribd2 <- sum(importanced2[wavelengths <= 700])
nir_contribd2 <- sum(importanced2[wavelengths > 700 & wavelengths <= 1300])
swir_contribd2 <- sum(importanced2[wavelengths > 1300])
total_contribd2 <- sum(importanced2)

cat("\nContribución por región:\n")
cat("VIS (≤700 nm):", round(visible_contribd2/total_contribd2*100, 1), "%\n")
cat("NIR (700-1300 nm):", round(nir_contribd2/total_contribd2*100, 1), "%\n")
cat("SWIR (>1300 nm):", round(swir_contribd2/total_contribd2*100, 1), "%\n")

# Consistencia con literatura
if(exists("prior_config")) {
  prior_expectationd2 <- prior_config$expected_inclusion
  literature_agreementd2 <- cor(importanced2, prior_expectationd2)
  cat("Correlación con priors:", round(literature_agreementd2, 3), "\n")
}

# ------------------------------------------------------------------------------día 3

fit_d3 = fit

B_0_samplesd3 <- rstan::extract(fit, "B_0")[[1]]
B_j_samplesd3 <- rstan::extract(fit, "B_j")[[1]]
pi_j_samplesd3 <- rstan::extract(fit, "pi_j")[[1]]

# === CALCULAR MÉTRICAS EN R  ===

# 1. Estadísticas posteriores
B_0_meand3 <- mean(B_0_samplesd3)
B_j_meand3 <- colMeans(B_j_samplesd3)
inclusion_probd3 <- colMeans(pi_j_samplesd3)

# 2. Coeficientes efectivos
effective_betad3 <- B_j_meand3 * inclusion_probd3

# 3. Métricas de importancia
importanced3 <- abs(effective_betad3)
selected_bandsd3 <- inclusion_probd3 > 0.5
num_selectedd3 <- sum(selected_bandsd3)

# 4. Predicciones
logit_predd3 <- B_0_meand3 + X_scaled %*% effective_betad3  # Usar X_scaled
y_pred_probd3 <- plogis(logit_predd3)  # inv_logit en R
y_predd3 <- as.numeric(y_pred_probd3 > 0.5)

# 5. Métricas de rendimiento
accuracyd3 <- mean(y_predd3 == y)
sensitivityd3 <- sum(y_predd3 == 1 & y == 1) / sum(y == 1)
specificityd3 <- sum(y_predd3 == 0 & y == 0) / sum(y == 0)

# === REPORTE (ACTUALIZADO) ===
cat("=== MÉTRICAS DEL MODELO - DÍA", dia_actual, "===\n")
cat("Precisión:", round(accuracyd3 * 100, 1), "%\n")
cat("Sensibilidad:", round(sensitivityd3 * 100, 1), "%\n")
cat("Especificidad:", round(specificityd3 * 100, 1), "%\n")
cat("Longitudes importantes (>50%):", num_selectedd3, "\n")
cat("Rango:", min(wavelengths[selected_bandsd3]), "-", max(wavelengths[selected_bandsd3]), "nm\n")


# Top 10 bandas más importantes (ACTUALIZADO)
top_bandsd3 <- order(importanced3, decreasing = TRUE)[1:10]
cat("\nTop 10 bandas:\n")
for(i in 1:10) {
  idx <- top_bandsd3[i]
  cat(sprintf("%.0f nm: π=%.3f, |β_eff|=%.4f\n", 
              wavelengths[idx], inclusion_probd3[idx], importanced3[idx]))
}

# === ANÁLISIS ADICIONAL ===
sparsity_ratiod3 <- mean(inclusion_probd3)
cat("\nSparsity ratio:", round(sparsity_ratiod3, 3), "\n")

# Incertidumbre en selección
inclusion_sdd3 <- apply(pi_j_samplesd3, 2, sd)
high_uncertaintyd3 <- sum(inclusion_sdd3 > 0.3)
cat("Bandas con alta incertidumbre (SD > 0.3):", high_uncertaintyd3, "\n")

# Contribución por región espectral
visible_contribd3 <- sum(importanced3[wavelengths <= 700])
nir_contribd3 <- sum(importanced3[wavelengths > 700 & wavelengths <= 1300])
swir_contribd3 <- sum(importanced3[wavelengths > 1300])
total_contribd3 <- sum(importanced3)

cat("\nContribución por región:\n")
cat("VIS (≤700 nm):", round(visible_contribd3/total_contribd3*100, 1), "%\n")
cat("NIR (700-1300 nm):", round(nir_contribd3/total_contribd3*100, 1), "%\n")
cat("SWIR (>1300 nm):", round(swir_contribd3/total_contribd3*100, 1), "%\n")

# Consistencia con literatura
if(exists("prior_config")) {
  prior_expectationd3 <- prior_config$expected_inclusion
  literature_agreementd3 <- cor(importanced3, prior_expectationd3)
  cat("Correlación con priors:", round(literature_agreementd3, 3), "\n")
}

# ------------------------------------------------------------------------------día 4

fit_d4 = fit

B_0_samplesd4 <- rstan::extract(fit, "B_0")[[1]]
B_j_samplesd4 <- rstan::extract(fit, "B_j")[[1]]
pi_j_samplesd4 <- rstan::extract(fit, "pi_j")[[1]]

# === CALCULAR MÉTRICAS EN R  ===

# 1. Estadísticas posteriores
B_0_meand4 <- mean(B_0_samplesd4)
B_j_meand4 <- colMeans(B_j_samplesd4)
inclusion_probd4 <- colMeans(pi_j_samplesd4)

# 2. Coeficientes efectivos
effective_betad4 <- B_j_meand4 * inclusion_probd4

# 3. Métricas de importancia
importanced4 <- abs(effective_betad4)
selected_bandsd4 <- inclusion_probd4 > 0.5
num_selectedd4 <- sum(selected_bandsd4)

# 4. Predicciones
logit_predd4 <- B_0_meand4 + X_scaled %*% effective_betad4  # Usar X_scaled
y_pred_probd4 <- plogis(logit_predd4)  # inv_logit en R
y_predd4 <- as.numeric(y_pred_probd4 > 0.5)

# 5. Métricas de rendimiento
accuracyd4 <- mean(y_predd4 == y)
sensitivityd4 <- sum(y_predd4 == 1 & y == 1) / sum(y == 1)
specificityd4 <- sum(y_predd4 == 0 & y == 0) / sum(y == 0)

# === REPORTE (ACTUALIZADO) ===
cat("=== MÉTRICAS DEL MODELO - DÍA", dia_actual, "===\n")
cat("Precisión:", round(accuracyd4 * 100, 1), "%\n")
cat("Sensibilidad:", round(sensitivityd4 * 100, 1), "%\n")
cat("Especificidad:", round(specificityd4 * 100, 1), "%\n")
cat("Longitudes importantes (>50%):", num_selectedd4, "\n")
cat("Rango:", min(wavelengths[selected_bandsd4]), "-", max(wavelengths[selected_bandsd4]), "nm\n")


# Top 10 bandas más importantes (ACTUALIZADO)
top_bandsd4 <- order(importanced4, decreasing = TRUE)[1:10]
cat("\nTop 10 bandas:\n")
for(i in 1:10) {
  idx <- top_bandsd4[i]
  cat(sprintf("%.0f nm: π=%.3f, |β_eff|=%.4f\n", 
              wavelengths[idx], inclusion_probd4[idx], importanced4[idx]))
}

# === ANÁLISIS ADICIONAL ===
sparsity_ratiod4 <- mean(inclusion_probd4)
cat("\nSparsity ratio:", round(sparsity_ratiod4, 3), "\n")

# Incertidumbre en selección
inclusion_sdd4 <- apply(pi_j_samplesd4, 2, sd)
high_uncertaintyd4 <- sum(inclusion_sdd4 > 0.3)
cat("Bandas con alta incertidumbre (SD > 0.3):", high_uncertaintyd4, "\n")

# Contribución por región espectral
visible_contribd4 <- sum(importanced4[wavelengths <= 700])
nir_contribd4 <- sum(importanced4[wavelengths > 700 & wavelengths <= 1300])
swir_contribd4 <- sum(importanced4[wavelengths > 1300])
total_contribd4 <- sum(importanced4)

cat("\nContribución por región:\n")
cat("VIS (≤700 nm):", round(visible_contribd4/total_contribd4*100, 1), "%\n")
cat("NIR (700-1300 nm):", round(nir_contribd4/total_contribd4*100, 1), "%\n")
cat("SWIR (>1300 nm):", round(swir_contribd4/total_contribd4*100, 1), "%\n")

# Consistencia con literatura
if(exists("prior_config")) {
  prior_expectationd4 <- prior_config$expected_inclusion
  literature_agreementd4 <- cor(importanced4, prior_expectationd4)
  cat("Correlación con priors:", round(literature_agreementd4, 3), "\n")
}


# ------------------------------------------------------------------------------día 5

fit_d5 = fit

B_0_samplesd5 <- rstan::extract(fit, "B_0")[[1]]
B_j_samplesd5 <- rstan::extract(fit, "B_j")[[1]]
pi_j_samplesd5 <- rstan::extract(fit, "pi_j")[[1]]

# === CALCULAR MÉTRICAS EN R  ===

# 1. Estadísticas posteriores
B_0_meand5 <- mean(B_0_samplesd5)
B_j_meand5 <- colMeans(B_j_samplesd5)
inclusion_probd5 <- colMeans(pi_j_samplesd5)

# 2. Coeficientes efectivos
effective_betad5 <- B_j_meand5 * inclusion_probd5

# 3. Métricas de importancia
importanced5 <- abs(effective_betad5)
selected_bandsd5 <- inclusion_probd5 > 0.5
num_selectedd5 <- sum(selected_bandsd5)

# 4. Predicciones
logit_predd5 <- B_0_meand5 + X_scaled %*% effective_betad5  # Usar X_scaled
y_pred_probd5 <- plogis(logit_predd5)  # inv_logit en R
y_predd5 <- as.numeric(y_pred_probd5 > 0.5)

# 5. Métricas de rendimiento
accuracyd5 <- mean(y_predd5 == y)
sensitivityd5 <- sum(y_predd5 == 1 & y == 1) / sum(y == 1)
specificityd5 <- sum(y_predd5 == 0 & y == 0) / sum(y == 0)

# === REPORTE ===
cat("=== MÉTRICAS DEL MODELO - DÍA", dia_actual, "===\n")
cat("Precisión:", round(accuracyd5 * 100, 1), "%\n")
cat("Sensibilidad:", round(sensitivityd5 * 100, 1), "%\n")
cat("Especificidad:", round(specificityd5 * 100, 1), "%\n")
cat("Longitudes importantes (>50%):", num_selectedd5, "\n")
cat("Rango:", min(wavelengths[selected_bandsd5]), "-", max(wavelengths[selected_bandsd5]), "nm\n")


# Top 10 bandas más importantes (ACTUALIZADO)
top_bandsd5 <- order(importanced5, decreasing = TRUE)[1:10]
cat("\nTop 10 bandas:\n")
for(i in 1:10) {
  idx <- top_bandsd5[i]
  cat(sprintf("%.0f nm: π=%.3f, |β_eff|=%.4f\n", 
              wavelengths[idx], inclusion_probd5[idx], importanced5[idx]))
}

# === ANÁLISIS ADICIONAL ===
sparsity_ratiod5 <- mean(inclusion_probd5)
cat("\nSparsity ratio:", round(sparsity_ratiod5, 3), "\n")

# Incertidumbre en selección
inclusion_sdd5 <- apply(pi_j_samplesd5, 2, sd)
high_uncertaintyd5 <- sum(inclusion_sdd5 > 0.3)
cat("Bandas con alta incertidumbre (SD > 0.3):", high_uncertaintyd5, "\n")

# Contribución por región espectral
visible_contribd5 <- sum(importanced5[wavelengths <= 700])
nir_contribd5 <- sum(importanced5[wavelengths > 700 & wavelengths <= 1300])
swir_contribd5 <- sum(importanced5[wavelengths > 1300])
total_contribd5 <- sum(importanced5)

cat("\nContribución por región:\n")
cat("VIS (≤700 nm):", round(visible_contribd5/total_contribd5*100, 1), "%\n")
cat("NIR (700-1300 nm):", round(nir_contribd5/total_contribd5*100, 1), "%\n")
cat("SWIR (>1300 nm):", round(swir_contribd5/total_contribd5*100, 1), "%\n")

# Consistencia con literatura
if(exists("prior_config")) {
  prior_expectationd5 <- prior_config$expected_inclusion
  literature_agreementd5 <- cor(importanced5, prior_expectationd5)
  cat("Correlación con priors:", round(literature_agreementd5, 3), "\n")
}

# ------------------------------------------------------------------------------día 6

fit_d6 = fit

B_0_samplesd6 <- rstan::extract(fit, "B_0")[[1]]
B_j_samplesd6 <- rstan::extract(fit, "B_j")[[1]]
pi_j_samplesd6 <- rstan::extract(fit, "pi_j")[[1]]

# === CALCULAR MÉTRICAS EN R  ===

# 1. Estadísticas posteriores
B_0_meand6 <- mean(B_0_samplesd6)
B_j_meand6 <- colMeans(B_j_samplesd6)
inclusion_probd6 <- colMeans(pi_j_samplesd6)

# 2. Coeficientes efectivos
effective_betad6 <- B_j_meand6 * inclusion_probd6

# 3. Métricas de importancia
importanced6 <- abs(effective_betad6)
selected_bandsd6 <- inclusion_probd6 > 0.5
num_selectedd6 <- sum(selected_bandsd6)

# 4. Predicciones
logit_predd6 <- B_0_meand6 + X_scaled %*% effective_betad6  # Usar X_scaled
y_pred_probd6 <- plogis(logit_predd6)  # inv_logit en R
y_predd6 <- as.numeric(y_pred_probd6 > 0.5)

# 5. Métricas de rendimiento
accuracyd6 <- mean(y_predd6 == y)
sensitivityd6 <- sum(y_predd6 == 1 & y == 1) / sum(y == 1)
specificityd6 <- sum(y_predd6 == 0 & y == 0) / sum(y == 0)

# === REPORTE (ACTUALIZADO) ===
cat("=== MÉTRICAS DEL MODELO - DÍA", dia_actual, "===\n")
cat("Precisión:", round(accuracyd6 * 100, 1), "%\n")
cat("Sensibilidad:", round(sensitivityd6 * 100, 1), "%\n")
cat("Especificidad:", round(specificityd6 * 100, 1), "%\n")
cat("Longitudes importantes (>50%):", num_selectedd6, "\n")
cat("Rango:", min(wavelengths[selected_bandsd6]), "-", max(wavelengths[selected_bandsd6]), "nm\n")


# Top 10 bandas más importantes (ACTUALIZADO)
top_bandsd6 <- order(importanced6, decreasing = TRUE)[1:10]
cat("\nTop 10 bandas:\n")
for(i in 1:10) {
  idx <- top_bandsd6[i]
  cat(sprintf("%.0f nm: π=%.3f, |β_eff|=%.4f\n", 
              wavelengths[idx], inclusion_probd6[idx], importanced6[idx]))
}

# === ANÁLISIS ADICIONAL ===
sparsity_ratiod6 <- mean(inclusion_probd6)
cat("\nSparsity ratio:", round(sparsity_ratiod6, 3), "\n")

# Incertidumbre en selección
inclusion_sdd6 <- apply(pi_j_samplesd6, 2, sd)
high_uncertaintyd6 <- sum(inclusion_sdd6 > 0.3)
cat("Bandas con alta incertidumbre (SD > 0.3):", high_uncertaintyd6, "\n")

# Contribución por región espectral
visible_contribd6 <- sum(importanced6[wavelengths <= 700])
nir_contribd6 <- sum(importanced6[wavelengths > 700 & wavelengths <= 1300])
swir_contribd6 <- sum(importanced6[wavelengths > 1300])
total_contribd6 <- sum(importanced6)

cat("\nContribución por región:\n")
cat("VIS (≤700 nm):", round(visible_contribd6/total_contribd6*100, 1), "%\n")
cat("NIR (700-1300 nm):", round(nir_contribd6/total_contribd6*100, 1), "%\n")
cat("SWIR (>1300 nm):", round(swir_contribd6/total_contribd6*100, 1), "%\n")

# Consistencia con literatura
if(exists("prior_config")) {
  prior_expectationd6 <- prior_config$expected_inclusion
  literature_agreementd6 <- cor(importanced6, prior_expectationd6)
  cat("Correlación con priors:", round(literature_agreementd6, 3), "\n")
}




# ------------------------------------------------------------------------------día 7

fit_d7 = fit

B_0_samplesd7 <- rstan::extract(fit, "B_0")[[1]]
B_j_samplesd7 <- rstan::extract(fit, "B_j")[[1]]
pi_j_samplesd7 <- rstan::extract(fit, "pi_j")[[1]]

# === CALCULAR MÉTRICAS EN R  ===

# 1. Estadísticas posteriores
B_0_meand7 <- mean(B_0_samplesd7)
B_j_meand7 <- colMeans(B_j_samplesd7)
inclusion_probd7 <- colMeans(pi_j_samplesd7)

# 2. Coeficientes efectivos
effective_betad7 <- B_j_meand7 * inclusion_probd7

# 3. Métricas de importancia
importanced7 <- abs(effective_betad7)
selected_bandsd7 <- inclusion_probd7 > 0.5
num_selectedd7 <- sum(selected_bandsd7)

# 4. Predicciones
logit_predd7 <- B_0_meand7 + X_scaled %*% effective_betad7  # Usar X_scaled
y_pred_probd7 <- plogis(logit_predd7)  # inv_logit en R
y_predd7 <- as.numeric(y_pred_probd7 > 0.5)

# 5. Métricas de rendimiento
accuracyd7 <- mean(y_predd7 == y)
sensitivityd7 <- sum(y_predd7 == 1 & y == 1) / sum(y == 1)
specificityd7 <- sum(y_predd7 == 0 & y == 0) / sum(y == 0)

# === REPORTE (ACTUALIZADO) ===
cat("=== MÉTRICAS DEL MODELO - DÍA", dia_actual, "===\n")
cat("Precisión:", round(accuracyd7 * 100, 1), "%\n")
cat("Sensibilidad:", round(sensitivityd7 * 100, 1), "%\n")
cat("Especificidad:", round(specificityd7 * 100, 1), "%\n")
cat("Longitudes importantes (>50%):", num_selectedd7, "\n")
cat("Rango:", min(wavelengths[selected_bandsd7]), "-", max(wavelengths[selected_bandsd7]), "nm\n")


# Top 10 bandas más importantes (ACTUALIZADO)
top_bandsd7 <- order(importanced7, decreasing = TRUE)[1:10]
cat("\nTop 10 bandas:\n")
for(i in 1:10) {
  idx <- top_bandsd7[i]
  cat(sprintf("%.0f nm: π=%.3f, |β_eff|=%.4f\n", 
              wavelengths[idx], inclusion_probd7[idx], importanced7[idx]))
}

# === ANÁLISIS ADICIONAL ===
sparsity_ratiod7 <- mean(inclusion_probd7)
cat("\nSparsity ratio:", round(sparsity_ratiod7, 3), "\n")

# Incertidumbre en selección
inclusion_sdd7 <- apply(pi_j_samplesd7, 2, sd)
high_uncertaintyd7 <- sum(inclusion_sdd7 > 0.3)
cat("Bandas con alta incertidumbre (SD > 0.3):", high_uncertaintyd7, "\n")

# Contribución por región espectral
visible_contribd7 <- sum(importanced7[wavelengths <= 700])
nir_contribd7 <- sum(importanced7[wavelengths > 700 & wavelengths <= 1300])
swir_contribd7 <- sum(importanced7[wavelengths > 1300])
total_contribd7 <- sum(importanced7)

cat("\nContribución por región:\n")
cat("VIS (≤700 nm):", round(visible_contribd7/total_contribd7*100, 1), "%\n")
cat("NIR (700-1300 nm):", round(nir_contribd7/total_contribd7*100, 1), "%\n")
cat("SWIR (>1300 nm):", round(swir_contribd7/total_contribd7*100, 1), "%\n")

# Consistencia con literatura
if(exists("prior_config")) {
  prior_expectationd7 <- prior_config$expected_inclusion
  literature_agreementd7 <- cor(importanced7, prior_expectationd7)
  cat("Correlación con priors:", round(literature_agreementd7, 3), "\n")
}


# ------------------------------------------------------------------------------día 8

fit_d8 = fit

B_0_samplesd8 <- rstan::extract(fit, "B_0")[[1]]
B_j_samplesd8 <- rstan::extract(fit, "B_j")[[1]]
pi_j_samplesd8 <- rstan::extract(fit, "pi_j")[[1]]

# === CALCULAR MÉTRICAS EN R  ===

# 1. Estadísticas posteriores
B_0_meand8 <- mean(B_0_samplesd8)
B_j_meand8 <- colMeans(B_j_samplesd8)
inclusion_probd8 <- colMeans(pi_j_samplesd8)

# 2. Coeficientes efectivos
effective_betad8 <- B_j_meand8 * inclusion_probd8

# 3. Métricas de importancia
importanced8 <- abs(effective_betad8)
selected_bandsd8 <- inclusion_probd8 > 0.5
num_selectedd8 <- sum(selected_bandsd8)

# 4. Predicciones
logit_predd8 <- B_0_meand8 + X_scaled %*% effective_betad8  # Usar X_scaled
y_pred_probd8 <- plogis(logit_predd8)  # inv_logit en R
y_predd8 <- as.numeric(y_pred_probd8 > 0.5)

# 5. Métricas de rendimiento
accuracyd8 <- mean(y_predd8 == y)
sensitivityd8 <- sum(y_predd8 == 1 & y == 1) / sum(y == 1)
specificityd8 <- sum(y_predd8 == 0 & y == 0) / sum(y == 0)

# === REPORTE (ACTUALIZADO) ===
cat("=== MÉTRICAS DEL MODELO - DÍA", dia_actual, "===\n")
cat("Precisión:", round(accuracyd8 * 100, 1), "%\n")
cat("Sensibilidad:", round(sensitivityd8 * 100, 1), "%\n")
cat("Especificidad:", round(specificityd8 * 100, 1), "%\n")
cat("Longitudes importantes (>50%):", num_selectedd8, "\n")
cat("Rango:", min(wavelengths[selected_bandsd8]), "-", max(wavelengths[selected_bandsd8]), "nm\n")


# Top 10 bandas más importantes (ACTUALIZADO)
top_bandsd8 <- order(importanced8, decreasing = TRUE)[1:10]
cat("\nTop 10 bandas:\n")
for(i in 1:10) {
  idx <- top_bandsd8[i]
  cat(sprintf("%.0f nm: π=%.3f, |β_eff|=%.4f\n", 
              wavelengths[idx], inclusion_probd8[idx], importanced8[idx]))
}

# === ANÁLISIS ADICIONAL ===
sparsity_ratiod8 <- mean(inclusion_probd8)
cat("\nSparsity ratio:", round(sparsity_ratiod8, 3), "\n")

# Incertidumbre en selección
inclusion_sdd8 <- apply(pi_j_samplesd8, 2, sd)
high_uncertaintyd8 <- sum(inclusion_sdd8 > 0.3)
cat("Bandas con alta incertidumbre (SD > 0.3):", high_uncertaintyd8, "\n")

# Contribución por región espectral
visible_contribd8 <- sum(importanced8[wavelengths <= 700])
nir_contribd8 <- sum(importanced8[wavelengths > 700 & wavelengths <= 1300])
swir_contribd8 <- sum(importanced8[wavelengths > 1300])
total_contribd8 <- sum(importanced8)

cat("\nContribución por región:\n")
cat("VIS (≤700 nm):", round(visible_contribd8/total_contribd8*100, 1), "%\n")
cat("NIR (700-1300 nm):", round(nir_contribd8/total_contribd8*100, 1), "%\n")
cat("SWIR (>1300 nm):", round(swir_contribd8/total_contribd8*100, 1), "%\n")

# Consistencia con literatura
if(exists("prior_config")) {
  prior_expectationd8 <- prior_config$expected_inclusion
  literature_agreementd8 <- cor(importanced8, prior_expectationd8)
  cat("Correlación con priors:", round(literature_agreementd8, 3), "\n")
}

# ------------------------------------------------------------------------------día 9

fit_d9 = fit

B_0_samplesd9 <- rstan::extract(fit, "B_0")[[1]]
B_j_samplesd9 <- rstan::extract(fit, "B_j")[[1]]
pi_j_samplesd9 <- rstan::extract(fit, "pi_j")[[1]]

# === CALCULAR MÉTRICAS EN R  ===

# 1. Estadísticas posteriores
B_0_meand9 <- mean(B_0_samplesd9)
B_j_meand9 <- colMeans(B_j_samplesd9)
inclusion_probd9 <- colMeans(pi_j_samplesd9)

# 2. Coeficientes efectivos
effective_betad9 <- B_j_meand9 * inclusion_probd9

# 3. Métricas de importancia
importanced9 <- abs(effective_betad9)
selected_bandsd9 <- inclusion_probd9 > 0.5
num_selectedd9 <- sum(selected_bandsd9)

# 4. Predicciones
logit_predd9 <- B_0_meand9 + X_scaled %*% effective_betad9  # Usar X_scaled
y_pred_probd9 <- plogis(logit_predd9)  # inv_logit en R
y_predd9 <- as.numeric(y_pred_probd9 > 0.5)

# 5. Métricas de rendimiento
accuracyd9 <- mean(y_predd9 == y)
sensitivityd9 <- sum(y_predd9 == 1 & y == 1) / sum(y == 1)
specificityd9 <- sum(y_predd9 == 0 & y == 0) / sum(y == 0)

# === REPORTE (ACTUALIZADO) ===
cat("=== MÉTRICAS DEL MODELO - DÍA", dia_actual, "===\n")
cat("Precisión:", round(accuracyd9 * 100, 1), "%\n")
cat("Sensibilidad:", round(sensitivityd9 * 100, 1), "%\n")
cat("Especificidad:", round(specificityd9 * 100, 1), "%\n")
cat("Longitudes importantes (>50%):", num_selectedd9, "\n")
cat("Rango:", min(wavelengths[selected_bandsd9]), "-", max(wavelengths[selected_bandsd9]), "nm\n")


# Top 10 bandas más importantes (ACTUALIZADO)
top_bandsd9 <- order(importanced9, decreasing = TRUE)[1:10]
cat("\nTop 10 bandas:\n")
for(i in 1:10) {
  idx <- top_bandsd9[i]
  cat(sprintf("%.0f nm: π=%.3f, |β_eff|=%.4f\n", 
              wavelengths[idx], inclusion_probd9[idx], importanced9[idx]))
}

# === ANÁLISIS ADICIONAL ===
sparsity_ratiod9 <- mean(inclusion_probd9)
cat("\nSparsity ratio:", round(sparsity_ratiod9, 3), "\n")

# Incertidumbre en selección
inclusion_sdd9 <- apply(pi_j_samplesd9, 2, sd)
high_uncertaintyd9 <- sum(inclusion_sdd9 > 0.3)
cat("Bandas con alta incertidumbre (SD > 0.3):", high_uncertaintyd9, "\n")

# Contribución por región espectral
visible_contribd9 <- sum(importanced9[wavelengths <= 700])
nir_contribd9 <- sum(importanced9[wavelengths > 700 & wavelengths <= 1300])
swir_contribd9 <- sum(importanced9[wavelengths > 1300])
total_contribd9 <- sum(importanced9)

cat("\nContribución por región:\n")
cat("VIS (≤700 nm):", round(visible_contribd9/total_contribd9*100, 1), "%\n")
cat("NIR (700-1300 nm):", round(nir_contribd9/total_contribd9*100, 1), "%\n")
cat("SWIR (>1300 nm):", round(swir_contribd9/total_contribd9*100, 1), "%\n")

# Consistencia con literatura
if(exists("prior_config")) {
  prior_expectationd9 <- prior_config$expected_inclusion
  literature_agreementd9 <- cor(importanced9, prior_expectationd9)
  cat("Correlación con priors:", round(literature_agreementd9, 3), "\n")
}

# ------------------------------------------------------------------------------día 10

fit_d10 = fit

B_0_samplesd10 <- rstan::extract(fit, "B_0")[[1]]
B_j_samplesd10 <- rstan::extract(fit, "B_j")[[1]]
pi_j_samplesd10 <- rstan::extract(fit, "pi_j")[[1]]

# === CALCULAR MÉTRICAS EN R  ===

# 1. Estadísticas posteriores
B_0_meand10 <- mean(B_0_samplesd10)
B_j_meand10 <- colMeans(B_j_samplesd10)
inclusion_probd10 <- colMeans(pi_j_samplesd10)

# 2. Coeficientes efectivos
effective_betad10 <- B_j_meand10 * inclusion_probd10

# 3. Métricas de importancia
importanced10 <- abs(effective_betad10)
selected_bandsd10 <- inclusion_probd10 > 0.5
num_selectedd10 <- sum(selected_bandsd10)

# 4. Predicciones
logit_predd10 <- B_0_meand10 + X_scaled %*% effective_betad10  # Usar X_scaled
y_pred_probd10 <- plogis(logit_predd10)  # inv_logit en R
y_predd10 <- as.numeric(y_pred_probd10 > 0.5)

# 5. Métricas de rendimiento
accuracyd10 <- mean(y_predd10 == y)
sensitivityd10 <- sum(y_predd10 == 1 & y == 1) / sum(y == 1)
specificityd10 <- sum(y_predd10 == 0 & y == 0) / sum(y == 0)

# === REPORTE (ACTUALIZADO) ===
cat("=== MÉTRICAS DEL MODELO - DÍA", dia_actual, "===\n")
cat("Precisión:", round(accuracyd10 * 100, 1), "%\n")
cat("Sensibilidad:", round(sensitivityd10 * 100, 1), "%\n")
cat("Especificidad:", round(specificityd10 * 100, 1), "%\n")
cat("Longitudes importantes (>50%):", num_selectedd10, "\n")
cat("Rango:", min(wavelengths[selected_bandsd10]), "-", max(wavelengths[selected_bandsd10]), "nm\n")


# Top 10 bandas más importantes (ACTUALIZADO)
top_bandsd10 <- order(importanced10, decreasing = TRUE)[1:10]
cat("\nTop 10 bandas:\n")
for(i in 1:10) {
  idx <- top_bandsd10[i]
  cat(sprintf("%.0f nm: π=%.3f, |β_eff|=%.4f\n", 
              wavelengths[idx], inclusion_probd10[idx], importanced10[idx]))
}

# === ANÁLISIS ADICIONAL ===
sparsity_ratiod10 <- mean(inclusion_probd10)
cat("\nSparsity ratio:", round(sparsity_ratiod10, 3), "\n")

# Incertidumbre en selección
inclusion_sdd10 <- apply(pi_j_samplesd10, 2, sd)
high_uncertaintyd10 <- sum(inclusion_sdd10 > 0.3)
cat("Bandas con alta incertidumbre (SD > 0.3):", high_uncertaintyd10, "\n")

# Contribución por región espectral
visible_contribd10 <- sum(importanced10[wavelengths <= 700])
nir_contribd10 <- sum(importanced10[wavelengths > 700 & wavelengths <= 1300])
swir_contribd10 <- sum(importanced10[wavelengths > 1300])
total_contribd10 <- sum(importanced10)

cat("\nContribución por región:\n")
cat("VIS (≤700 nm):", round(visible_contribd10/total_contribd10*100, 1), "%\n")
cat("NIR (700-1300 nm):", round(nir_contribd10/total_contribd10*100, 1), "%\n")
cat("SWIR (>1300 nm):", round(swir_contribd10/total_contribd10*100, 1), "%\n")

# Consistencia con literatura
if(exists("prior_config")) {
  prior_expectationd10 <- prior_config$expected_inclusion
  literature_agreementd10 <- cor(importanced10, prior_expectationd10)
  cat("Correlación con priors:", round(literature_agreementd10, 3), "\n")
}

# ------------------------------------------------------------------------------día 11

fit_d11 = fit

B_0_samplesd11 <- rstan::extract(fit, "B_0")[[1]]
B_j_samplesd11 <- rstan::extract(fit, "B_j")[[1]]
pi_j_samplesd11 <- rstan::extract(fit, "pi_j")[[1]]

# === CALCULAR MÉTRICAS EN R  ===

# 1. Estadísticas posteriores
B_0_meand11 <- mean(B_0_samplesd11)
B_j_meand11 <- colMeans(B_j_samplesd11)
inclusion_probd11 <- colMeans(pi_j_samplesd11)

# 2. Coeficientes efectivos
effective_betad11 <- B_j_meand11 * inclusion_probd11

# 3. Métricas de importancia
importanced11 <- abs(effective_betad11)
selected_bandsd11 <- inclusion_probd11 > 0.5
num_selectedd11 <- sum(selected_bandsd11)

# 4. Predicciones
logit_predd11 <- B_0_meand11 + X_scaled %*% effective_betad11  # Usar X_scaled
y_pred_probd11 <- plogis(logit_predd11)  # inv_logit en R
y_predd11 <- as.numeric(y_pred_probd11 > 0.5)

# 5. Métricas de rendimiento
accuracyd11 <- mean(y_predd11 == y)
sensitivityd11 <- sum(y_predd11 == 1 & y == 1) / sum(y == 1)
specificityd11 <- sum(y_predd11 == 0 & y == 0) / sum(y == 0)

# === REPORTE (ACTUALIZADO) ===
cat("=== MÉTRICAS DEL MODELO - DÍA", dia_actual, "===\n")
cat("Precisión:", round(accuracyd11 * 100, 1), "%\n")
cat("Sensibilidad:", round(sensitivityd11 * 100, 1), "%\n")
cat("Especificidad:", round(specificityd11 * 100, 1), "%\n")
cat("Longitudes importantes (>50%):", num_selectedd11, "\n")
cat("Rango:", min(wavelengths[selected_bandsd11]), "-", max(wavelengths[selected_bandsd11]), "nm\n")


# Top 10 bandas más importantes (ACTUALIZADO)
top_bandsd11 <- order(importanced11, decreasing = TRUE)[1:10]
cat("\nTop 10 bandas:\n")
for(i in 1:10) {
  idx <- top_bandsd11[i]
  cat(sprintf("%.0f nm: π=%.3f, |β_eff|=%.4f\n", 
              wavelengths[idx], inclusion_probd11[idx], importanced11[idx]))
}

# === ANÁLISIS ADICIONAL ===
sparsity_ratiod11 <- mean(inclusion_probd11)
cat("\nSparsity ratio:", round(sparsity_ratiod11, 3), "\n")

# Incertidumbre en selección
inclusion_sdd11 <- apply(pi_j_samplesd11, 2, sd)
high_uncertaintyd11 <- sum(inclusion_sdd11 > 0.3)
cat("Bandas con alta incertidumbre (SD > 0.3):", high_uncertaintyd11, "\n")

# Contribución por región espectral
visible_contribd11 <- sum(importanced11[wavelengths <= 700])
nir_contribd11 <- sum(importanced11[wavelengths > 700 & wavelengths <= 1300])
swir_contribd11 <- sum(importanced11[wavelengths > 1300])
total_contribd11 <- sum(importanced11)

cat("\nContribución por región:\n")
cat("VIS (≤700 nm):", round(visible_contribd11/total_contribd11*100, 1), "%\n")
cat("NIR (700-1300 nm):", round(nir_contribd11/total_contribd11*100, 1), "%\n")
cat("SWIR (>1300 nm):", round(swir_contribd11/total_contribd11*100, 1), "%\n")

# Consistencia con literatura
if(exists("prior_config")) {
  prior_expectationd11 <- prior_config$expected_inclusion
  literature_agreementd11 <- cor(importanced11, prior_expectationd11)
  cat("Correlación con priors:", round(literature_agreementd11, 3), "\n")
}

# ------------------------------------------------------------------------------día 12

fit_d12 = fit

B_0_samplesd12 <- rstan::extract(fit, "B_0")[[1]]
B_j_samplesd12 <- rstan::extract(fit, "B_j")[[1]]
pi_j_samplesd12 <- rstan::extract(fit, "pi_j")[[1]]

# === CALCULAR MÉTRICAS EN R  ===

# 1. Estadísticas posteriores
B_0_meand12 <- mean(B_0_samplesd12)
B_j_meand12 <- colMeans(B_j_samplesd12)
inclusion_probd12 <- colMeans(pi_j_samplesd12)

# 2. Coeficientes efectivos
effective_betad12 <- B_j_meand12 * inclusion_probd12

# 3. Métricas de importancia
importanced12 <- abs(effective_betad12)
selected_bandsd12 <- inclusion_probd12 > 0.5
num_selectedd12 <- sum(selected_bandsd12)

# 4. Predicciones
logit_predd12 <- B_0_meand12 + X_scaled %*% effective_betad12  # Usar X_scaled
y_pred_probd12 <- plogis(logit_predd12)  # inv_logit en R
y_predd12 <- as.numeric(y_pred_probd12 > 0.5)

# 5. Métricas de rendimiento
accuracyd12 <- mean(y_predd12 == y)
sensitivityd12 <- sum(y_predd12 == 1 & y == 1) / sum(y == 1)
specificityd12 <- sum(y_predd12 == 0 & y == 0) / sum(y == 0)

# === REPORTE (ACTUALIZADO) ===
cat("=== MÉTRICAS DEL MODELO - DÍA", dia_actual, "===\n")
cat("Precisión:", round(accuracyd12 * 100, 1), "%\n")
cat("Sensibilidad:", round(sensitivityd12 * 100, 1), "%\n")
cat("Especificidad:", round(specificityd12 * 100, 1), "%\n")
cat("Longitudes importantes (>50%):", num_selectedd12, "\n")
cat("Rango:", min(wavelengths[selected_bandsd12]), "-", max(wavelengths[selected_bandsd12]), "nm\n")


# Top 10 bandas más importantes (ACTUALIZADO)
top_bandsd12 <- order(importanced12, decreasing = TRUE)[1:10]
cat("\nTop 10 bandas:\n")
for(i in 1:10) {
  idx <- top_bandsd12[i]
  cat(sprintf("%.0f nm: π=%.3f, |β_eff|=%.4f\n", 
              wavelengths[idx], inclusion_probd12[idx], importanced12[idx]))
}

# === ANÁLISIS ADICIONAL ===
sparsity_ratiod12 <- mean(inclusion_probd12)
cat("\nSparsity ratio:", round(sparsity_ratiod12, 3), "\n")

# Incertidumbre en selección
inclusion_sdd12 <- apply(pi_j_samplesd12, 2, sd)
high_uncertaintyd12 <- sum(inclusion_sdd12 > 0.3)
cat("Bandas con alta incertidumbre (SD > 0.3):", high_uncertaintyd12, "\n")

# Contribución por región espectral
visible_contribd12 <- sum(importanced12[wavelengths <= 700])
nir_contribd12 <- sum(importanced12[wavelengths > 700 & wavelengths <= 1300])
swir_contribd12 <- sum(importanced12[wavelengths > 1300])
total_contribd12 <- sum(importanced12)

cat("\nContribución por región:\n")
cat("VIS (≤700 nm):", round(visible_contribd12/total_contribd12*100, 1), "%\n")
cat("NIR (700-1300 nm):", round(nir_contribd12/total_contribd12*100, 1), "%\n")
cat("SWIR (>1300 nm):", round(swir_contribd12/total_contribd12*100, 1), "%\n")

# Consistencia con literatura
if(exists("prior_config")) {
  prior_expectationd12 <- prior_config$expected_inclusion
  literature_agreementd12 <- cor(importanced12, prior_expectationd12)
  cat("Correlación con priors:", round(literature_agreementd12, 3), "\n")
}

# ------------------------------------------------------------------------------día 13

fit_d13 = fit

B_0_samplesd13 <- rstan::extract(fit, "B_0")[[1]]
B_j_samplesd13 <- rstan::extract(fit, "B_j")[[1]]
pi_j_samplesd13 <- rstan::extract(fit, "pi_j")[[1]]

# === CALCULAR MÉTRICAS EN R  ===

# 1. Estadísticas posteriores
B_0_meand13 <- mean(B_0_samplesd13)
B_j_meand13 <- colMeans(B_j_samplesd13)
inclusion_probd13 <- colMeans(pi_j_samplesd13)

# 2. Coeficientes efectivos
effective_betad13 <- B_j_meand13 * inclusion_probd13

# 3. Métricas de importancia
importanced13 <- abs(effective_betad13)
selected_bandsd13 <- inclusion_probd13 > 0.5
num_selectedd13 <- sum(selected_bandsd13)

# 4. Predicciones
logit_predd13 <- B_0_meand13 + X_scaled %*% effective_betad13  # Usar X_scaled
y_pred_probd13 <- plogis(logit_predd13)  # inv_logit en R
y_predd13 <- as.numeric(y_pred_probd13 > 0.5)

# 5. Métricas de rendimiento
accuracyd13 <- mean(y_predd13 == y)
sensitivityd13 <- sum(y_predd13 == 1 & y == 1) / sum(y == 1)
specificityd13 <- sum(y_predd13 == 0 & y == 0) / sum(y == 0)

# === REPORTE (ACTUALIZADO) ===
cat("=== MÉTRICAS DEL MODELO - DÍA", dia_actual, "===\n")
cat("Precisión:", round(accuracyd13 * 100, 1), "%\n")
cat("Sensibilidad:", round(sensitivityd13 * 100, 1), "%\n")
cat("Especificidad:", round(specificityd13 * 100, 1), "%\n")
cat("Longitudes importantes (>50%):", num_selectedd13, "\n")
cat("Rango:", min(wavelengths[selected_bandsd13]), "-", max(wavelengths[selected_bandsd13]), "nm\n")


# Top 10 bandas más importantes (ACTUALIZADO)
top_bandsd13 <- order(importanced13, decreasing = TRUE)[1:10]
cat("\nTop 10 bandas:\n")
for(i in 1:10) {
  idx <- top_bandsd13[i]
  cat(sprintf("%.0f nm: π=%.3f, |β_eff|=%.4f\n", 
              wavelengths[idx], inclusion_probd13[idx], importanced13[idx]))
}

# === ANÁLISIS ADICIONAL ===
sparsity_ratiod13 <- mean(inclusion_probd13)
cat("\nSparsity ratio:", round(sparsity_ratiod13, 3), "\n")

# Incertidumbre en selección
inclusion_sdd13 <- apply(pi_j_samplesd13, 2, sd)
high_uncertaintyd13 <- sum(inclusion_sdd13 > 0.3)
cat("Bandas con alta incertidumbre (SD > 0.3):", high_uncertaintyd13, "\n")

# Contribución por región espectral
visible_contribd13 <- sum(importanced13[wavelengths <= 700])
nir_contribd13 <- sum(importanced13[wavelengths > 700 & wavelengths <= 1300])
swir_contribd13 <- sum(importanced13[wavelengths > 1300])
total_contribd13 <- sum(importanced13)

cat("\nContribución por región:\n")
cat("VIS (≤700 nm):", round(visible_contribd13/total_contribd13*100, 1), "%\n")
cat("NIR (700-1300 nm):", round(nir_contribd13/total_contribd13*100, 1), "%\n")
cat("SWIR (>1300 nm):", round(swir_contribd13/total_contribd13*100, 1), "%\n")

# Consistencia con literatura
if(exists("prior_config")) {
  prior_expectationd13 <- prior_config$expected_inclusion
  literature_agreementd13 <- cor(importanced13, prior_expectationd13)
  cat("Correlación con priors:", round(literature_agreementd13, 3), "\n")
}


# ------------------------------------------------------------------------------día 14


fit_d14 = fit 


B_0_samplesd14 <- rstan::extract(fit, "B_0")[[1]]
B_j_samplesd14 <- rstan::extract(fit, "B_j")[[1]]
pi_j_samplesd14 <- rstan::extract(fit, "pi_j")[[1]]

# === CALCULAR MÉTRICAS EN R  ===

# 1. Estadísticas posteriores
B_0_meand14 <- mean(B_0_samplesd14)
B_j_meand14 <- colMeans(B_j_samplesd14)
inclusion_probd14 <- colMeans(pi_j_samplesd14)

# 2. Coeficientes efectivos
effective_betad14 <- B_j_meand14 * inclusion_probd14

# 3. Métricas de importancia
importanced14 <- abs(effective_betad14)
selected_bandsd14 <- inclusion_probd14 > 0.5
num_selectedd14 <- sum(selected_bandsd14)

# 4. Predicciones
logit_predd14 <- B_0_meand14 + X_scaled %*% effective_betad14  # Usar X_scaled
y_pred_probd14 <- plogis(logit_predd14)  # inv_logit en R
y_predd14 <- as.numeric(y_pred_probd14 > 0.5)

# 5. Métricas de rendimiento
accuracyd14 <- mean(y_predd14 == y)
sensitivityd14 <- sum(y_predd14 == 1 & y == 1) / sum(y == 1)
specificityd14 <- sum(y_predd14 == 0 & y == 0) / sum(y == 0)

# === REPORTE (ACTUALIZADO) ===
cat("=== MÉTRICAS DEL MODELO - DÍA", dia_actual, "===\n")
cat("Precisión:", round(accuracyd14 * 100, 1), "%\n")
cat("Sensibilidad:", round(sensitivityd14 * 100, 1), "%\n")
cat("Especificidad:", round(specificityd14 * 100, 1), "%\n")
cat("Longitudes importantes (>50%):", num_selectedd14, "\n")
cat("Rango:", min(wavelengths[selected_bandsd14]), "-", max(wavelengths[selected_bandsd14]), "nm\n")


# Top 10 bandas más importantes (ACTUALIZADO)
top_bandsd14 <- order(importanced14, decreasing = TRUE)[1:10]
cat("\nTop 10 bandas:\n")
for(i in 1:10) {
  idx <- top_bandsd14[i]
  cat(sprintf("%.0f nm: π=%.3f, |β_eff|=%.4f\n", 
              wavelengths[idx], inclusion_probd14[idx], importanced14[idx]))
}

# === ANÁLISIS ADICIONAL ===
sparsity_ratiod14 <- mean(inclusion_probd14)
cat("\nSparsity ratio:", round(sparsity_ratiod14, 3), "\n")

# Incertidumbre en selección
inclusion_sdd14 <- apply(pi_j_samplesd14, 2, sd)
high_uncertaintyd14 <- sum(inclusion_sdd14 > 0.3)
cat("Bandas con alta incertidumbre (SD > 0.3):", high_uncertaintyd14, "\n")

# Contribución por región espectral
visible_contribd14 <- sum(importanced14[wavelengths <= 700])
nir_contribd14 <- sum(importanced14[wavelengths > 700 & wavelengths <= 1300])
swir_contribd14 <- sum(importanced14[wavelengths > 1300])
total_contribd14 <- sum(importanced14)

cat("\nContribución por región:\n")
cat("VIS (≤700 nm):", round(visible_contribd14/total_contribd14*100, 1), "%\n")
cat("NIR (700-1300 nm):", round(nir_contribd14/total_contribd14*100, 1), "%\n")
cat("SWIR (>1300 nm):", round(swir_contribd14/total_contribd14*100, 1), "%\n")

# Consistencia con literatura
if(exists("prior_config")) {
  prior_expectationd14 <- prior_config$expected_inclusion
  literature_agreementd14 <- cor(importanced14, prior_expectationd14)
  cat("Correlación con priors:", round(literature_agreementd14, 3), "\n")
}

# ------------------------------------------------------------------------------día 15


fit_d15 = fit

B_0_samplesd15 <- rstan::extract(fit, "B_0")[[1]]
B_j_samplesd15 <- rstan::extract(fit, "B_j")[[1]]
pi_j_samplesd15 <- rstan::extract(fit, "pi_j")[[1]]

# === CALCULAR MÉTRICAS EN R  ===

# 1. Estadísticas posteriores
B_0_meand15 <- mean(B_0_samplesd15)
B_j_meand15 <- colMeans(B_j_samplesd15)
inclusion_probd15 <- colMeans(pi_j_samplesd15)

# 2. Coeficientes efectivos
effective_betad15 <- B_j_meand15 * inclusion_probd15

# 3. Métricas de importancia
importanced15 <- abs(effective_betad15)
selected_bandsd15 <- inclusion_probd15 > 0.5
num_selectedd15 <- sum(selected_bandsd15)

# 4. Predicciones
logit_predd15 <- B_0_meand15 + X_scaled %*% effective_betad15  # Usar X_scaled
y_pred_probd15 <- plogis(logit_predd15)  # inv_logit en R
y_predd15 <- as.numeric(y_pred_probd15 > 0.5)

# 5. Métricas de rendimiento
accuracyd15 <- mean(y_predd15 == y)
sensitivityd15 <- sum(y_predd15 == 1 & y == 1) / sum(y == 1)
specificityd15 <- sum(y_predd15 == 0 & y == 0) / sum(y == 0)

# === REPORTE (ACTUALIZADO) ===
cat("=== MÉTRICAS DEL MODELO - DÍA", dia_actual, "===\n")
cat("Precisión:", round(accuracyd15 * 100, 1), "%\n")
cat("Sensibilidad:", round(sensitivityd15 * 100, 1), "%\n")
cat("Especificidad:", round(specificityd15 * 100, 1), "%\n")
cat("Longitudes importantes (>50%):", num_selectedd15, "\n")
cat("Rango:", min(wavelengths[selected_bandsd15]), "-", max(wavelengths[selected_bandsd15]), "nm\n")


# Top 10 bandas más importantes (ACTUALIZADO)
top_bandsd15 <- order(importanced15, decreasing = TRUE)[1:10]
cat("\nTop 10 bandas:\n")
for(i in 1:10) {
  idx <- top_bandsd15[i]
  cat(sprintf("%.0f nm: π=%.3f, |β_eff|=%.4f\n", 
              wavelengths[idx], inclusion_probd15[idx], importanced15[idx]))
}

# === ANÁLISIS ADICIONAL ===
sparsity_ratiod15 <- mean(inclusion_probd15)
cat("\nSparsity ratio:", round(sparsity_ratiod15, 3), "\n")

# Incertidumbre en selección
inclusion_sdd15 <- apply(pi_j_samplesd15, 2, sd)
high_uncertaintyd15 <- sum(inclusion_sdd15 > 0.3)
cat("Bandas con alta incertidumbre (SD > 0.3):", high_uncertaintyd15, "\n")

# Contribución por región espectral
visible_contribd15 <- sum(importanced15[wavelengths <= 700])
nir_contribd15 <- sum(importanced15[wavelengths > 700 & wavelengths <= 1300])
swir_contribd15 <- sum(importanced15[wavelengths > 1300])
total_contribd15 <- sum(importanced15)

cat("\nContribución por región:\n")
cat("VIS (≤700 nm):", round(visible_contribd15/total_contribd15*100, 1), "%\n")
cat("NIR (700-1300 nm):", round(nir_contribd15/total_contribd15*100, 1), "%\n")
cat("SWIR (>1300 nm):", round(swir_contribd15/total_contribd15*100, 1), "%\n")

# Consistencia con literatura
if(exists("prior_config")) {
  prior_expectationd15 <- prior_config$expected_inclusion
  literature_agreementd15 <- cor(importanced15, prior_expectationd15)
  cat("Correlación con priors:", round(literature_agreementd15, 3), "\n")
}


#====================================================================================
#====================================================================================

#días 10 - etapa tardía

fit_d10t = fit

B_0_samplesd10t <- rstan::extract(fit, "B_0")[[1]]
B_j_samplesd10t <- rstan::extract(fit, "B_j")[[1]]
pi_j_samplesd10t <- rstan::extract(fit, "pi_j")[[1]]

# === CALCULAR MÉTRICAS EN R  ===

# 1. Estadísticas posteriores
B_0_meand10t <- mean(B_0_samplesd10t)
B_j_meand10t <- colMeans(B_j_samplesd10t)
inclusion_probd10t <- colMeans(pi_j_samplesd10t)

# 2. Coeficientes efectivos
effective_betad10t <- B_j_meand10t * inclusion_probd10t

# 3. Métricas de importancia
importanced10t <- abs(effective_betad10t)
selected_bandsd10t <- inclusion_probd10t > 0.5
num_selectedd10t <- sum(selected_bandsd10t)

# 4. Predicciones
logit_predd10t <- B_0_meand10t + X_scaled %*% effective_betad10t  # Usar X_scaled
y_pred_probd10t <- plogis(logit_predd10t)  # inv_logit en R
y_predd10t <- as.numeric(y_pred_probd10t > 0.5)

# 5. Métricas de rendimiento
accuracyd10t <- mean(y_predd10t == y)
sensitivityd10t <- sum(y_predd10t == 1 & y == 1) / sum(y == 1)
specificityd10t <- sum(y_predd10t == 0 & y == 0) / sum(y == 0)

# === REPORTE (ACTUALIZADO) ===
cat("=== MÉTRICAS DEL MODELO - DÍA", dia_actual, "===\n", "Precisión:", round(accuracyd10t * 100, 1), "%\n",
"Sensibilidad:", round(sensitivityd10t * 100, 1), "%\n",
"Especificidad:", round(specificityd10t * 100, 1), "%\n",
"Longitudes importantes (>50%):", num_selectedd10t, "\n")


