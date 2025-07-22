# MODELO XGBOOST
# Comparación: Datos Espectrales vs Datos Espectrales + Info Bayesiana

library(xgboost)
library(caret)
library(dplyr)


# PREPARACIÓN DEL DATASET 

create_unified_dataset <- function() {
  
  # Días disponibles
  dias_disponibles <- 1:15
  
  # Listas
  all_spectral_data <- list()
  all_bayesian_data <- list()
  all_labels <- list()
  all_days <- list()
  
  for(dia in dias_disponibles) {
    
    tryCatch({
      spectral_day <- X_scaled
      
      # Variables bayesianas del día
      var_inclusion <- paste0("inclusion_probd", dia)
      var_effective <- paste0("effective_betad", dia)
      var_importance <- paste0("importanced", dia)
      
      if(exists(var_inclusion) && exists(var_effective) && exists(var_importance)) {
        
        inclusion_prob <- get(var_inclusion)
        effective_beta <- get(var_effective)
        importance <- get(var_importance)
        
        # Información bayesiana resumida para agregar como features
        bayesian_summary <- data.frame(
          # Métricas de inclusión
          inclusion_mean = mean(inclusion_prob, na.rm = TRUE),
          inclusion_max = max(inclusion_prob, na.rm = TRUE),
          inclusion_median = median(inclusion_prob, na.rm = TRUE),
          inclusion_q75 = quantile(inclusion_prob, 0.75, na.rm = TRUE),
          
          # Métricas de coeficientes efectivos
          effective_mean = mean(abs(effective_beta), na.rm = TRUE),
          effective_max = max(abs(effective_beta), na.rm = TRUE),
          effective_std = sd(effective_beta, na.rm = TRUE),
          
          # Métricas de importancia
          importance_mean = mean(importance, na.rm = TRUE),
          importance_max = max(importance, na.rm = TRUE),
          importance_std = sd(importance, na.rm = TRUE),
          
          # Métricas de selección
          sparsity_ratio = mean(inclusion_prob > 0.5, na.rm = TRUE),
          n_selected = sum(inclusion_prob > 0.5, na.rm = TRUE),
          
          # Información del día
          day_post_inoculation = dia,
          day_sin = sin(2 * pi * dia / 15),  # Componente cíclica
          day_cos = cos(2 * pi * dia / 15)   # Componente cíclica
        )
        
        # Replicar información bayesiana para cada muestra del día
        n_samples <- nrow(spectral_day)
        bayesian_features <- bayesian_summary[rep(1, n_samples), ]
        
        # Almacenar datos
        all_spectral_data[[length(all_spectral_data) + 1]] <- spectral_day
        all_bayesian_data[[length(all_bayesian_data) + 1]] <- bayesian_features
        all_labels[[length(all_labels) + 1]] <- y  # Etiquetas Fusarium/No-Fusarium
        all_days[[length(all_days) + 1]] <- rep(dia, n_samples)
        
        cat("Día", dia, "procesado -", n_samples, "muestras\n")
      }
      
    }, error = function(e) {
      cat("Error procesando día", dia, ":", e$message, "\n")
    })
  }
  
  # Combinar todos los datos
  if(length(all_spectral_data) > 0) {
    
    # Dataset solo con datos espectrales
    spectral_matrix <- do.call(rbind, all_spectral_data)
    days_vector <- do.call(c, all_days)
    labels_vector <- do.call(c, all_labels)
    
    # Agregar información del día como feature
    dataset_spectral <- data.frame(
      spectral_matrix,
      day_post_inoculation = days_vector,
      day_sin = sin(2 * pi * days_vector / 15),
      day_cos = cos(2 * pi * days_vector / 15),
      fusarium = labels_vector
    )
    
    # Dataset con información bayesiana
    bayesian_matrix <- do.call(rbind, all_bayesian_data)
    
    dataset_bayesian <- data.frame(
      spectral_matrix,
      bayesian_matrix,
      fusarium = labels_vector
    )
    
    cat("- Muestras totales:", nrow(dataset_spectral), "\n")
    cat("- Features espectrales:", ncol(spectral_matrix), "\n")
    cat("- Features bayesianas:", ncol(bayesian_matrix) - 3, "\n")  # -3 por day info duplicada
    cat("- Casos Fusarium:", sum(labels_vector), "\n")
    cat("- Casos No-Fusarium:", sum(1 - labels_vector), "\n")
    
    return(list(
      spectral = dataset_spectral,
      bayesian = dataset_bayesian,
      success = TRUE
    ))
    
  } else {
    cat("No se pudieron crear los datasets\n")
    return(list(success = FALSE))
  }
}


# FUNCIÓN DE ENTRENAMIENTO Y EVALUACIÓN

train_robust_xgboost <- function(data, model_name, use_cv = TRUE) {
  
  cat("\nEntrenando", model_name, "...\n")
  
  # Separar features y target
  X <- as.matrix(data[, !names(data) %in% "fusarium"])
  y <- data$fusarium
  

  set.seed(42)
  train_idx <- createDataPartition(y, p = 0.75, list = FALSE)
  
  X_train <- X[train_idx, ]
  y_train <- y[train_idx]
  X_test <- X[-train_idx, ]
  y_test <- y[-train_idx]
  
  # Parámetros optimizados
  xgb_params <- list(
    objective = "binary:logistic",
    eval_metric = "auc",
    max_depth = 4,              
    eta = 0.05,                 
    subsample = 0.8,
    colsample_bytree = 0.6,     
    min_child_weight = 3,       
    reg_alpha = 0.1,           
    reg_lambda = 0.1,           
    scale_pos_weight = 1        
  )
  
  # Validación cruzada
  if(use_cv) {
    cv_result <- xgb.cv(
      data = X_train,
      label = y_train,
      params = xgb_params,
      nrounds = 500,
      nfold = 5,
      early_stopping_rounds = 20,
      verbose = 0,
      stratified = TRUE
    )
    
    best_nrounds <- cv_result$best_iteration
    cat("Mejor número de rounds:", best_nrounds, "\n")
  } else {
    best_nrounds <- 100
  }
  
  # Entrenar modelo final
  model <- xgboost(
    data = X_train,
    label = y_train,
    params = xgb_params,
    nrounds = best_nrounds,
    verbose = 0
  )
  
  # Predicciones
  pred_prob <- predict(model, X_test)
  pred_class <- ifelse(pred_prob > 0.5, 1, 0)
  
  # Métricas de evaluación
  accuracy <- mean(pred_class == y_test)
  
  # Sensitivity y Specificity (manejar casos edge)
  tp <- sum(pred_class == 1 & y_test == 1)
  fp <- sum(pred_class == 1 & y_test == 0)
  tn <- sum(pred_class == 0 & y_test == 0)
  fn <- sum(pred_class == 0 & y_test == 1)
  
  sensitivity <- ifelse(tp + fn > 0, tp / (tp + fn), NA)
  specificity <- ifelse(tn + fp > 0, tn / (tn + fp), NA)
  precision <- ifelse(tp + fp > 0, tp / (tp + fp), NA)
  f1_score <- ifelse(!is.na(sensitivity) & !is.na(precision) & 
                       (sensitivity + precision) > 0, 
                     2 * (sensitivity * precision) / (sensitivity + precision), NA)
  
  # AUC
  auc_value <- tryCatch({
    if(length(unique(y_test)) > 1) {
      roc_obj <- pROC::roc(y_test, pred_prob, quiet = TRUE)
      as.numeric(pROC::auc(roc_obj))
    } else {
      NA
    }
  }, error = function(e) NA)
  
  # Importancia de features (top 20)
  importance_matrix <- xgb.importance(model = model)
  top_features <- head(importance_matrix$Feature, 20)
  
  return(list(
    model_name = model_name,
    accuracy = accuracy,
    sensitivity = sensitivity,
    specificity = specificity,
    precision = precision,
    f1_score = f1_score,
    auc = auc_value,
    n_features = ncol(X),
    n_train = nrow(X_train),
    n_test = nrow(X_test),
    best_nrounds = best_nrounds,
    top_features = top_features,
    model = model  # Guardar modelo para análisis posterior
  ))
}


#############################################
datasets <- create_unified_dataset()

if(datasets$success) {
  
  # Entrenar modelo solo con datos espectrales
  cat("\n" %+% "="*60)
  result_spectral <- train_robust_xgboost(
    datasets$spectral, 
    "XGBoost_Spectral"
  )
  
  # Entrenar modelo con datos espectrales + información bayesiana
  cat("\n" %+% "="*60)
  result_bayesian <- train_robust_xgboost(
    datasets$bayesian, 
    "XGBoost_Spectral_Bayesian"
  )

  # COMPARACIÓN DE RESULTADOS

  
  # Crear tabla comparativa
  comparison <- data.frame(
    Metric = c("Accuracy", "Sensitivity", "Specificity", "Precision", "F1-Score", "AUC", "N_Features"),
    Spectral_Only = c(
      round(result_spectral$accuracy, 4),
      round(result_spectral$sensitivity, 4),
      round(result_spectral$specificity, 4),
      round(result_spectral$precision, 4),
      round(result_spectral$f1_score, 4),
      round(result_spectral$auc, 4),
      result_spectral$n_features
    ),
    Spectral_Bayesian = c(
      round(result_bayesian$accuracy, 4),
      round(result_bayesian$sensitivity, 4),
      round(result_bayesian$specificity, 4),
      round(result_bayesian$precision, 4),
      round(result_bayesian$f1_score, 4),
      round(result_bayesian$auc, 4),
      result_bayesian$n_features
    )
  )
  
  # Calcular mejoras
  comparison$Improvement <- round(
    comparison$Spectral_Bayesian - comparison$Spectral_Only, 4
  )
  comparison$Improvement_Pct <- round(
    (comparison$Improvement / comparison$Spectral_Only) * 100, 2
  )
  
  print(comparison)
}