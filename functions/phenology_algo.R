
#PHENOLOGY ALGORITHMS

#--------------------------------------
#         Simple GDD model
#--------------------------------------
#Growing degree day model optimization
gdd_model <- function(temp, par) {
  # split out parameters from a simple
  # vector of parameter values
  temp_threshold <- par[1]
  gdd_crit <- par[2]
  
  # accumulate growing degree days for
  # temperature data
  gdd <- cumsum(ifelse(temp > temp_threshold, temp - temp_threshold, 0))
  
  # figure out when the number of growing
  # degree days exceeds the minimum value
  # required for leaf development, only
  # return the first value
  doy <- unlist(which(gdd >= gdd_crit)[1])
  
  return(doy)
}

#--------------------------------------
#   RMSE for simple GDD model
#--------------------------------------
#Phenology model calibration
# run model and compare to true values
# returns the RMSE
rmse_gdd <- function(par, data) {
  
  # split out data
  drivers <- data$drivers
  validation <- data$validation
  
  # calculate phenology predictions
  # and put in a data frame
  predictions <- drivers |>
    group_by(year) |>
    summarise(
      predictions = gdd_model(
        temp = tmean,
        par = par
      )
    )
  
  predictions <- left_join(predictions, validation, by = "year")
  
  rmse <- predictions |>
    summarise(
      rmse = sqrt(mean((predictions - doy)^2, na.rm = TRUE))
    ) |>
    pull(rmse)
  
  # return rmse value
  return(rmse)
}

#--------------------------------------
# Optimization: model "chilling + GDD"
#--------------------------------------
# par = c(chill_threshold, required_chill, gdd_threshold, gdd_crit)
gdd_chill_model <- function(temp, par) {
  
  chill_threshold <- par[1]
  required_chill  <- par[2]
  gdd_threshold   <- par[3]
  gdd_crit        <- par[4]
  
  # chilling phase : we add one unit when tmean < seuil
  chill_units <- cumsum(ifelse(temp < chill_threshold, 1, 0))
  chill_reached <- which(chill_units >= required_chill)[1]
  
  # if the chilling is not reached, it returns NA
  if (is.na(chill_reached)) return(NA_integer_)
  
  # GDD phase begin from the day when the chilling is reached
  # recomuputation of the sum from this day
  temp_after_chill <- temp[chill_reached:length(temp)]
  gdd <- cumsum(ifelse(temp_after_chill > gdd_threshold,
                       temp_after_chill - gdd_threshold, 0))
  
  doy_offset <- which(gdd >= gdd_crit)[1]
  if (is.na(doy_offset)) return(NA_integer_)
  
  doy <- chill_reached + doy_offset - 1
  return(doy)
}

# --------------------------------------
#  RMSE for the model chilling + GDD
# --------------------------------------
rmse_gdd_chill <- function(par, data) {
  
  drivers <- data$drivers
  validation <- data$validation
  
  predictions <- drivers |>
    dplyr::group_by(year) |>
    dplyr::summarise(
      predictions = gdd_chill_model(
        temp = tmean,
        par = par
      )
    )
  
  predictions <- dplyr::left_join(predictions, validation, by = "year")
  
  rmse <- predictions |>
    dplyr::summarise(
      rmse = sqrt(mean((predictions - doy)^2, na.rm = TRUE))
    ) |>
    dplyr::pull(rmse)
  
  return(rmse)
}

# --------------------------------------
#       Photo-thermal GDD model
#par = c(temp_threshold, gdd_crit, photoperiod_crit)
# --------------------------------------
gdd_photo_model <- function(temp, photoperiod, par) {
  
  temp_threshold   <- par[1]
  gdd_crit         <- par[2]
  photoperiod_crit <- par[3]
  
  # accumulate GDD
  gdd <- cumsum(ifelse(temp > temp_threshold, temp - temp_threshold, 0))
  
  # find first day where both conditions are met
  idx <- which(gdd >= gdd_crit & photoperiod >= photoperiod_crit)[1]
  
  if (is.na(idx)) return(NA_integer_)
  
  return(idx)
}

# --------------------------------------
#    RMSE for photo-thermal model
# --------------------------------------
rmse_gdd_photo_fast <- function(par, data) {
  
  drivers    <- data$drivers
  validation <- data$validation
  
  years <- sort(unique(drivers$year))
  n_yrs <- length(years)
  preds <- numeric(n_yrs)
  obs   <- numeric(n_yrs)
  
  for (i in seq_along(years)) {
    yr <- years[i]
    idx <- drivers$year == yr
    
    preds[i] <- gdd_photo_model(
      temp        = drivers$tmean[idx],
      photoperiod = drivers$photoperiod[idx],
      par         = par
    )
    
    obs[i] <- validation$doy[validation$year == yr]
  }
  
  valid <- !is.na(preds) & !is.na(obs)
  sqrt(mean((preds[valid] - obs[valid])^2))
}


# --------------------------------------
# RMSE for photo-thermal model (slower)
# --------------------------------------
gdd_photo_model <- function(temp, photoperiod, par) {
  
  temp_threshold   <- par[1]
  gdd_crit         <- par[2]
  photoperiod_crit <- par[3]
  
  # GDD cumulés
  gdd <- cumsum(ifelse(temp > temp_threshold, temp - temp_threshold, 0))
  
  # premier jour où on a à la fois assez de GDD et assez de photopériode
  idx <- which(gdd >= gdd_crit & photoperiod >= photoperiod_crit)[1]
  
  if (is.na(idx)) return(NA_integer_)
  return(idx)
}

# ----------------------------------------
# RMSE pour le modèle photo-thermal
# version "rapide" (sans dplyr dans la boucle)
# ----------------------------------------
rmse_gdd_photo <- function(par, data) {
  
  drivers    <- data$drivers
  validation <- data$validation
  
  years <- sort(unique(drivers$year))
  n_yrs <- length(years)
  preds <- numeric(n_yrs)
  obs   <- numeric(n_yrs)
  
  for (i in seq_along(years)) {
    yr  <- years[i]
    idx <- drivers$year == yr
    
    preds[i] <- gdd_photo_model(
      temp        = drivers$tmean[idx],
      photoperiod = drivers$photoperiod[idx],
      par         = par
    )
    
    obs[i] <- validation$doy[validation$year == yr]
  }
  
  valid <- !is.na(preds) & !is.na(obs)
  sqrt(mean((preds[valid] - obs[valid])^2))
}

# ----------------------------------------
# Photo-thermal GDD model
# par = c(temp_threshold, gdd_crit, photoperiod_crit)
# ----------------------------------------
gdd_photo_model <- function(temp, photoperiod, par) {
  
  temp_threshold   <- par[1]
  gdd_crit         <- par[2]
  photoperiod_crit <- par[3]
  
  # GDD cumulés (vectorisé)
  gdd <- cumsum(pmax(temp - temp_threshold, 0))
  
  # premier jour où on a à la fois assez de GDD et assez de photopériode
  idx <- which(gdd >= gdd_crit & photoperiod >= photoperiod_crit)[1]
  
  if (is.na(idx)) return(NA_integer_)
  idx
}

# ----------------------------------------
# RMSE pour le modèle photo-thermal
# version "fast" (sans dplyr)
# ----------------------------------------
rmse_gdd_photo <- function(par, data) {
  
  drivers    <- data$drivers
  validation <- data$validation
  
  years <- sort(unique(drivers$year))
  n_yrs <- length(years)
  preds <- rep(NA_real_, n_yrs)
  obs   <- rep(NA_real_, n_yrs)
  
  for (i in seq_along(years)) {
    yr  <- years[i]
    idx <- drivers$year == yr
    
    # si, pour une raison quelconque, pas de données -> on saute
    if (!any(idx)) next
    
    # 1) prédiction du modèle pour cette année
    p <- try(
      gdd_photo_model(
        temp        = drivers$tmean[idx],
        photoperiod = drivers$photoperiod[idx],
        par         = par
      ),
      silent = TRUE
    )
    
    # si le modèle plante ou renvoie une erreur -> on ignore cette année
    if (inherits(p, "try-error")) next
    
    # si le modèle renvoie plusieurs valeurs, on prend la moyenne
    if (length(p) > 1) p <- mean(p, na.rm = TRUE)
    
    preds[i] <- p
    
    # 2) observation pour cette année
    o <- validation$doy[validation$year == yr]
    
    if (length(o) == 0) {
      # pas d'observation -> obs reste NA
      next
    } else if (length(o) > 1) {
      # plusieurs obs pour la même année -> on prend la moyenne
      o <- mean(o, na.rm = TRUE)
    }
    
    obs[i] <- o
  }
  
  # années où on a à la fois une prédiction et une obs
  valid <- is.finite(preds) & is.finite(obs)
  
  # si aucune année valide -> grosse pénalité
  if (!any(valid)) return(1e12)
  
  rmse <- sqrt(mean((preds[valid] - obs[valid])^2))
  
  # par sécurité, si encore NaN / Inf -> pénalité
  if (!is.finite(rmse) || is.na(rmse)) rmse <- 1e12
  
  rmse
}
# -------------------------------------------------
# Leave-one-year-out cross-validation for GDD model
# -------------------------------------------------
loocv_rmse_gdd <- function(par, data) {
  
  drivers    <- data$drivers
  validation <- data$validation
  
  years <- sort(unique(validation$year))
  n_yrs <- length(years)
  
  preds_cv <- numeric(n_yrs)
  obs_cv   <- numeric(n_yrs)
  
  for (i in seq_along(years)) {
    yr_test <- years[i]
    
    # "jeu de test" = une seule année
    test_idx  <- validation$year == yr_test
    obs_cv[i] <- validation$doy[test_idx]
    
    # prédiction pour cette année avec les paramètres globaux "par"
    temp_year <- drivers$tmean[drivers$year == yr_test]
    
    preds_cv[i] <- gdd_model(
      temp = temp_year,
      par  = par
    )
  }
  
  valid <- !is.na(preds_cv) & !is.na(obs_cv)
  sqrt(mean((preds_cv[valid] - obs_cv[valid])^2))
}
