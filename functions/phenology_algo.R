
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
  
  yrs   <- sort(unique(drivers$year))
  n_yrs <- length(yrs)
  
  preds <- numeric(n_yrs)
  obs   <- numeric(n_yrs)
  
  for (i in seq_along(yrs)) {
    yr  <- yrs[i]
    idx <- drivers$year == yr
    
    preds[i] <- gdd_photo_model(
      temp        = drivers$tmean[idx],
      photoperiod = drivers$photoperiod[idx],
      par         = par
    )
    
    # DOY observé pour cette année
    obs[i] <- validation$doy[validation$year == yr]
  }
  
  valid <- !is.na(preds) & !is.na(obs)
  sqrt(mean((preds[valid] - obs[valid])^2))
}
