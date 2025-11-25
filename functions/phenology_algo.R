
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

# --------------------------------------
# Triangular temperature response model
# par = c(Tmin, Topt, Tmax, Hcrit)
# --------------------------------------
gdd_tri_model <- function(temp, par) {
  
  Tmin  <- par[1]
  Topt  <- par[2]
  Tmax  <- par[3]
  Hcrit <- par[4]
  
  # contraintes de base : si non respectées -> pas de solution
  if (Tmin >= Topt || Topt >= Tmax) {
    return(NA_integer_)
  }
  
  # taux de développement journalier (0–1)
  dev <- numeric(length(temp))
  
  # zone croissante Tmin–Topt
  idx1 <- temp > Tmin & temp < Topt
  dev[idx1] <- (temp[idx1] - Tmin) / (Topt - Tmin)
  
  # zone décroissante Topt–Tmax
  idx2 <- temp >= Topt & temp < Tmax
  dev[idx2] <- (Tmax - temp[idx2]) / (Tmax - Topt)
  
  # dev = 0 en dehors [Tmin, Tmax]
  
  # somme cumulée
  dev_cum <- cumsum(dev)
  
  # premier jour où on atteint Hcrit
  doy <- which(dev_cum >= Hcrit)[1]
  
  if (is.na(doy)) return(NA_integer_)
  return(doy)
}

# --------------------------------------
# RMSE for triangular model
# --------------------------------------
rmse_gdd_tri <- function(par, data) {
  
  drivers    <- data$drivers
  validation <- data$validation
  
  years <- sort(unique(drivers$year))
  n_yrs <- length(years)
  preds <- rep(NA_real_, n_yrs)
  obs   <- rep(NA_real_, n_yrs)
  
  for (i in seq_along(years)) {
    yr  <- years[i]
    idx <- drivers$year == yr
    
    # si pas de données pour cette année -> on saute
    if (!any(idx)) next
    
    # prédiction triangulaire
    p <- try(
      gdd_tri_model(
        temp = drivers$tmean[idx],
        par  = par
      ),
      silent = TRUE
    )
    
    if (inherits(p, "try-error") || length(p) == 0) next
    preds[i] <- p
    
    # obs PhenoCam pour cette année
    o <- validation$doy[validation$year == yr]
    if (length(o) == 0) next
    if (length(o) > 1) o <- mean(o, na.rm = TRUE)
    obs[i] <- o
  }
  
  valid <- is.finite(preds) & is.finite(obs)
  if (!any(valid)) return(1e12)
  
  rmse <- sqrt(mean((preds[valid] - obs[valid])^2))
  if (!is.finite(rmse) || is.na(rmse)) rmse <- 1e12
  
  rmse
}