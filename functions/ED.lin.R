
control_variance = function(model, constant_variance, training_dose){
  res = residuals(model)
  if(constant_variance){
    ctrl_variance = mean(res^2) # MSE - use residuals of all samples
  }else{
    ctrl_variance = mean(res[training_dose==0]^2) # RMSE - use residuals of all controls samples (dose=0)
  }
  return(ctrl_variance)
}


#' This function estimate effective doses for linear moels
#' @param lmObject fitted object
#' @params dataframe dataframe containing experimental data (dose, exp)
#' @param respLev is a number specifying the response level
#' @param ci confidence interval
#' @param dose dose value
#' @importFrom stats predict
#' @importFrom stats approx
#' @importFrom stats coef
#' @export

"ED.lin" <- function(lmObject, dataframe=dataframe,constant_variance, ci = 0.95, 
                     rl, training_dose, strictly_monotonic = TRUE )
{
  
  monotonic_behaviour = monotonicity(fittedModel = lmObject, dose = training_dose, 
                                     dataframe = dataframe, 
                                     strictly_monotonic= strictly_monotonic)
  
  starting_point = predict(lmObject, data.frame(dose = 0))
  respLev = rl * sqrt(control_variance(lmObject, constant_variance = constant_variance, training_dose = training_dose))
  doseRange = seq(min(training_dose), max(training_dose), length.out = 1000)
  conf_interval <- stats::predict(lmObject, newdata=data.frame(dose=doseRange), interval="confidence",level = ci)
  
  if(monotonic_behaviour == -1){
    decreasing = TRUE
    diff_resp = starting_point - respLev
    
    bmd <-  stats::approx(x = conf_interval[,"fit"], y = doseRange, xout = diff_resp)$y
    bmdl <- stats::approx(x = conf_interval[,"lwr"], y = doseRange, xout = diff_resp)$y
    bmdu <- stats::approx(x = conf_interval[,"upr"], y = doseRange, xout = diff_resp)$y
    
    half_resp = max(conf_interval[,"fit"]) - (abs(max(conf_interval[,"fit"]) - min(conf_interval[,"fit"]))/2)
    half_con <-  stats::approx(x = conf_interval[,"fit"], y = doseRange, xout = half_resp)$y
    
  }else if(monotonic_behaviour == 1){
    decreasing = FALSE
    diff_resp = starting_point + respLev
    
    bmd <-  stats::approx(x = conf_interval[,"fit"], y = doseRange, xout = diff_resp)$y
    bmdl <- stats::approx(x = conf_interval[,"upr"], y = doseRange, xout = diff_resp)$y
    bmdu <- stats::approx(x = conf_interval[,"lwr"], y = doseRange, xout = diff_resp)$y
    
    half_resp = max(conf_interval[,"fit"]) - (abs(max(conf_interval[,"fit"]) - min(conf_interval[,"fit"]))/2)
    half_con <-  stats::approx(x = conf_interval[,"fit"], y = doseRange, xout = half_resp)$y
    
  }else{
    decreasing = NA
    return(NULL)
  }
  
  # if(is.na(bmdl)){
  #   bmdl = min(training_dose)
  # }
  # if(is.na(bmdu)){
  #   bmdu = max(training_dose)
  # }
  
  mm = matrix(data = c(bmd, bmdl, bmdu,half_con,decreasing),nrow = 1,ncol = 5)
  colnames(mm) = c("BMD","BMDL","BMDU","IC50/EC50","Decreasing")
  return(mm)
}

# "ED.lin" <- function(lmObject, respLev, ci = 0.8, dose )
# {
#   yVal <- lmObject$"model"[, 1]
#   xVal <- dose#lmObject$"model"[, 2]
#   parCoef <- coef(lmObject)
#   lparco <- length(parCoef)
#   decreasing <- ((lparco == 2) && (parCoef[lparco] < 0)) || ((lparco == 3) && (parCoef[lparco] > 0))
#   
#   starting_point = mean(yVal[xVal == 0])
#   
#   if(decreasing){
#     diff_resp = starting_point - respLev
#     conf_interval <- predict(lmObject, newdata=data.frame(x=xVal), interval="confidence",level = ci)
#     bmd <-  approx(x = lmObject$fitted.values, y = dose, xout = diff_resp)$y
#     bmdl <- approx(x = conf_interval[,"lwr"], y = dose, xout = diff_resp)$y
#     bmdu <- approx(x = conf_interval[,"upr"], y = dose, xout = diff_resp)$y
#     ic50 = max(yVal) - (abs(max(yVal) - min(yVal))/2)
#     half_con <-  approx(x = lmObject$fitted.values, y = dose, xout = ic50)$y
#   }else{
#     diff_resp = starting_point + respLev
#     conf_interval <- predict(lmObject, newdata=data.frame(x=xVal), interval="confidence",level = ci)
#     bmd <-  approx(x = lmObject$fitted.values, y = dose, xout = diff_resp)$y
#     bmdl <- approx(x = conf_interval[,"upr"], y = dose, xout = diff_resp)$y
#     bmdu <- approx(x = conf_interval[,"lwr"], y = dose, xout = diff_resp)$y
#     ic50 = max(yVal) - (abs(max(yVal) - min(yVal))/2)
#     half_con <-  approx(x = lmObject$fitted.values, y = dose, xout = ic50)$y
#     
#   }
#   
#   #library(chemCal)
#   
#   # plot(xVal, yVal, xlab="x", ylab="y", main="Regression")
#   # #abline(lmObject, col="lightblue")
#   # lines(xVal, conf_interval[,1],col="lightblue")
#   # lines(xVal, conf_interval[,2], col="blue", lty=2)
#   # lines(xVal, conf_interval[,3], col="blue", lty=2)
#   # abline(h=diff_resp)
#   
#   if(is.na(bmdl)){
#     bmdl = min(dose)
#   }
#   if(is.na(bmdu)){
#     bmdu = max(dose)
#   }
#   
#   mm = matrix(data = c(bmd, bmdl, bmdu,half_con,decreasing),nrow = 1,ncol = 5)
#   colnames(mm) = c("BMD","BMDL","BMDU","IC50/EC50","Decreasing")
#   return(mm)
# }


#' #' This function estimate effective doses for linear moels
#' #' @param lmObject fitted object
#' #' @param respLev is a number specifying the response level
#' #' @param ci confidence interval
#' #' @param dose dose value
#' #' @importFrom stats predict
#' #' @importFrom stats approx
#' #' @importFrom stats coef
#' #' @export
#' 
#' EDLin <- function(lmObject, respLev, ci = 0.8, dose )
#' {
#'   # print(dose)
#'   yVal <- lmObject$"model"[, 1]
#'   # print(yVal)
#'   xVal <- dose#lmObject$"model"[, 2]
#'   # print(xVal)
#'   #parCoef <- stats::coef(lmObject)
#'   # print(parCoef)
#'   
#'   #lparco <- length(parCoef)
#'   # print(lparco)
#'   
#'   monotonic_behaviour = monotonicity(lmObject, xVal)
#'   #decreasing <- ((lparco == 2) && (parCoef[lparco] < 0)) || ((lparco == 3) && (parCoef[lparco] > 0))
#'   # print(decreasing)
#'   
#'   starting_point = mean(yVal[xVal == 0])
#'   # print("starting_point")
#'   # print(starting_point)
#'   
#'   if(monotonic_behaviour == -1){
#'     decreasing = TRUE
#'     # print("decreasing")
#'     diff_resp = starting_point - respLev
#'     
#'     
#'     # conf_interval <- stats::predict(lmObject, newdata=data.frame(x=xVal), interval="confidence",level = ci)
#'     conf_interval <- stats::predict(lmObject, newdata=data.frame(dose = xVal), interval="confidence",level = ci)
#'     
#'     # print("conf_interval decreasing")
#'     # print(conf_interval)
#'     
#'     bmd <-  stats::approx(x = lmObject$fitted.values, y = dose, xout = diff_resp)$y
#'     # print("bmd decreasing")
#'     # print(bmd)
#'     
#'     bmdl <- stats::approx(x = conf_interval[,"lwr"], y = dose, xout = diff_resp)$y
#'     # print("bmdl decreasing")
#'     # print(bmdl)
#'     
#'     bmdu <- stats::approx(x = conf_interval[,"upr"], y = dose, xout = diff_resp)$y
#'     # print("bmdl decreasing")
#'     # print(bmdu)
#'     
#'     ic50 = max(yVal) - (abs(max(yVal) - min(yVal))/2)
#'     # print("ic50 decreasing")
#'     # print(ic50)
#'     
#'     half_con <-  stats::approx(x = lmObject$fitted.values, y = dose, xout = ic50)$y
#'     # print("half_con decreasing")
#'     # print(half_con)
#'   }else if(monotonic_behaviour == 1){
#'     decreasing = FALSE
#'     # print("increasing")
#'     diff_resp = starting_point + respLev
#'     conf_interval <- stats::predict(lmObject, newdata=data.frame(dose=xVal), interval="confidence",level = ci)
#'     
#'     bmd <-  stats::approx(x = lmObject$fitted.values, y = dose, xout = diff_resp)$y
#'     bmdl <- stats::approx(x = conf_interval[,"upr"], y = dose, xout = diff_resp)$y
#'     bmdu <- stats::approx(x = conf_interval[,"lwr"], y = dose, xout = diff_resp)$y
#'     ic50 = max(yVal) - (abs(max(yVal) - min(yVal))/2)
#'     half_con <-  stats::approx(x = lmObject$fitted.values, y = dose, xout = ic50)$y
#'     
#'   }else{
#'     decreasing = NA
#'     return(NULL)
#'   }
#'   
#'   #library(chemCal)
#'   
#'   # plot(xVal, yVal, xlab="x", ylab="y", main="Regression")
#'   # #abline(lmObject, col="lightblue")
#'   # lines(xVal, conf_interval[,1],col="lightblue")
#'   # lines(xVal, conf_interval[,2], col="blue", lty=2)
#'   # lines(xVal, conf_interval[,3], col="blue", lty=2)
#'   # abline(h=diff_resp)
#'   
#'   if(is.na(bmdl)){
#'     bmdl = min(dose)
#'   }
#'   if(is.na(bmdu)){
#'     bmdu = max(dose)
#'   }
#'   
#'   mm = matrix(data = c(bmd, bmdl, bmdu,half_con,decreasing),nrow = 1,ncol = 5)
#'   colnames(mm) = c("BMD","BMDL","BMDU","IC50/EC50","Decreasing")
#'   return(mm)
#' }
#' 
