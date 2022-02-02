#' This function fit the linear, polynonial, exponential and hill models
#' @param dataframe dataframe containing the data for the fitting. It has to be a two column dataframe where the first column is the doses and the second one the expression values
#' @param mNames vector of character with the name of the models to be fitted
#' @param sorted character string determining according to which criterion the model fits are ranked.
#' @param linreg logical indicating whether or not additionally polynomial regression models (linear, quadratic, and cubic models) should be fitted (they could be useful for a kind of informal lack-of-test consideration for the models specified, capturing unexpected departures).
#' @param powreg logical indicating whether or not additionally power regression models should be fitted (they could be useful for a kind of informal lack-of-test consideration for the models specified, capturing unexpected departures).
#' @param pow power to be used in the power regression model
#' @param expreg logical indicating whether or not additionally exponential regression models should be fitted (they could be useful for a kind of informal lack-of-test consideration for the models specified, capturing unexpected departures).
#' @param hillreg logical indicating whether or not additionally hill regression models should be fitted (they could be useful for a kind of informal lack-of-test consideration for the models specified, capturing unexpected departures).
#' @param hillN  N value for hill model
#' @param Kd Kd value for hill model
#' @param icfct function for supplying the information criterion to be used. AIC and BIC are two options.
#' @param contData boolean value. TRUE for continuous data
#' @return a list containing a matrix with evaluation metrics and a list with the fitted models
#' @importFrom stats AIC
#' @importFrom stats logLik
#' @importFrom stats anova
#' @importFrom stats lm
#' @importFrom stats pf
#' @importFrom stats update
#'
#' @export
#'
fit_models_lm = function(dataframe = dataframe, mNames=mNames,
                         pow = pow, hillN = hillN, Kd = Kd,
                         sorted = c("IC", "Res var", "Lack of fit", "no"), 
                         contData = TRUE,icfct = stats::AIC){
  
  sorted <- match.arg(sorted, c("IC", "Res var", "Lack of fit", "no"))
  
  retMat = c()
  
  cnames <- c("logLik", "IC", "Lack of fit")
  if(contData){
    cnames <- c(cnames, "Res var")
  }
  
  FittedModels = list()
  
  linreg = any(mNames %in% c("Linear", "Quadratic", "Cubic"))
  powreg = any(mNames %in% c( "Power2","Power3","Power4"))
  expreg = any(mNames %in% c("Exponential"))
  hillreg = any(mNames %in% c("Hill05","Hill1","Hill2","Hill3","Hill4","Hill5"))
  
  if (linreg) {
    
    linMod = c("Linear", "Quadratic", "Cubic")[c("Linear", "Quadratic", "Cubic") %in% mNames]
    
    linFitList <- list("Linear" = stats::lm(expr ~ dose, data = dataframe),
                       "Quadratic" = stats::lm(expr ~ dose + I(dose * dose), data = dataframe),
                       "Cubic" = stats::lm(expr ~ dose + I(dose * dose) + I(dose * dose *
                                                                              dose), data = dataframe))
    FittedModels = c(FittedModels, linFitList[linMod])
    
    linModMat = do.call(rbind,lapply(linFitList[linMod], function(listObj) {
        c(stats::logLik(listObj), icfct(listObj), NA, (summary(listObj)$sigma)^2)
      }))
    
    pureErr = c()
    for(i in rownames(linModMat)){
      if(i == "Linear") pureErr = c(pureErr, alr3::pureErrorAnova(linFitList[[1]])[[5]][3])
      if(i == "Quadratic") pureErr = c(pureErr,alr3::pureErrorAnova(linFitList[[2]])[[5]][4])
      if(i == "Cubic") pureErr = c(pureErr,alr3::pureErrorAnova(linFitList[[3]])[[5]][5])
      
    }
    
    linModMat[,3]=  pureErr
    colnames(linModMat) <- cnames[1:4]
    
    retMat <- rbind(retMat, linModMat)
  }
  
  if(powreg){ #fit power model
    
    powerList = list()
    for(pp in pow){
      if(pp ==2) mod_pow = stats::lm(expr ~ I(dose^2), data = dataframe)
      if(pp ==3) mod_pow = stats::lm(expr ~ I(dose^3), data = dataframe)
      if(pp ==4) mod_pow = stats::lm(expr ~ I(dose^4), data = dataframe)
      powerList[[paste("Power",pp,sep="")]] = mod_pow
      powModMat = matrix(c(stats::logLik(mod_pow), icfct(mod_pow), alr3::pureErrorAnova(mod_pow)[[5]][3], (summary(mod_pow)$sigma)^2), nrow = 1, ncol = 4)
      rownames(powModMat) = paste("Power",pp,sep="")
      colnames(powModMat) <- cnames[1:4]
      retMat <- rbind(retMat, powModMat)
    }
    FittedModels = c(FittedModels, powerList)
    
  }
  if(expreg){ #fit power model
    #save(dataframe, file = "expreg.RData")
    exp_x = exp(dataframe$dose)
    if(sum(is.infinite(exp_x))==0){
      exp_mod = stats::lm(expr ~ I(exp(dose)), data = dataframe)
      sm = summary(exp_mod)
      lof = 1 - stats::pf(sm$fstatistic[1], sm$fstatistic[2], sm$fstatistic[3])
      expModMat = matrix(c(stats::logLik(exp_mod), icfct(exp_mod), lof, (summary(exp_mod)$sigma)^2), nrow = 1, ncol = 4) # instead of lof pureErrorAnova(exp_mod)[[5]][3]
      rownames(expModMat) = "Exponential"
      colnames(expModMat) <- cnames[1:4]
      retMat <- rbind(retMat, expModMat)
      FittedModels = c(FittedModels, list("Exponential"=exp_mod))
      
    }else{
      expModMat = matrix(NA,1,4)
      rownames(expModMat) = "Exponential"
      retMat <- rbind(retMat, expModMat)
      FittedModels = c(FittedModels, list("Exponential"=NA))
      
    }
    
  }
  if(hillreg){
    hillFitted = list()
    
    hill_names = c("Hill05","Hill1","Hill2","Hill3","Hill4","Hill5")
    hill_params = c(0.5,1,2,3,4,5)
    
    for(hN in hillN){
      
      if(hN == 0.5)	hill_mod <- stats::lm(expr ~ I(dose^0.5 / (Kd + dose^0.5)),data = dataframe)
      if(hN == 1)	hill_mod <- stats::lm(expr ~ I(dose^1 / (Kd + dose^1)),data = dataframe)
      if(hN == 2)	hill_mod <- stats::lm(expr ~ I(dose^2 / (Kd + dose^2)),data = dataframe)
      if(hN == 3)	hill_mod <- stats::lm(expr ~ I(dose^3 / (Kd + dose^3)),data = dataframe)
      if(hN == 4)	hill_mod <- stats::lm(expr ~ I(dose^4 / (Kd + dose^4)),data = dataframe)
      if(hN == 5)	hill_mod <- stats::lm(expr ~ I(dose^5 / (Kd + dose^5)),data = dataframe)
      
      hillFitted[[hill_names[which(hill_params %in% hN)]]] = hill_mod
      #effect_plot(y.nls, pred = dose, interval = TRUE, plot.points = TRUE)
      hillModMat = matrix(c(stats::logLik(hill_mod), icfct(hill_mod), alr3::pureErrorAnova(hill_mod)[[5]][3], (summary(hill_mod)$sigma)^2), nrow = 1, ncol = 4)
      if(hN == 0.5) hN = "05"
      rownames(hillModMat) = paste("Hill",hN,sep="")
      colnames(hillModMat) <- cnames[1:4]
      retMat <- rbind(retMat, hillModMat)
    }
    FittedModels = c(FittedModels, hillFitted)
    
  }
  
  
  colnames(retMat) <- cnames
  
  if (sorted != "no" & nrow(retMat)>1) {
    retMat = retMat[order(retMat[, sorted]), ]
  }
  
  return(list(retMat=retMat, FittedModels=FittedModels))
  
}
