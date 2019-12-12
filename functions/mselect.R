mselect2 = function (object, fctList = NULL, nested = FALSE, 
                     sorted = c("IC", "Res var", "Lack of fit", "no"), 
                     linreg = FALSE, powreg = TRUE, pow = 2, expreg = TRUE,hillreg = TRUE,
                     hillN = 2, Kd = 10,
                     icfct = AIC) 
{
  sorted <- match.arg(sorted)
  if (!is.logical(nested)) {
    stop("'nested' argument takes only the values: FALSE, TRUE")
  }
  contData <- identical(object$type, "continuous")
  nestedInd <- 3 + contData + nested
  mc <- match.call()
  lenFL <- length(fctList)
  retMat <- matrix(0, lenFL + 1, 3 + contData + nested)
  retMat[1, 1] <- logLik(object)
  retMat[1, 2] <- icfct(object)
  retMat[1, 3] <- modelFit(object)[2, 5]
  if (contData) {
    tryRV <- try(summary(object)$resVar, silent = TRUE)
    if (!inherits("tryRV", "try-error")) {
      retMat[1, 4] <- tryRV
    }
    else {
      retMat[1, 4] <- NA
    }
  }
  if (nested) {
    retMat[1, nestedInd] <- NA
  }
  fctList2 <- rep("", lenFL + 1)
  fctList2[1] <- object$fct$name
  if (!is.null(fctList)) {
    prevObj <- object
    for (i in 1:lenFL) {
      tempObj <- try(update(object, fct = fctList[[i]]), 
                     silent = TRUE)
      fctList2[i + 1] <- fctList[[i]]$name
      if (!inherits(tempObj, "try-error")) {
        retMat[i + 1, 1] <- logLik(tempObj)
        retMat[i + 1, 2] <- icfct(tempObj)
        retMat[i + 1, 3] <- modelFit(tempObj)[2, 5]
        if (contData) {
          tryRV2 <- try(summary(tempObj)$resVar, silent = TRUE)
          if (!inherits("tryRV2", "try-error")) {
            retMat[i + 1, 4] <- tryRV2
          }
          else {
            retMat[i + 1, 4] <- NA
          }
        }
        if (nested) {
          retMat[i + 1, nestedInd] <- anova(prevObj, 
                                            tempObj, details = FALSE)[2, 5]
        }
      }
      else {
        retMat[i + 1, ] <- NA
      }
      prevObj <- tempObj
    }
  }
  rownames(retMat) <- as.vector(unlist(fctList2))
  cnames <- c("logLik", "IC", "Lack of fit")
  if (contData) {
    cnames <- c(cnames, "Res var")
  }
  if (nested) {
    cnames <- c(cnames, "Nested F test")
  }
  colnames(retMat) <- cnames
  if (linreg) {
    drcData <- as.data.frame(object$data[, c(2, 1)])
    names(drcData) <- c("yVec", "xVec")
    linFitList <- list(lm(yVec ~ xVec, data = drcData), 
                       lm(yVec ~ xVec + I(xVec * xVec), data = drcData), 
                       lm(yVec ~ xVec + I(xVec * xVec) + I(xVec * xVec * 
                                                             xVec), data = drcData))
    linModMat <- matrix(unlist(lapply(linFitList, function(listObj) {
      c(logLik(listObj), icfct(listObj), NA, (summary(listObj)$sigma)^2)
    })), 3, 4, byrow = TRUE)
    linModMat[1,3] = alr3::pureErrorAnova(linFitList[[1]])[[5]][3] #1
    linModMat[2,3] = alr3::pureErrorAnova(linFitList[[2]])[[5]][4] #2
    linModMat[3,3] = alr3::pureErrorAnova(linFitList[[3]])[[5]][5] #3
    
    rownames(linModMat) <- c("Linear", "Quadratic", "Cubic")
    colnames(linModMat) <- cnames[1:4]
    if (nested) {
      retMat <- retMat[, 1:4]
    }
    retMat <- rbind(retMat, linModMat)
  }
  if(powreg){ #fit power model
    for(pp in pow){
      drcData <- as.data.frame(object$data[, c(2, 1)])
      names(drcData) <- c("yVec", "xVec")
      mod_pow = lm(yVec ~ I(xVec^pp), data = drcData)
      powModMat = matrix(c(logLik(mod_pow), icfct(mod_pow), alr3::pureErrorAnova(mod_pow)[[5]][3], (summary(mod_pow)$sigma)^2), nrow = 1, ncol = 4)
      rownames(powModMat) = paste("Power",pp,sep="")  
      colnames(powModMat) <- cnames[1:4]
      retMat <- rbind(retMat, powModMat) 
    }
  }
  if(expreg){ #fit power model
    drcData <- as.data.frame(object$data[, c(2, 1)])
    names(drcData) <- c("yVec", "xVec")
    
    #save(drcData, file = "expreg.RData")
    exp_x = exp(drcData$xVec)
    if(sum(is.infinite(exp_x))==0){
      exp_mod = lm(yVec ~ I(exp(xVec)), data = drcData)
      sm = summary(exp_mod)
      lof = 1 - pf(sm$fstatistic[1], sm$fstatistic[2], sm$fstatistic[3]) 
      expModMat = matrix(c(logLik(exp_mod), icfct(exp_mod), lof, (summary(exp_mod)$sigma)^2), nrow = 1, ncol = 4) # instead of lof pureErrorAnova(exp_mod)[[5]][3]
      rownames(expModMat) = "Exponential"  
      colnames(expModMat) <- cnames[1:4]
      retMat <- rbind(retMat, expModMat)
    }else{
      expModMat = matrix(NA,1,4)
      rownames(expModMat) = "Exponential"  
      retMat <- rbind(retMat, expModMat)
    }
    
  }
  if(hillreg){
    for(hN in hillN){
      drcData <- as.data.frame(object$data[, c(2, 1)])
      names(drcData) <- c("yVec", "xVec")
      ill_mod <- lm(yVec ~ I(xVec^hillN / (Kd + xVec^hillN)),data = drcData)
      #effect_plot(y.nls, pred = dose, interval = TRUE, plot.points = TRUE)
      illModMat = matrix(c(logLik(ill_mod), icfct(ill_mod), alr3::pureErrorAnova(ill_mod)[[5]][3], (summary(ill_mod)$sigma)^2), nrow = 1, ncol = 4)
      if(hN == 0.5) hN = "05"
      rownames(illModMat) = paste("Hill",hN,sep="")   
      colnames(illModMat) <- cnames[1:4]
      retMat <- rbind(retMat, illModMat)
    }
  }
  if (sorted != "no") {
    return(retMat[order(retMat[, sorted]), ])
  }
  else {
    return(retMat)
  }
}

# library(nlme)
# library(alr3)
# library(mixtox)
# 
# Richard<-function(x,Em,n,b,s){Em/((1+10^(n*(b-x)))^s)}
# Gompertz<-function (x,Em,n,i){Em*(exp(-exp(ln(10)*n*(i-x))))}
# HillModif<-function(x,Em,b,p){Em/((1+10^(b-x))^p)}
# 
# Hill1<-function(x, Em, n, D){Em/(1 + 10^(n * (D - x)))}
# 
# EstH.Pop<-function(DataSet){
#   InitVal<-function(DataSet){
#     xy<-sortedXyData(DataSet$dose,DataSet$expr,DataSet)
#     Em<-max(xy[c(2)])
#     D<-NLSstClosestX(xy,Em/2)
#     n<-1
#     value<-c(Em,n,D)
#     value
#   }
#   
#   hill = Hill1(dose,Em,n,D)
#   expval = DataSet$expr
#   df = data.frame(expval = expval, hill = hill)
#   DataSet.nlme<-nls(formula = expval~hill,data = df,start=c(Em=InitVal(DataSet)[c(1)],n=InitVal(DataSet)[c(2)],D=InitVal(DataSet)[c(3)]))
# }
# 
# EstH.Pop(DataSet)

