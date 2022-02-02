
#' this function select the the optimal model between the one fitted
#' @param formula formula to be fitted. Default value expr~dose
#' @param dataframe dataframe containing the data for the fitting
#' @param sel_mod_list vector of integers to specify the models to be fitted. Fossible models are "LL.2","LL.3","LL.3u","LL.4","LL.5","W1.2","W1.3","W1.4","W2.2","W2.3","W2.4","BC.4","BC.5","LL2.2","LL2.3","LL2.4","LL2.5","AR.2","AR.3","MM.2","MM.3","Linear", "Quadratic", "Cubic","Power2","Power3","Power4","Exponential","Hill05","Hill1","Hill2","Hill3","Hill4","Hill5"
#' @param Kd Kd parameter for hill models
#' @return the optimal fitted model
#' @importFrom stats lm
#' @export
fit_models_mselect2 = function(formula=expr~dose, dataframe,sel_mod_list=sel_mod_list, Kd = 10){#, hillN=2, pow = 2){
  f_list = list(drc::LL.2(),drc::LL.3(),drc::LL.3u(),drc::LL.4(),drc::LL.5(),
                drc::W1.2(),drc::W1.3(),drc::W1.4(),drc::W2.2(),drc::W2.3(),drc::W2.4(),
                drc::BC.4(),drc::BC.5(),
                drc::LL2.2(),drc::LL2.3(),drc::LL2.4(),drc::LL2.5(),
                drc::AR.2(),drc::AR.3(),
                drc::MM.2(),drc::MM.3())
  f_names = c("LL.2()","LL.3()","LL.3u()","LL.4()","LL.5()",
              "W1.2()","W1.3()","W1.4()","W2.2()","W2.3()","W2.4()",
              "BC.4()","BC.5()",
              "LL2.2()","LL2.3()","LL2.4()","LL2.5()",
              "AR.2()","AR.3()",
              "MM.2()","MM.3()")
  
  mNames_all = c("LL.2","LL.3","LL.3u","LL.4","LL.5",
                 "W1.2","W1.3","W1.4","W2.2","W2.3","W2.4",
                 "BC.4","BC.5",
                 "LL2.2","LL2.3","LL2.4","LL2.5",
                 "AR.2","AR.3",
                 "MM.2","MM.3","Linear", "Quadratic", "Cubic",
                 "Power2","Power3","Power4","Exponential",
                 "Hill05","Hill1","Hill2","Hill3","Hill4","Hill5")
  
  
  for(i in 1:length(f_list)){
    f_list[[i]]$name = mNames_all[i]
  }
  names(f_list) = f_names
  
  hillN = c(0.5,1,2,3,4,5)[which(c("Hill05","Hill1","Hill2","Hill3","Hill4","Hill5") %in% mNames_all[sel_mod_list])]
  pow = c(2,3,4)[which(c("Power2","Power3","Power4") %in% mNames_all[sel_mod_list])]
  Kd = 10
  
  drc_selected_models = sel_mod_list[sel_mod_list <= 21] # models selected that can be fitted with the drc package
  lm_selected_models = sel_mod_list[sel_mod_list > 21] # models selected that can be fitted with the drc package
  
  X_drc = NULL
  X_lm = NULL
  lm_mod_list = NULL
  drc_mod_list = NULL
  
  if(length(drc_selected_models)>0){ # if the user select any of the drc models
    selected_models = drc_selected_models
    f_list = f_list[selected_models]
    f_names = f_names[selected_models]
    mNames = mNames_all[selected_models]
    
    # gives the optimal fitted model between the drc models selectetd from the user ordered by lowest AIC
    mod.internal <- fit_models(formula = formula,dataframe = dataframe, f_list=f_list,f_names=f_names,mNames=mNames)
    
    if(is.null(mod.internal)==FALSE){
      # models to compare with the best model based on AIC
      models_to_compare = f_names[-which(f_names %in% mod.internal$mod_name)]
      if (length(models_to_compare) > 0) fctList = f_list[models_to_compare] else fctList = NULL
      
      # get AIC, lack of fit pvalue and logLik for all the models selected by the user
      # object = mod.internal$opt_mod
      # save(object,fctList, formula, dataframe, file = "tmp_mselect2.RData")
      X_drc = mselect2(object = mod.internal$opt_mod, 
                       fctList = fctList,sorted = "IC", 
                       dataframe = dataframe)

      if(is.null(fctList)){
        X_drc = matrix(X_drc, nrow = 1,dimnames = list(mod.internal$mod_name, names(X_drc)))
        rownames(X_drc) = gsub(pattern = "\\(\\)",replacement = "",x =rownames(X_drc))
      }
      
      mname = gsub(pattern = "\\(\\)",replacement = "",x = mod.internal$mod_name)
      X_drc[which(rownames(X_drc) %in% mname)[1],3] = mod.internal$loof_test$`p value`[2]
      
      toSelect = rownames(X_drc) %in% mNames
      if (sum(toSelect) == 1) X_drc = matrix(data = X_drc[toSelect,],nrow = 1,dimnames = list(rownames(X_drc)[rownames(X_drc) %in% mNames],names(X_drc[toSelect,]))) else X_drc = X_drc[toSelect,]
      
      toRem = which(is.na(X_drc[,3]))
      
      if(length(toRem) == nrow(X_drc)){ #all the lack of fit are NA
        X_drc = NULL
      } else{
        if(length(toRem)>0){
          #		if(length(toRem)< nrow(X_drc)){
          X_drc = X_drc[-toRem,,drop=FALSE]
          #		}
        }
        
        mod_name = rownames(X_drc)[1]
        opt_drc_mod <- mod.internal$opt_mod#drm(formula = formula, data=dataframe,type="continuous", fct=f_list[[which(f_names %in% paste(mod_name,"()",sep=""))]])
        drc_mod_list = mod.internal$mod_list[rownames(X_drc)]
      }
    }
    
  }
  
  if(length(lm_selected_models)>0){ # if the user select any of the lm models
    selected_models = sel_mod_list[sel_mod_list>21]
    mNames = mNames_all[selected_models]
    
    fitt_res = fit_models_lm(dataframe = dataframe, mNames=mNames, pow = pow, hillN = hillN, Kd = Kd,sorted = "IC", contData = TRUE)
    if(is.null(fitt_res)==FALSE){
      X_lm = fitt_res$retMat
      
      toRem = which(is.na(X_lm[,3]))
      
      if(length(toRem) == nrow(X_lm)){ #all the lack of fit are NA
        X_lm = NULL
      } else{
        if(length(toRem)>0){
          #		if(length(toRem)< nrow(X_drc)){
          X_lm = X_lm[-toRem,,drop=FALSE]
          #		}
        }
        
        mod_name = rownames(X_lm)[1]
        opt_linear_mod = fitt_res$FittedModels[[rownames(X_lm)[1]]]
        lm_mod_list = fitt_res$FittedModels[rownames(X_lm)]
      }
      
    }
    
  }
  
  if(is.null(X_drc) & is.null(X_lm)){
    # print("no mod fit")
    return(NULL)
  }
  
  X = rbind(X_drc, X_lm)
  if(nrow(X)>1)  X = X[order(X[, "IC"]), ]
  
  mod_name = rownames(X)[1] # name of optimal model by IC
  if(which(mNames_all %in% mod_name)<=21) mod2 = opt_drc_mod else mod2 = opt_linear_mod # assign the final model to return
  
  # list with all the fitted models ordered by AIC value
  fitted_models = c(lm_mod_list,drc_mod_list)
  fitted_models = fitted_models[rownames(X)]
  
  return(list(opt_mod = mod2, mod_name = mod_name,loof_test = X[1,3],aic = X[1,2], X = X,fitted_models=fitted_models))
  
  
}
