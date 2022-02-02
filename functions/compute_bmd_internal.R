#' This function compute the BMD values for all the experiments by using the compute_bmd function.
#' @param mod object from fit_models_mselect2
#' @param dataframe  pheno data tables
#' @param constantVar boolean specifying if constant or non constant variance is assumed
#' @param dose vector of doses
#' @param conf_interval confidence interval value
#' @param min_dose min dose
#' @param max_dose max dose
#' @param max_low_dos_perc_allowd percentage of variability of minumum dose allowed
#' @param max_max_dos_perc_allowd percentage of variablity of maximum dose allowed
#' @param interval_type  string specifying the type of interval to be used. Default ="delta"
#' @param loofth integer specifying the lack of fit pvalue
#' @param ratio_filter boolean specifying if filtering is applied on bmd, bmdl and bmdu ratio values
#' @param bmd_bmdl_th threshold for the bmd/bmdl ratio. Default = 20 meaning that genes whose bmd/bmdl > 20 are removed
#' @param bmdu_bmd_th threshold for the bmdu/bmd ratio. Default = 20 meaning that genes whose bmdu/bmd > 20 are removed
#' @param bmdu_bmdl_th threshold for the bmdu/bmdl ratio. Default = 40 meaning that genes whose bmu/bmdl > 40 are removed
#' @param first_only boolean value. If true only the best model (based on AIC criterion) is used to computet BMD/BMDL/BMDU/IC50 values. Otherwise, the models will be screended and the first one with lowest AIC and predicted alues satisfying all the filtering criteria will be selected
#' @param filter_bounds_bmdl_bmdu boolean specifyinig if models with bmdl=0 or bmdu = max dose should be removed
#' @return a list containinig fitted models and resume values
#' @export
compute_bmd_internal = function(mod, rl=1.349,dataframe, strictly_monotonic = TRUE,
                                first_only = TRUE,
                                constantVar, dose, conf_interval, min_dose,max_dose,
                                interval_type, loofth = 0.1,
                                max_low_dos_perc_allowd = 0, max_max_dos_perc_allowd = 0,
                                ratio_filter = FALSE, bmd_bmdl_th = 20,
                                bmdu_bmd_th = 20, bmdu_bmdl_th = 40,
                                filter_bounds_bmdl_bmdu = FALSE){
  
  exp = dataframe$exp
  dose = dataframe$dose
  
  n_models = nrow(mod$X)
  if(first_only) max_index = 1 else max_index  = n_models
  
  for(i in 1:max_index){
    opt_mod = mod$fitted_models[[i]]
    mod_name = names(mod$fitted_models)[i]
    # mod$X[mod_name,]
    
    bmd_val = tryCatch({
      if(mod_name %in% c("Linear","Quadratic","Cubic","Power2","Power3","Power4","Exponential",
                         "Hill05","Hill1","Hill2","Hill3","Hill4","Hill5")){
        
        # bmd_val = EDLin(lmObject = opt_mod,respLev = sd_level, dose=dose, ci = conf_interval)
        bmd_val = ED.lin(lmObject = opt_mod,dataframe=dataframe,constant_variance = constantVar,
                         rl = rl,ci=conf_interval,training_dose = dose, strictly_monotonic=strictly_monotonic)
        
        if(is.null(bmd_val) == FALSE){
          bmd_val = cbind(bmd_val, matrix(c(mod_name, mod$X[mod_name,3]),1,2))
          colnames(bmd_val) = c("BMD", "BMDL","BMDU","IC50/EC50","Decreasing","MOD_NAME","LOFPVal")#,"ANOVAPVal")
        }else{
          bmd_val = NULL
        }
      }else{
        monotonic_behaviour = monotonicity(fittedModel=opt_mod,dose = dose,
                                           range.length = 1000,dataframe = dataframe,
                                           strictly_monotonic = strictly_monotonic)
        
        starting_point = predict(opt_mod, data.frame(dose = 0))
        sd_level = rl *  sqrt(control_variance(model = opt_mod, constant_variance = constantVar, training_dose =dose ))
        
        
        if(monotonic_behaviour == 0){
          bmd_val = NULL
        }else{
          if(monotonic_behaviour==1){
            response_level = starting_point + sd_level
            decreasing = FALSE
          }else if (monotonic_behaviour==-1){
            response_level = starting_point - sd_level
            decreasing = TRUE
          }
          
          bmd_val = drc::ED(object = opt_mod, respLev = response_level, interval = interval_type, level = conf_interval, type = "absolute", display=FALSE)[, c(1,3,4), drop = FALSE]
          
          doseRange = seq(min(dose), max(dose), length.out = 1000)
          yFit <- stats::predict(opt_mod, newdata=data.frame(dose=doseRange))
          
          ic50 = max(yFit) - ((max(yFit) - min(yFit))/2)
          icval = drc::ED(object = opt_mod, respLev = ic50, interval = interval_type, level = conf_interval, type = "absolute")[, 1, drop = FALSE]
          
          
          bmd_val = cbind(bmd_val,icval, decreasing)
          bmd_val = cbind(bmd_val, matrix(c(mod_name, mod$X[mod_name,3]),1,2))
          colnames(bmd_val) = c("BMD", "BMDL","BMDU","IC50/EC50","Decreasing","MOD_NAME","LOFPVal")#,"ANOVAPVal")
          
        }
      }
      bmd_val
    }, error = function(e) {
      print(e)
      bmd_val = NULL
    })
    
    
    if(!is.null(bmd_val)){
      # print(bmd_val)
      
      pass_filter = inner.check.model(bmd = as.numeric(bmd_val[1,"BMD"]),
                                      bmdl= as.numeric(bmd_val[1,"BMDL"]),
                                      bmdu= as.numeric(bmd_val[1,"BMDU"]),
                                      ic50= as.numeric(bmd_val[1,"IC50/EC50"]),
                                      fitting.pval= as.numeric(bmd_val[1,"LOFPVal"]),
                                      min_dose=min_dose, max_dose=max_dose, loofth = loofth)
                                      # max_low_dos_perc_allowd = max_low_dos_perc_allowd,
                                      # max_max_dos_perc_allowd = max_max_dos_perc_allowd,
                                      # ratio_filter = ratio_filter, bmd_bmdl_th = bmd_bmdl_th,
                                      # bmdu_bmd_th = bmdu_bmd_th, bmdu_bmdl_th = bmdu_bmdl_th,
                                      # filter_bounds_bmdl_bmdu = filter_bounds_bmdl_bmdu)
      
      # print(paste("MODEL:", mod_name, "PASS FILTER:", pass_filter))
      
      if(pass_filter){
        toRetList = list("bmd_val" = bmd_val, "mod" = opt_mod)
        return(toRetList)
      }
    }
    
    
  }
  
  return(NULL)
}

inner.check.model = function(bmd, bmdl, bmdu, ic50,
                             fitting.pval,
                             min_dose,
                             max_dose,
                             loofth = 0.1,
                             max_low_dos_perc_allowd = 0,
                             max_max_dos_perc_allowd = 0,
                             ratio_filter = FALSE,
                             bmd_bmdl_th = 20,
                             bmdu_bmd_th = 20,
                             bmdu_bmdl_th = 40,
                             filter_bounds_bmdl_bmdu = FALSE){
  
  # none of the predicted values are NA
  values_are_not_na = (is.na(bmd)==FALSE) & (is.na(bmdl)==FALSE) & (is.na(bmdu)==FALSE) & (is.na(ic50)==FALSE)
  
  #values are not negative
  values_are_not_negative = (bmd>=0) & (bmdl>=0) & (bmdu>=0) & (ic50>=0)
  
  # tolerannce on the bmd minimum and maximum value
  # bmd_pass_low_dose_check = (bmd > (min_dose - (min_dose*max_low_dos_perc_allowd))) #& (bmdl > (min_dose - (min_dose*max_low_dos_perc_allowd))) & (bmdu > (min_dose - (min_dose*max_low_dos_perc_allowd)))
  # bmd_pass_max_dose_check = (bmd < (max_dose - (max_dose*max_max_dos_perc_allowd))) #& bmdl < (max_dose - (max_dose*max_max_dos_perc_allowd)) & bmdu < (max_dose - (max_dose*max_max_dos_perc_allowd))
  
  # lack of fit pvalue
  pass_loof_test = fitting.pval >= loofth
  
  # if pass is true it will be an accepted model
  pass = values_are_not_na & values_are_not_negative & pass_loof_test # & bmd_pass_low_dose_check & bmd_pass_max_dose_check
  
  # if(ratio_filter){
  #   pass_ratio_filter = (bmd/bmdl > bmd_bmdl_th) & (bmdu/bmd > bmdu_bmd_th) & (bmdu/bmdl > bmdu_bmdl_th)
  #   pass = pass & pass_ratio_filter
  # }
  # 
  # if(filter_bounds_bmdl_bmdu){
  #   bounds = bmdl>0 & bmdu< max_dose
  #   pass = pass & bounds
  #   
  # }
  
  return(pass)
  
}

#' Function to compute numerical derivative of the model's predicted values
#'
#' @param fittedModel is a model fitted to the gene expression data
#' @param dose is the numeric vecotr of doses
#' @return range.length number of interpolated points between the minimum and maximum doses. Default value is 1000
#' @keywords internal
#' @export
# return -1 if decreasing; 0 if not monotonic and 1 if increasing
monotonicity = function(fittedModel,dose,range.length = 1000, dataframe, strictly_monotonic){
  
  if(strictly_monotonic){
    step = (max(dose)-min(dose))/range.length
    x = seq(min(dose),max(dose),length.out = range.length)
    
    f1  = predict(fittedModel, newdata = data.frame(dose = x[1:(length(x)-1)]))
    f2  = predict(fittedModel, newdata = data.frame(dose = x[2:length(x)]))
    
    deriv = (f2-f1)/step
    
    if(all(deriv>0)) return(1)
    if(all(deriv<0)) return(-1)
    return(0)
  }else{
    coef_mono = coef(lm(formula = expr~dose, data = dataframe))["dose"]
    return(sign(coef_mono))
  }
  

}
