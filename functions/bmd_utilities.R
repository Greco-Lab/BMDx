
#' this function plots the optimal fitted model for a specific gene
#' @param model_list list of models obtained with the compute_bmd or BMD_filters function
#' @param gene_name gene name. it has to be one of the model in the list
#' @param doses numeric vector of doses

#' @return list containing a dataframe for every file sheet
#' @export

plot_model = function(model_list, gene_name, doses=c(0,5,10,20)){
  plot(model_list[[gene_name]]$opt_mod, xaxt="n")
  axis(side=1, at=doses, labels = TRUE)
}


#' this function reads all sheets of an excel file and return it as a list of dataframe
#' @param filename path of the file
#' @param tibble boolean specifying if the file has to be treated as a tibble
#' @return list containing a dataframe for every file sheet
#' @export

read_excel_allsheets <- function(filename, tibble = FALSE) {
  sheets <- readxl::excel_sheets(filename)
  x <- lapply(sheets, function(X){
    
    y = readxl::read_excel(filename, sheet = X)
    y = as.data.frame(y)
    #y[,1] = as.character(as.vector(y[[1]]))
    rownames(y) = y[,1]
    y = y[,-1]
    y
  })
  
  if(!tibble) x <- lapply(x, as.data.frame)
  names(x) <- sheets
  return(x)
  
}

#' this function run the anova analysis for every experimetn
#' @param pheno_data_list list of pheno data table
#' @param expression_list list of expression matrices
#' @param time_point_index column number of time point in the pheno data table
#' @param dose_index column number of dose in the pheno data table
#' @param sample_index column number of sample id in the pheno data table
#' @param adj.pval is a boolean parameter. If true the anova pvalue will be adjusted by fdr correction
#' @param p.th is the threshold for anova pvalues
#' @param time_point_id is a string specifying one of the time points in the dataset, or the string "All" meaning that the analysis will be run for all the time points.
#' @param nCores is the number of cores used to run the experiments. Parallelization is applied at gene level
#' @return a list containing multiple items
#' \item{list_of_filtered_expression_values}{filtered matrix of expression values containing only the genes that survive the anova test}
#' \item{list_of_variable_genes}{genes that survive the anova test}
#' \item{list_of_non_variable_genes}{genes that do not pass the anova test}
#' \item{list_of_matrices_with_anova_pvalue}{list of matrices containing the pvalue of the anova }
#' @export
#'
compute_anova_multiple_experiments = function(pheno_data_list,expression_list, time_point_index = 4, dose_index = 2, time_point_id = "All", sample_index = 1, adj.pval = TRUE, p.th=0.01, nCores = 1){
  list_of_filtered_expression_values = list()
  list_of_variable_genes = list()
  list_of_non_variable_genes = list()
  list_of_matrices_with_anova_pvalue = list()
  
  for(i in 1:length(pheno_data_list)){
    
    pTable = pheno_data_list[[i]]
    timep = unique(pTable[,time_point_index])
    
    if(time_point_id == "All"){
      for(tp in  timep){
        print(tp)
        
        pvalues_genes = compute_anova(exp_data = expression_list[[i]],
                                      pheno_data=pTable,
                                      time_t=tp,
                                      tpc = time_point_index,
                                      dc = dose_index,
                                      sc = sample_index, nCores = nCores)
        
        if(is.null(pvalues_genes)){
          print("Cannot perform anova with only one dose-level. check dose and time point columns")
          return(NULL)
        }else{
          print("Anova running...")
        }
        
        c(filtered_expression_values, variable_genes, non_variable_genes,matrix_with_anova_pvalue) %<-% filter_anova(exp_data=expression_list[[i]],
                                                                                                                     pvalues_genes=pvalues_genes,
                                                                                                                     adj.pval = adj.pval,
                                                                                                                     p.th=p.th)
        
        list_of_filtered_expression_values[[names(pheno_data_list)[i]]][[as.character(tp)]] = filtered_expression_values
        list_of_variable_genes[[names(pheno_data_list)[i]]][[as.character(tp)]] = variable_genes
        list_of_non_variable_genes[[names(pheno_data_list)[i]]][[as.character(tp)]] = non_variable_genes
        list_of_matrices_with_anova_pvalue[[names(pheno_data_list)[i]]][[as.character(tp)]] = matrix_with_anova_pvalue
      }
    }else{
      pvalues_genes = compute_anova(exp_data = expression_list[[i]],
                                    pheno_data=pTable,
                                    time_t=tp,
                                    tpc = time_point_index,
                                    dc = dose_index,
                                    sc = sample_index)
      
      c(filtered_expression_values, variable_genes, non_variable_genes,matrix_with_anova_pvalue) %<-% filter_anova(exp_data=expression_list[[i]],
                                                                                                                   pvalues_genes=pvalues_genes,
                                                                                                                   adj.pval = adj.pval,
                                                                                                                   p.th=p.th)
      
      list_of_filtered_expression_values[[names(pheno_data_list)[i]]][[as.character(input$time_point_id)]] = filtered_expression_values
      list_of_variable_genes[[names(pheno_data_list)[i]]][[as.character(input$time_point_id)]] = variable_genes
      list_of_non_variable_genes[[names(pheno_data_list)[i]]][[as.character(input$time_point_id)]] = non_variable_genes
      list_of_matrices_with_anova_pvalue[[names(pheno_data_list)[i]]][[as.character(input$time_point_id)]] = matrix_with_anova_pvalue
    }
    
    
  }
  
  
  gVars = list()
  gVars$list_of_filtered_expression_values = list_of_filtered_expression_values # EXP_FIL
  gVars$list_of_variable_genes = list_of_variable_genes
  gVars$list_of_non_variable_genes = list_of_non_variable_genes #not_var_genes
  gVars$list_of_matrices_with_anova_pvalue= list_of_matrices_with_anova_pvalue #PValMat
  
  return(gVars)
}

#' this function perform anova for every gene with respect to time and dose. If a gene does not survive to the filtering, it is not used in the BMD analysis
#' @param exp_data expression matrix
#' @param pheno_data pheno data table
#' @param time_t time point at which perform the analyis. It has to be one of the values reported in the pheno data table
#' @param tpc column number of time point in the pheno data table
#' @param dc column number of dose in the pheno data table
#' @param sc column number of sample id in the pheno data table
#' @param adj.pval is a boolean parameter. If true the anova pvalue will be adjusted by fdr correction
#' @param p.th is the threshold for anova pvalues
#' @param nCores is the number of cores used to run the experiments. Parallelization is applied at gene level
#' @return an list containing three items
#' \item{filt_exp}{filtered matrix of expression values containing only the genes that survive the anova test}
#' \item{not_var_genes}{genes that do not survive the anova test}
#' \item{var_genes}{genes passing the anova test}
#' @export

compute_anova = function(exp_data, pheno_data, time_t=24,tpc = 4, dc = 2, sc = 1, adj.pval = TRUE, p.th=0.01, nCores = 1){
  if(tpc > ncol(pheno_data))
    stop("'time point column bigger than the size of the pheno data table!")
  
  if(dc > ncol(pheno_data))
    stop("'dose column bigger than the size of the pheno data table!")
  
  if(sc > ncol(pheno_data))
    stop("'sample id column bigger than the size of the pheno data table!")
  
  if(!time_t %in% pheno_data[,tpc])
    stop("'time point not available in the pheno data table!")
  
  print("Sample ID column")
  print(sc)
  
  print("samples id --->>>>")
  print(as.character(pheno_data[,sc]))
  
  print(sum(as.character(pheno_data[,sc]) %in% colnames(exp_data)))
  
  exp_data = exp_data[,as.character(pheno_data[,sc])]
  exp_data = as.matrix(exp_data)
  
  df_timei = pheno_data[which(pheno_data[,tpc] %in% time_t),]
  
  #for each gene in the dataset compute the ANOVA across the different doses
  #pvalues_genes = c()
  
  cl <- makeCluster(nCores)
  registerDoParallel(cl)
  
  pvalues_genes = foreach(i = 1:nrow(exp_data), .combine=c) %dopar% {
    #take the value of the genes across the samples at time 1d
    exp = exp_data[i,as.character(df_timei[,sc])]
    
    #the datafame contains the expression values for genes i for samples at time 1d and their doses
    anova_df = data.frame(exp = as.numeric(exp),dose=df_timei[,dc])
    
    unique_ds = length(unique(anova_df$dose))
    
    if(unique_ds==1){
      return(NULL)
    }
    
    anova_res = aov(exp~dose, anova_df)
    pvalue = unlist(summary(anova_res))["Pr(>F)1"]
    #pvalues_genes = c(pvalues_genes,pvalue)
    pvalue
    # Increment the progress bar, and update the detail text.
    #incProgress(1/nrow(exp_data), detail = paste("Anova Gene", i))
    
  }
  
  stopCluster(cl)
  
  
  
  if(is.null(pvalues_genes)){
    return(NULL)
  }else{
    names(pvalues_genes) = rownames(exp_data)
    return(pvalues_genes)
  }
  
  
}

filter_anova = function(exp_data, pvalues_genes, adj.pval = TRUE, p.th=0.01){
  
  if(adj.pval){
    pvalues_genes = p.adjust(pvalues_genes,method = "fdr")
  }
  
  idx = which(pvalues_genes<=p.th)
  #print(length(idx))
  
  if(length(idx)==0){
    print("No genes passed the anova filtering")
    return(list(filt_exp=NULL, not_var_genes=NULL,var_genes = NULL,PValMat=NULL))
  }
  
  var_genes = pvalues_genes[idx]
  not_var_genes = pvalues_genes[-idx]
  
  if(length(idx)==1){ # if its only one gene filt_exp still have to be a matrix
    filt_exp = matrix(0, nrow = 1, ncol = ncol(exp_data))
    filt_exp[1,] = as.numeric(exp_data[names(var_genes),])
    rownames(filt_exp) = names(var_genes)
    colnames(filt_exp) = colnames(exp_data)
    #filt_exp = exp_data[names(var_genes),]
    
  }else{
    filt_exp = exp_data[names(var_genes),]
  }
  
  
  PValMat = cbind(names(pvalues_genes), pvalues_genes)
  rownames(PValMat) = NULL
  PValMat = as.data.frame(PValMat)
  PValMat[,2] = as.numeric(as.vector(PValMat[,2]))
  
  if(adj.pval)
    colnames(PValMat) = c("Gene","adj.pvalue")
  else
    colnames(PValMat) = c("Gene","pvalue")
  
  PValMat[,2] = round(PValMat[,2],4)
  toRet = list(filt_exp=filt_exp, not_var_genes=not_var_genes,var_genes = var_genes,PValMat=PValMat)
  
}


#' This function compute the BMD values for all the experiments by using the compute_bmd function.
#' @param filtered_expression_data_list list of expression matrices
#' @param pheno_data_list list of pheno data tables
#' @param time_point_index column number of time point in the pheno data table
#' @param dose_index column number of dose in the pheno data table
#' @param sample_index column number of sample id in the pheno data table
#' @param adj.pval is a boolean parameter. If true the anova pvalue will be adjusted by fdr correction
#' @param p.th is the threshold for anova pvalues

#' @return an list containing three items
#' \item{filt_exp}{filtered matrix of expression values containing only the genes that survive the anova test}
#' \item{not_var_genes}{genes that do not survive the anova test}
#' \item{var_genes}{genes passing the anova test}

#' @export
run_bmd_multiple_experiment = function(filtered_expression_data_list,
                                       pheno_data_list,
                                       interval_type = "delta",
                                       time_point_index = 4,
                                       dose_index = 2, sample_index = 1,
                                       sel_mod_list = c(19,21,22,23,25,27), rl = 1.349, lack_of_fit_pvalue= 0.1, nCores = 2){
  
  list_of_bmd_fitted_models_and_resume_table = list()
  list_of_bmd_fitted_models_and_resume_table_filtered = list()
  list_of_bmd_resume_table_filtered = list()
  
  for(j in 1:length(pheno_data_list)){
    
    pTable = pheno_data_list[[j]]
    timep = as.numeric(names(filtered_expression_data_list[[j]]))
    
    maxDose = max(as.numeric(unique(pheno_data_list[[j]][,dose_index])))
    minDose = min(as.numeric(unique(pheno_data_list[[j]][,dose_index])))
    
    for(i in timep){
      
      list_of_bmd_fitted_models_and_resume_table[[names(filtered_expression_data_list)[j]]][[as.character(i)]]  = compute_bmd(exp_data=filtered_expression_data_list[[names(pheno_data_list)[j]]][[as.character(i)]],
                                                                                                                              pheno_data=pheno_data_list[[names(pheno_data_list)[j]]],
                                                                                                                              time_t=as.character(i), #interval_type = input$Interval,
                                                                                                                              tpc = time_point_index,
                                                                                                                              dc = dose_index,
                                                                                                                              sc = sample_index, sel_mod_list = sel_mod_list,rl = rl,nCores=nCores)
      
      list_of_bmd_fitted_models_and_resume_table_filtered[[names(pheno_data_list)[j]]][[as.character(i)]]   = BMD_filters(BMDRes = list_of_bmd_fitted_models_and_resume_table[[names(filtered_expression_data_list)[j]]][[as.character(i)]],max_dose = as.numeric(maxDose), min_dose = as.numeric(minDose),loofth = lack_of_fit_pvalue)
      list_of_bmd_resume_table_filtered[[names(pheno_data_list)[j]]][[as.character(i)]]  =  list_of_bmd_fitted_models_and_resume_table_filtered[[names(pheno_data_list)[j]]][[as.character(i)]]$BMDValues_filtered
    }
  }
  
  gVars = list()
  gVars$list_of_bmd_fitted_models = list_of_bmd_fitted_models_and_resume_table
  gVars$list_of_bmd_fitted_models_and_resume_table_filtered = list_of_bmd_fitted_models_and_resume_table_filtered
  gVars$list_of_bmd_resume_table_filtered = list_of_bmd_resume_table_filtered
  return(gVars)
  
}

#' Function to combine the results of the parallel random split method
#'
#' Given the results of the  random splits method this function computes the active size for each lambda
#' @param x is the output iniziatilization, if present in the .init parameter of foreach
#' @param ... values returned by the foreach function
#' @return a list of combined objects
#' @keywords internal
#' @examples
#'

listComb = function(x, ...) {
  # create a list with all the first results of the foreach iteration
  y1 = lapply(list(...), function(y) y[[1]])
  
  print("in comb")
  print(y1)
  
  # create a matrix from this list
  x1 = matrix(unlist(y1),byrow = T,ncol = length(y1[[1]]))
  #print(x1)
  
  if(length(y1)>1){
    y2 = lapply(list(...), function(y) y[[2]])
    print(y2)
    list(x1,y2)
  }else{
    list(x1,c())
  }
  
  
}


#' This function compute the BMD value for all the genes that have an ANOVA pvalue (across the doses) less than 0.05.
#' It fits 21 different models and identify the optimal one by using the AIC creterion. The optimal one is the one with minimum AIC
#' The BMD and BMDL values associated to this model are reported along with the Lack of fitness pvalue calculated with the modelFit function
#' @param exp_data expression matrix
#' @param pheno_data pheno data table
#' @param time_t time point at which perform the analyis. It has to be one of the values reported in the pheno data table
#' @param tpc column number of time point in the pheno data table
#' @param dc column number of dose in the pheno data table
#' @param sc column number of sample id in the pheno data table
#' @param sel_mod_list
#' @param rl
#' @param nCores the number of cores used to run the analysis
#' @return an list containing three items
#' \item{filt_exp}{filtered matrix of expression values containing only the genes that survive the anova test}
#' \item{not_var_genes}{genes that do not survive the anova test}
#' \item{var_genes}{genes passing the anova test}

#' @export

compute_bmd = function(exp_data,pheno_data,time_t=4,interval_type = "delta",tpc = 4, dc = 2, sc = 1,
                       sel_mod_list = c(19,21,22,23,25,27), rl = 1.349, constantVar = FALSE, nCores=2){
  #,Kd = 10, hillN=2, pow = 2){
  if(tpc > ncol(pheno_data))
    stop("'time point column bigger than the size of the pheno data table!")
  
  if(dc > ncol(pheno_data))
    stop("'dose column bigger than the size of the pheno data table!")
  
  if(sc > ncol(pheno_data))
    stop("'sample id column bigger than the size of the pheno data table!")
  
  if(!time_t %in% pheno_data[,tpc])
    stop("'time point not available in the pheno data table!")
  
  df_timei = pheno_data[which(pheno_data[,tpc] %in% time_t),]
  
  #withProgress(message = 'Computing BMD', detail = paste("BMD Time Point: ",time_t, " Gene:  0/",nrow(exp_data),sep=""), value = 1, {
  # BMDValues = c()
  # opt_models_list = list()
  cl <- makeCluster(nCores)
  registerDoParallel(cl)
  
  #listComb,.multicombine=TRUE, .init=list(c(), c())
  res <- foreach(i=1:nrow(exp_data),.combine="c", 
                 .export = c("fit_models_mselect","fit_models","modelFit","mselect2","ED.lin"),
                 .errorhandling="stop") %dopar% {
                   
                   #for(i in 1:nrow(exp_data)){
                   
                   #BMD needs to be evaluated on the whole data, the models differs based on how they compute the mean
                   exp = exp_data[i,as.character(df_timei[,sc])]
                   exp = as.numeric(exp)
                   dose = df_timei[,dc]
                   dose = as.numeric(as.vector(dose))
                   
                   if(constantVar){
                     sd_level = sd(exp) * rl
                   }else{
                     sd_level = sd(exp[which(dose==0)]) * rl
                   }
                   
                   # print("Response level, mean of the controls * 1.349 ---------------->>>>>>>>>>")
                   # print(sd_level)
                   
                   df_gi = data.frame(dose=dose,expr=exp)
                   df_gi$dose = as.numeric(as.vector(df_gi$dose))
                   
                   mod = fit_models_mselect(formula=expr~dose, dataframe=df_gi,sel_mod_list=sel_mod_list,Kd = 10)#, hillN=hillN, pow = pow)
                   
                   if(is.null(mod)){
                     print("no mod fit this gene")
                   }else{
                     #opt_models_list[[rownames(exp_data)[i]]] = mod
                     print("Optimal Model -------->>>>>>")
                     print(mod$mod_name)
                     bmd_val = tryCatch({
                       if(mod$mod_name %in% c("Linear","Quadratic","Cubic","Power2","Power3","Power4","Exponential",
                                              "Hill05","Hill1","Hill2","Hill3","Hill4","Hill5")){
                         bmd_val = ED.lin(lmObject = mod$opt_mod,respLev = sd_level, dose=df_gi$dose)
                       }else{
                         mmin = min(mod$opt_mod$predres[which(dose==min(dose)),1])
                         mmax = max(mod$opt_mod$predres[which(dose==max(dose)),1])
                         if(mmax>mmin){
                           response_level = mean(exp[which(dose==0)]) + sd_level
                           decreasing = FALSE
                         }else{
                           response_level = mean(exp[which(dose==0)]) - sd_level
                           decreasing = TRUE
                         }
                         bmd_val = ED(object = mod$opt_mod, respLev = response_level, interval = interval_type, level = 0.05, type = "relative")[, c(1,3,4), drop = FALSE]
                         
                         ic50 = max(df_gi$expr) - ((max(df_gi$expr) - min(df_gi$expr))/2)
                         icval = ED(object = mod$opt_mod, respLev = ic50, interval = interval_type, level = 0.05, type = "relative")[, 1, drop = FALSE]
                         bmd_val = cbind(bmd_val,icval, decreasing)
                         
                       }
                       colnames(bmd_val) = c("BMD","BMDL", "BMDU","IC50/EC50","Decreasing")
                       bmd_val
                     }, error = function(e) {
                       cat("Error\n")
                       return(NULL)
                     })
                     
                     #if(!is.null(bmd_val))BMDValues = rbind(BMDValues,c(rownames(exp_data)[i],bmd_val,mod$mod_name,mod$loof_test))
                     if(!is.null(bmd_val)) {
                       bmd_values = c(rownames(exp_data)[i],bmd_val,mod$mod_name,mod$loof_test)
                       #return(list(bmd_values = bmd_values, mod = mod))
                       toRet = list()
                       toRet[[paste(rownames(exp_data)[i],"BMDValues",sep="_")]] = bmd_values
                       toRet[[paste(rownames(exp_data)[i],"Mod",sep="_")]] = mod
                       return(toRet)
                     }
                   }
                   
                   #incProgress(1/nrow(exp_data), detail = paste("BMD Time Point: ",time_t," Gene: ",i ,"/",nrow(exp_data),sep=""))
                 }
  #})
  stopCluster(cl)
  
  print("xxx")
  if(length(res)==0){
    return(list(BMDValues = NULL,opt_models_list=NULL))
  }
  
  toRem = which(names(res)=="")
  if(length(toRem)>0){
    res = res[-toRem]
  }
  
  BMDValues = c()
  opt_models_list = list()
  idx = seq(from = 1,to = length(res),by = 2)
  # indici = 1
  # 
  # for(ii in idx){
  #   BMDValues = rbind(BMDValues, res[[ii]])
  #   opt_models_list[[indici]] = res[[ii+1]]
  #   indici = indici + 1
  # }
  
  values_idx = grep(x = names(res),pattern = "BMDValues")
  mod_idx = grep(x = names(res),pattern = "Mod")
  
  genigsub = gsub(pattern = "_BMDValues",replacement = "",x = names(res))
  
  for(ii in idx){
    print(ii)
    BMDValues = rbind(BMDValues, res[[ii]])
    opt_models_list[[res[[ii]][1]]] = res[[ii+1]]
  }
  
  # close(pb)
  #BMDValues = res[[1]]
  #opt_models_list = res[[2]]
  BMDValues = cbind(BMDValues,rownames(exp_data)[BMDValues[,1]])
  
  if(length(BMDValues)==0){
    return(list(BMDValues = NULL,opt_models_list=NULL))
  }
  
  
  colnames(BMDValues) = c("Gene","BMD", "BMDL","BMDU","IC50/EC50","Decreasing","MOD_NAME","LOFPVal","ANOVAPVal")
  BMDValues = as.data.frame(BMDValues)
  BMDValues$BMD = round(as.numeric(as.vector(BMDValues$BMD)),4)
  BMDValues$BMDL = round(as.numeric(as.vector(BMDValues$BMDL)),4)
  BMDValues$LOFPVal = round(as.numeric(as.vector(BMDValues$LOFPVal)),4)
  BMDValues$ANOVAPVal = round(as.numeric(as.vector(BMDValues$ANOVAPVal)),4)
  BMDValues = BMDValues[, 1:8]
  
  #names(opt_models_list) = BMDValues[,1]
  
  return(list(BMDValues = BMDValues,opt_models_list=opt_models_list))
}



#this function select the the optimal model between the one fitted
library(ggplot2)
fit_models_mselect = function(formula=expr~dose, dataframe=df_gi,sel_mod_list=sel_mod_list, Kd = 10){#, hillN=2, pow = 2){
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
  
  mNames = c("LL.2","LL.3","LL.3u","LL.4","LL.5",
             "W1.2","W1.3","W1.4","W2.2","W2.3","W2.4",
             "BC.4","BC.5",
             "LL2.2","LL2.3","LL2.4","LL2.5",
             "AR.2","AR.3",
             "MM.2","MM.3","Linear", "Quadratic", "Cubic",
             "Power2","Power3","Power4","Exponential",
             "Hill05","Hill1","Hill2","Hill3","Hill4","Hill5")
  
  hillN = c(0.5,1,2,3,4,5)[which(c("Hill05","Hill1","Hill2","Hill3","Hill4","Hill5") %in% mNames[sel_mod_list])]
  pow = c(2,3,4)[which(c("Power2","Power3","Power4") %in% mNames[sel_mod_list])]
  
  selected_models = sel_mod_list[sel_mod_list<=21] # models selected that can be fitted with the drc package
  
  if(length(selected_models)==0){ # if the user do not select any model from the drc package i fit one of them so I can run the mselect2 function
    selected_models = 19
  }
  
  f_list = f_list[selected_models]
  f_names = f_names[selected_models]
  mNames = mNames[sel_mod_list]
  mod.internal <- fit_models(formula = formula,dataframe = dataframe, f_list=f_list,f_names=f_names)
  
  if(is.null(mod.internal)){
    cat("Mod is null\n")
    return(NULL)
  }else{
    
    X = mselect2(object = mod.internal$opt_mod, fctList = f_list,sorted = "IC",
                 linreg= TRUE, powreg = TRUE, pow = pow, expreg = TRUE, hillreg = TRUE, hillN = hillN, Kd = Kd)
    
    mname = gsub(pattern = "\\(\\)",replacement = "",x = mod.internal$mod_name)
    X[which(rownames(X) %in% mname)[1],3] = mod.internal$loof_test$`p value`[2]
    
    toSelect = rownames(X) %in% mNames
    if(sum(toSelect) == 1){
      ss = matrix(data = X[toSelect,],nrow = 1,dimnames = list(rownames(X)[rownames(X) %in% mNames],names(X[toSelect,])))
      X = ss
    }else{
      X = X[toSelect,]
    }
    
    
    toRem = which(is.na(X[,3]))
    
    if(length(toRem) == nrow(X)){ #all the lack of fit are NA
      return(NULL)
      #toRem = c()
      
    }
    
    if(length(toRem)>0){
      if(length(toRem)< nrow(X)){
        X = X[-toRem,,drop=FALSE]
      }
    }
    
    mod_name = rownames(X)[1]
    
    # print("OPtimal model -------------------->")
    # print(mod_name)
    if(mod_name %in% "Linear"){
      mod2 = lm(formula,dataframe)
      #effect_plot(mod2, pred = dose, interval = TRUE, plot.points = TRUE)
    }
    else{
      if(mod_name %in% "Quadratic"){
        mod2 = lm(expr ~ dose + I(dose * dose), data = dataframe)
        #effect_plot(mod2, pred = dose, interval = TRUE, plot.points = TRUE)
        
      }else{
        if(mod_name %in% "Cubic"){
          mod2 = lm(expr~dose + I(dose * dose) +  I(dose * dose * dose), data = dataframe)
        }else{
          if(mod_name %in% c("Power2","Power3","Power4")){
            sp = as.numeric(gsub(pattern = "Power",replacement = "",x = mod_name))
            print("Power mod -------------->>>>>>>>>>>>>>>>>>>>")
            print(sp)
            mod2 = lm(expr ~ I(dose^sp), data = dataframe)
            #effect_plot(mod2, pred = dose, interval = TRUE, plot.points = TRUE)
            
          }else{
            if(mod_name %in% "Exponential"){
              mod2 = lm(expr ~ I(exp(dose)), data = dataframe)
              #effect_plot(mod2, pred = dose, interval = TRUE, plot.points = TRUE)
              
            }else{
              if(mod_name %in% c("Hill05","Hill1","Hill2","Hill3","Hill4","Hill5")){
                hlN = gsub(pattern = "Hill",replacement = "",x = mod_name)
                if(hlN == "05") hlN = 0.05
                hlN = as.numeric(hlN)
                print("Hill mod -------------->>>>>>>>>>>>>>>>>>>>")
                print(hlN)
                mod2 <- lm(expr ~ I(dose^hlN / (Kd + dose^hlN)),data = dataframe)
                #mod2 <- lm(expr ~ I(dose^n / (Kd + dose^n)),data = dataframe)
                
                #effect_plot(mod2, pred = dose, interval = TRUE, plot.points = TRUE)
                #c(logLik(y.nls), icfct(y.nls), pureErrorAnova(y.nls)[[5]][3], (summary(y.nls)$sigma)^2)
                #plot(y.nls)
                
              }else{
                
                print("OPtimal model for drm -------------------->")
                print(paste(mod_name,"()",sep=""))
                mod2 <- mod.internal$opt_mod#drm(formula = formula, data=dataframe,type="continuous", fct=f_list[[which(f_names %in% paste(mod_name,"()",sep=""))]])
              }
            }
          }
        }
      }
    }
    
    return(list(opt_mod = mod2, mod_name = mod_name,loof_test = X[1,3],aic = X[1,2], X = X))
    
  }
  
}

#This function fit the models available as the fct parameter of function drm in the drc package
fit_models = function(formula=expr~dose, dataframe=df_gi,f_list,f_names){
  
  # f_list = list(LL.2(),LL.3(),LL.3u(),LL.4(),LL.5(),
  #               W1.2(),W1.3(),W1.4(),W2.2(),W2.3(),W2.4(),
  #               BC.4(),BC.5(),
  #               LL2.2(),LL2.3(),LL2.4(),LL2.5(),
  #               AR.2(),AR.3(),
  #               MM.2(),MM.3())
  # f_names = c("LL.2()","LL.3()","LL.3u()","LL.4()","LL.5()",
  #             "W1.2()","W1.3()","W1.4()","W2.2()","W2.3()","W2.4()",
  #             "BC.4()","BC.5()",
  #             "LL2.2()","LL2.3()","LL2.4()","LL2.5()",
  #             "AR.2()","AR.3()",
  #             "MM.2()","MM.3()")
  #
  mod_list = list()
  for(i in 1:length(f_list)){
    mod_list[[i]] = tryCatch({
      mod.internal.internal <- drc::drm(formula = formula, data=dataframe,type="continuous", fct=f_list[[i]])
    }, error = function(e) {
      cat("Error\n")
      return(NULL)
    })
  }
  
  if(length(mod_list)==0){
    print("No mod did the fitting\n")
    return(NULL)
  }
  
  # NB: lowest AIC are preferred! remember, they can be negative! https://stats.stackexchange.com/questions/84076/negative-values-for-aic-in-general-mixed-model
  AIC_val = c()
  for(i in 1:length(mod_list)){
    mod = mod_list[[i]]
    if(is.null(mod)){
      AIC_val = c(AIC_val,NA)
    }else{
      AIC_val = c(AIC_val,AIC(mod))
    }
  }
  
  opt_mod = mod_list[[which.min(AIC_val)]]
  #We compute the lack of fit test. In this case we can filter model with a pvalue greater that 0.05
  x = drc::modelFit(opt_mod)
  
  return(list(opt_mod = opt_mod,mod_name = f_names[which.min(AIC_val)],loof_test = x))
}


# This function filters the gene resulting from the compute_bmd function.
# In particular the genes with BMD value greter than the maximum dose and the genes with a pvalue smaller that 0.1 are removed

BMD_filters = function(BMDRes,max_dose = 20, min_dose, max_low_dos_perc_allowd = 0, max_max_dos_perc_allowd = 0, loofth = 0.1){
  BMDValues = BMDRes$BMDValues
  BMDModels = BMDRes$opt_models_list
  
  bmd_to_rem_low_dose = which(as.numeric(BMDValues[,"BMD"])<= (min_dose - (min_dose*0.1)))
  
  bmd_to_rem = which(as.numeric(BMDValues[,"BMD"])>=max_dose)
  loof_to_rem = which(as.numeric(BMDValues[,"LOFPVal"])<=loofth)
  
  bmd_na = union(which(is.na(BMDValues[,"BMD"])),
                 union(which(is.na(BMDValues[,"BMDL"])),
                       which(is.na(BMDValues[,"BMDU"]))))
  
  bmd_neg = union(which(as.numeric(BMDValues[,"BMD"])<0),
                  union(which(as.numeric(BMDValues[,"BMDL"])<0),
                        which(as.numeric(BMDValues[,"BMDU"])<0)))
  
  bounds = union(which(as.numeric(BMDValues[,"BMDL"])<=0), which(as.numeric(BMDValues[,"BMDU"])>= max_dose))
  
  # print("BMD to rem  -->> ")
  # print(bmd_to_rem)
  # print("LOOF to rem  -->> ")
  # print(loof_to_rem)
  
  to_rem = union(union(bmd_to_rem,loof_to_rem), union(bmd_na, bmd_neg))
  to_rem = union(to_rem, bmd_to_rem_low_dose)
  to_rem = union(to_rem, bounds)
  # print("object to rem: ---------------->>>>>>>>")
  # print(to_rem)
  
  if(length(to_rem)>0){
    BMDValues_filtered = BMDValues[-to_rem,]
    BMDModels_filtered = BMDModels[-to_rem]
  }else{
    BMDValues_filtered = BMDValues
    BMDModels_filtered = BMDModels
  }
  
  #BMDModels_filtered = BMDModels[-c(BMDValues[toRem,1])]
  # print("Check that for every gene i have a model")
  # print(sum(BMDValues_filtered[,1] %in% names(BMDModels_filtered)))
  # print("dim(BMDValues -->")
  # print(dim(BMDValues_filtered))
  # print("length(BMDModels)")
  # print(length(BMDModels_filtered))
  
  return(list(BMDValues_filtered=BMDValues_filtered,BMDModels_filtered=BMDModels_filtered))
}


#
# pathway_enrichment_kegga = function(X){
#   library(org.Mm.eg.db)
#   xx <- as.list(org.Mm.egALIAS2EG)
#   # Remove pathway identifiers that do not map to any entrez gene id
#   xx <- xx[!is.na(xx)]
#   entrez_g = unlist(xx[X])
#
#   library(clusterProfiler)
#   kk = kegga(de=entrez_g,  species = "Mm")
#
#   kk = as.data.frame(kk)
#   return(kk)
# }
#
# pathway_enrichment = function(X){
#   library(org.Mm.eg.db)
#   xx <- as.list(org.Mm.egALIAS2EG)
#   # Remove pathway identifiers that do not map to any entrez gene id
#   xx <- xx[!is.na(xx)]
#   entrez_g = unlist(xx[X])
#
#   library(clusterProfiler)
#   kk <- enrichKEGG(gene          = entrez_g,
#                    pAdjustMethod ="none",
#                    organism      = 'mmu',
#                    qvalueCutoff  = 0.1,
#                    pvalueCutoff  = 0.05)
#   kk = as.data.frame(kk)
#   return(kk)
# }
#
# GO_enrichment_goana = function(X){
#   library(org.Mm.eg.db)
#   xx <- as.list(org.Mm.egALIAS2EG)
#   # Remove pathway identifiers that do not map to any entrez gene id
#   xx <- xx[!is.na(xx)]
#   entrez_g = unlist(xx[X])
#
#   library(clusterProfiler)
#   kk <- goana(de         = entrez_g,
#               species = "Mm")
#   kk = as.data.frame(kk)
#   return(kk)
# }
#
# GO_enrichment = function(X){
#   library(org.Mm.eg.db)
#   xx <- as.list(org.Mm.egALIAS2EG)
#   # Remove pathway identifiers that do not map to any entrez gene id
#   xx <- xx[!is.na(xx)]
#   entrez_g = unlist(xx[X])
#
#   library(clusterProfiler)
#   kk <- enrichGO(gene         = entrez_g,
#                  pAdjustMethod="fdr",
#                  OrgDb =  'org.Mm.eg.db',
#                  pvalueCutoff = 0.05)
#   kk = as.data.frame(kk)
#   return(kk)
# }
#
#
# venn_diagram_plot = function(genes1d,genes3d,genes28d){
#   library(VennDiagram)
#
#   genes1d3d = intersect(genes1d,genes3d)
#   genes1d28d = intersect(genes1d,genes28d)
#   genes3d28d = intersect(genes3d,genes28d)
#   genes1d3d28d = intersect(genes1d,intersect(genes3d,genes28d))
#
#   plot(0,xaxt='n',yaxt='n',bty='n',pch='',ylab='',xlab='')
#   draw.triple.venn(area1 = length(genes1d),area2 = length(genes3d),area3 = length(genes28d),
#                    n12 = length(genes1d3d),n23 = length(genes3d28d),n13 = length(genes1d28d),n123 = length(genes1d3d28d),
#                    category = c("1d", "3d", "28d"),col=rainbow(3))
#
# }
#
# barplot_days = function(BMDValues_1d_filtered,BMDValues_3d_filtered,BMDValues_28d_filtered){
#   df = data.frame(nGenes = c(nrow(BMDValues_1d_filtered),nrow(BMDValues_3d_filtered),nrow(BMDValues_28d_filtered)),
#                   Day = c("1d","3d","28d"))
#   library(ggplot2)
#   p<-ggplot(data=df, aes(x=Day, y=nGenes)) +
#     geom_bar(stat="identity")
#
#   p
#
# }
#
# plot_boxplot = function(BMDValues_1d_filtered,BMDValues_3d_filtered,BMDValues_28d_filtered){
#
#   par(mfrow=c(1,2))
#   boxplot(as.numeric(BMDValues_1d_filtered[,"BMD"]),
#           as.numeric(BMDValues_3d_filtered[,"BMD"]),
#           as.numeric(BMDValues_28d_filtered[,"BMD"]),ylab = "BMDValues")
#
#   boxplot(as.numeric(BMDValues_1d_filtered[,"BMDL"]),
#           as.numeric(BMDValues_3d_filtered[,"BMDL"]),
#           as.numeric(BMDValues_28d_filtered[,"BMDL"]),ylab = "BMDLValues")
#
# }
#
# find_genes_associated_to_time_and_cell = function(exp_data,pheno_data){
#   cell_list = list()
#   for(cell_name in c("Macro","Lympho","Neutro","Eosino","Epithel")){
#     pvalues_genes = c()
#     pb = txtProgressBar(min = 1,max=nrow(exp_data),style=3)
#     for(i in 1:nrow(exp_data)){
#       #take the value of the genes across the samples at time 1d
#       exp = as.numeric(as.vector(exp_data[i,]))
#       #the datafame contains the expression values for genes i for samples at time 1d and their doses
#
#       anova_df = data.frame(exp = exp,balCount=pheno_data[,cell_name], time = pheno_data[,"Day"],dose=pheno_data[,"Dose"])
#       anova_df$time = gsub(pattern = "d",replacement = "",x = anova_df$time)
#       anova_df$time = as.numeric(as.vector(anova_df$time))
#       anova_df$dose = gsub(pattern = "ug",replacement = "",x = anova_df$dose)
#       anova_df$dose = as.character(anova_df$dose)
#       anova_df$dose[anova_df$dose %in% "none"] = 0
#       anova_df$dose = as.numeric(anova_df$dose)
#       # fitting regression between gene expression, balcounts and time,
#       # having both balcount and time as continuous variables
#       #anova_res = lm(exp ~ balCount * time * dose, anova_df)
#       anova_res = lm(exp ~ balCount : time : dose, anova_df)
#       xx = summary(anova_res)$coefficients
#       pvalues = xx[,"Pr(>|t|)"]
#
#       pvalues_genes = rbind(pvalues_genes,pvalues)
#       setTxtProgressBar(pb,i)
#     }
#     close(pb)
#     rownames(pvalues_genes) = rownames(exp_data)
#     #  pvalues_genes = apply(pvalues_genes,2,p.adjust,method="fdr")
#     cell_list[[cell_name]]= pvalues_genes
#   }
#   return(cell_list)
# }
#
# find_genes_associated_to_time_and_cell_limma_each_time_each_dose = function(exp_data,pheno_data,cell_type = c("Macro","Eosino","Lympho","Neutro","Epithel")){
#   cell_list = list()
#   for(cell_name in cell_type){
#     print(cell_name)
#     pvalues_genes = c()
#     pb = txtProgressBar(min = 1,max=nrow(exp_data),style=3)
#     for(i in 1:nrow(exp_data)){
#       #take the value of the genes across the samples at time 1d
#       exp = as.numeric(as.vector(exp_data[i,]))
#       #the datafame contains the expression values for genes i for samples at time 1d and their doses
#
#       anova_df = data.frame(exp = exp,balCount=pheno_data[,cell_name])#, time = pheno_data[,"Day"],dose=pheno_data[,"Dose"])
#       #anova_df$time = gsub(pattern = "d",replacement = "",x = anova_df$time) #ulla data
#       #anova_df$time = gsub(pattern = "Day ",replacement = "",x = anova_df$time) #carole data
#
#       #anova_df$time = as.numeric(as.vector(anova_df$time))
#       #anova_df$dose = gsub(pattern = "ug",replacement = "",x = anova_df$dose) #ulla data
#       #anova_df$dose = as.character(anova_df$dose) #ulla data
#       #anova_df$dose[anova_df$dose %in% "none"] = 0 #ulla data
#
#       #anova_df$dose = as.numeric(anova_df$dose)
#       # fitting regression between gene expression, balcounts and time,
#       # having both balcount and time as continuous variables
#       #anova_res = lm(exp ~ balCount * time * dose, anova_df)
#       #scatterplot(exp ~  balCount + dose + time , anova_df)
#       if(cell_name %in% c("Eosino","Lympho")){k = 1} else{ k=3}
#       library(splines)
#       X_neutro =anova_df$balCount#ns(anova_df$balCount, df=k)
#       X_neutro[is.na(X_neutro)] = 0
#       design_neutro = model.matrix(~X_neutro, anova_df)
#       #design_neutro = cbind(rep(1,nrow(X_neutro)),X_neutro,anova_df$time,anova_df$dose)
#       fit_neutro = lmFit(anova_df$exp, design_neutro)
#       fit2_neutro = eBayes(fit_neutro)
#       sig.genes.neutro = topTable(fit2_neutro)
#       #sig.genes.neutro = topTable(fit2_neutro, number=100000, p.value=.05, coef="Neutro")
#
#       pvalues = sig.genes.neutro#[,"adj.P.Val"]
#
#       pvalues_genes = rbind(pvalues_genes,pvalues)
#       setTxtProgressBar(pb,i)
#     }
#     close(pb)
#     rownames(pvalues_genes) = rownames(exp_data)
#     #  pvalues_genes = apply(pvalues_genes,2,p.adjust,method="fdr")
#
#     SMat = sign(pvalues_genes[,1:k])
#     if(k==1){
#       Sign = SMat
#     }else{
#       Sign = apply(SMat,MARGIN = 1,FUN = function(row){
#         c(-1,1)[which.max(c(sum(row==-1) ,sum(row==1)))]
#       })
#     }
#
#     pvalues_genes = cbind(pvalues_genes,Sign)
#     cell_list[[cell_name]]= pvalues_genes
#   }
#   return(cell_list)
# }
#
#
# find_genes_associated_to_time_and_cell_limma = function(exp_data,pheno_data,cell_type = c("Macro","Eosino","Lympho","Neutro","Epithel")){
#   cell_list = list()
#   for(cell_name in cell_type){
#     print(cell_name)
#     pvalues_genes = c()
#     pb = txtProgressBar(min = 1,max=nrow(exp_data),style=3)
#     for(i in 1:nrow(exp_data)){
#       #take the value of the genes across the samples at time 1d
#       exp = as.numeric(as.vector(exp_data[i,]))
#       #the datafame contains the expression values for genes i for samples at time 1d and their doses
#
#       anova_df = data.frame(exp = exp,balCount=pheno_data[,cell_name], time = pheno_data[,"Day"],dose=pheno_data[,"Dose"])
#       anova_df$time = gsub(pattern = "d",replacement = "",x = anova_df$time) #ulla data
#       #anova_df$time = gsub(pattern = "Day ",replacement = "",x = anova_df$time) #carole data
#
#       #anova_df$time = as.numeric(as.vector(anova_df$time))
#       anova_df$dose = gsub(pattern = "ug",replacement = "",x = anova_df$dose) #ulla data
#       anova_df$dose = as.character(anova_df$dose) #ulla data
#       anova_df$dose[anova_df$dose %in% "none"] = 0 #ulla data
#
#       anova_df$dose = as.numeric(anova_df$dose)
#       # fitting regression between gene expression, balcounts and time,
#       # having both balcount and time as continuous variables
#       #anova_res = lm(exp ~ balCount * time * dose, anova_df)
#       #scatterplot(exp ~  balCount + dose + time , anova_df)
#       if(cell_name %in% c("Eosino","Lympho")){k = 1} else{ k=3}
#       library(splines)
#       X_neutro <- ns(anova_df$balCount, df=k)
#       X_neutro[is.na(X_neutro)] = 0
#       design_neutro = model.matrix(~X_neutro + time + dose, anova_df)
#       #design_neutro = cbind(rep(1,nrow(X_neutro)),X_neutro,anova_df$time,anova_df$dose)
#       fit_neutro = lmFit(anova_df$exp, design_neutro)
#       fit2_neutro = eBayes(fit_neutro)
#       sig.genes.neutro = topTable(fit2_neutro)
#       #sig.genes.neutro = topTable(fit2_neutro, number=100000, p.value=.05, coef="Neutro")
#
#       pvalues = sig.genes.neutro#[,"adj.P.Val"]
#
#       pvalues_genes = rbind(pvalues_genes,pvalues)
#       setTxtProgressBar(pb,i)
#     }
#     close(pb)
#     rownames(pvalues_genes) = rownames(exp_data)
#     #  pvalues_genes = apply(pvalues_genes,2,p.adjust,method="fdr")
#
#     SMat = sign(pvalues_genes[,1:k])
#     if(k==1){
#       Sign = SMat
#     }else{
#       Sign = apply(SMat,MARGIN = 1,FUN = function(row){
#         c(-1,1)[which.max(c(sum(row==-1) ,sum(row==1)))]
#       })
#     }
#
#     pvalues_genes = cbind(pvalues_genes,Sign)
#     cell_list[[cell_name]]= pvalues_genes
#   }
#   return(cell_list)
# }
#
# # compute the pathways for the genes coming from the limma model with spline decomposition of the bal counts
#
# find_pathways_from_correlated_genes_only_interactions = function(DAT){
#   Path_res = c()
#   Gene_list = list()
#   #pb = txtProgressBar(min = 2,max = length(DAT),style=3)
#   for(i in 1:length(DAT)){
#
#     M = DAT[[i]]
#     colnames(M) = "ballCount:dose:time"
#     for(j in 1:ncol(M)){
#       infoMat = cbind(rownames(M)[M[,j]<0.05],round(M[M[,j]<0.05,j],4))
#       colnames(infoMat) = c("Gene","PValue")
#
#       p_adj = p.adjust(as.numeric(infoMat[,2]),method = "fdr")
#       idX = which(p_adj<0.05)
#       if(length(idX)>0){
#         infoMat = infoMat[idX,]
#       }
#       gene_list = infoMat
#       pat = pathway_enrichment(X = rownames(gene_list))
#       if(nrow(pat)>0){
#         pat = cbind(pat, names(DAT)[i],colnames(M)[j])
#         Path_res = rbind(Path_res,pat)
#       }
#     }
#     Gene_list[[names(DAT)[i]]] = gene_list
#     #setTxtProgressBar(pb,i)
#   }
#   #close(pb)
#   colnames(Path_res)[10:11] = c("Cell","Type")
#   return(list(Path_res=Path_res,Gene_list=Gene_list))
# }
#
#
# # compute the pathways for the genes coming from the model with * in the formula
# # anova_res = lm(exp ~ balCount * time * dose, anova_df)
# find_pathways_from_correlated_genes = function(DAT){
#   Path_res = c()
#   Gene_list = list()
#   pb = txtProgressBar(min = 2,max = length(DAT),style=3)
#   for(i in 2:length(DAT)){
#     gl = list()
#     M = DAT[[i]]
#     for(j in 2:8){
#       infoMat = cbind(rownames(M)[M[,j]<0.05],round(M[M[,j]<0.05,j],4))
#       colnames(infoMat) = c("Gene","PValue")
#
#       p_adj = p.adjust(as.numeric(infoMat[,2]),method = "fdr")
#       idX = which(p_adj<0.05)
#       if(length(idX)>0){
#         infoMat = infoMat[idX,]
#       }
#       gene_list = infoMat
#       gl[[colnames(M)[j]]] = gene_list
#       pat = pathway_enrichment(X = gene_list)
#       if(nrow(pat)>0){
#         pat = cbind(pat, names(DAT)[i],colnames(M)[j])
#         Path_res = rbind(Path_res,pat)
#       }
#     }
#     Gene_list[[names(DAT)[i]]] = gl
#     setTxtProgressBar(pb,i)
#   }
#   close(pb)
#   colnames(Path_res)[10:11] = c("Cell","Type")
#   return(list(Path_res=Path_res,Gene_list=Gene_list))
# }
#
# find_correlated_genes_and_dose_dep_genes = function(DAT,BMDVal=BMDValues_1d_filtered){
#   genes = rownames(BMDVal)
#   FM = c()
#   pb = txtProgressBar(min = 2,max = length(DAT),style=3)
#   for(i in 2:length(DAT)){
#     M = DAT[[i]]
#     for(j in 1:8){
#       X = rownames(M)[M[,j]<0.05]
#       common = intersect(genes,X)
#       pat = pathway_enrichment(X =common)
#
#       if(nrow(pat)>0){
#         pat = cbind(names(DAT)[i],colnames(M)[j],length(genes),length(common),pat)
#       }
#       FM = rbind(FM,pat)
#     }
#     setTxtProgressBar(pb,i)
#   }
#   close(pb)
#   colnames(FM)[1:4] = c("Cell","Type","NGenes","NGenesDoseResp")
#   return(FM)
# }
#
# plot_countings = function(pheno_data,main_title,leg_pos){
#   idx1 = which(pheno_data$Day == "1d")
#   idx3 = which(pheno_data$Day == "3d")
#   idx28 = which(pheno_data$Day == "28d")
#
#   m1d = colMeans(pheno_data[idx1,10:14])
#   m3d = colMeans(pheno_data[idx3,10:14])
#   m28d = colMeans(pheno_data[idx28,10:14])
#
#
#
#   M = rbind(m1d,m3d,m28d)
#   matplot(M,type="l",col=rainbow(5),lwd = 3,xaxt = "n",main = main_title,ylab ="Cell counts")
#   axis(1, at=c(1, 2, 3), labels=c("1d", "3d", "28d"))
#   legend(x = leg_pos,legend = colnames(M),
#          fill = rainbow(5),cex = 0.7,ncol = 2)
#
# }
#
# #Returns the jaccard index of the
# intersection_gene_between_countings_all = function(DAT){
#   Mat = c()
#   col_names = c()
#   for(i in 2:6){
#     for(k in 2:8){
#       Mat = cbind(Mat,as.numeric(DAT[[i]][,k]<0.05))
#       col_names = c(col_names,paste(names(DAT)[i],colnames(DAT[[i]])[k],sep="_"))
#     }
#   }
#   rownames(Mat) = rownames(DAT[[1]])
#   colnames(Mat) = col_names
#   return(Mat)
# }
#
# #Returns the jaccard index of the genes shared by different experimental conditions
# intersection_gene_between_countings = function(DAT){
#   Mat = c()
#   Mat2 = c()
#   names_g = c()
#   for(i in 2:5){
#     for(j in (i+1):6){
#       int_genes = c()
#       for(k in 2:8){
#         gi = rownames(DAT[[i]])[DAT[[i]][,k]<0.05]
#         gj = rownames(DAT[[j]])[DAT[[j]][,k]<0.05]
#         int_genes = c(int_genes, round(100 * length(intersect(gi,gj))/length(union(gi,gj)),2))
#         int_genes2 = c(int_genes,paste(length(intersect(gi,gj)),"/",length(union(gi,gj)),sep=""))
#
#       }
#       names(int_genes) = c("balCount","time","dose","balCount * time", "balCount * dose", "time * dose", "balCount * time * dose")
#       names(int_genes2) = c("balCount","time","dose","balCount * time", "balCount * dose", "time * dose", "balCount * time * dose")
#
#       names_g = c(names_g,paste(names(DAT)[i],names(DAT)[j],sep="_"))
#
#       Mat = rbind(Mat,int_genes)
#       Mat2 = rbind(Mat2,int_genes2)
#
#     }
#   }
#   rownames(Mat) = names_g
#   rownames(Mat2) = names_g
#
#   return(list(Mat=Mat,Mat2=Mat2))
# }
#
# perc_rel_genes = function(DAT){
#   Mat = c()
#   names_g = c()
#   for(i in 2:6){
#     M = DAT[[i]]
#     pgenes = c()
#     for(j in 2:4){
#       pgenes = c(pgenes, round(100 * (sum(M[,j]<0.05)/nrow(M)),2))
#     }
#     Mat = rbind(Mat,c(names(DAT)[i],pgenes))
#   }
#
#   rownames(Mat) = Mat[,1]
#   Mat = Mat[,-1]
#   colnames(Mat) = c(c("balCount","time","balCount * time"))
#   library(reshape)
#   df = melt(Mat)
#   colnames(df) = c("Cell","Type","Perc")
#   gp = ggplot(data=df, aes(x = Type, y = Perc,fill=Cell)) +  geom_bar(stat="identity", position=position_dodge())
#
#   return(list(Mat=Mat,gp=gp))
# }
#
#
# barplot_ggplot = function(Mat){
#   library(reshape)
#   library(ggplot2)
#   df = melt(Mat)
#   colnames(df) = varnames = c("Cell","Type","JaccardIndex")
#   gp = ggplot(data=df, aes(x = Type, y = JaccardIndex,fill=Cell)) +  geom_bar(stat="identity", position=position_dodge())
#   gp
#   return(gp)
# }
#
# create_upset_path = function(PAT){
#   pathways = unique(PAT$Description)
#   cells = unique(PAT$Cell)
#   Mat = matrix(0,nrow=length(pathways),ncol=length(cells))
#   rownames(Mat) = pathways
#   colnames(Mat) = as.character(cells)
#
#   for(celi in cells){
#     idx = which(PAT$Cell %in% celi)
#     Mat[PAT$Description[idx],celi] = 1
#   }
#
#   return(Mat)
# }
#
# create_genes_tables = function(GSE,gene_list){
#   library(data.table)
#   library(xlsx)
#
#   for(i in 1:length(gene_list)){
#     max_sixe = max(unlist(lapply(gene_list[[i]],FUN =function(mat)dim(mat)[1])))
#     MAT = matrix("",max_sixe,2*length(names(gene_list[[i]])))
#     names_MAT = c()
#     start = 1
#     stop = 2
#     for(j in 1:length(gene_list[[i]])){
#       names_MAT = c(names_MAT,paste("Genes",names(gene_list[[i]])[j],sep="_"),paste("PVal",names(gene_list[[i]])[j],sep="_"))
#       MAT[1:nrow(gene_list[[i]][[j]]),start:stop] = gene_list[[i]][[j]]
#       start = start +2
#       stop = stop + 2
#     }
#     colnames(MAT) = names_MAT
#     write.table(MAT,sep="\t",file = paste(GSE,"_",names(gene_list)[i],".txt",sep=""),quote = FALSE,row.names = FALSE)
#     write.xlsx(MAT, paste(GSE,"_",names(gene_list)[i],".xlsx",sep=""))
#   }
#   #Mat = matrix("",nrow = max_size)
#
#   #Macro_genes = t(rbindlist(lapply(gene_list[[i]], function(x) data.table(t(x))),
#   #                           fill = TRUE))
#   #colnames(Macro_genes) =names(PAT_res$Gene_list[[1]])
#
#
# }
#
# gene_neutro_across_nano = function(exp_data,pheno_data){
#   cell_list = list()
#   for(cell_name in c("Dead.cells","Macro","Lympho","Neutro","Eosino","Epithel")){
#     pvalues_genes = c()
#     pb = txtProgressBar(min = 1,max=nrow(exp_data),style=3)
#     for(i in 1:nrow(exp_data)){
#       #take the value of the genes across the samples at time 1d
#       exp = as.numeric(as.vector(exp_data[i,]))
#       #the datafame contains the expression values for genes i for samples at time 1d and their doses
#
#       anova_df = data.frame(exp = exp,balCount=pheno_data[,cell_name], time = pheno_data[,"Day"],dose=pheno_data[,"Dose"])
#       anova_df$time = gsub(pattern = "d",replacement = "",x = anova_df$time)
#       anova_df$time = as.numeric(as.vector(anova_df$time))
#       anova_df$dose = gsub(pattern = "ug",replacement = "",x = anova_df$dose)
#       anova_df$dose = as.character(anova_df$dose)
#       anova_df$dose[anova_df$dose %in% "none"] = 0
#       anova_df$dose = as.numeric(anova_df$dose)
#       # fitting regression between gene expression, balcounts and time,
#       # having both balcount and time as continuous variables
#       anova_res = lm(exp ~ balCount * time * dose, anova_df)
#       xx = summary(anova_res)$coefficients
#       pvalues = xx[,"Pr(>|t|)"]
#
#       pvalues_genes = rbind(pvalues_genes,pvalues)
#       setTxtProgressBar(pb,i)
#     }
#     close(pb)
#     rownames(pvalues_genes) = rownames(exp_data)
#     #  pvalues_genes = apply(pvalues_genes,2,p.adjust,method="fdr")
#     cell_list[[cell_name]]= pvalues_genes
#   }
#   return(cell_list)
# }
#
#
# find_genes_associated_to_time_cel_nano = function(exp_data,pheno_data){
#   cell_list = list()
#   for(cell_name in "Neutro"){ #c("Dead.cells","Macro","Lympho","Neutro","Eosino","Epithel")
#     pvalues_genes = c()
#     pb = txtProgressBar(min = 1,max=nrow(exp_data),style=3)
#     for(i in 1:nrow(exp_data)){
#       #take the value of the genes across the samples at time 1d
#       exp = as.numeric(as.vector(exp_data[i,]))
#       #the datafame contains the expression values for genes i for samples at time 1d and their doses
#
#       anova_df = data.frame(exp = exp,balCount=pheno_data[,cell_name], time = pheno_data[,"Day"],dose=pheno_data[,"Dose"], treatment = pheno_data[,"CT"])
#       #  anova_df$time = gsub(pattern = "d",replacement = "",x = anova_df$time)
#       #  anova_df$time = as.numeric(as.vector(anova_df$time))
#       #  anova_df$dose = gsub(pattern = "ug",replacement = "",x = anova_df$dose)
#       #  anova_df$dose = as.character(anova_df$dose)
#       #  anova_df$dose[anova_df$dose %in% "none"] = 0
#       #  anova_df$dose = as.numeric(anova_df$dose)
#       # fitting regression between gene expression, balcounts and time,
#       # having both balcount and time as continuous variables
#       #anova_res1 = lm(exp ~ balCount, anova_df)
#       anova_res = lm(exp ~ balCount : dose : time : treatment, anova_df)
#       xx = summary(anova_res)$coefficients
#       pvalues = xx[,"Pr(>|t|)"]
#       tscore = xx[,"t value"]
#
#       pvalues_genes = rbind(pvalues_genes,c(pvalues[2:7],pvalues[2:7]<0.05,tscore[2:7],sum(pvalues[2:7]<0.05)))
#       setTxtProgressBar(pb,i)
#     }
#     close(pb)
#     rownames(pvalues_genes) = rownames(exp_data)
#     pvalues_genes=pvalues_genes[order(pvalues_genes[,19],decreasing = T),]
#
#     for(i in 1:6){
#       pvalues_genes[,i] = p.adjust(pvalues_genes[,i],method="fdr")
#     }
#
#     for(i in 1:nrow(pvalues_genes)){
#       pvalues_genes[i,8:13] = pvalues_genes[i,1:6]<0.05
#       pvalues_genes[i,19] = sum(pvalues_genes[i,1:7]<0.05)
#     }
#
#
#     #  pvalues_genes = apply(pvalues_genes,2,p.adjust,method="fdr")
#     cell_list[[cell_name]]= pvalues_genes
#   }
#   return(cell_list)
# }
#
#
