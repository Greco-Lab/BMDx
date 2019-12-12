
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
trend_test = function(pheno_data_list,expression_list, time_point_index = 4, dose_index = 2, time_point_id = "All", sample_index = 1, adj.pval = TRUE, p.th=0.01, nCores = 1){
  list_of_filtered_expression_values = list()
  list_of_variable_genes = list()
  list_of_non_variable_genes = list()
  list_of_matrices_with_anova_pvalue = list()
  
  withProgress(message = 'Running TREND TEST', detail = paste("Experiment: ",1 ,"/",length(pheno_data_list),sep=""), value = 1, {
    
  for(i in 1:length(pheno_data_list)){
    
    pTable = pheno_data_list[[i]]
    timep = unique(pTable[,time_point_index])
    
    if(time_point_id == "All"){
      for(tp in  timep){
        print(tp)
        
        pvalues_genes = compute_trend_test(exp_data = expression_list[[i]],
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
        
        c(filtered_expression_values, variable_genes, non_variable_genes,matrix_with_anova_pvalue) %<-% filter_pvalues(exp_data=expression_list[[i]],
                                                                                                                       pvalues_genes=pvalues_genes,
                                                                                                                       adj.pval = adj.pval,
                                                                                                                       p.th=p.th)
        
        list_of_filtered_expression_values[[names(pheno_data_list)[i]]][[as.character(tp)]] = filtered_expression_values
        list_of_variable_genes[[names(pheno_data_list)[i]]][[as.character(tp)]] = variable_genes
        list_of_non_variable_genes[[names(pheno_data_list)[i]]][[as.character(tp)]] = non_variable_genes
        list_of_matrices_with_anova_pvalue[[names(pheno_data_list)[i]]][[as.character(tp)]] = matrix_with_anova_pvalue
      }
    }else{
      pvalues_genes = compute_trend_test(exp_data = expression_list[[i]],
                                         pheno_data=pTable,
                                         time_t=tp,
                                         tpc = time_point_index,
                                         dc = dose_index,
                                         sc = sample_index)
      
      c(filtered_expression_values, variable_genes, non_variable_genes,matrix_with_anova_pvalue) %<-% filter_pvalues(exp_data=expression_list[[i]],
                                                                                                                     pvalues_genes=pvalues_genes,
                                                                                                                     adj.pval = adj.pval,
                                                                                                                     p.th=p.th)
      
      list_of_filtered_expression_values[[names(pheno_data_list)[i]]][[as.character(input$time_point_id)]] = filtered_expression_values
      list_of_variable_genes[[names(pheno_data_list)[i]]][[as.character(input$time_point_id)]] = variable_genes
      list_of_non_variable_genes[[names(pheno_data_list)[i]]][[as.character(input$time_point_id)]] = non_variable_genes
      list_of_matrices_with_anova_pvalue[[names(pheno_data_list)[i]]][[as.character(input$time_point_id)]] = matrix_with_anova_pvalue
    }
    
    
  }
  })
  
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

compute_trend_test = function(exp_data, pheno_data, time_t=24,tpc = 4, dc = 2, sc = 1, adj.pval = TRUE, p.th=0.01, nCores = 1){
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
    
    anova_df2 = aggregate(anova_df[, 1], list(anova_df$dose), median)
    colnames(anova_df2) = c("dose","exp")
    
    trend_test = trend::mk.test(anova_df$exp, alternative="two.sided", continuity=FALSE)
    pvalue = trend_test$p.value
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

filter_pvalues = function(exp_data, pvalues_genes, adj.pval = TRUE, p.th=0.01){
  
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