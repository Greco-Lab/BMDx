options(java.parameters = c("-XX:+UseConcMarkSweepGC", "-Xmx8192m"))

library(drc)
library(bmd)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(org.Rn.eg.db)
library(gProfileR)
library(ggplot2)
library(plotly)
library(shiny)
library(shinyjs)
library(shinyBS)
library(ggplotify)
library(RColorBrewer)
library(reshape)
library(GO.db)
library(GOSim)
library(igraph)
library(gtools)
library(alr3)
library(jtools)
library(plotly)
library(tibble)
library(readxl)
library(DT)
library(randomcoloR)
library(shinycssloaders)
library(zeallot)
library(gplots)
library(xlsx)
library(tidyverse)
library(gridExtra)
library(grid)
library(gtable)
library(shinycssloaders)
library(UpSetR)
library(XLConnect)
library(doParallel)
library(trend)

set.seed(12)
load(file ="data/updated_kegg_hierarhcy.RData")
load("data/mm_reactome_hierarchy.RData")
load("data/hm_reactome_hierarchy.RData")
load("data/rat_reactome_hierarchy.RData")
load("data/rat_kegg_hyerarchy.RData")
reduced_kegg_hierarchy = kegg_hierarchy
source("functions/bmd_utilities.R")
source("functions/path_grid_plot.R")
source("functions/enrich_functions.R")
source("functions/goHier.R")
source("functions/kegg_mats.R")
source("functions/mselect.R")
source("functions/ED.lin.R")
source("functions/bmdx_utilities.R")
source("functions/trend_test.R")

# this instruction increase the maximum size of the file that can be uploaded in a shiny application
options(shiny.maxRequestSize=300*1024^2) 

celTable <- NULL
phTable <- NULL
extractedList <- NULL

gxTable <- NULL
gxCorMat <- NULL
tmpMethod <- NULL

#Function to hide bsCollapsePanel
hideBSCollapsePanel <- function(session, panel.name)
{
  #print("In hideBSCollapsePanel")
  session$sendCustomMessage(type='hideBSCollapsePanelMessage', message=list(name=panel.name))
  #print("Exiting hideBSCollapsePanel")
}


#Function to get names of columns with datatype character
factorize_cols <- function(phTable, idx){
  for(i in idx){
    phTable[,i] <- factor(phTable[,i])
  }
  return(phTable)
}


shinyServer(function(input, output, session) {

  # Global variable list
  gVars <- shiny::reactiveValues(
    phTable=NULL,             # setting pheno data matrix
    sampleColID = FALSE,      # setting sample column number
    doseColID = NULL,         # setting dose column number
    TPColID=NULL,             # setting TP column number
    inputGx=NULL,             # setting gene expression matrix
    
    
    KEGG_MAT=NULL,
    lev1_h = NULL,
    lev2_h = NULL,
    lev3_h = NULL,
    n_groups = NULL,
    hierarchy = NULL,
    GList = NULL,
    pheno = NULL, #original pheno, used as backup
    exp_ann = NULL,
    reduced_kegg_hierarchy = NULL,
    toPlot = NULL,
    toPlotMap = NULL,
    clust_mat = NULL,
    samplesID = NULL,
    nSamples = NULL,
    nPath = NULL,
    DATA = NULL,
    
    rgList=NULL,
    celDir=NULL,
    totalSamples=NULL,
    filteredSamples=NULL,
    removedSamples=NULL,
    norm.data=NULL,
    pcChoices=NULL,
    comb.data=NULL,
    agg.data=NULL,
    comps=list(),
    loadedRaw=FALSE,
    QC_passed=FALSE,
    filtered=FALSE,
    normalized=FALSE,
    corrected=FALSE,
    gxTable=NULL, 
    dgxTable=NULL,
    clustered = 1
  )
  
  gVars$sepChoices <- c("TAB", ",", ";", "SPACE", "OTHER")
  gVars$quoteChoices <- c(NA, "SINGLE", "DOUBLE")
  gVars$pvAdjChoices <- c("Holm"="holm", "Hochberg"="hochberg", "Hommel"="hommel", "Bonferroni"="bonferroni", "Benjamini & Hochberg"="BH", "Benjamini & Yekutieli"="BY", "False Detection Rate"="fdr", "None"="none")
  gVars$normChoices <- c("Between Arrays"="BA", "Quantile"="quantile", "Variance Stabilizing"="vsn", "Cyclic Loess"="cl")
  gVars$baChoices <- c("None"="none", "Scale"="scale", "Quantile"="quantile", "Cyclic Loess"="cyclicloess")
  
  
  #gVars$inputGx <- 
  observeEvent(input$upload_gx_submit, {
    gxFile <- input$gx
    if (is.null(gxFile))
      return(NULL)
    
    print("Reading loaded gx File...")
    
    gx = read_excel_allsheets(filename = gxFile$datapath,tibble = FALSE)
    for(i in 1:length(gx)){
      gi = gx[[i]]
      gi = as.matrix(gi)
      
      gii = matrix(0, nrow = nrow(gi), ncol = ncol(gi))
      rownames(gii) = rownames(gi)
      colnames(gii) = colnames(gi)
      for(j in 1:nrow(gi)){
        gii[j,] = as.numeric(gi[j,])
      }
      gii = gii[, gVars$phTable[[i]][,gVars$sampleColID]]
      gx[[i]] = gii
    }
    
    gVars$inputGx <- gx
    

    shinyBS::toggleModal(session, "importGxModal", toggle="close")
    shinyBS::updateButton(session, "import_expr_submit", style="success", icon=icon("check-circle"))
    shinyBS::updateCollapse(session, "bsSidebar1", open="GENE FILTERING", style=list("LOAD EXPRESSION MATRIX"="success","GENE FILTERING"="warning"))
    print("open anova")
    shiny::updateTabsetPanel(session, "display",selected = "gExpTab")
    
  })
  
  gVars$gxTable <- shiny::reactive({
    if (is.null(gVars$inputGx) || is.null(gVars$dgxTable))
      return(NULL)
    
    gx <- gVars$inputGx
    gx <- gx[rownames(gVars$dgxTable),]
    return(gx)
  })
  
  gVars$inputDgx <- eventReactive(input$load_dgx_submit, {
    if(is.null(input$dgx))
      return(NULL)
    
    dgxFile <- input$dgx
    sepS <- input$sepS
    sepT <- input$sepT
    sepChar=NULL
    if(sepS=="OTHER"){
      sepChar <- sepT
    }else{
      if(sepS=="TAB"){
        sepChar="\t"
      }else if(sepS=="SPACE"){
        sepChar=" "
      }else{
        sepChar=sepS
      }
    }
    print(sepChar)
    
    quote <- input$quote
    print(quote)
    print(is.na(quote))
    if(is.na(quote) || quote=="NA"){
      quote <- ""
    }else if(quote=="SINGLE"){
      quote <- "'"
    }else if(quote=="DOUBLE"){
      quote <- '"'
    }
    
    rowNames <- NULL
    con <- file(dgxFile$datapath, "r", blocking=FALSE)
    fileLines <- readLines(con)
    fileLines <- gsub("\t\t", "\tNA\t", fileLines)
    fileLines <- gsub("\t$", "\tNA", fileLines)
    close(con)
    colNumByRowDist <- table(sapply(fileLines, function(x) {length(strsplit(x, sepChar)[[1]])}, USE.NAMES=FALSE))
    if(length(colNumByRowDist) > 1){
      fileLines2 <- fileLines[-1]
      colNumByRowDist2 <- table(sapply(fileLines2, function(x) {length(strsplit(x, sepChar)[[1]])}, USE.NAMES=FALSE))
      if(length(colNumByRowDist2) > 1){
        shinyjs::info(paste0("Separating character '", sepChar, "', results in inconsistent number of columns!\n\nPlease check the input file format and select the correct separating character!"))
        gVars$dgxLoaded <- NULL
        return(NULL)
      }
      rowNames <- 1
    }
    
    gVars$dgxLoaded <- 1
    gVars$dgxFileName <- dgxFile$name
    dgxTable <- read.csv(dgxFile$datapath, row.names=1, header=TRUE, sep=sepChar, stringsAsFactors=FALSE, quote=quote, as.is=TRUE, strip.white=TRUE, check.names=FALSE)
    
    print("str(dgxTable) -- after:")
    print(str(dgxTable))
    return(dgxTable)
  })
  
  gVars$dgxColChoices <- reactive({
    if(is.null(gVars$dgxLoaded))
      return(c("NA"))
    
    choicesVec <- seq(1,ncol(gVars$inputDgx()))
    choicesNames <- paste0("Column ", choicesVec)
    names(choicesVec) <- choicesNames
    return(choicesVec)
  })
  
  output$dgxDT <- DT::renderDataTable({
    if(is.null(gVars$inputDgx))
      return(NULL)
    
    dgxTable <- gVars$inputDgx()
    colnames(dgxTable) <- paste0(colnames(dgxTable), " [", c(1:ncol(dgxTable)), "]")
    DT::datatable(dgxTable, filter="none", 
                  options = list(
                    ch = list(regex=TRUE, caseInsensitive=FALSE), 
                    scrollX=TRUE, 
                    pageLength=2,
                    lengthMenu=c(1,2,3),
                    ordering=FALSE
                  )
    )
  },server=TRUE)
  
  
  #upload an excell file with N sheets such as the number of experiments that we are considering
  gVars$inputPh <- eventReactive(input$load_pheno_submit, {
    if(is.null(input$fPheno))
      return(NULL)
    
    phFile <- input$fPheno
    print("separator -->")
    phTable = read_excel_allsheets(filename = phFile$datapath,tibble = FALSE)
    
    phTable = phTable
    
    gVars$phLoaded <- 1
    print("str(phTable) -- after:")
    print(str(phTable))
    
    return(phTable)
  })
  
  gVars$phColTypes <- reactive({
    if(is.null(gVars$inputPh()))
      return(NULL)
    
    phTable <- gVars$inputPh()[[1]]
    colNames <- list("c"=NULL,"n"=NULL)
    coltypes <- unlist(lapply(phTable, class))
    coltypes.charOnly.idx <- which(coltypes=="character")
    coltypes.nonChar.idx <- which(!coltypes=="character")
    coltypes.charOnly.len <- length(coltypes.charOnly.idx)
    coltypes.nonChar.len <- length(coltypes.nonChar.idx)
    if(coltypes.charOnly.len>0){
      colNames[["c"]] <- colnames(phTable)[coltypes.charOnly.idx]
    }
    if(coltypes.nonChar.len>0){
      colNames[["n"]] <- colnames(phTable)[coltypes.nonChar.idx]
    }
    print(str(colNames))
    return(colNames)
  })
  
  output$phRowsText <- renderText({
    if(is.null(gVars$inputPh())){
      nRow <- "NA"
    }else{
      nRow <- nrow(gVars$inputPh()[[1]])
    }
    return(paste0("Samples: ", nRow))
  })
  
  output$phColsText <- renderText({
    if(is.null(gVars$inputPh())){
      nCol <- "NA"
    }else{
      nCol <- ncol(gVars$inputPh()[[1]])
    }
    return(paste0("Variables: ", nCol))
  })
  
  gVars$phColChoices <- reactive({
    if(is.null(gVars$phLoaded))
      return(c("NA"))
    phenoList = gVars$inputPh()
    choicesVec <- seq(1,ncol(phenoList[[1]]))
    choicesNames <- paste0("Variable ", choicesVec)
    names(choicesVec) <- choicesNames
    return(choicesVec)
  })
  
  output$phenoDT <- DT::renderDataTable({
    # if(is.null(gVars$phTable))
    #   return(NULL)
    
    
    phenoList = gVars$inputPh()
    phTable <- phenoList[[1]]
    colnames(phTable) <- paste0(colnames(phTable), " [", c(1:ncol(phTable)), "]")
    DT::datatable(phTable, filter="none",
                  options = list(
                    ch = list(regex=TRUE, caseInsensitive=FALSE),
                    scrollX=TRUE,
                    pageLength=2,
                    lengthMenu=c(1,2,3),
                    ordering=F
                  )
    )
  },server=TRUE)
  
  output$phenoTypesRH <- rhandsontable::renderRHandsontable({
    shiny::validate(
      need(!is.null(gVars$inputPh()), "No phenotype file!")
    )
    phenoList = gVars$inputPh()
    phTable <- phenoList[[1]]
    colNames <- paste0(colnames(phTable), " [", c(1:ncol(phTable)), "]")
    colTypesList <- gVars$phColTypes()
    print(str(colTypesList))
    colClass <- unlist(lapply(phTable, class))
    print(str(colClass))

    colTypesDF <- data.frame(Variable=colNames, Type="-", Class=colClass, t(head(phTable)), stringsAsFactors=FALSE)
    colnames(colTypesDF)[c(4:ncol(colTypesDF))] <- paste0("Sample", c(1:(ncol(colTypesDF)-3)))
    if(!is.null(colTypesList[["c"]])){
      print("In colTypes C...")
      cIdx <- which(colnames(phTable) %in% colTypesList[["c"]])
      print(cIdx)
      if(length(cIdx>0)){
        colTypesDF[cIdx,"Type"] <- "factor"
      }
    }
    if(!is.null(colTypesList[["n"]])){
      print("In colTypes N...")
      nIdx <- which(colnames(phTable) %in% colTypesList[["n"]])
      print(nIdx)
      if(length(nIdx>0)){
        colTypesDF[nIdx,"Type"] <- "vector"
      }
    }
    rhandsontable::rhandsontable(colTypesDF, rowHeaders=NULL, readOnly=TRUE, contextMenu=FALSE) %>%
      hot_col("Type", readOnly=FALSE, type="dropdown", source=c("factor","vector"),
              renderer = "function (instance, td, row, col, prop, value, cellProperties) {
              Handsontable.renderers.NumericRenderer.apply(this, arguments);
              if (value == 'factor') {
              td.style.background = 'lightblue';
              } else if (value == 'vector') {
              td.style.background = 'lightgreen';
              }
              }"
              )
})
  
  observeEvent(input$skip_anova_filtering_button, {
      shiny::validate(
        need(!is.null(gVars$phTable), "No Phenotype File Provided!")
      )
      
      shiny::validate(
        need(!is.null(gVars$inputGx), "No Expression File Provided!")
      )
      
      print("skip anova")
      
      #timep = unique(gVars$phTable[,gVars$TPColID])
      EXP_FIL_List = list()
      VG_List = list()
      NVG_List = list()
      PValMat_List = list()
      
      for(i in 1:length(gVars$phTable)){
        
        pTable = gVars$phTable[[i]]
        timep = unique(pTable[,gVars$TPColID])
        
        for(tp in timep){
          EXP_FIL_List[[names(gVars$phTable)[i]]][[as.character(tp)]] = gVars$inputGx[[i]]
          VG_List[[names(gVars$phTable)[i]]][[as.character(tp)]] = rownames(gVars$inputGx[[i]])
          NVG_List[[names(gVars$phTable)[i]]][[as.character(tp)]] = c()
          PValMat_List[[names(gVars$phTable)[i]]][[as.character(tp)]] = NULL
        }
        
      }
      
      gVars$EXP_FIL = EXP_FIL_List # EXP_FIL
      gVars$var_genes = VG_List #var_genes
      gVars$not_var_genes = NVG_List #not_var_genes
      gVars$PValMat= PValMat_List #PValMat
      
      print("no anova")
      
      shinyBS::updateButton(session, "anova_filtering_button", style="success", icon=icon("check-circle"))
      shinyBS::updateButton(session, "skip_anova_filtering_button", style="success", icon=icon("check-circle"))
      
      shinyBS::updateCollapse(session, "bsSidebar1", open="COMPUTE BMD", style=list("GENE FILTERING"="success","COMPUTE BMD"="warning"))
      

  })
  
  observeEvent(input$trend_analysis,{
    shiny::validate(
      need(!is.null(gVars$phTable), "No Phenotype File Provided!")
    )
    
    shiny::validate(
      need(!is.null(gVars$inputGx), "No Phenotype File Provided!")
    )
    
  print("in trend test")
    
    trend_test_res = trend_test(pheno_data_list=gVars$phTable,
                                expression_list=gVars$inputGx,
                                time_point_index = gVars$TPColID,
                                dose_index = gVars$doseColID,
                                time_point_id = input$time_point_id,
                                sample_index = gVars$sampleColID,
                                adj.pval = input$TrendAdjBool,
                                p.th=as.numeric(input$trendPvalTh),
                                nCores = as.numeric(input$trendNCores))

    if(is.null(trend_test_res)==FALSE){
    #if(trend_completed){
      gVars$EXP_FIL = trend_test_res$list_of_filtered_expression_values # EXP_FIL
      gVars$var_genes = trend_test_res$list_of_variable_genes #var_genes
      gVars$not_var_genes =trend_test_res$list_of_non_variable_genes #not_var_genes
      gVars$PValMat= trend_test_res$list_of_matrices_with_anova_pvalue #PValMat
      
      shinyBS::toggleModal(session, "computeTrend", toggle="close")
      shinyBS::updateButton(session, "trend_filtering_button", style="success", icon=icon("check-circle"))
      shinyBS::updateButton(session, "skip_anova_filtering_button", style="success", icon=icon("check-circle"))
      shinyBS::updateCollapse(session, "bsSidebar1", open="COMPUTE BMD", style=list("GENE FILTERING"="success","COMPUTE BMD"="warning"))
      shiny::updateTabsetPanel(session, "display",selected = "AnovaTab")
    }
      
  })
  
  observeEvent(input$anova_analysis, {
    shiny::validate(
      need(!is.null(gVars$phTable), "No Phenotype File Provided!")
    )
    
    shiny::validate(
      need(!is.null(gVars$inputGx), "No Phenotype File Provided!")
    )
    
    EXP_FIL_List = list()
    VG_List = list()
    NVG_List = list()
    PValMat_List = list()
    anova_completed = TRUE
    
    withProgress(message = 'Running ANOVA', detail = paste("Experiment: ",1 ,"/",length(gVars$phTable),sep=""), value = 1, {
      
    
    for(i in 1:length(gVars$phTable)){
      
      pTable = gVars$phTable[[i]]
      timep = unique(pTable[,gVars$TPColID])
      
      if(input$time_point_id == "All"){
        for(tp in  timep){
          print(tp)
          
          print("ciao")
          pvalues_genes = compute_anova(exp_data = gVars$inputGx[[i]], #[runif(n = 500,min = 1,max = nrow(gVars$inputGx)),], # TOGLIERE LA SELEZIONE RANDOM DEI GENI
                                        pheno_data=gVars$phFactor[[i]], 
                                        time_t=tp,
                                        tpc = gVars$TPColID, 
                                        dc = gVars$doseColID, 
                                        sc = gVars$sampleColID,nCores = as.numeric(input$anovaNCores))
          
          if(is.null(pvalues_genes)){
            
            output$anova_update <- renderText({ 
              "Cannot perform anova with only one dose-level. check dose and time point columns"
            })
            
            #shinyalert("Oops!", "Cannot perform anova with only one dose-level. check dose and time point columns", type = "error")
            anova_completed=FALSE
            return(NULL)
          }else{
            output$anova_update <- renderText({ 
              "Anova running..."
            })
          }
          
          c(EXP_FIL, var_genes, not_var_genes,PValMat) %<-% filter_anova(exp_data=gVars$inputGx[[i]], 
                                                                         pvalues_genes=pvalues_genes, 
                                                                         adj.pval = input$adjBool, 
                                                                         p.th=as.numeric(input$anovaPvalTh))
          
          EXP_FIL_List[[names(gVars$phTable)[i]]][[as.character(tp)]] = EXP_FIL
          VG_List[[names(gVars$phTable)[i]]][[as.character(tp)]] = var_genes
          NVG_List[[names(gVars$phTable)[i]]][[as.character(tp)]] = not_var_genes
          PValMat_List[[names(gVars$phTable)[i]]][[as.character(tp)]] = PValMat
        }
      }else{
        pvalues_genes = compute_anova(exp_data = gVars$inputGx, #[runif(n = 500,min = 1,max = nrow(gVars$inputGx)),], # TOGLIERE LA SELEZIONE RANDOM DEI GENI
                                      pheno_data=gVars$phTable, 
                                      time_t=input$time_point_id,
                                      tpc = gVars$TPColID, 
                                      dc = gVars$doseColID, 
                                      sc = gVars$sampleColID)
        
        c(EXP_FIL, var_genes, not_var_genes,PValMat) %<-% filter_anova(exp_data=gVars$inputGx, 
                                                                       pvalues_genes=pvalues_genes, 
                                                                       adj.pval = input$adjBool, 
                                                                       p.th=as.numeric(input$anovaPvalTh))
        EXP_FIL_List[[names(gVars$phTable)[i]]][[as.character(input$time_point_id)]] = EXP_FIL
        VG_List[[names(gVars$phTable)[i]]][[as.character(input$time_point_id)]] = var_genes
        NVG_List[[names(gVars$phTable)[i]]][[as.character(input$time_point_id)]] = not_var_genes
        PValMat_List[[names(gVars$phTable)[i]]][[as.character(input$time_point_id)]] = PValMat
      }
      
      incProgress(1/length(gVars$phTable), detail = paste("ANOVA Experiment: ",i+1," Experiment: ",i+1 ,"/",length(gVars$phTable),sep=""))
      
    }
    
})
    
    if(anova_completed){
      gVars$EXP_FIL = EXP_FIL_List # EXP_FIL
      gVars$var_genes = VG_List #var_genes
      gVars$not_var_genes = NVG_List #not_var_genes
      gVars$PValMat= PValMat_List #PValMat
      
      shinyBS::toggleModal(session, "computeAnova", toggle="close")
      shinyBS::updateButton(session, "anova_filtering_button", style="success", icon=icon("check-circle"))
      shinyBS::updateButton(session, "skip_anova_filtering_button", style="success", icon=icon("check-circle"))
      shinyBS::updateCollapse(session, "bsSidebar1", open="COMPUTE BMD", style=list("GENE FILTERING"="success","COMPUTE BMD"="warning"))
      shiny::updateTabsetPanel(session, "display",selected = "AnovaTab")
    }
    
  })
  

  observeEvent(input$bmd_analysis, {
    shiny::validate(
      need(!is.null(gVars$phTable), "No Phenotype File Provided!")
    )
    
    shiny::validate(
      need(!is.null(gVars$EXP_FIL), "No Phenotype File Provided!")
    )
    
    print("BMD")
    
    withProgress(message = 'Running BMD', detail = paste("Experiment: ",1 ,"/",length(gVars$phTable),sep=""), value = 1, {
    MQ_BMDList = list()
    MQ_BMDListFiltered = list()
    MQ_BMDListFilteredValues = list()
    
    for(j in 1:length(gVars$phTable)){
      print("this is the j ----> ")
      print(j)
      
      pTable = gVars$phTable[[j]]
      timep = as.numeric(names(gVars$EXP_FIL[[j]]))#unique(pTable[,gVars$TPColID])
      print(timep)
      
      maxDose = max(as.numeric(unique(gVars$phTable[[j]][,gVars$doseColID])))
      minDose = min(as.numeric(unique(gVars$phTable[[j]][,gVars$doseColID])))
      
      for(i in timep){

        print("xxx")
        
        MQ_BMDList[[names(gVars$EXP_FIL)[j]]][[as.character(i)]]  = compute_bmd(exp_data=gVars$EXP_FIL[[names(gVars$phTable)[j]]][[as.character(i)]], pheno_data=gVars$phTable[[names(gVars$phTable)[j]]],
                              time_t=as.character(i), #interval_type = input$Interval,
                              tpc = gVars$TPColID, 
                              dc = gVars$doseColID, 
                              sc = gVars$sampleColID, sel_mod_list = as.numeric(input$ModGroup),rl = as.numeric(input$RespLev),constantVar = input$constantVar, nCores = as.numeric(input$BMDNCores))

        MQ_BMDListFiltered[[names(gVars$phTable)[j]]][[as.character(i)]]   = BMD_filters(BMDRes = MQ_BMDList[[names(gVars$EXP_FIL)[j]]][[as.character(i)]],max_dose = as.numeric(maxDose),min_dose = as.numeric(minDose),max_low_dos_perc_allowd = as.numeric(input$min_dose_perc),max_max_dos_perc_allowd = as.numeric(input$max_dose_perc)  ,loofth = as.numeric(input$LOOF))
        MQ_BMDListFilteredValues[[names(gVars$phTable)[j]]][[as.character(i)]]  =  MQ_BMDListFiltered[[names(gVars$phTable)[j]]][[as.character(i)]]$BMDValues_filtered
      }
      incProgress(1/length(gVars$phTable), detail = paste("BMD Experiment: ",j+1," Experiment: ",j+1 ,"/",length(gVars$phTable),sep=""))
      
    }
    })

    gVars$MQ_BMD = MQ_BMDList
    gVars$MQ_BMD_filtered = MQ_BMDListFiltered
    gVars$MQ_BMDListFilteredValues = MQ_BMDListFilteredValues
    gVars$MQ_BMD_filtered2 = gVars$MQ_BMD_filtered
    
    shinyBS::toggleModal(session, "computeBMD", toggle="close")
    shinyBS::updateButton(session, "bmd_button", style="success", icon=icon("check-circle"))
    shinyBS::updateCollapse(session, "bsSidebar1", open="PATHWAY ENRICHMENT", style=list("COMPUTE BMD"="success","PATHWAY ENRICHMENT"="warning"))
    shiny::updateTabsetPanel(session, "display",selected = "BMDTab")
  })
  
  output$upsetplot = renderPlot({
    print(input$range)
    print(class(input$range[1]))
    
    # if(is.null(gVars$GList) || is.null(gVars$pheno)){
    #   print('no input')
    #   return(NULL)
    # }
    shiny::validate(
      need(!is.null(gVars$MQ_BMD_filtered), "No BMD performed!")
    )
    
    GList = list()
    index = 1
    
    for(j in names(gVars$EXP_FIL)){ # for each experiment
      
      pTable = gVars$phTable[[j]]
      timep = unique(pTable[,gVars$TPColID])
      print(timep)
      
      for(i in timep){
        BMDFilMat = gVars$MQ_BMD_filtered[[j]][[as.character(i)]]$BMDValues_filtered
        
        genelist = BMDFilMat[,c("Gene","BMD")]
        GList[[paste(j,as.character(i),sep="_")]] = genelist
        index = index + 1
        
        print(i)
      }   
    }
    
    #print(input$range)
    
    th1 = input$range[1]/100 #0
    th2 = input$range[2]/100 #10
     
    listUpset = list()
    for(i in 1:length(GList)){
      
      mini = quantile(GList[[i]][,2],probs = th1)
      maxi = quantile(GList[[i]][,2],probs = th2)
      idx = which(GList[[i]][,2]>= mini & GList[[i]][,2]<= maxi)
      listUpset[[names(GList)[i]]] = as.character(GList[[i]][idx,1])
    }
    
    # drugs = unlist(lapply(X = strsplit(names(listUpset),split = "_"), FUN = function(elem)elem[1]))
    # TP = unlist(lapply(X = strsplit(names(listUpset),split = "_"), FUN = function(elem)elem[2]))
    # 
    # drcolor = rainbow(length(unique(drugs)))
    # names(drcolor) = unique(drugs)
    upset(fromList(listUpset),nsets = length(listUpset),nintersects = input$maxInt, group.by = input$groupby,order.by = input$orderby)
   # upset(fromList(listUpset),nsets = length(listUpset),nintersects = 50,group.by = "sets",order.by = "freq")
    
  })
  
  output$downloadFunmapponeInput <- downloadHandler(
    
    filename = function() {
      paste("input_funmappone.xlsx", sep = "")
    },
    content = function(file) {
      
      if(is.null(gVars$GList) || is.null(gVars$pheno)){
        print('no input')
        
        return(NULL)
      }
      print('Thank you for clicking')
      
      shinyjs::html(id="loadingText", "Saving tables")
      shinyjs::show(id="loading-content")
      on.exit({
        print("inside on exit")
        shinyjs::hide(id="loading-content", anim=TRUE, animType="fade")    
      })
      
      print("I'm saving anova tables")
      
      write.xlsx(gVars$GList[[1]], file, sheetName = names(gVars$GList)[1], row.names = FALSE) 
      
      if(length(gVars$GList)>1){
        for(i in 2:length(gVars$GList)){
          write.xlsx(gVars$GList[[i]], file, sheetName = names(gVars$GList)[i], append = TRUE, row.names = FALSE) 
        }
      }
      
      write.xlsx(gVars$pheno, file, sheetName = "Groups", append = TRUE, row.names = FALSE) 
      
      
      print("Funmappone input stored!")
    }
  )
  

  observeEvent(input$enrichment_analysis, {
    shiny::validate(
      need(!is.null(gVars$MQ_BMD_filtered), "No BMD performed!")
    )
    shinyjs::html(id="loadingText", "COMPUTING ENRICHMENT")
    shinyjs::show(id="loading-content")
    
    # thrs = input$filterBMD
    # print(thrs)
    
    print("Enrichment Analysis")
    nExp = length(gVars$EXP_FIL)
    
    EnrichRes = list()
    
    tryCatch(
      {
      GList = list()
      grouping_experiment = c()
      grouping_TP = c()
      grouping_unique = c()
      index = 1
      
      for(j in names(gVars$EXP_FIL)){ # for each experiment
        
        pTable = gVars$phTable[[j]]
        timep = as.numeric(names(gVars$EXP_FIL[[j]]))#unique(pTable[,gVars$TPColID])
        
        #timep = unique(pTable[,gVars$TPColID])
        print(timep)
        
        for(i in timep){
          BMDFilMat = gVars$MQ_BMD_filtered[[j]][[as.character(i)]]$BMDValues_filtered
          
          # tths = quantile(BMDFilMat[,"BMD"],probs = thrs/100)
          # iidx = which(BMDFilMat[,"BMD"]>=tths[1] & BMDFilMat[,"BMD"]<=tths[2])
          
          genelist = BMDFilMat[,c("Gene","BMD")]
          GList[[paste(j,as.character(i),sep="_")]] = genelist
          grouping_experiment = c(grouping_experiment,j)
          grouping_TP = c(grouping_TP,as.character(i))
          grouping_unique = c(grouping_unique,index)
          index = index + 1
          
          print(i)
        }   
      }
  
      organism = input$organism
      annType = input$idtype

        GList = convert_genes(organism = input$organism, GList=GList, annType = input$idtype)
        Mp = cbind(names(GList),1:length(GList),grouping_experiment,grouping_TP)  #DF[[length(DF)]]
        #Mp = cbind(names(GList),grouping)  #DF[[length(DF)]]
        
        colnames(Mp) = c("Chemical","Grouping", "Experiment","Timepoint")
        pheno = cbind(Mp[,2],Mp[,1])
        
        #save(GList, pheno, file = "funmapponeInput.RData")
        
        gVars$GList = GList
        gVars$pheno = pheno
        gVars$exp_ann = gVars$pheno

        DTa = matrix("",ncol = length(GList),nrow = max(unlist(lapply(GList, FUN = nrow))))
        for(i in 1:(length(GList))){
          gli = as.character(GList[[i]][,1])
          if(length(gli)>0) DTa[1:nrow(GList[[i]]),i]= gli
        }
        
        colnames(DTa) = names(GList)
        
        if(is.null(DTa)){
          gVars$DATA = NULL
        }else{
          gVars$DATA = DTa
        }

        gVars$toPlot <- NULL # refresh plot map in the Plot Maps tab
        gVars$clust_mat <- NULL
        gVars$KEGG_MAT <- NULL
        
        DAT = gVars$DATA   
        if(input$organism == "Mouse"){          
          org = "mmu"
          reactome_hierarchy = mm_reactome_hierarchy
          reactome_hierarchy$ID = unlist(reactome_hierarchy$Pathway)
          org_enrich = "mmusculus"
          reactome_hierarchy[,1] = mouse_map[reactome_hierarchy[,1],2]
          reactome_hierarchy[,2] = mouse_map[unlist(reactome_hierarchy[,2]),2]
          reactome_hierarchy[,3] = mouse_map[unlist(reactome_hierarchy[,3]),2]
        }
        if(input$organism == "Human"){
          org = "hsa"
          reactome_hierarchy = hm_reactome_hierarchy
          reactome_hierarchy$ID = unlist(reactome_hierarchy$Pathway)
          reactome_hierarchy[,1] = human_map[reactome_hierarchy[,1],2]
          reactome_hierarchy[,2] = human_map[unlist(reactome_hierarchy[,2]),2]
          reactome_hierarchy[,3] = human_map[unlist(reactome_hierarchy[,3]),2]
          org_enrich = "hsapiens"
        }
        if(input$organism == "Rat"){
          org = "Rat"
          reactome_hierarchy = rat_reactome_hierarchy
          reactome_hierarchy$ID = unlist(reactome_hierarchy$Pathway)
          reactome_hierarchy[,1] = rat_map[reactome_hierarchy[,1],2]
          reactome_hierarchy[,2] = rat_map[unlist(reactome_hierarchy[,2]),2]
          reactome_hierarchy[,3] = rat_map[unlist(reactome_hierarchy[,3]),2]
          org_enrich = "rnorvegicus"
          
        }
        
        ####### remove reactome duplicates
        reactome_hierarchy=unique(reactome_hierarchy)
        
        # Compute enrichment
        annType = input$EnrichType
        GOType = input$GOType
        
        if(annType=="KEGG"){
          gVars$hierarchy = kegg_hierarchy
          type_enrich = "KEGG"
        }
        if(annType=="REACTOME"){
          type_enrich = "REAC"
          gVars$hierarchy = reactome_hierarchy
        }
        if(annType=="GO"){
          if(GOType == "BP"){
            #create geograph object
            makeGOGraph(ont = "bp") -> geograph
            #convert graphNEL into igraph
            igraph.from.graphNEL(geograph) -> igraphgeo
            #make igraph object undirected
            igraphgeo = as.undirected(igraphgeo)
            #set root as BP root term
            root = "GO:0008150"
            type_enrich = "GO:BP"
          }
          if(GOType == "CC"){
            #EnrichDatList = all_GO_CC
            makeGOGraph(ont = "cc") -> geograph
            igraph.from.graphNEL(geograph) -> igraphgeo
            igraphgeo = as.undirected(igraphgeo)
            root="GO:0005575"
            type_enrich = "GO:CC"
          }
          if(GOType == "MF"){
            #EnrichDatList = all_GO_MF
            makeGOGraph(ont = "mf") -> geograph
            igraph.from.graphNEL(geograph) -> igraphgeo
            igraphgeo = as.undirected(igraphgeo)
            root="GO:0003674"
            type_enrich = "GO:MF"
          }
        }
        
        GL = gVars$GList
        pval = as.numeric(input$pvalueTh)
        adjust_method = input$pcorrection
        org = org_enrich
        type = type_enrich
        print("Before enrichment")

        if(input$pcorrection == "none"){ 
          print("Nominal PValue")
          #as.numeric(input$pvalueTh)
          EnrichDatList = lapply(gVars$GList,enrich,type_enrich,org_enrich,1,"bonferroni", sig = FALSE, mis = as.numeric(input$min_intersection), only_annotated=input$only_annotated)
          for(i in 1:length(EnrichDatList)){
            ERi = EnrichDatList[[i]]
            ERi$pValueAdj = ERi$pValueAdj / length(ERi$pValueAdj)
            ERi$pValue = ERi$pValue / length(ERi$pValue)
            EnrichDatList[[i]] = ERi
          }
        }else{
          EnrichDatList = lapply(gVars$GList,enrich,type_enrich,org_enrich,as.numeric(input$pvalueTh),input$pcorrection, sig = TRUE, mis = as.numeric(input$min_intersection),only_annotated=input$only_annotated)
        }
        
        gVars$EnrichDatList = EnrichDatList
        
        print("After enrichment")
        
        if(input$EnrichType == "GO"){
          #find the list of term into  the enriched list
          #res2 = filterGO(EnrichDatList,go_type = GOType)
          
          res2 = filterGO(EnrichDatList,go_type = input$GOType)
          EnrichDatList = res2$EnrichDatList
          go_terms = res2$goTerm
          print("compute GO hierarchy")
          
          #Compute all the shortest path between the root to the go_terms
          asp = all_shortest_paths(graph = igraphgeo,from = root, to = go_terms)
          
          #reduce the shortest path to length 3 to build a 3 level hierarchy
          go_hierarchy = matrix("",nrow = length(asp$res),ncol = 3)
          
          for(i in 1:length(asp$res)){
            nn = names(asp$res[[i]])
            nn = nn[2:length(nn)]
            if(length(nn)<3){
              nn2 = rep(nn[length(nn)],3)
              nn2[1:length(nn)] = nn
            }else{
              nn2 = c(nn[1:2],nn[length(nn)])
            }
            go_hierarchy[i,] =nn2   
          }
          
          go_hierarchy = unique(go_hierarchy)
          colnames(go_hierarchy) = c("level1","level2","level3")
          go_hierarchy = as.data.frame(go_hierarchy)
          go_hierarchy$ID = go_hierarchy[,3]
          
          idx = which(is.na(go_hierarchy[,1]))
          if(length(idx)>0){
            go_hierarchy = go_hierarchy[-idx,]
          }
          go_hierarchy[,1] = Term(GOTERM[as.character(go_hierarchy[,1])])
          go_hierarchy[,2] = Term(GOTERM[as.character(go_hierarchy[,2])])
          go_hierarchy[,3] = Term(GOTERM[as.character(go_hierarchy[,3])])
          
          # colnames(hier_names) = c("Level1","Level2","Pathway")
          gVars$hierarchy = go_hierarchy
          print(gVars$hierarchy)
          
        }
        
        print(head(gVars$hierarchy))
        
        NCol = ncol(gVars$GList[[1]])
        print("NCol -------- >>>>> ")
        print(NCol)
        
        #if(input$fileType %in% "GenesOnly"){
        if(input$MapValueType == "PVAL"){
          print("only genes")
          M = kegg_mat_p(EnrichDatList,hierarchy = gVars$hierarchy)
        }else{
          
          if(NCol==1){
            shinyjs::info("No Modification provided! Enrichment will be performed by using only pvalues")
            M = kegg_mat_p(EnrichDatList,hierarchy = gVars$hierarchy)
            
          }else{
            M1 = kegg_mat_p(EnrichDatList,hierarchy = gVars$hierarchy)
            M2 = kegg_mat_fc(EnrichDatList = EnrichDatList,hierarchy = gVars$hierarchy,GList = gVars$GList, summ_fun=get(input$aggregation))
            
            print("LOGGGGG")
            
            # if(input$MapValueType == "PVAL"){
            #   M = M1
            # }
            if(input$MapValueType == "FC"){
              M = M2
            }
            if(input$MapValueType == "FCPV"){
              M = M2 * -log(M1)
            }
            rownames(M) = rownames(M1)
            colnames(M) = colnames(M1)
          }
          
        }
        
        gVars$KEGG_MAT = M
        
        hierarchy <- collapse_paths(kegg_hierarchy =gVars$hierarchy,kegg_mat_cell = gVars$KEGG_MAT, collapse_level = 3)
        mat <- hierarchy[[1]]
        hier <- hierarchy[[2]]
        
        print(hier)
        
        gVars$lev1_h = as.list(c("All",unique(as.character(hier[,1]))))
        
        gVars$lev2_h = list("All" = "All")
        for(i in unique(hier[,1])){
          gVars$lev2_h[[i]] = as.character(unique(hier[which(hier[,1] %in% i),2]))
        }
        
        gVars$lev3_h = list("All" = "All")
        for(i in unique(hier[,2])){
          gVars$lev3_h[[i]] = as.character(unique(hier[which(hier[,2] %in% i),3]))
        }
        
        # code added here
        #Compute genes matrix
        Mgenes <- kegg_mat_genes(EnrichDatList=EnrichDatList, hierarchy=gVars$hierarchy)
        print("dim(Mgenes)")
        print(dim(Mgenes))
        print("str(Mgenes)")
        print(str(Mgenes))
        
        gVars$KEGG_MAT_GENES <- Mgenes
        
        
        on.exit({
          print("inside on exit")
          shinyjs::hide(id="loading-content", anim=TRUE, animType="fade")    
        })
        
        output$updatedPat <- renderText({
          DF = DATA()
          shiny::validate(need(expr = !is.null(DF),message = "") )
          
          x = paste("Pathway computed! Number of pathways for each sample: ")
          
          return(x)
        })
        
        output$colSumsPat <- DT::renderDataTable({
          DF = gVars$KEGG_MAT
          shiny::validate(need(expr = !is.null(DF),message = "") )
          
          M = matrix("",1,ncol = nrow(DF))
          for(i in 1:nrow(DF)){
            M[1,i]= sum(is.na(DF[i,])==FALSE)
          }
          colnames(M) = rownames(DF)
          #print(M)
          DT::datatable(M, filter = "none", options = list(scrollX = TRUE))
          
          #if(sum(M) ==0) gVars$noEnrich = TRUE
          
        })
        
        updateTabsetPanel(session = session,inputId = "page_id", selected = "PlotMaps")
        
      },
      error = function(e) {
        # return a safeError if a parsing error occurs
        shinyjs::info(e$message)
        shinyjs::hide(id="loading-content", anim=TRUE, animType="fade")    
       })
    
    shinyBS::toggleModal(session, "enrichPathways", toggle="close")
    shinyBS::updateButton(session, "enrich_button", style="success", icon=icon("check-circle"))
    shinyBS::updateCollapse(session, "bsSidebar1", open="COMPARE TIMEPOINTS", style=list("PATHWAY ENRICHMENT"="success","TIMEPOINTS"="warnings"))

    shiny::updateTabsetPanel(session, "display",selected = "enrTab")
    
  })
  
  output$enrich_table <- DT::renderDataTable({
    gVars$EnrichRes
    
    if(is.null(gVars$EnrichRes))
      return(NULL)
    
    enr_tab <- gVars$EnrichRes[[input$time_point_id_visual3]]
    print(head(enr_tab))
 
    DT::datatable(enr_tab, filter="none",
                  options = list(
                    search = list(regex=TRUE, caseInsensitive=FALSE),
                    scrollX=TRUE,
                    ordering=F
                  )
    )
    
  })

  output$downloadAnovaData <- downloadHandler(
    filename = function() {
      paste("anova_table.xlsx", sep = "")
    },
    content = function(file) {
      if(length(gVars$PValMat)==0){ 
        print("No anova table to save!")
        return(NULL)
      }
      
      shinyjs::html(id="loadingText", "Saving tables")
      shinyjs::show(id="loading-content")
      on.exit({
        print("inside on exit")
        shinyjs::hide(id="loading-content", anim=TRUE, animType="fade")    
      })
      
      print("I'm saving anova tables")
      
      firstCall = TRUE
      
      for(exp in names(gVars$PValMat)){
        print(exp)
        times = names(gVars$PValMat[[exp]])
        for(tp in times){
          print(tp)
          if(firstCall){
            if(nrow(gVars$PValMat[[exp]][[tp]])>0){
              write.xlsx(gVars$PValMat[[exp]][[tp]], file, sheetName = paste(exp,tp,sep = "_"), append = FALSE)
              xlcFreeMemory()
              firstCall = FALSE
              print("first")
            }
          }else{
            if(nrow(gVars$PValMat[[exp]][[tp]])>0){
              write.xlsx(gVars$PValMat[[exp]][[tp]], file, sheetName = paste(exp,tp,sep = "_"), append = TRUE)
              xlcFreeMemory()
              print("appending")
            }
          }
          
        }
      }
      
      # write.xlsx(gVars$PValMat[[1]], file, sheetName = names(gVars$PValMat)[1]) 
      # 
      # if(length(gVars$PValMat)>1){
      #   for(i in 2:length(gVars$PValMat)){
      #     write.xlsx(gVars$PValMat[[i]], file, sheetName = names(gVars$PValMat)[i], append = TRUE) 
      #   }
      # }
      
      print("Anova table stored!")
    }
  )
  
  output$downloadBMDData <- downloadHandler(
    filename = function() {
      paste("bmd_table.xlsx", sep = "")
    },
    content = function(file) {
      if(length(gVars$MQ_BMD_filtered)==0){ 
        print("No bmd table to save!")
        return(NULL)
      }
      shinyjs::html(id="loadingText", "Saving tables")
      shinyjs::show(id="loading-content")
      on.exit({
        print("inside on exit")
        shinyjs::hide(id="loading-content", anim=TRUE, animType="fade")    
      })
      
      print("I'm saving bmd tables")
      
      # write.xlsx(gVars$MQ_BMD_filtered[[1]]$BMDValues_filtered, file, sheetName = names(gVars$MQ_BMD_filtered)[1]) 
      # if(length(gVars$MQ_BMD_filtered)>1){
      #   for(i in 2:length(gVars$MQ_BMD_filtered)){
      #     if(nrow(gVars$MQ_BMD_filtered[[i]]$BMDValues_filtered)>0){
      #       write.xlsx(gVars$MQ_BMD_filtered[[i]]$BMDValues_filtered, file, sheetName =  names(gVars$MQ_BMD_filtered)[i], append = TRUE) 
      #     }
      #   }
      # }
      
      
      firstCall = TRUE

      for(exp in names(gVars$MQ_BMD_filtered)){
        print(exp)
        times = names(gVars$MQ_BMD_filtered[[exp]])
        for(tp in times){
          print(tp)
          if(firstCall){
            if(nrow(gVars$MQ_BMD_filtered[[exp]][[tp]]$BMDValues_filtered)>0){
              write.xlsx(gVars$MQ_BMD_filtered[[exp]][[tp]]$BMDValues_filtered, file, sheetName = paste(exp,tp,sep = "_"), append = FALSE)
              xlcFreeMemory()
              firstCall = FALSE
              print("first")
            }
          }else{
            if(nrow(gVars$MQ_BMD_filtered[[exp]][[tp]]$BMDValues_filtered)>0){
              write.xlsx(gVars$MQ_BMD_filtered[[exp]][[tp]]$BMDValues_filtered, file, sheetName = paste(exp,tp,sep = "_"), append = TRUE)
              xlcFreeMemory()
              print("appending")
            }
          }

        }
      }
      
      
      print("BMD table stored!")
    }
  )
  
  # output$downloadBMDData <- downloadHandler(
  #   # if(is.null(gVars$BMD_tab)){
  #   #   return(NULL)
  #   # }
  #   
  #   filename = function() {
  #     paste("bmd_table.csv", sep = "")
  #   },
  #   content = function(file) {
  #     write.csv(gVars$BMD_tab, file, row.names = FALSE)
  #   }
  # )
  
  output$Anova_table <- DT::renderDataTable({
    if(is.null(gVars$PValMat))
      return(NULL)
    
    if(is.null(input$time_point_id_visual))
      return(NULL)
    
    print("TP da analizzare:")
    print(input$time_point_id_visual)
    
    Anova_tab <- gVars$PValMat[[input$experiment_anova]][[input$time_point_id_visual]]
    print(head(Anova_tab))
    print(head(Anova_tab))
    
    print(head(Anova_tab))
    gVars$Anova_tab=Anova_tab
    DT::datatable(Anova_tab, filter="top",
                  options = list(
                    search = list(regex=TRUE, caseInsensitive=FALSE),
                    scrollX=TRUE,
                    ordering=T,
                    rowCallback = JS('
                                                  function(nRow, aData, iDisplayIndex, iDisplayIndexFull) {
                                                  // Bold and green cells for conditions
                                                  if (parseFloat(aData[2]) <=0.05)
                                                  $("td:eq(3)", nRow).css("font-weight", "bold");
                                                  }')
                  )
    )
  },server=TRUE) #,sanitize.text.function=function(x){x}
  
  
  # output$AnovaVenn = renderPlot({
  #   if(is.null(gVars$PValMat))
  #     return(NULL)
  #   
  #   print(gVars$var_genes)
  #   gplots::venn(gVars$var_genes,show.plot = TRUE)
  # })
  
  output$anovaPlot <- renderPlot({
    if(is.null(gVars$PValMat))
      return(NULL)
    
    if(length(gVars$PValMat)==0)
      return(NULL)
    

    #Anova_tab <-gVars$PValMat[[1]]
    #Anova_tab <- gVars$PValMat[[input$time_point_id_visual]]
    Anova_tab <- gVars$PValMat[[input$experiment_anova]][[input$time_point_id_visual]]
    
    print("Nomi time points: ")
    print(names(gVars$PValMat))
    
    x <-  c(sum(as.numeric(Anova_tab[,2])<=0.05), sum(as.numeric(Anova_tab[,2])>0.05))
    labels <-  c("p.val <= 0.05","p.val>0.05")
    
    piepercent<- round(100*x/sum(x), 1)
    
    # Plot the chart.
    pie(x, labels = piepercent, main = "Filtering Result",col = rainbow(length(x)))
    legend("topright", c("Significant Genes","Non Significant Genes"), cex = 0.8,
           fill = rainbow(length(x)))
    
  })
  
  
  output$BMDVenn <- renderPlot({
    if(is.null(gVars$MQ_BMD_filtered)){ 
      print("Null BMD")
      return(NULL)
    }
    
    listFV = list()
    for(i in 1:length(gVars$MQ_BMD_filtered)){
      listFV[[names(gVars$MQ_BMD_filtered)[i]]] = gVars$MQ_BMD_filtered[[i]]$BMDValues_filtered[,1]
    }
    venn(listFV) 
  }) #server = TRUE

  output$meanBMD_timepoint = renderPlotly({
    if(is.null(gVars$EnrichDatList)){ 
      print("Null enrichment")
      return(NULL)
    }
    
    PD = c()
    
    ER <- gVars$EnrichDatList # result of enrichment
    Hier = gVars$hierarchy
    
    for(i in 1:length(ER)){
      print(i)
      
      ERi = ER[[i]]
      if(nrow(ERi)==0) next
      
      print(ERi)
      elemn = names(ER)[i]
      exper = strsplit(x = elemn, split = "_")[[1]][1]
      timep = strsplit(x = elemn, split = "_")[[1]][2]
      
      BMDFilMat = gVars$MQ_BMD_filtered[[exper]][[as.character(timep)]]$BMDValues_filtered #filtered BMD genes
      
      if(!is.null(BMDFilMat) & nrow(ERi)>0){
        
        for(npat in 1:nrow(ER[[i]])){
          gi = unlist(strsplit(x = as.character(ERi[npat,2]),split = ","))
          
          goID = ERi[npat,"annID"]
          PatName = ERi[npat,"Description"]
          PatPval = ERi[npat,"pValueAdj"]
          NGenes = length(gi)
          # PatNGenes = ERi[npat,"gID"]
          # PatSize = ERi[npat,]
          
          idx = which(tolower(BMDFilMat[,1]) %in% tolower(gi))
          BMD = as.numeric(BMDFilMat[idx,"BMD"])
          names(BMD) = BMDFilMat[idx,"Gene"]
          bmd = mean(BMD)
          
          levidx = which(Hier[,"ID"] %in% goID)
          levName = as.character(Hier[levidx[1],1])
          
          #PD = rbind(PD, cbind(names(gVars$MQ_BMD_filtered)[i],PatName,PatPval, bmd, NGenes, levName))
          PD = rbind(PD, cbind(exper, timep, PatName,PatPval, bmd, NGenes, levName))
          
        }
      }
    }
    
    colnames(PD) = c("Experiment","TimePoint","FunCategory", "Adj.pval","MeanBMD","NGenes", "LevName")
    PD = as.data.frame(PD)
    PD$TimePoint = factor(PD$TimePoint, level=sort(as.numeric(as.vector(unique(PD$TimePoint)))))
    PD$Adj.pval = as.numeric(as.vector(PD$Adj.pval))
    PD$MeanBMD = as.numeric(as.vector(PD$MeanBMD))
    PD$logPval = -log(PD$Adj.pval)
    PD$logPval = -log(PD$Adj.pval)
    PD$NGenes = as.numeric(as.vector(PD$NGenes))
    gVars$MeanBMDTimePoint = PD
    ggp = ggplot(data=PD, aes(MeanBMD, fill=TimePoint)) + 
      geom_histogram() + facet_grid(. ~ Experiment)
    
    ggplotly(ggp)
    
  })
  
  
  output$pathway_bubble = renderPlotly({
    if(is.null(gVars$EnrichDatList)){ 
      print("Null enrichment")
      return(NULL)
    }
    
    PD = c()
    
    ER <- gVars$EnrichDatList # result of enrichment
    Hier = gVars$hierarchy
    
    for(i in 1:length(ER)){
      print(i)
      
      ERi = ER[[i]]
      if(nrow(ERi)==0) next
      
      print(ERi)
      elemn = names(ER)[i]
      exper = strsplit(x = elemn, split = "_")[[1]][1]
      timep = strsplit(x = elemn, split = "_")[[1]][2]
      
      BMDFilMat = gVars$MQ_BMD_filtered[[exper]][[as.character(timep)]]$BMDValues_filtered #filtered BMD genes
      if(!is.null(BMDFilMat) & nrow(ERi)>0){
   
          for(npat in 1:nrow(ER[[i]])){
          gi = unlist(strsplit(x = as.character(ERi[npat,2]),split = ","))
          
          goID = ERi[npat,"annID"]
          PatName = ERi[npat,"Description"]
          PatPval = ERi[npat,"pValueAdj"]
          NGenes = length(gi)

          idx = which(tolower(BMDFilMat[,1])%in% tolower(gi))
          BMD = as.numeric(BMDFilMat[idx,"BMD"])
          names(BMD) = BMDFilMat[idx,"Gene"]
          bmd = mean(BMD)
          
          levidx = which(Hier[,"ID"] %in% goID)
          levName = as.character(Hier[levidx[1],1])
          
          PD = rbind(PD, cbind(exper, timep,PatName,PatPval, bmd, NGenes, levName))
        }
    }
    }
    
    colnames(PD) = c("Experiment","TimePoint","FunCategory", "Adj.pval","MeanBMD","NGenes", "LevName")
    PD = as.data.frame(PD)
    PD$TimePoint = factor(PD$TimePoint, level=sort(as.numeric(as.vector(unique(PD$TimePoint)))))
    PD$Adj.pval = as.numeric(as.vector(PD$Adj.pval))
    PD$MeanBMD = as.numeric(as.vector(PD$MeanBMD))
    PD$logPval = -log(PD$Adj.pval)
    PD$logPval = -log(PD$Adj.pval)
    PD$NGenes = as.numeric(as.vector(PD$NGenes))
    gVars$MeanBMDTimePoint = PD
    #save(PD, file = "../pathway_bubble.RData")
    
    p <- PD %>%
     ggplot(aes(MeanBMD,logPval, label = FunCategory, color = LevName, size = NGenes)) + # color=FunCategory
      geom_point() +
      scale_x_log10() +
      theme_bw() + labs(x = "Mean BMD", y = "-log(P.Value)") +
      facet_grid(Experiment ~ TimePoint)
    
    ggplotly(p)
  })
  
  
  
  output$path_bmd_dist = renderPlotly({
    # shiny::validate(
    #   need(is.null(input$BMD_table_rows_selected), "Select at least one row of the matrix!")
    # )
    
    if(length(input$PAT_table_rows_selected) == 0){
      return(NULL)
    }
    
    print("xyz")
    
    selectedrowindex = input$PAT_table_rows_selected[length(input$PAT_table_rows_selected)]
    selectedrowindex = as.numeric(selectedrowindex)
    
    ER = gVars$EnrichDatList[[input$time_point_id_visualPat]]
    PatName = ER[selectedrowindex,"annID"]
    PatGenes = ER[selectedrowindex,"gID"]
    
    exp_tp = unlist(strsplit(input$time_point_id_visualPat,"_"))
    
    BMDFilMat = gVars$MQ_BMD_filtered[[exp_tp[1]]][[exp_tp[2]]]$BMDValues_filtered
    
    
    #BMDFilMat = gVars$MQ_BMD_filtered[[input$time_point_id_visualPat]]$BMDValues_filtered
    gi = unlist(strsplit(x = ER[input$PAT_table_rows_selected,2],split = ","))
    idx = which(tolower(BMDFilMat[,1])%in% tolower(gi))
    BMD = as.numeric(BMDFilMat[idx,"BMD"])
    BMDL = as.numeric(BMDFilMat[idx,"BMDL"])
    BMDU = as.numeric(as.vector(BMDFilMat[idx,"BMDU"]))
    
    names(BMD) = BMDFilMat[idx,"Gene"]
    names(BMDL) = names(BMDU) = names(BMD)
    print("plotto un pathwayyyyyyyy - >>>>>>>>>> ")
    print(BMD)
    #save(selectedrowindex,ER,BMDFilMat, gi, idx, BMD, BMDL, BMDU, file = "../path_plot.RData")
    
    #BMD = sort(BMD, decreasing = F)
    
    BMD = data.frame(gene = names(BMD), bmd = BMD, bmdl = BMDL, bmdu = BMDU)
    BMD = BMD[order(BMD$bmd),]
    BMD$gene = factor(x = BMD$gene, levels = BMD$gene)
    
    p = ggplot(data=BMD, aes(x=gene, y=bmd, group=1, label1 = bmdl, label2 = bmdu)) +
      geom_line()+
      geom_point() +
      geom_ribbon(aes(ymin=BMD$bmdl, ymax=BMD$bmdu), linetype=2, alpha=0.1) +
      labs(y = "BMDL - BMD - BMDU", x = "Gene")
      
    ggplotly(p)
    #plot(BMD, type = "l")
    
  })
  
  output$PAT_table = DT::renderDataTable({
    # if(is.null(gVars$MQ_BMD_filtered)){ 
    #   print("Null BMD")
    #   return(NULL)
    # }
    print("ciao")
    if((input$time_point_id_visualPat %in% names(gVars$EnrichDatList))==FALSE){
      print("No enrichment for this TP")
      return(NULL)
    }

    PATWAY_tab <- gVars$EnrichDatList[[input$time_point_id_visualPat]]
    PATWAY_tab = PATWAY_tab[,c(1,3,4,5)]
    DT::datatable(PATWAY_tab, filter="top",
                  options = list(
                    search = list(regex=TRUE, caseInsensitive=FALSE),
                    scrollX=TRUE,
                    ordering=T
                  )
    )
  })
  
  

  output$BMD_BMDL = renderPlotly({
    if(is.null(gVars$MQ_BMD_filtered)){ 
      print("Null BMD")
      return(NULL)
    }
    
    BMD_tab <- gVars$MQ_BMD_filtered[[input$experiment_bmd]]#[[input$time_point_id_visual2]]$BMDValues_filtered
    
    #BMD_tab <- gVars$MQ_BMD_filtered#[[input$time_point_id_visual2]]$BMDValues_filtered
    #save(BMD_tab, file="BMD_table.RData")
    
    DF = c()
    for(i in 1:length(BMD_tab)){ # for each time point
      print(BMD_tab[[i]]$BMDValues_filtered)
      if(!is.null(BMD_tab[[i]]$BMDValues_filtered)){
        if(nrow(BMD_tab[[i]]$BMDValues_filtered)>=1){
          DF = rbind(DF, cbind(names(BMD_tab)[i],
                           BMD_tab[[i]]$BMDValues_filtered[,"BMD"],
                           BMD_tab[[i]]$BMDValues_filtered[,"BMDL"],
                           as.character(BMD_tab[[i]]$BMDValues_filtered[,"MOD_NAME"])))
        }
      }
    }
    colnames(DF) = c("TimePoint", "BMD","BMDL","MOD_NAME")
    print(head(DF))
    DF = as.data.frame(DF)
    
    DF$TimePoint = factor(DF$TimePoint, level=sort(as.numeric(as.vector(unique(DF$TimePoint)))))
    DF$MOD_NAME = as.factor( DF$MOD_NAME)
    DF$BMD = as.numeric(as.vector(DF$BMD))
    DF$BMDL = as.numeric(as.vector(DF$BMDL))
    
    #save(DF, file = "../BMD_BMDL.RData")
    
    p = ggplot(DF, aes(x=BMDL, y=BMD, color = MOD_NAME)) +
      geom_point(shape=1) +  facet_grid(. ~ TimePoint)
    
    ggplotly(p)
    # ggplot(data=DF, aes(Model)) + 
    #   geom_histogram()+  facet_grid(. ~ TimePoint)
    
  })
  
  output$BMD_dist_models = renderPlot({
    if(is.null(gVars$MQ_BMD_filtered)){ 
      print("Null BMD")
      return(NULL)
    }
    
    BMD_tab <- gVars$MQ_BMD_filtered[[input$experiment_bmd]]#[[input$time_point_id_visual2]]$BMDValues_filtered
    
    DF = c()
    for(i in 1:length(BMD_tab)){ # for each time point
      print(BMD_tab[[i]]$BMDValues_filtered)
      if(!is.null(BMD_tab[[i]]$BMDValues_filtered)){
        if(nrow((BMD_tab[[i]]$BMDValues_filtered))>=1){
          xx = cbind(names(BMD_tab)[i],as.character(BMD_tab[[i]]$BMDValues_filtered[,"MOD_NAME"]))
          print(xx)
          DF = rbind(DF, xx)
        }

      }
    }
    colnames(DF) = c("TimePoint", "Model")
    DF = as.data.frame(DF)
    DF$TimePoint = factor(DF$TimePoint, level=sort(as.numeric(as.vector(unique(DF$TimePoint)))))
    DF = melt(table(DF$TimePoint, DF$Model))
    colnames(DF) = c("TimePoint", "Model", "Freq")
    DF$TimePoint = factor(DF$TimePoint, level=sort(as.numeric(as.vector(unique(DF$TimePoint)))))
    
    mm = unique(DF$TimePoint)
    for(i in mm){
      iid = which(DF$TimePoint %in% i)
      DF$Freq[iid] = DF$Freq[iid]/sum(DF$Freq[iid]) * 100
    }
    
    #save(DF, file = "../mod_hist.RData")
    
    print(head(DF))
    ggplot(data=DF, aes(x = "", y = Freq, fill = Model)) + geom_bar(width = 1, stat = "identity") + 
      coord_polar("y", start=0)    +  facet_grid(. ~ TimePoint)
      #geom_histogram()+  facet_grid(. ~ TimePoint)
    
  })
  
  output$BMD_BMDL_BMDU_by_model = renderPlotly({
    if(is.null(gVars$MQ_BMD_filtered)){ 
      print("Null BMD")
      return(NULL)
    }
    
    
    BMD_tab <- gVars$MQ_BMD_filtered[[input$experiment_bmd]]#[[input$time_point_id_visual2]]$BMDValues_filtered

    DF = c()
    for(i in 1:length(BMD_tab)){ # for each time point
      print(BMD_tab[[i]]$BMDValues_filtered)
      if(!is.null(BMD_tab[[i]]$BMDValues_filtered)){
        if(nrow(BMD_tab[[i]]$BMDValues_filtered)>=1){
          
          mi = rbind(cbind(names(BMD_tab)[i],BMD_tab[[i]]$BMDValues_filtered[,"BMD"],   as.character(BMD_tab[[i]]$BMDValues_filtered[,"MOD_NAME"]),"BMD"),
                     cbind(names(BMD_tab)[i],BMD_tab[[i]]$BMDValues_filtered[,"BMDL"],  as.character(BMD_tab[[i]]$BMDValues_filtered[,"MOD_NAME"]),"BMDL"),
                     cbind(names(BMD_tab)[i],BMD_tab[[i]]$BMDValues_filtered[,"BMDU"],  as.character(BMD_tab[[i]]$BMDValues_filtered[,"MOD_NAME"]), "BMDU"))
          
          DF = rbind(DF, mi)
        }
      }
    }
    colnames(DF) = c("TimePoint", "Value","Model","ValueType")
    print(head(DF))
    DF = as.data.frame(DF)
    
    DF$TimePoint = factor(DF$TimePoint, level=sort(as.numeric(as.vector(unique(DF$TimePoint)))))
    DF$ValueType = factor(DF$ValueType, level=c("BMDL","BMD","BMDU"))
    
    DF$Value = as.numeric(as.vector(DF$Value))

    p = ggplot(DF, aes(x=Model, y=Value, fill = ValueType)) + geom_boxplot(width=0.5) +facet_grid(. ~ TimePoint)
    ggplotly(p) %>%layout(boxmode = "group")
    
  })
  
  output$BMD_dist_TP = renderPlotly({
    if(is.null(gVars$MQ_BMD_filtered)){ 
      print("Null BMD")
      return(NULL)
    }
    print("in BMD_dist_TP --------------------->>>>>>>>>>>>>>>>>>>>>")
    BMD_tab <- gVars$MQ_BMD_filtered[[input$experiment_bmd]]#[[input$time_point_id_visual2]]$BMDValues_filtered

    DF = c()
    for(i in 1:length(BMD_tab)){
      print(BMD_tab[[i]]$BMDValues_filtered)
      
      if(!is.null(BMD_tab[[i]]$BMDValues_filtered)){
        if(nrow(BMD_tab[[i]]$BMDValues_filtered)>=1){
          xx = cbind(names(BMD_tab)[i],BMD_tab[[i]]$BMDValues_filtered[,"BMD"])
          print(xx)
          DF = rbind(DF, xx)
        }
      }
    }
    
    print(head(DF))
    colnames(DF) = c("TimePoint", "BMD")
    DF = as.data.frame(DF)
    DF$TimePoint = factor(DF$TimePoint, level=sort(as.numeric(as.vector(unique(DF$TimePoint)))))
    DF$BMD = as.numeric(as.vector(DF$BMD))
    # toRem = union(which(DF$BMD< as.numeric(input$bmd_th_min)), which(DF$BMD > as.numeric(input$bmd_th_max)))
    # 
    # if(length(toRem)>0){
    #   DF = DF[-toRem,]
    # }
    
    DF$BMD = as.numeric(as.vector(DF$BMD))
    #save(DF, file = "../BMD_hist.RData")
    print("in DF --------------------->>>>>>>>>>>>>>>>>>>>>")
    print(head(DF))
    
    p = ggplot(data=DF, aes(x = BMD)) + 
        geom_histogram(aes(y = ..density..)) + 
        geom_density(alpha = .2, fill="#FF6655") +
        facet_grid(. ~ TimePoint)

    ggplotly(p)
  })
  
  output$BMD_pval_fitting = renderPlotly({
    if(is.null(gVars$MQ_BMD_filtered)){ 
      print("Null BMD")
      return(NULL)
    }
    BMD_tab <- gVars$MQ_BMD_filtered[[input$experiment_bmd]]#[[input$time_point_id_visual2]]$BMDValues_filtered
    
    DF = c()
    for(i in 1:length(BMD_tab)){
      if(!is.null(BMD_tab[[i]]$BMDValues_filtered)){
        if(nrow(BMD_tab[[i]]$BMDValues_filtered)>=1){
          DF = rbind(DF, cbind(names(BMD_tab)[i],BMD_tab[[i]]$BMDValues_filtered[,"LOFPVal"]))
        }
      }
    }
    colnames(DF) = c("TimePoint", "LOFPVal")
    DF = as.data.frame(DF)
    DF$TimePoint = factor(DF$TimePoint, level=sort(as.numeric(as.vector(unique(DF$TimePoint)))))
    DF$LOFPVal = as.numeric(as.vector(DF$LOFPVal))
    
    p = ggplot(data=DF, aes(x = LOFPVal)) + 
      geom_histogram(aes(y = ..density..)) + 
      geom_density(alpha = .2, fill="#FF6655") +
      facet_grid(. ~ TimePoint)
    
    ggplotly(p)
  })
    
  output$NGTime = renderPlotly({
    if(is.null(gVars$MQ_BMD_filtered)){ 
      print("Null BMD")
      return(NULL)
    }
    
    BMD_tab <- gVars$MQ_BMD_filtered[[input$experiment_bmd]]
    
    GL = list()
    NG = c()
    for(i in 1:length(BMD_tab)){
      BMD_tab2 <- BMD_tab[[i]]$BMDValues_filtered
      if(!is.null(BMD_tab2) & nrow(BMD_tab2)>0){
        GL[[names(BMD_tab)[i]]] =  BMD_tab2[,"Gene"]
        NG = rbind(NG, c(names(BMD_tab)[i], nrow(BMD_tab2)))
      }else{
        NG = rbind(NG, c(names(BMD_tab)[i], 0))
      }
    }
    
    colnames(NG) = c("TimePoint","NGenes")
    NG = as.data.frame(NG)
    NG$NGenes = as.numeric(as.vector(NG$NGenes))
    NG$TimePoint = factor(as.character(NG$TimePoint),levels = sort(as.numeric(as.vector(NG$TimePoint))))
    
    p<-ggplot(data=NG, aes(x=TimePoint, y=NGenes)) +
      geom_bar(stat="identity")
    
    ggplotly(p)
    
  })
  
  
  output$gene_bmd_plot = renderPlotly({
    if(is.null(gVars$MQ_BMD_filtered)){ return(NULL) }
    if(is.null(input$geneClust)) {return(NULL)}
    if(is.null(input$intersectionName)) {return(NULL)}
    
    print("gene_bmd_plot")
    BMD_tab <- gVars$MQ_BMD_filtered[[input$experiment_bmd]]

    if(length(BMD_tab)>1){
      
      c(GL, ItemsList) %<-%  venn_diagram_bmd_genes_across_time_point(BMD_tab)
      print("xxx")
      c(XX, hls, BMD_gene) %<-% create_gene_bmd_dataframe_and_cluster_genes_by_bmd(bmd_list=BMD_tab,
                                                                                   ItemsList=ItemsList,
                                                                                   hmethod=input$geneClustMeth, 
                                                                                   nclust=as.numeric(input$geneClust),
                                                                                   intersectionName=input$intersectionName)
  
      if(unique(XX$Cluster==0)){
        p = ggplot(data=XX, aes(x=BMD, y=Gene, color = BMD)) +
          geom_point() 
      }else{
        p = ggplot(data=XX, aes(x=TimePoint, y=BMD, group=Gene, color = Gene)) +
          geom_line(linetype = "dashed")+
          geom_point() + facet_grid(~Cluster)
      }
      
      ggplotly(p)
    }
  })
  
  
  output$nclustvenn <- renderUI({
    if(is.null(gVars$MQ_BMD_filtered)){ return(NULL) }
    print("nclustvenn")
    
    BMD_tab <- gVars$MQ_BMD_filtered[[input$experiment_bmd]]
    
    
    if(length(BMD_tab)>1){
      c(GL, ItemsList) %<-%  venn_diagram_bmd_genes_across_time_point(BMD_tab)
      
      if(is.null(input$input$intersectionName)){
        ng = length(ItemsList[[length(ItemsList)]])
      }else{
        ng = length(ItemsList[[input$input$intersectionName]])
      }
      #
      selectInput(inputId = "geneClust", label = "Number of clusters",choices = as.list(1:ng), selected = 1)
      
    }
    
  })
  
  output$intersectionNameUI <- renderUI({
    if(is.null(gVars$MQ_BMD_filtered)){ return(NULL) }
    print("intersectionNameUI")
    
    BMD_tab <- gVars$MQ_BMD_filtered[[input$experiment_bmd]]
    
    
    if(length(BMD_tab)>1){
      c(GL, ItemsList) %<-%  venn_diagram_bmd_genes_across_time_point(BMD_tab)
      
      ng = length(ItemsList[[length(ItemsList)]])
      selectInput(inputId = "intersectionName", label = "Select Gene Set",choices = as.list(names(ItemsList)), selected = names(ItemsList)[[length(ItemsList)]])
      
    }
})
  

  
  
  output$VennDF = DT::renderDataTable({
    if(is.null(gVars$MQ_BMD_filtered)){ return(NULL) }
    
    BMD_tab <- gVars$MQ_BMD_filtered[[input$experiment_bmd]]
    
    
    if(length(BMD_tab)>1){
      c(GL, ItemsList) %<-%  venn_diagram_bmd_genes_across_time_point(BMD_tab)
      gVars$ItemsList = ItemsList
      
      df = build_dataframe_from_genes_in_venn_diagrams(ItemsList)
      gVars$venn_genes_df = df
      
      DT::datatable(df, filter="top",
                    options = list(
                      search = list(regex=TRUE, caseInsensitive=FALSE),
                      scrollX=TRUE,
                      ordering=T,
                      lengthMenu = c(5, 10, 15, 20, 40, 60, 100)
                    )
      )
    }

  })

  output$NGVenn = renderPlot({
    if(is.null(gVars$MQ_BMD_filtered)){ return(NULL) }
    print("NGVenn")
    
    BMD_tab <- gVars$MQ_BMD_filtered[[input$experiment_bmd]]
    
    
    if(length(BMD_tab)>1){
      c(GL, ItemsList) %<-%  venn_diagram_bmd_genes_across_time_point(BMD_tab)
      gVars$ItemsList = ItemsList
      gVars$GL = GL
      par(mar = c(0,0,0,0))
      venn(GL)
    }
    
  })
  
  output$timePointSelTab <- renderUI({
    shiny::validate(
      need(!is.null(gVars$EnrichDatList), "No Enrichment!")
    )
    
    print("Rending UI: ")
    print(names(gVars$EnrichDatList))
    print(min(as.numeric(names(gVars$EnrichDatList))))
    
    selectInput("time_point_id_table", "Time Points", choices=names(gVars$EnrichDatList), selected = names(gVars$EnrichDatList)[1])#c("All",unique(gVars$phTable[,gVars$TPColID])))
  })
  
  output$PatTable = DT::renderDataTable({
    shiny::validate(
      need(!is.null(gVars$EnrichDatList), "No Enrichment!")
    )
    
    ER <- gVars$EnrichDatList[[input$time_point_id_table]]
    ER = ER[,c(5,1:4)]
    DT::datatable(ER, filter="top",
                  options = list(
                    search = list(regex=TRUE, caseInsensitive=FALSE),
                    scrollX=TRUE,
                    ordering=T
                  )
    )
    
  })
  

  
  output$downloadEnrichedPathwayTables = downloadHandler(
    filename = function() {
      paste("terms_enrichment_table.xlsx", sep = "")
    },
    content = function(file) {
      if(length(gVars$EnrichDatList)==0){ 
        print("No enrichment tables to save!")
        return(NULL)
      }
      
      shinyjs::html(id="loadingText", "Saving tables")
      shinyjs::show(id="loading-content")
      on.exit({
        print("inside on exit")
        shinyjs::hide(id="loading-content", anim=TRUE, animType="fade")    
      })
      
      print("I'm saving enrichment tables")
      for(i in 1:length(gVars$EnrichDatList)){
        if(nrow(gVars$EnrichDatList[[i]])>0){
          write.xlsx(gVars$EnrichDatList[[1]], file, sheetName = names(gVars$EnrichDatList)[1]) 
          startI = i
          break
        }
      }
      
      if(length(gVars$PValMat)>1){
        for(i in (startI+1):length(gVars$EnrichDatList)){
          if(nrow(gVars$EnrichDatList[[i]])>0){
            write.xlsx(gVars$EnrichDatList[[i]], file, sheetName =  names(gVars$EnrichDatList)[i], append = TRUE) 
          }
        }
      }
      print("Enrichment table stored!")
    }
  )
  
  output$ExptimeSelBMD <- renderUI({
    if(is.null(gVars$EXP_FIL)){
      return(NULL)
    }
    print("Rending UI: ")
    print(names(gVars$EXP_FIL))
    print(min(as.numeric(names(gVars$EXP_FIL))))
    
    selectInput("experiment_bmd", "Experiment", choices=names(gVars$EXP_FIL), selected = names(gVars$EXP_FIL)[1])#c("All",unique(gVars$phTable[,gVars$TPColID])))
  })
  
  output$timePointSel2 <- renderUI({
    if(is.null(input$experiment_bmd)){
      return(NULL)
    }
    
    print("Rending UI: ")
    print(input$experiment_bmd)
    #print(min(as.numeric(names(gVars$EXP_FIL))))
    
    expList = gVars$EXP_FIL[[input$experiment_bmd]]
    
    selectInput("time_point_id_visual2", "Time Points", choices=names(expList), selected = min(as.numeric(names(expList))))#c("All",unique(gVars$phTable[,gVars$TPColID])))
  })
  
  
  
  output$BMD_table <- DT::renderDataTable({
    if(is.null(gVars$MQ_BMD_filtered)){ 
      print("Null BMD")
      return(NULL)
    }
    print("BMD_tab")
 
    gVars$MQ_BMD_filtered = gVars$MQ_BMD_filtered2
    
    
    BMD_tab <- gVars$MQ_BMD_filtered[[input$experiment_bmd]][[input$time_point_id_visual2]]$BMDValues_filtered
  
    BMD_tab = BMD_tab#[,1:5]
    print("printing BMD table")
    print(head(BMD_tab))
    gVars$BMD_tab = BMD_tab
    
    BMD_tab[,"BMD"]  =  as.numeric(as.vector(BMD_tab[,"BMD"]))
    BMD_tab[,"BMDL"] =  as.numeric(as.vector(BMD_tab[,"BMDL"]))
    BMD_tab[,"BMDU"] =  as.numeric(as.vector(BMD_tab[,"BMDU"]))
        
    #save(BMD_tab, file="../BMD_tab.RData")
    
    #DT::datatable(BMD_tab,options = list(pageLength = 5))
    
    BMD_tab
    
  },selection = "single", server = TRUE)
  
  # remember to remove
  #load(file="../BMD_tab.RData")
  # mtcars2 = mtcars[, 1:8]
  # output$x3 = DT::renderDataTable( DT::datatable(BMD_tab), server = TRUE, selection = "single")
  # 
  # # print the selected indices
  # output$x4 = renderPrint({
  #   s = input$x3_rows_selected
  #   if (length(s)) {
  #     cat('These rows were selected:\n\n')
  #     cat(s, sep = ', ')
  #   }
  # })
  
  
  output$bmd_fitting = renderPlot({
    if(length(input$BMD_table_rows_selected) == 0){
      return(NULL)
    }
    
    print("test")
    print(input$BMD_table_rows_selected)
    
    selectedrowindex = input$BMD_table_rows_selected[length(input$BMD_table_rows_selected)]
    print(length(input$BMD_table_rows_selected))
    selectedrowindex = as.numeric(selectedrowindex)
    print(selectedrowindex)
    
    BMDVF = gVars$MQ_BMD_filtered[[input$experiment_bmd]][[input$time_point_id_visual2]]$BMDValues_filtered
    
    geneName = as.character(gVars$MQ_BMD_filtered[[input$experiment_bmd]][[input$time_point_id_visual2]]$BMDValues_filtered[selectedrowindex,"Gene"])
    geneBMD = as.numeric(gVars$MQ_BMD_filtered[[input$experiment_bmd]][[input$time_point_id_visual2]]$BMDValues_filtered[selectedrowindex,"BMD"])
    geneBMDL = as.numeric(gVars$MQ_BMD_filtered[[input$experiment_bmd]][[input$time_point_id_visual2]]$BMDValues_filtered[selectedrowindex,"BMDL"])
    geneBMDU = as.numeric(as.vector(gVars$MQ_BMD_filtered[[input$experiment_bmd]][[input$time_point_id_visual2]]$BMDValues_filtered[selectedrowindex,"BMDU"]))
    modName = as.character(gVars$MQ_BMD_filtered[[input$experiment_bmd]][[input$time_point_id_visual2]]$BMDValues_filtered[selectedrowindex,"MOD_NAME"])
    icval = as.numeric(as.vector(gVars$MQ_BMD_filtered[[input$experiment_bmd]][[input$time_point_id_visual2]]$BMDValues_filtered[selectedrowindex,"IC50/EC50"]))
    decreasing = as.numeric(as.vector(gVars$MQ_BMD_filtered[[input$experiment_bmd]][[input$time_point_id_visual2]]$BMDValues_filtered[selectedrowindex,"Decreasing"]))
    
    print(paste("Selected rows --> ", selectedrowindex))
    print(paste("Selected rows gene name--> ", geneName))
    
    #NB: in all three lines I was using 1 before, instead of input$time_point_id_visual
    print(names(gVars$MQ_BMD[[input$experiment_bmd]][[input$time_point_id_visual2]]$opt_models_list))
    print(geneName %in% names(gVars$MQ_BMD[[input$experiment_bmd]][[input$time_point_id_visual2]]$opt_models_list))
    
    #Select mod to plot
    modToPlot = gVars$MQ_BMD[[input$experiment_bmd]][[input$time_point_id_visual2]]$opt_models_list[[geneName]]$opt_mod
   # modToPlot = fm$CD163$opt_mod
    print("Class optimal model ------------------>>>>>>>>>>>>>>> ")
    class(modToPlot)
    
    print("Saving plot file")
    print("END Saving plot file")

    print(geneBMD)
    print(geneBMDL)
    print(geneBMDU)
    

    if(modName %in% c("Linear", "Quadratic", "Cubic", 
                      "Power2","Power3","Power4","Exponential",
                      "Hill05","Hill1","Hill2","Hill3","Hill4","Hill5")){
      
        #save(modToPlot,geneBMD, geneBMDL, geneBMDU, file = "../linear_mod_to_plot.RData")
        ggp = effect_plot(model = modToPlot, pred = dose, interval = TRUE, plot.points = TRUE) +
        scale_y_continuous(limits = c( min(modToPlot$model[,1]),  max(modToPlot$model[,1]))) +
        geom_segment(aes(x = geneBMD, y = min(modToPlot$model[,1]), 
                         xend = geneBMD, yend = predict(modToPlot, data.frame(dose = geneBMD))),color = "blue",linetype="dashed") +
        geom_point(aes(x=geneBMD, y=predict(modToPlot, data.frame(dose = geneBMD))), size = 2,colour="blue") +
        geom_segment(aes(x = geneBMDL, y = min(modToPlot$model[,1]), 
                         xend = geneBMDL, yend = predict(modToPlot, data.frame(dose = geneBMDL))),color = "red",linetype="dashed") +
        geom_point(aes(x=geneBMDL, y=predict(modToPlot, data.frame(dose = geneBMDL))), size = 2,colour="red") +
        geom_segment(aes(x = geneBMDU, y = min(modToPlot$model[,1]), 
                         xend = geneBMDU, yend = predict(modToPlot, data.frame(dose = geneBMDU))),color = "green",linetype="dashed") +
        geom_point(aes(x=geneBMDU, y=predict(modToPlot, data.frame(dose = geneBMDU))), size = 2,colour="green") +
        geom_segment(aes(x = icval, y = min(modToPlot$model[,1]), 
                           xend = icval, yend = predict(modToPlot, data.frame(dose = icval))),color = "black",linetype="dashed") +
        geom_point(aes(x=icval, y=predict(modToPlot, data.frame(dose = icval))), size = 2,colour="black") +
        geom_segment(aes(x = 0, predict(modToPlot, data.frame(dose = icval)), 
                           xend = icval, yend = predict(modToPlot, data.frame(dose = icval))),color = "black",linetype="dashed") +
        geom_segment(aes(x = 0, predict(modToPlot, data.frame(dose = geneBMD)), 
                         xend = geneBMD, yend = predict(modToPlot, data.frame(dose = geneBMD))),color = "blue",linetype="dashed") +
        geom_segment(aes(x = 0, predict(modToPlot, data.frame(dose = geneBMDL)), 
                         xend = geneBMDL, yend = predict(modToPlot, data.frame(dose = geneBMDL))),color = "red",linetype="dashed") +
        geom_segment(aes(x = 0, predict(modToPlot, data.frame(dose = geneBMDU)), 
                         xend = geneBMDU, yend = predict(modToPlot, data.frame(dose = geneBMDU))),color = "green",linetype="dashed") +
        scale_colour_manual(name="Values",values=c(BMDL="red", BMD="blue", BMDU="green"))
        
     }else{
      
      demo.fits <- expand.grid(conc=seq(min(modToPlot$data[,1]), max(modToPlot$data[,1]), length=100))
      pm <- predict(modToPlot, newdata=demo.fits, interval="confidence") 
      demo.fits$p <- pm[,1]
      demo.fits$pmin <- pm[,2]
      demo.fits$pmax <- pm[,3]
      
      ggp = ggplot(modToPlot$data[,1:2], aes(x = dose, y = expr)) +
        geom_point() +
        geom_ribbon(data=demo.fits, aes(x=conc,y=p,ymin=pmin, ymax=pmax), alpha=0.2) +
        geom_line(data=demo.fits, aes(x=conc, y=p)) +
        geom_segment(aes(x = geneBMD, y = min(modToPlot$data$expr), 
                         xend = geneBMD, yend = modToPlot$curve[[1]](geneBMD)),color = "blue",linetype="dashed") +
        geom_point(aes(x=geneBMD, y=modToPlot$curve[[1]](geneBMD)), size = 2,colour="blue") +
        geom_segment(aes(x = 0, modToPlot$curve[[1]](geneBMD), 
                         xend = geneBMD, yend = modToPlot$curve[[1]](geneBMD)),color = "blue",linetype="dashed") +
        geom_segment(aes(x = geneBMDL, y = min(modToPlot$data$expr), 
                         xend = geneBMDL, yend = modToPlot$curve[[1]](geneBMDL)),color = "red",linetype="dashed") +
        geom_point(aes(x=geneBMDL, y=modToPlot$curve[[1]](geneBMDL)), size = 2,colour="red") +
        geom_segment(aes(x = 0, modToPlot$curve[[1]](geneBMDL), 
                         xend = geneBMDL, yend = modToPlot$curve[[1]](geneBMDL)),color = "red",linetype="dashed") +
        geom_segment(aes(x = geneBMDU, y = min(modToPlot$data$expr), 
                         xend = geneBMDU, yend = modToPlot$curve[[1]](geneBMDU)),color = "green",linetype="dashed") +
        geom_point(aes(x=geneBMDU, y=modToPlot$curve[[1]](geneBMDU)), size = 2,colour="green") +
        geom_segment(aes(x = 0, modToPlot$curve[[1]](geneBMDU), 
                         xend = geneBMDU, yend = modToPlot$curve[[1]](geneBMDU)),color = "green",linetype="dashed")+
        geom_segment(aes(x = icval, y = min(modToPlot$model[,1]), 
                         xend = geneBMD, yend = predict(modToPlot, data.frame(dose = icval))),color = "black",linetype="dashed") +
        geom_point(aes(x=icval, y=predict(modToPlot, data.frame(dose = icval))), size = 2,colour="black") +
        geom_segment(aes(x = 0, predict(modToPlot, data.frame(dose = icval)), 
                         xend = icval, yend = predict(modToPlot, data.frame(dose = icval))),color = "black",linetype="dashed") +
        theme_bw()+ theme(panel.border = element_blank()) +
        theme(panel.grid.major = element_line(colour = "black",linetype="dashed",size=0.1)) +
        theme(panel.grid.minor = element_line(colour = "black",linetype="dashed",size=0.1)) +
        scale_colour_manual(name="Values",
                            values=c(BMDL="red", BMD="blue", BMDU="green"))
     }
    ggp

  })
  
  output$bmd_fitting_legend = renderPlot({
    if(length(input$BMD_table_rows_selected) == 0){
      return(NULL)
    }
    selectedrowindex = input$BMD_table_rows_selected[length(input$BMD_table_rows_selected)]
    selectedrowindex = as.numeric(selectedrowindex)

    par(mar = c(0,0,0,0))
    plot(1, type="n", xlab="", ylab="", xlim=c(0, 2), ylim=c(0, 2), xaxt = "n", yaxt = "n", frame.plot = FALSE)
    legend("center",  legend=c("BMD","BMDL","BMDU","IC50/EC50"), fill = c("blue","red","green","black"), ncol=4, bty = "n")
  })

  
  observeEvent(input$upload_pheno_submit, {
    shiny::validate(
      need(!is.null(gVars$inputPh()), "No Phenotype File Provided!")
    )
    
    print("upload")
    
    phTable <- gVars$inputPh()
    doseColID <- NULL
    TPColID <- NULL
    sampleColID <- as.integer(input$sampleIDCol)
    gVars$sampleColID = sampleColID # setting sample ID column in the variable list
    
    doseColID <- as.integer(input$doseCol)
    gVars$doseColID = doseColID # setting dose ID column in the variable list
    
    TPColID <- as.integer(input$TPCol)
    gVars$TPColID = TPColID # setting time point ID column in the variable list
    
    for(i in 1:length(phTable)){
      phTable[[i]] = phTable[[i]][order(phTable[[i]][,doseColID]),]
    }
    

    #fileNameColID <- as.integer(input$fileNameCol)
    doseColName <- NULL
    TPColName <- NULL
    sampleColName <- lapply(phTable, FUN = function(elem){colnames(elem)[as.integer(input$sampleIDCol)]})
    sampleIDs = lapply(1:length(phTable), FUN = function(i)phTable[[i]][sampleColName[[i]]])
    
    for(i in 1:length(phTable)){
      if(any(duplicated(sampleIDs[[i]]))){
        shinyjs::info(paste0("Sample IDs are not unique!\n\nPlease check the phenotype data and ensure that the selected Sample ID column has unique values!"))
        return(NULL)
      }else if(any(sampleIDs[[i]]=="") || any(sampleIDs[[i]]==" ") || any(is.na(sampleIDs[[i]]))){
        shinyjs::info(paste0("Sample IDs contain BLANK and/or NA values!\n\nPlease check the phenotype data and ensure that the selected Sample ID column has complete information!"))
        return(NULL)
      }
    }
    
    for(i in 1:length(phTable)){
      #Make sample ID column character
      if(is.integer(phTable[[i]][,sampleColName[[i]]])){
        print("Updating sample ID column from integer to character!")
        phTable[[i]][,sampleColName[[i]]] <- as.character(phTable[[i]][,sampleColName[[i]]])
      }
    }

    
    #Create factorized phTable
    print("Create factorized phTable")
    print(is.null(input$phenoTypesRH))
    print(str(input$phenoTypesRH))
    phRH <- rhandsontable::hot_to_r(input$phenoTypesRH)
    print(str(phRH))
    factorIdx <- which(phRH$Type=="factor")
    print(factorIdx)
    factorCols <- NULL
    if(length(factorIdx)>0){
      factorCols <- lapply(phTable, FUN = function(elem) colnames(elem)[factorIdx])
      phFactor <- lapply(1:length(phTable), FUN = function(i) factorize_cols(phTable[[i]], factorIdx))
    }else{
      phFactor <- phTable
    }
    print(str(phFactor))
    print(str(phTable))
    
    gVars$phTable <- phTable
    gVars$phFactor <- phFactor
    gVars$factorCols <- factorCols
    gVars$doseColID <- doseColID
    gVars$TPColID <- TPColID
    #gVars$fileNameColID <- fileNameColID
    gVars$sampleColID <- sampleColID
    gVars$doseColName <- doseColName
    gVars$TPColName <- TPColName
    #gVars$fileNameColName <- fileNameColName
    gVars$sampleColName <- sampleColName
    gVars$totalSamples <- nrow(phTable)
    gVars$filteredSamples <- nrow(phTable)
    gVars$removedSamples <- 0
    gVars$removedSamplesInfo <- NULL
    
    shinyBS::toggleModal(session, "importPhenoModal", toggle="close")
    shinyBS::updateButton(session, "import_pheno_submit", style="success", icon=icon("check-circle"))
    shinyBS::updateCollapse(session, "bsSidebar1", open="LOAD EXPRESSION MATRIX", style=list("LOAD PHENOTYPE DATA"="success","LOAD EXPRESSION MATRIX"="warning"))
    print("open expMat")

    
  })
  
  
  # output$kpi_summary_box_1 <- renderValueBox({
  #   valueBox(
  #     value = sprintf("%s", compress(245923)),
  #     subtitle = sprintf("KPI 1 (%.1f%%)", 8.9202),
  #     icon = icon("arrow-up"),
  #     color = "green"
  #   )
  # })
  
  gVars$pcChoices <- reactive({
    phTable <- gVars$phTable
    if (is.null(phTable))
      return(c("NA"))
    
    colNames <- colnames(phTable)[apply(phTable, 2, function(x){res<-any(is.na(x));!res})]
    if(length(colNames)==0){
      shinyjs::info("No available phenotype variables.\nAll variables contain NAs\nCheck 'Keep NA Variables' to include variables with NA values")
      return(c("NA"))
    }
    return(colNames)
  })
  
  gVars$pcChoicesSV <- reactive({
    phTable <- gVars$phTable
    if (is.null(phTable))
      return(c("NA"))
    
    if(!is.null(gVars$svaSV)){
      svaSV <- as.data.frame(apply(gVars$svaSV, 2, factor))
      #svaSVc <- as.data.frame(apply(gVars$svaSVc, 2, factor))
      svaSVc <- gVars$svaSVc
      phTable <- cbind(phTable, svaSV, svaSVc)
    }
    
    colNames <- colnames(phTable)[apply(phTable, 2, function(x){res<-any(is.na(x));!res})]
    if(length(colNames)==0){
      shinyjs::info("No available phenotype variables.\nAll variables contain NAs\nCheck 'Keep NA Variables' to include variables with NA values")
      return(c("NA"))
    }
    return(colNames)
  })
  
  gVars$condChoices <- reactive({
    if(is.null(gVars$phTable))
      return(c("NA"))
    
    varI <- input$varILimma
    if (is.null(varI) || is.na(varI))
      return(c("NA"))
    
    phTable <- gVars$phTable
    conds <- levels(factor(phTable[,varI]))
    return(conds)
  })
  
  observeEvent(input$filterPh, {
    if(is.null(gVars$phTable) || is.null(input$filtered_rows_selected))
      return(NULL)
    
    doseColID <- gVars$doseColID
    TPColID <- gVars$TPColID
    fileNameColID <- gVars$fileNameColID
    sampleColID <- gVars$sampleColID
    doseColName <- gVars$doseColName
    TPColName <- gVars$TPColName
    fileNameColName <- gVars$fileNameColName
    sampleColName <- gVars$sampleColName
    phTable <- gVars$phTable
    phRowsSel <- as.integer(input$filtered_rows_selected)
    print("Sel Rows:")
    print(phRowsSel)
    phRows2Remove <- NULL
    removedSampleIDs <- NULL
    if (length(phRowsSel)>0){
      #phRows2Remove <- which(rownames(phTable) %in% phRowsSel)
      phRows2Remove <- phRowsSel
      print("Rows to Remove:")
      print(phRows2Remove)
      removedSampleIDs <- phTable[phRows2Remove,sampleColName]
      removedSamplesInfo <- phTable[phRows2Remove,c(fileNameColName, doseColName, TPColName)]
      phTable <- phTable[-phRows2Remove,]
      cat("Removed Sample IDs :", paste0(removedSampleIDs, collapse=","), "\n")
    }
    
    #Update expression objects
    if(!is.null(gVars$norm.data)){
      print("Subsetting gVars$norm.data")
      rmIdx <- which(colnames(gVars$norm.data) %in% removedSampleIDs)
      #if(length(rmIdx)>1){
      if(length(rmIdx)>0){
        cat("Removing index ", rmIdx, " from gVars$norm.data matrix!\n")
        gVars$norm.data <- gVars$norm.data[,-rmIdx]
      }else{
        print("Nothing was removed!")
      }
    }
    if(!is.null(gVars$rgList)){
      print("Subsetting gVars$rgList")
      #removedArrays <- phTable[phRows2Remove,fileNameColID]
      removedArrays <- phTable[phRows2Remove,fileNameColName]
      print("removedArrays: ")
      print(removedArrays)
      #removedArrays <- phTable[phRows2Remove,fileNameColName, drop=FALSE]
      removedArrays <- gsub("\\.[a-zA-Z]+$", "", removedArrays)
      print("removedArrays: ")
      print(removedArrays)
      print("table(removedArrays): ")
      print(table(removedArrays))
      array2Remove <- names(which(table(removedArrays)>1))
      #if(length(array2Remove)>1){
      if(length(array2Remove)>0){
        gVars$rgList <- gVars$rgList[,-array2Remove]
      }else{
        print("Nothing was removed!")
      }
    }
    if(!is.null(gVars$rgList.norm)){
      print("Subsetting gVars$rgList.norm")
      rIdx <- which(colnames(gVars$rgList.norm$R) %in% removedSampleIDs)
      rmCheck <- 0
      if(length(rIdx)>0){
        cat("Removing index ", rIdx, " from R matrix!\n")
        gVars$rgList.norm$R <- gVars$rgList.norm$R[,-rIdx]
        rmCheck <- 1
      }
      gIdx <- which(colnames(gVars$rgList.norm$G) %in% removedSampleIDs)
      if(length(gIdx)>0){
        cat("Removing index ", gIdx, " from G matrix!\n")
        gVars$rgList.norm$G <- gVars$rgList.norm$G[,-gIdx]
        rmCheck <- 1
      }
      if(rmCheck==0){
        print("Nothing was removed!")
      }
    }
    if(!is.null(gVars$RGset)){
      print("Subsetting gVars$RGset")
      #removedArrays <- phTable[phRows2Remove,fileNameColID]
      removedArrays <- phTable[phRows2Remove,fileNameColName]
      rmIdx <- which(sampleNames(gVars$RGset) %in% removedArrays)
      #if(length(rmIdx)>1){
      if(length(rmIdx)>0){
        cat("Removing index ", rmIdx, " from gVars$RGset matrix!\n")
        gVars$RGset <- gVars$RGset[,-rmIdx]
      }else{
        print("Nothing was removed!")
      }
    }
    if(!is.null(gVars$Mset)){
      print("Subsetting gVars$Mset")
      #removedArrays <- phTable[phRows2Remove,fileNameColID]
      removedArrays <- phTable[phRows2Remove,fileNameColName]
      rmIdx <- which(sampleNames(gVars$Mset) %in% removedArrays)
      #if(length(rmIdx)>1){
      if(length(rmIdx)>0){
        cat("Removing index ", rmIdx, " from gVars$Mset matrix!\n")
        gVars$Mset <- gVars$Mset[,-rmIdx]
      }else{
        print("Nothing was removed!")
      }
    }
    if(!is.null(gVars$expr.data)){
      print("Subsetting gVars$expr.data")
      rmIdx <- which(colnames(gVars$expr.data) %in% removedSampleIDs)
      #if(length(rmIdx)>1){
      if(length(rmIdx)>0){
        cat("Removing index ", rmIdx, " from gVars$expr.data matrix!\n")
        gVars$expr.data <- gVars$expr.data[,-rmIdx]
      }else{
        print("Nothing was removed!")
      }
    }
    if(!is.null(gVars$nc.data)){
      print("Subsetting gVars$nc.data")
      rmIdx <- which(colnames(gVars$nc.data) %in% removedSampleIDs)
      #if(length(rmIdx)>1){
      if(length(rmIdx)>0){
        cat("Removing index ", rmIdx, " from gVars$nc.data matrix!\n")
        gVars$nc.data <- gVars$nc.data[,-rmIdx]
      }else{
        print("Nothing was removed!")
      }
    }
    if(!is.null(gVars$c.data)){
      print("Subsetting gVars$c.data")
      rmIdx <- which(colnames(gVars$c.data) %in% removedSampleIDs)
      #if(length(rmIdx)>1){
      if(length(rmIdx)>0){
        cat("Removing index ", rmIdx, " from gVars$c.data matrix!\n")
        gVars$c.data <- gVars$c.data[,-rmIdx]
      }else{
        print("Nothing was removed!")
      }
    }
    if(!is.null(gVars$comb.data)){
      print("Subsetting gVars$comb.data")
      rmIdx <- which(colnames(gVars$comb.data) %in% removedSampleIDs)
      #if(length(rmIdx)>1){
      if(length(rmIdx)>0){
        cat("Removing index ", rmIdx, " from gVars$comb.data matrix!\n")
        gVars$comb.data <- gVars$comb.data[,-rmIdx]
      }else{
        print("Nothing was removed!")
      }
    }
    if(!is.null(gVars$agg.data)){
      print("Subsetting gVars$agg.data")
      rmIdx <- which(colnames(gVars$agg.data) %in% removedSampleIDs)
      #if(length(rmIdx)>1){
      if(length(rmIdx)>0){
        cat("Removing index ", rmIdx, " from gVars$agg.data matrix!\n")
        gVars$agg.data <- gVars$agg.data[,-rmIdx]
      }else{
        print("Nothing was removed!")
      }
    }
    if(!is.null(gVars$comb.sva.data)){
      print("Subsetting gVars$comb.sva.data")
      rmIdx <- which(colnames(gVars$comb.sva.data) %in% removedSampleIDs)
      #if(length(rmIdx)>1){
      if(length(rmIdx)>0){
        cat("Removing index ", rmIdx, " from gVars$comb.sva.data matrix!\n")
        gVars$comb.sva.data <- gVars$comb.sva.data[,-rmIdx]
      }else{
        print("Nothing was removed!")
      }
    }
    
    #Remove columns with single level data
    nrlevels <- apply(phTable, 2, function(x){length(levels(factor(x)))})
    nrlevels.singular <- which(nrlevels==1)
    
    remInfo <- 0
    remStr <- ""
    if(length(nrlevels.singular)>0){
      remInfo <- 1
      remStr <- paste0(remStr, "Following columns are removed because they contain only single repeated value:\n[",paste0(names(nrlevels.singular), collapse=", "), "]")
      if(length(nrlevels.singular)==ncol(phTable)){
        remStr <- paste0(remStr, "\n\nNo column survived filtering!!! Please define phenotype data columns with singular data type.")
        shinyjs::info(remStr)
        return(NULL)
      }
      col2rem <- which(colnames(phTable) %in% names(nrlevels.singular))
      phTable <- phTable[,-col2rem, drop=F]
    }
    
    #Inform user with the columns removed from the data frame
    if(remInfo==1){
      shinyjs::info(remStr)
    }
    
    factorCols <- gVars$factorCols
    if(is.null(factorCols)){
      phFactor <- phTable
    }else{
      tmpIdx <- which(factorCols %in% colnames(phTable))
      if(length(tmpIdx)>0){
        factorCols <- factorCols[tmpIdx]
        factorIdx <- which(colnames(phTable) %in% factorCols)
        phFactor <- factorize_cols(phTable, factorIdx)
      }else{
        phFactor <- phTable
      }
    }
    gVars$phTable <- phTable
    gVars$phFactor <- phFactor
    gVars$filteredSamples <- nrow(phTable)
    gVars$removedSamples <- gVars$totalSamples - nrow(phTable)
    gVars$removedSamplesInfo <- rbind(gVars$removedSamplesInfo, removedSamplesInfo)
    print("Removed Samples:")
    print(gVars$removedSamplesInfo)
  })
  
  output$filtered <- DT::renderDataTable({
    # shiny::validate(
    #   need(!is.null(gVars$phTable), "No Phenodata File Provided!")
    # )
    if(is.null(gVars$phTable))
      return(NULL)
    
    phTable <- gVars$phTable[[1]]
    DT::datatable(phTable, filter="none",
                  options = list(
                    search = list(regex=TRUE, caseInsensitive=FALSE),
                    scrollX=TRUE,
                    ordering=F
                  )
    )
  },server=TRUE)

  output$gExpMat <- DT::renderDataTable({
    
    if(is.null(gVars$inputGx))
      return(NULL)
    
    # if((class(gVars$inputGx) != "matrix") || (class(gVars$inputGx) != "data.frame"))
    #   return(NULL)
    
    print(gVars$inputGx)

    inputGx <- gVars$inputGx[[1]]
    DT::datatable(inputGx, filter="none",
                  options = list(
                    search = list(regex=TRUE, caseInsensitive=FALSE),
                    scrollX=TRUE,
                    ordering=F
                  )
    )
  },server=TRUE)
  
  observe({
    disable_ann <- FALSE
    disable_de <- FALSE
    disable_sva <- FALSE
    disable_combat <- FALSE
    disable_corr_skip <- FALSE
    disable_limma <- FALSE
    
    if(is.null(input$fPheno)){
      shinyjs::disable("load_pheno_submit")
    }else{
      shinyjs::enable("load_pheno_submit")
    }
    
    if(is.null(gVars$phTable)){
      shinyjs::disable("dirButton")
    }else{
      shinyjs::enable("dirButton")
    }
    
    if(is.null(gVars$phLoaded)){
      shinyjs::hide("phenoPreviewDiv")
    }else{
      shinyjs::show("phenoPreviewDiv")
    }
  })
  
  output$loading_text <- renderText({
    #loadingText <- gVars$loadingText
    loadingText <- "LOADING"
    return(loadingText)
  })
  
  output$selDoseCol <- renderUI({
    selectInput("doseCol", "Dose Variable", choices=gVars$phColChoices())
  })
  
  output$selTPCol <- renderUI({
    selectInput("TPCol", "Time Point Variable", choices=gVars$phColChoices())
  })
  
  output$selSampleIDCol <- renderUI({
    selectInput("sampleIDCol", "Sample ID Variable", choices=gVars$phColChoices())
  })
  
  
  output$selQuote <- renderUI({
    selectInput("quote", "Quotes", choices=gVars$quoteChoices, selected=gVars$quoteChoices[1])
  })
  
  output$selSep <- renderUI({
    #selectInput("sepS", "Column Seperator", choices=gVars$sepChoices, selected=gVars$sepChoices[1])
    selectInput("sepS", "Field Seperator", choices=gVars$sepChoices, selected=gVars$sepChoices[1])
  })
  
  
  output$selGxSep <- renderUI({
    selectInput("gxSepS", "Column Seperator", choices=gVars$sepChoices, selected=gVars$sepChoices[1])
  })
  
  output$selGxQuote <- renderUI({
    selectInput("gxQuote", "Quotes", choices=gVars$quoteChoices, selected=gVars$quoteChoices[1])
  })
  
  output$timePoint <- renderUI({
    shiny::validate(
      need(!is.null(gVars$phTable), "No Pheno Data file!")
    )
    print("Ciao")
    if(length(gVars$phTable)>1){
      selectInput("time_point_id", "Time Points", choices=c("All"))
      
    }else{
      selectInput("time_point_id", "Time Points", choices=c("All",unique(gVars$phTable[[1]][,gVars$TPColID])))
    }
    
  })
  
  output$trendTimePoint <- renderUI({
    shiny::validate(
      need(!is.null(gVars$phTable), "No Pheno Data file!")
    )
    print("Ciao")
    if(length(gVars$phTable)>1){
      selectInput("time_point_id", "Time Points", choices=c("All"))
      
    }else{
      selectInput("time_point_id", "Time Points", choices=c("All",unique(gVars$phTable[[1]][,gVars$TPColID])))
    }
    
  })
  

  output$bmd_checkbox = renderUI({
    #shiny::validate(need(!is.null(gVars$phTable), "No Pheno Data file!"))
    shiny::validate(need(is.null(gVars$BMDSettings), "No BMD Settings!"))
    
    if(input$BMDSettings == "All"){
      selected = c(1,2,4,5,6,7,8,12,13,18:34) 
    }
    if(input$BMDSettings == "Regulatory"){
      selected = 22:34
    }
    if(input$BMDSettings == "Custom"){
      selected = c(22,23,25,29,30,31)#c(19,21,22,23,25,27)
    }
    if(input$BMDSettings == "Degree of Freedom"){
      print("degree of freedom")
      if(is.null(gVars$phTable)){
        nDose = 0
      }else{
        nDoses = c()
        for(i in 1:length(gVars$phTable)){
          nDoses = c(nDoses,length(unique(gVars$phTable[[i]][,gVars$doseColID])))
        }
      }
      if(length(nDoses)>0){
        nDose = min(nDoses)-1
        DF =  nDose - 1
        modDF = c(2,3,Inf, 4,5,2,3,4,Inf,Inf,Inf,4,5,Inf,Inf,Inf,Inf,2,3,2,3,1,2,3,2,3,4,1,1,1,2,3,4,5)
        selected = which(modDF<=DF)
      }else{
        selected = c()
      }
    }
    
    tags$div(align = 'left',class = 'multicol', 
             checkboxGroupInput("ModGroup", label = "Models available", 
                                choices = list("Linear" = 22, "Quadratic" = 23, "Cubic" = 24,
                                               "Power2" = 25, "Power3" = 26,"Power4" = 27,
                                               "Exponential" = 28,
                                               "Hill05"= 29, "Hill1" = 30,
                                               "Hill2" = 31, "Hill3" = 32, 
                                               "Hill4" = 33, "Hill5" = 34,
                                               "AR.2"=18,"AR.3"= 19,
                                               "MM.2"=20,"MM.3"= 21
                                               #"LL.2"=1,"LL.3"=2,
                                               #"LL.4"=4,"LL.5"=5,
                                               #"W1.2"=6,"W1.3"=7,"W1.4"=8,
                                               #"BC.4"=12,"BC.5"=13
                                               #"LL2.2"=14,"LL2.3"=15,"LL2.4"=16,"LL2.5"=17,"LL.3u"=3,"W2.2"=9,"W2.3"=10,"W2.4"=11
                                ),
                                selected = selected,inline = FALSE)
    )
  })
  
  # output$MaxDose <- renderUI({
  #   shiny::validate(
  #     need(!is.null(gVars$phTable), "No Pheno Data file!")
  #   )
  #   selectInput("max_dose_input", "Maximum dose", choices=unique(gVars$phTable[[input$experiment_BMD]][,gVars$doseColID]), selected = max(as.numeric(unique(gVars$phTable[[input$experiment_BMD]][,gVars$doseColID]))))
  # })


  
  
  output$ExptimeSelAnova <- renderUI({
    if(is.null(gVars$EXP_FIL)){
      return(NULL)
    }
    print("Rending UI: ")
    print(names(gVars$EXP_FIL))

    selectInput("experiment_anova", "Experiment", choices=names(gVars$EXP_FIL), selected = names(gVars$EXP_FIL)[1])#c("All",unique(gVars$phTable[,gVars$TPColID])))
  })
  
  output$timePointSel <- renderUI({
    if(is.null(gVars$EXP_FIL)){
      return(NULL)
    }
    print("Rending UI: ")
    print(names(gVars$EXP_FIL))
    print(min(as.numeric(names(gVars$EXP_FIL))))
    
    expList = gVars$EXP_FIL[[input$experiment_anova]]
    
    selectInput("time_point_id_visual", "Time Points", choices=names(expList), 
                selected = min(as.numeric(names(expList))))
  })
  

  
  output$timePointSelPat <- renderUI({
    if(is.null(gVars$EnrichDatList)){
      return(NULL)
    }
    
    print("Rending UI: ")
    print(names(gVars$EXP_FIL))
    print(min(as.numeric(names(gVars$EXP_FIL))))
    
    #selectInput("time_point_id_visualPat", "Time Points", choices=names(gVars$EnrichDatList), selected = min(as.numeric(names(gVars$EnrichDatList))))#c("All",unique(gVars$phTable[,gVars$TPColID])))
    selectInput("time_point_id_visualPat", "Time Points", choices=names(gVars$EnrichDatList), selected = names(gVars$EnrichDatList)[1])#c("All",unique(gVars$phTable[,gVars$TPColID])))
    
  })
  

  # output$timePointSel3 <- renderUI({
  #   if(is.null(gVars$EXP_FIL)){
  #     return(NULL)
  #   }
  #   
  #   print("Rending UI: ")
  #   print(names(gVars$EXP_FIL))
  #   print(min(as.numeric(names(gVars$EXP_FIL))))
  #   
  #   selectInput("time_point_id_visua32", "Time Points", choices=names(gVars$EXP_FIL), selected = min(as.numeric(names(gVars$EXP_FIL))))#c("All",unique(gVars$phTable[,gVars$TPColID])))
  # })
  
  output$chose_lev1 <- renderUI({
    print("Render gui")
    selectInput("lev1", "Level 1", gVars$lev1_h,multiple = TRUE,selected = "All")
  })
  
  output$nClust <- renderUI({
    print("prova")
    if(is.null(dim(gVars$KEGG_MAT))){
      opt = c("N/A")
    }else{
      opt = 1:nrow(gVars$KEGG_MAT)
    }
    selectInput("nc", "Number of clusters", opt,multiple = FALSE,selected = 1)
  })
  
  output$chose_lev2 <- renderUI({
    # print(input$lev1)
    print("inside chose lev 2")
    
    need(input$lev1 != "", "Please select a level 1 object")
    
    if(is.null(input$lev1)){
      selectInput("lev2", "Level 2", gVars$lev2_h,multiple = TRUE,selected = "All")
    }else{
      if("All" %in% input$lev1){
        print("all in lev1")
        
        selectInput("lev2", "Level 2", gVars$lev2_h,multiple = TRUE,selected = "All")
      }else{
        if(length(input$lev1)==1){
          print("length is 1")
          #selectInput("lev2", "Level 2", as.list(lev2_h[[input$lev1]]),multiple = TRUE)
          selectInput("lev2", "Level 2", as.list(c("All",gVars$lev2_h[[input$lev1]])),multiple = TRUE,selected = "All")
          
        }else{
          print("length is greater than 1")
          #selectInput("lev2", "Level 2", lev2_h[input$lev1],multiple = TRUE)
          selectInput("lev2", "Level 2", c("All",gVars$lev2_h[input$lev1]),multiple = TRUE,selected = "All")
          
        }
      }      
    }
  })
  
  output$chose_lev3 <- renderUI({
    print("inside chose lev 2")
    need(input$lev1 != "", "Please select a level 1 object")
    need(input$lev2 != "", "Please select a level 2 object")
    
    if(is.null(input$lev2)){
      selectInput("lev3", "Level 3", gVars$lev3_h,multiple = TRUE,selected = "All")
    }else{
      if("All" %in% input$lev2){
        print("all in lev1")
        
        selectInput("lev3", "Level 3", gVars$lev3_h,multiple = TRUE,selected = "All")
      }else{
        if(length(input$lev2)==1){
          print("length is 1")
          selectInput("lev3", "Level 3", as.list(c("All",gVars$lev3_h[[input$lev2]])),multiple = TRUE,selected = "All")
        }else{
          print("length is greater than 1")
          selectInput("lev3", "Level 3", c("All",gVars$lev3_h[input$lev2]),multiple = TRUE,selected = "All")
        }
      }      
    }
  })
  
  output$selectColumn <- renderUI({
    selectInput("colID","Select samples",c("All",rownames(gVars$KEGG_MAT)),multiple=TRUE,selected = "All")
  })

  output$contents <- DT::renderDataTable({ #renderTable
    print("Inside contents")
    
    DF = DATA()
    shiny::validate(need(expr = !is.null(DF),message = "Waiting for input file!") )
    
    DT::datatable(DF, filter = "none", options = list(scrollX = TRUE))

  })
  
  output$updatedTable <- renderText({
    print("updata table")
    DF = DATA()
    shiny::validate(need(expr = !is.null(DF),message = "") )
    
    x = paste("Number of genes for each sample: ",sep="")
    return(x)
  })
  
  output$colSums <- DT::renderDataTable({
    print("Inside colSums")
    
    
    DF = DATA()
    shiny::validate(need(expr = !is.null(DF),message = "") )
    
    M = matrix("",1,ncol = ncol(DF))
    for(i in 1:ncol(DF)){
      M[1,i]= sum(DF[,i]!="")
    }
    colnames(M) = colnames(DF)
    DT::datatable(M, filter = "none", options = list(scrollX = TRUE))
    
    #return(M)
  })
  
  
  # observeEvent(input$computePathways,{
  # })
  # 
  output$heatmap <- renderPlot({
    shiny::validate(need(expr = !is.null(gVars$toPlot),message = "No data to plot"))
    print(class(gVars$toPlot))
    print(as.ggplot(gVars$toPlot))
    
  }, width = function(){
    if(!is.null(gVars$nSamples)){
      #
      mwidth = (gVars$length_path * 15) + (gVars$nSamples * 20)
      print(paste("MY sample IS ---->", gVars$nSamples))
      print(paste("MY path len IS ---->", gVars$length_path))
      print(paste("MY WIDTH IS ---->", mwidth))
      
      mwidth = max(600, mwidth)
      
      return(mwidth)
    }else{
      return("auto")
    }
  },  height = function(){
    
    if(!is.null(gVars$nPath)){
      mysize = (gVars$nPath* 20 ) + 10 * max(sapply(gVars$exp_ann,nchar))
      print(paste("MY HEIGHT IS ---->", mysize))
      
      print(paste("MY path IS ---->", gVars$nPath))
      mysize = min(max(600, mysize),30e3)
      
      return(mysize)
    }else{
      return("auto")
    }
  })
  
  observeEvent(input$do, {
    
    shinyjs::html(id="loadingText", "Rendering Map")
    shinyjs::show(id="loading-content")
    on.exit({
      print("inside on exit")
      shinyjs::hide(id="loading-content", anim=TRUE, animType="fade")    
    })
    # print("input lev1 is the following ---->")
    # print(input$lev1)
    need(is.null(input$lev1), "Please select a level 1 object")
    need(is.null(input$lev2), "Please select a level 2 object")
    need(is.null(input$lev3), "Please select a level 3 object")
    need(is.null(input$lev1), "Please select a level 1 object")
    need(is.null(input$lev2), "Please select a level 2 object")
    need(is.null(input$lev3), "Please select a level 3 object")
    
    print("INSIDE RENDER PLOT ----->>>>")
    print("Inside object event input$do")
    
    l1 = input$lev1
    l2 = input$lev2
    l3 = input$lev3
    
    #    print(head(gVars$hierarchy))
    gVars$reduced_kegg_hierarchy = update_hierarchy(kegg_hierarchy = gVars$hierarchy ,l1,l2,l3)
    #    print("Reduced Kegg hierarchy -->")
    #    print(head(gVars$reduced_kegg_hierarchy))
    
    gVars$samplesID = input$colID
    
    #    print(gVars$samplesID)
    
    if(!is.null(gVars$clust_mat)){
      # print("I did cluster")
      if("All" %in% gVars$samplesID){
        #  print("I WANT ALL SAMPLES")
        gVars$toPlotMap =gVars$clust_mat
      }else{
        gVars$toPlotMap = gVars$clust_mat[rownames(gVars$clust_mat) %in% gVars$samplesID,]
      }
    }else{
      if("All" %in% gVars$samplesID){
        #  print("I WANT ALL SAMPLES")
        gVars$toPlotMap = gVars$KEGG_MAT
      }else{
        gVars$toPlotMap = gVars$KEGG_MAT[rownames(gVars$KEGG_MAT) %in% gVars$samplesID,]
      }
      
    }
    
    print("before collapse")
    #controlla anche la selezione di righe e colonne
    #gVars$toPlot = plot_function(gcl = gVars$clust_mat,kegg_h = gVars$reduced_kegg_hierarchy,plt_mat = gVars$toPlotMap,input_n=as.numeric(input$level))
    
    kegg_nano_1 <- collapse_paths(kegg_hierarchy = gVars$reduced_kegg_hierarchy,kegg_mat_cell = gVars$toPlotMap, collapse_level = as.numeric(input$level))
    #extract collapsed matrix and collapsed hierarachy
    mat <- kegg_nano_1[[1]]
    if(sum(dim(gVars$toPlotMap) == c(1,1))==2){
      rownames(mat) = rownames(gVars$toPlotMap) # toPlotMap e' 1x123
    }
    hier <- kegg_nano_1[[2]]
    
    shiny::validate(need(expr = ncol(mat)>0, message = "No result for the enrichment or the filters are too restrictive. Please enlarge your selection"))
    
    if(is.null(gVars$clust_mat)){
      gVars$exp_ann = gVars$pheno#cbind(c(rep(1:5,5),1),rownames(mat))
      if("All" %in% gVars$samplesID == FALSE){
        gVars$exp_ann = gVars$exp_ann[gVars$exp_ann[,2] %in% gVars$samplesID,]
      }
    }
    
    xxx = gVars$exp_ann
    #save(mat, hier,xxx,file="demo/demo.RData")
    
    ############################   DISCRETIZE MAT if user chose "discrete and there are negative and positve values"
    ############################
    mat_to_Plot=mat
    print("mat_to_Plot dim ------------------->")
    print(dim(mat_to_Plot))
    gVars$nSamples = nrow(mat_to_Plot)
    gVars$nPath = ncol(mat_to_Plot)
    
    MLP = max(unlist(lapply(X = colnames(mat_to_Plot),FUN = nchar)))
    gVars$length_path = MLP
    
    print("maximum character length ------------------->")
    print(MLP)
    
    
    if (input$continuous=="discrete"){
      mat_to_Plot[mat_to_Plot<0]=-1
      mat_to_Plot[mat_to_Plot>0]=1
      isDiscrete = T
    }else{
      isDiscrete = F
    }
    ########################################################
    #gVars$toPlot = plot_grid(path_mat = mat_to_Plot,path_hier = hier,experiment_ann =  gVars$exp_ann,discrete =  isDiscrete,level_col = as.numeric(input$level),square_colors=c(),color_leg=c(),path_text_size = 12,treat_text_size = 12)
    if(input$doGrouping){
      print("grouping selected")
      
      print(dim(mat_to_Plot))
      
      print(gVars$exp_ann)
      
      if(nrow(mat_to_Plot)==1){
        rownames(mat_to_Plot) = gVars$exp_ann[1,2]
      }
      
      print("xxx")      
      gVars$toPlot = plot_grid(path_mat = mat_to_Plot,path_hier = hier,experiment_ann =  gVars$exp_ann,discrete =  isDiscrete,level_col = max(1,as.numeric(input$level)-1),square_colors=c(),color_leg=c(),path_text_size = 12,treat_text_size = 12, asRatio=(input$aspectRatio))
    }else{
      print("grouping NOT selected")
      level_n = max(1,as.numeric(input$level)-1)
      fake_hier = hier
      fake_hier[,level_n] = rep("",nrow(fake_hier))
      
      gVars$toPlot = plot_grid(path_mat = mat_to_Plot,path_hier = fake_hier,experiment_ann = gVars$exp_ann ,discrete =  isDiscrete,level_col =level_n,square_colors=c(),color_leg=c(),path_text_size = 12,treat_text_size = 12, asRatio=(input$aspectRatio))
      
    }
    
    #check for display size, pop up message if too large
    if(!is.null(gVars$nPath) && ((gVars$nPath* 20 ) + 10 * max(sapply(gVars$exp_ann,nchar))) > 30e3){
      print("exceeding dimensions")
      shinyjs::info("Warning: too many functional categories, the map might not be readable. Download the PDF for better resolution image.")
    }
    
    print("afterPLOTSPLOTSPLOTSPLOTSPLOTSPLOTSPLOTS")
  })
  
  output$downloadData <- downloadHandler(
    filename = paste(tempdir(),"/maps.pdf",sep=""),
    
    content = function(file) {
      wd = as.numeric(input$img_width)
      hi = as.numeric(input$img_height)
      
      pdf("www/map.pdf",width = wd,height = hi)
      plot(gVars$toPlot)
      dev.off()
      
      file.copy("www/map.pdf", file)
    }
  )
  
  observeEvent(input$resetCluster,{
    
    shinyjs::html(id="loadingText", "Rendering Map")
    shinyjs::show(id="loading-content")
    on.exit({
      print("inside on exit")
      shinyjs::hide(id="loading-content", anim=TRUE, animType="fade")    
    })
    gVars$exp_ann = gVars$pheno#cbind(c(rep(1:5,5),1),rownames(gVars$KEGG_MAT))
    gVars$clust_mat = NULL
    
    if("All" %in% gVars$samplesID){
      #  print("I WANT ALL SAMPLES")
      gVars$toPlotMap = gVars$KEGG_MAT
    }else{
      gVars$toPlotMap = gVars$KEGG_MAT[rownames(gVars$KEGG_MAT) %in% gVars$samplesID,]
    }
    
    
    #gVars$toPlotMap = gVars$toPlotMap
    #gVars$toPlot = plot_function(gcl = gVars$clust_mat,kegg_h = gVars$reduced_kegg_hierarchy,plt_mat = gVars$toPlotMap,input_n=as.numeric(input$level))
    kegg_nano_1 <- collapse_paths(kegg_hierarchy = gVars$reduced_kegg_hierarchy,kegg_mat_cell = gVars$toPlotMap, collapse_level = as.numeric(input$level))
    #extract collapsed matrix and collapsed hierarachy
    mat <- kegg_nano_1[[1]]
    hier <- kegg_nano_1[[2]]
    
    shiny::validate(need(expr = ncol(mat)>0, message = "No result for the enrichment or the filters are too restrictive. Please enlarge your selection"))
    
    print("Matrix dimension -->")
    print(dim(mat))
    print("Hierarchy dimension -->")
    print(dim(hier))
    #plot the collapsed matrix
    
    if(is.null(gVars$clust_mat)){
      gVars$exp_ann = gVars$pheno#cbind(c(rep(1:5,5),1),rownames(mat))
      if("All" %in% gVars$samplesID == FALSE){
        gVars$exp_ann = gVars$exp_ann[gVars$exp_ann[,2] %in% gVars$samplesID,]
      }
    }
    
    xxx = gVars$exp_ann
    #save(mat, hier,xxx,file="demo/demo.RData")
    
    ############################   DISCRETIZE MAT if user chose "discrete and there are negative and positve values"
    ############################
    mat_to_Plot=mat
    
    if (input$continuous=="discrete"){
      mat_to_Plot[mat_to_Plot<0]=-1
      mat_to_Plot[mat_to_Plot>0]=1
      isDiscrete = T
    }else{
      isDiscrete = F
    }
    ########################################################
    
    if(input$doGrouping){
      print("grouping selected")
      gVars$toPlot = plot_grid(path_mat = mat_to_Plot,path_hier = hier,experiment_ann =  gVars$exp_ann,discrete =  isDiscrete,level_col = max(1,as.numeric(input$level)-1),square_colors=c(),color_leg=c(),path_text_size = 12,treat_text_size = 12, asRatio=(input$aspectRatio))
    }else{
      print("grouping NOT selected")
      level_n = max(1,as.numeric(input$level)-1)
      fake_hier = hier
      fake_hier[,level_n] = rep("",nrow(fake_hier))
      
      gVars$toPlot = plot_grid(path_mat = mat_to_Plot,path_hier = fake_hier,experiment_ann = gVars$exp_ann ,discrete =  isDiscrete,level_col =level_n,square_colors=c(),color_leg=c(),path_text_size = 12,treat_text_size = 12, asRatio=(input$aspectRatio))
      
    }
    print("after PLOTS in reset cluster")  
    gVars$hls <- NULL
    gVars$cls <- NULL
    gVars$clustered <- 1
  })
  

  
  observeEvent(input$doCluster,{
    
    shinyjs::html(id="loadingText", "Rendering Map")
    shinyjs::show(id="loading-content")
    on.exit({
      print("inside on exit")
      shinyjs::hide(id="loading-content", anim=TRUE, animType="fade")    
    })
    
    print("Clustering columns")
    M=gVars$KEGG_MAT
    #Jaccard index distance
    
    D = matrix(data = 0,nrow = nrow(M),ncol = nrow(M))
    rownames(D) = colnames(D) = rownames(M)
    for(i in 1:(nrow(M)-1)){
      pi= colnames(M)[!is.na(M[i,])]
      for(j in (i+1):nrow(M)){
        pj= colnames(M)[!is.na(M[j,])]
        D[i,j] = D[j,i] = length(intersect(pi,pj))/length(union(pi,pj))    
      }
    }
    idx= which(rowSums(!is.na(M)) == 0)
    print(idx)
    
    if(length(idx)>0){
      D[idx,] = 0
      D[,idx] = 0
      D[idx,idx] = 1
    }
    diag(D) = 1
    
    D1 = matrix(data = 0,nrow = nrow(M),ncol = nrow(M))
    rownames(D1) = colnames(D1) = rownames(M)
    
    for(i in 1:(nrow(M)-1)){
      idx1 = which(!is.na(M[i,]))
      print(idx1)
      for(j in (i+1):nrow(M)){
        idx2 = which(!is.na(M[j,]))
        print(idx2)
        idx = intersect(idx1,idx2)
        print(idx)
        
        if(length(idx)>0){
          D1[i,j] = D1[j,i] = as.numeric(dist(t(cbind(M[i,idx],M[j,idx]))))
        }
        
      }
    }
    
    D = 1-D
    
    #View(D)
    
    print(class(D1))
    print(dim(D1))
    
    range01 <- function(x){(x-min(x))/(max(x)-min(x))}
    D1 = range01(D1)
    print(class(D1))
    print(dim(D1))
    
    #View(D1)
    
    if(input$Distance %in% "euclidean"){
      DD = D1
    }else{
      if(input$Distance %in% "jaccard"){
        DD = D
      }else{
        DD = (D + D1)/2
      }
    }
    
    print(class(DD))
    print(dim(DD))
    #View(DD)
    
    hls = hclust(as.dist(DD),method = input$ClusterMethod)
    gVars$hls <- hls
    
    #plot(hls)
    
    output$hclust_plot = renderPlot({
      plot(hls,xlab="", sub="",hang = -1)
      if(as.numeric(input$nc)>1){
        rect.hclust(tree = hls,k = as.numeric(input$nc))
      }
    })
    
    print(rownames(gVars$KEGG_MAT))
    print(hls$order)
    gVars$clust_mat = gVars$KEGG_MAT[hls$order,]
    
    
    cat("CHECK ROWNAMES -->")
    print(rownames(gVars$clust_mat))
    nClust = as.numeric(input$nc)
    
    cls = cutree(hls,k=nClust)
    cls = cls[hls$order]
    gVars$cls <- cls
    
    cat("CLUSTERING RESULTS -->")
    print(cls)
    
    gVars$exp_ann = cbind(cls,names(cls))#nrow(gVars$KEGG_MAT)
    
    if("All" %in% gVars$samplesID == FALSE){
      gVars$exp_ann = gVars$exp_ann[gVars$exp_ann[,2] %in% gVars$samplesID,]
    }
    
    print(gVars$exp_ann)
    gVars$toPlotMap = gVars$clust_mat
    #gVars$toPlot = plot_function(kegg_h = gVars$reduced_kegg_hierarchy,plt_mat = gVars$clust_mat,input_n=as.numeric(input$level))
    kegg_nano_1 <- collapse_paths(kegg_hierarchy = gVars$reduced_kegg_hierarchy,kegg_mat_cell = gVars$toPlotMap, collapse_level = as.numeric(input$level))
    #extract collapsed matrix and collapsed hierarachy
    mat <- kegg_nano_1[[1]]
    hier <- kegg_nano_1[[2]]
    
    shiny::validate(need(expr = ncol(mat)>0, message = "No result for the enrichment or the filters are too restrictive. Please enlarge your selection"))
    
    print("Matrix dimension -->")
    print(dim(mat))
    print("Hierarchy dimension -->")
    print(dim(hier))

    if(is.null(gVars$clust_mat)){
      gVars$exp_ann = gVars$pheno#cbind(c(rep(1:5,5),1),rownames(mat))
      if("All" %in% gVars$samplesID == FALSE){
        gVars$exp_ann = gVars$exp_ann[gVars$exp_ann[,2] %in% gVars$samplesID,]
      }
    }
    
    ############################   DISCRETIZE MAT if user chose "discrete and there are negative and positve values"
    ############################
    mat_to_Plot=mat
    
    if (input$continuous=="discrete"){
      mat_to_Plot[mat_to_Plot<0]=-1
      mat_to_Plot[mat_to_Plot>0]=1
      isDiscrete = T
    }else{
      isDiscrete = F
    }
    ########################################################
    
    if(input$doGrouping){
      print("grouping selected")
      gVars$toPlot = plot_grid(path_mat = mat_to_Plot,path_hier = hier,experiment_ann =  gVars$exp_ann,discrete =  isDiscrete,level_col = max(1,as.numeric(input$level)-1),square_colors=c(),color_leg=c(),path_text_size = 12,treat_text_size = 12, asRatio=(input$aspectRatio))
    }else{
      print("grouping NOT selected")
      level_n = max(1,as.numeric(input$level)-1)
      fake_hier = hier
      fake_hier[,level_n] = rep("",nrow(fake_hier))
      
      gVars$toPlot = plot_grid(path_mat = mat_to_Plot,path_hier = fake_hier,experiment_ann = gVars$exp_ann ,discrete =  isDiscrete,level_col =level_n,square_colors=c(),color_leg=c(),path_text_size = 12,treat_text_size = 12, asRatio=(input$aspectRatio))
      
    }
    
    print("after PLOTS in clustering")
    gVars$clustered <- gVars$clustered + 1
    
  })

  # GENES HEATMAP CODE
  
  output$heatmapGenes <- renderPlot({
    shiny::validate(need(expr = !is.null(gVars$toPlotGenes), message = "No data to plot"))
    
    #Singular gene map
    print("class(gVars$toPlotGenes)")
    print(class(gVars$toPlotGenes))
    #print(gVars$toPlotGenes)
    #print(as.ggplot(gVars$toPlotGenes))
    grid::grid.draw(gVars$toPlotGenes)
    
    # #Get from list of grobs
    # print("class(gVars$toPlotGenes$gplotGene)")
    # print(class(gVars$toPlotGenes$gplotGene))
    # #print(gVars$toPlotGenes)
    # print(as.ggplot(gVars$toPlotGenes$gplotGene))
  }, width = function(){
    if(!is.null(gVars$nSamples)){
      mwidth = (gVars$length_gene * 15) + (gVars$nSamplesGene * 20)
      print(paste("MY sample IS ---->", gVars$nSamplesGene))
      print(paste("MY path len IS ---->", gVars$length_gene))
      print(paste("MY WIDTH IS ---->", mwidth))
      
      mwidth = max(600, mwidth)
      return(mwidth)
    }else{
      return("auto")
    }
  }, height = function(){
    if(!is.null(gVars$nGenes)){
      exp_ann <- gVars$exp_ann
      print("exp_ann")
      print(exp_ann)
      if(is.null(exp_ann) || nrow(exp_ann)==0){
        exp_ann <- gVars$pheno
      }
      mysize = (gVars$nGenes* 20 ) + 10 * max(sapply(exp_ann, nchar))
      print(paste("MY HEIGHT IS ---->", mysize))
      print(paste("MY path IS ---->", gVars$nGenes))
      
      mysize = min(max(600, mysize),30e3)
      return(mysize)
    }else{
      return("auto")
    }
  })
  
  output$choosePath <- renderUI({
    print("Render gui")
    if(input$levelGene==1){
      levLocal <- gVars$lev1_h
    }else if(input$levelGene==2){
      levLocal <- gVars$lev2_h
    }else if(input$levelGene==3){
      levLocal <- gVars$lev3_h
    }else{
      print("Incorrect level provided!")
      levLocal <- c("NONE")
    }
    
    idxAll <- which(tolower(levLocal)=="all")
    if(length(idxAll)>0){
      print("!!!!!!!! REPLACING 'ALL' !!!!!!!")
      levLocal[idxAll] <- "NONE"
      if(length(names(levLocal))>0){
        names(levLocal)[idxAll] <- "NONE"
      }
    }else{
      print("!!!!!!!! ALL NOT FOUND !!!!!!!")
    }
    
    # print("Getting genes matrix...")
    # Mgenes = kegg_mat_genes(EnrichDatList=gVars$EnrichDatList, hierarchy=gVars$hierarchy)
    
    selectInput("selPath", "Choose Term", levLocal, multiple=FALSE, selected="NONE")
  })
  
  #Generate gene heatmap
  observeEvent(c(
    input$doGeneHeatMap,
    gVars$clustered
  ), {
    if(is.null(gVars$KEGG_MAT_GENES) || nrow(gVars$KEGG_MAT_GENES)==0){
      shinyjs::info("Enrichment matrix is NULL or Empty!!!")
      return(NULL)
    }
    
    if(is.null(input$selPath) || input$selPath=="NONE"){
      shinyjs::info("Please select a Term!")
      return(NULL)
    }
    
    shinyjs::html(id="loadingText", "Rendering Map")
    shinyjs::show(id="loading-content")
    on.exit({
      print("inside on exit")
      shinyjs::hide(id="loading-content", anim=TRUE, animType="fade")
    })
    # print("input lev1 is the following ---->")
    # print(input$lev1)
    
    
    print("INSIDE RENDER GENE HEATMAP ----->>>>")
    print("Inside object event input$doGeneHeatMap")
    
    selLevel <- as.integer(input$levelGene)
    selTerm <- input$selPath
    kegg_hierarchy <- gVars$hierarchy
    
    print("selLevel")
    print(selLevel)
    
    print("selTerm")
    print(selTerm)
    
    print("dim(kegg_hierarchy)")
    print(dim(kegg_hierarchy))
    
    print("str(kegg_hierarchy)")
    print(str(kegg_hierarchy))
    
    #Filter hierarchy by selected level and term
    kegg_hierarchy_sel <- dplyr::filter(kegg_hierarchy, kegg_hierarchy[,selLevel]==selTerm)
    kegg_hierarchy_sel_pathIDs <- unique(as.character(kegg_hierarchy_sel$ID))
    print("kegg_hierarchy_sel_pathIDs")
    print(kegg_hierarchy_sel_pathIDs)
    
    ##Get enriched gene matrix
    KEGG_MAT_GENES <- gVars$KEGG_MAT_GENES
    print("str(KEGG_MAT_GENES)")
    print(str(KEGG_MAT_GENES))
    
    pathIdMapIdx <- which(kegg_hierarchy_sel_pathIDs %in% colnames(KEGG_MAT_GENES))
    
    if(length(pathIdMapIdx)==0){
      shinyjs::info("No enriched information for the selected TERM!")
    }else{
      kegg_hierarchy_sel_pathIDs <- kegg_hierarchy_sel_pathIDs[pathIdMapIdx]
    }
    
    KEGG_MAT_GENES_sel <- KEGG_MAT_GENES[, kegg_hierarchy_sel_pathIDs]
    
    print("dim(KEGG_MAT_GENES_sel)")
    print(dim(KEGG_MAT_GENES_sel))
    
    print("str(KEGG_MAT_GENES_sel)")
    print(str(KEGG_MAT_GENES_sel))
    
    KEGG_MAT_GENES_enriched_genes <- unique(unlist(KEGG_MAT_GENES_sel))
    
    print("length(KEGG_MAT_GENES_enriched_genes)")
    print(length(KEGG_MAT_GENES_enriched_genes))
    
    print("KEGG_MAT_GENES_enriched_genes")
    print(KEGG_MAT_GENES_enriched_genes)
    
    #    print(head(gVars$hierarchy))
    #gVars$reduced_kegg_hierarchy = update_hierarchy(kegg_hierarchy = gVars$hierarchy ,l1,l2,l3)
    #    print("Reduced Kegg hierarchy -->")
    reduced_kegg_hierarchy <- gVars$reduced_kegg_hierarchy
    
    #gVars$samplesID = input$colID
    #    print(gVars$samplesID)
    samplesID <- gVars$samplesID
    
    # if(!is.null(gVars$clust_mat)){
    #   # print("I did cluster")
    #   if("All" %in% gVars$samplesID){
    #     #  print("I WANT ALL SAMPLES")
    #     gVars$toPlotMap =gVars$clust_mat
    #   }else{
    #     gVars$toPlotMap = gVars$clust_mat[rownames(gVars$clust_mat) %in% gVars$samplesID,]
    #   }
    # }else{
    #   if("All" %in% gVars$samplesID){
    #     #  print("I WANT ALL SAMPLES")
    #     gVars$toPlotMap = gVars$KEGG_MAT
    #   }else{
    #     gVars$toPlotMap = gVars$KEGG_MAT[rownames(gVars$KEGG_MAT) %in% gVars$samplesID,]
    #   }
    #
    # }
    
    #Gene tables uploaded by user
    GList <- gVars$GList
    egGeneTable <- GList[[1]]
    print("dim(egGeneTable)")
    print(dim(egGeneTable))
    
    #Score type selected by user
    selScoreType <- input$selScoreType
    
    if(selScoreType %in% c("comb", "pval") && ncol(egGeneTable)<3){
      shinyjs::info("Missing required scores!")
      return(NULL)
    }
    
    GListMod <- NULL
    if(selScoreType=="comb"){
      GListMod <- lapply(GList,
                         function(x){
                           x <- mutate(x, dscore=x[2] * -log10(x[3]))
                           x <- x[,c(1,4)]
                           return(x)
                         }
      )
    }else{
      valMapVec <- c("lfc"=2, "pval"=3)
      selColIDx <- unname(valMapVec[selScoreType])
      print("selColIDx")
      print(selColIDx)
      GListMod <- lapply(GList,
                         function(x){
                           x <- x[,c(1,selColIDx)]
                           return(x)
                         }
      )
    }
    
    if(is.null(GListMod)){
      shinyjs::info("Unexpected response from score selection!")
      return(NULL)
    }
    
    ##Get the reduced data frame from combining all samples
    geneTableColNames <- colnames(egGeneTable)
    geneIDCol <- geneTableColNames[1]
    print("geneIDCol")
    print(geneIDCol)
    reducedGeneTable <- purrr::reduce(GListMod, dplyr::full_join, by=geneIDCol)
    print("str(reducedGeneTable)")
    print(str(reducedGeneTable))
    
    print("str(gVars$toPlotMap)")
    print(str(gVars$toPlotMap))
    
    reducedGeneTableSel <- dplyr::filter(reducedGeneTable, tolower(reducedGeneTable[,geneIDCol]) %in% tolower(KEGG_MAT_GENES_enriched_genes))
    print("dim(reducedGeneTableSel)")
    print(dim(reducedGeneTableSel))
    
    print("str(reducedGeneTableSel)")
    print(str(reducedGeneTableSel))
    
    reducedGeneTableSel <- tibble::column_to_rownames(reducedGeneTableSel, geneIDCol)
    
    print("rownames(reducedGeneTableSel)")
    print(rownames(reducedGeneTableSel))
    
    print("names(GList)")
    print(names(GList))
    
    colnames(reducedGeneTableSel) <- names(GList)
    print("dim(reducedGeneTableSel)")
    print(dim(reducedGeneTableSel))
    
    print("str(reducedGeneTableSel)")
    print(str(reducedGeneTableSel))
    
    if(ncol(reducedGeneTableSel)==0){
      shinyjs::info("No result for the enrichment or the filters are too restrictive. Please enlarge your selection")
      return(NULL)
    }
    
    mat_to_Plot <- t(as.matrix(reducedGeneTableSel))
    
    KEGG_MAT <- gVars$KEGG_MAT
    pathIdMapIdx <- which(kegg_hierarchy_sel_pathIDs %in% colnames(KEGG_MAT))
    gridPlotTypeGoup <- NULL
    if(length(pathIdMapIdx)==0){
      shinyjs::info("Did not find TERM enrichment information in the pathway matrix!!!")
    }else{
      print("############ Adding path information to gene plot ############")
      kegg_hierarchy_sel_pathIDs <- kegg_hierarchy_sel_pathIDs[pathIdMapIdx]
      KEGG_MAT_sel <- KEGG_MAT[, kegg_hierarchy_sel_pathIDs, drop=FALSE]
      
      print("str(KEGG_MAT_sel)")
      print(str(KEGG_MAT_sel))
      
      print("head(KEGG_MAT_sel)")
      print(head(KEGG_MAT_sel))
      
      res_collapsed <- collapse_paths(kegg_hierarchy=kegg_hierarchy_sel, kegg_mat_cell=KEGG_MAT_sel, collapse_level=as.numeric(selLevel))
      #extract collapsed matrix and collapsed hierarachy
      KEGG_MAT_sel_collaped <- res_collapsed[[1]]
      KEGG_MAT_sel_collaped_hier <- res_collapsed[[2]]
      
      print("str(KEGG_MAT_sel_collaped)")
      print(str(KEGG_MAT_sel_collaped))
      
      print("head(KEGG_MAT_sel_collaped)")
      print(head(KEGG_MAT_sel_collaped))
      
      print("str(KEGG_MAT_sel_collaped_hier)")
      print(str(KEGG_MAT_sel_collaped_hier))
      
      #return(NULL)
      
      gridPlotTypeGoup <- factor(
        c(
          rep("Pathway", ncol(KEGG_MAT_sel_collaped)*nrow(mat_to_Plot)),
          rep("Genes", ncol(mat_to_Plot)*nrow(mat_to_Plot))
        ),
        levels=c("Pathway", "Genes")
      )
      mat_to_Plot <- cbind(KEGG_MAT_sel_collaped, mat_to_Plot)
    }
    
    #xxx = gVars$exp_ann
    
    ############################   DISCRETIZE MAT if user chose "discrete and there are negative and positve values"
    ############################
    
    print("mat_to_Plot dim ------------------->")
    print(dim(mat_to_Plot))
    gVars$nSamplesGene = nrow(mat_to_Plot)
    gVars$nGenes = ncol(mat_to_Plot)
    
    #Get maximum character length for the columns
    MLG = max(unlist(lapply(X = colnames(mat_to_Plot),FUN = nchar)))
    gVars$length_gene = MLG
    print("maximum character length ------------------->")
    print(MLG)
    
    #Performing discretization
    if (input$continuous=="discrete"){
      mat_to_Plot[mat_to_Plot<0]=-1
      mat_to_Plot[mat_to_Plot>0]=1
      isDiscrete = T
    }else{
      isDiscrete = F
    }
    
    ########################################################
    
    ##No categories to display for genes
    # if(input$doGrouping){
    #   print("grouping selected")
    #   print(dim(mat_to_Plot))
    #   print(gVars$exp_ann)
    #
    #   gVars$toPlot = plot_grid(path_mat = mat_to_Plot,path_hier = hier,experiment_ann =  gVars$exp_ann,discrete =  isDiscrete,level_col = max(1,as.numeric(input$level)-1),square_colors=c(),color_leg=c(),path_text_size = 12,treat_text_size = 12, asRatio=(input$aspectRatio))
    # }else{
    print("grouping NOT selected")
    
    exp_ann <- gVars$exp_ann
    print("exp_ann")
    print(exp_ann)
    if(is.null(exp_ann) || nrow(exp_ann)==0){
      exp_ann <- gVars$pheno
    }
    print("exp_ann")
    print(exp_ann)
    
    #Get clustering variables
    hls <- gVars$hls
    cls <- gVars$cls
    
    if(!is.null(hls)){
      mat_to_Plot = mat_to_Plot[hls$order,]
    }
    
    if(!is.null(cls)){
      exp_ann = cbind(cls,names(cls))
    }
    print("exp_ann")
    print(exp_ann)
    
    gVars$toPlotGenes = plot_grid_genes(path_mat=mat_to_Plot, experiment_ann=exp_ann, gene_group=gridPlotTypeGoup, discrete=isDiscrete, square_colors=c(), color_leg=c(), treat_text_size=12, asRatio=(input$aspectRatio))
    # }
    
    #check for display size, pop up message if too large
    if(!is.null(gVars$nGenes) && ((gVars$nGenes* 20 ) + 10 * max(sapply(exp_ann, nchar))) > 30e3){
      print("Exceeding dimensions")
      shinyjs::info("Warning: too many functional categories, the map might not be readable. Download the PDF for better resolution image.")
    }
    print("Exiting doGeneHeatMap!")
  }, ignoreInit = TRUE)
  
  observe({
    M = gVars$DATA
    if(is.null(M)){
      shinyjs::disable("computePathways")
    }else{
      shinyjs::enable("computePathways")
    }
    
    if(is.null(gVars$KEGG_MAT)){
      shinyjs::disable("do")
      shinyjs::disable("downloadData")
      shinyjs::disable("doCluster")
      shinyjs::disable("resetCluster")
    }else{
      shinyjs::enable("do")
      shinyjs::enable("downloadData")
      shinyjs::enable("doCluster")
      shinyjs::enable("resetCluster")
    }
  })
  
  
  
  
})
