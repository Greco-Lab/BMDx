suppressMessages(library(shiny))
suppressMessages(library(shinyjs))
suppressMessages(library(shinyBS))
suppressMessages(library(shinydashboard))
suppressMessages(library(shinyFiles))
suppressMessages(library(DT))
suppressMessages(library(rhandsontable))
suppressMessages(library(shinycssloaders))
library(plotly)
library(xtable)

appCSS <- "
//.main-header { z-index: 100000; }
.main-sidebar { background-color: white !important;}
.sidebar { color: black; max-height: 900px; overflow-y: scroll; }
.sidebar a { color: blue; margin: 20px !important;}
.content-wrapper { margin-left: 15%; }
//.panel { background-color: #222d32; }
.panel-title a { font-weight: bold; color: white !important; }
.panel-warning .panel-heading { background-color: #00c0ef; }
.panel-warning { border-color: #8de9ff; }
.panel-sample .panel-heading { background-color: #2C6C99; }
.panel-sample { border-color: #2C6C99; }
.panel-danger .panel-heading { background-color: #dd4b39; }
.panel-success .panel-heading { background-color: #00a65a; }
.panel-info .panel-heading { background-color: #7e46ff; }
.panel-primary .panel-heading { background-color: #3079ae; }
.multicol { 
height: 150px;
-webkit-column-count: 5; /* Chrome, Safari, Opera */ 
-moz-column-count: 5;    /* Firefox */ 
column-count: 5; 
-moz-column-fill: auto;
-column-fill: auto;
} 
#loading-content {
position: absolute;
background: white !important;
opacity: 0.8;
z-index: 1000000;
left: 0;
right: 0;
top: 0;
bottom: 0;
font-size: 50px;
text-align: center;
color: #black;
}
#loading-gif { 
opacity: 0.8; 
display: block;
margin-left: auto;
margin-right: auto;
vertical-align: middle;
z-index: 1000000;
}
"

jsCode <- "
callback = 'function(table) {
table.on('click.dt', 'tr', function() {
table.$('tr.selected').removeClass('selected');
$(this).toggleClass('selected');            
Shiny.onInputChange('rows',
table.rows('.selected').data()[0][0]);
});
}'
"

fluidPage(
  useShinyjs(),
  inlineCSS(appCSS),
  hidden(div(id="loading-content",
             img(id="loading-gif", src="screen-loading.gif"),
             p(id="loadingText", "WORKING"),
             p("...")
  )),
  dashboardPage(
    dashboardHeader(title="BMDx: dose response analysis for expression data", titleWidth="35%"),
    dashboardSidebar(disable=FALSE,
                     bsCollapse(id="bsSidebar1", open="LOAD PHENOTYPE DATA",
                                bsCollapsePanel("LOAD PHENOTYPE DATA", style="warning",
                                                fluidRow(
                                                  column(12, align="center",
                                                         shinyBS::bsButton("import_pheno_submit", label="Import Phenotype Data", style="danger",icon=icon("exclamation-circle")),
                                                         shinyBS::bsTooltip("import_pheno_submit", "Launch a graphical window, to configure import of phenotype data from a file!", placement="bottom")
                                                  )
                                                )
                                ),
                                bsCollapsePanel("LOAD EXPRESSION MATRIX", style="danger",
                                                fluidRow(
                                                  column(12, align="center",
                                                         shinyBS::bsButton("import_expr_submit", label="Import Expression Matrix", style="danger", icon=icon("exclamation-circle")),
                                                         shinyBS::bsTooltip("import_expr_submit", "Launch a graphical window, to configure import of expression matrix from a file!", placement="bottom")
                                                         
                                                  )
                                                )
                                ),
                                bsCollapsePanel("GENE FILTERING", style="danger",
                                                fluidRow(
                                                  column(12, align="center",
                                                         shinyBS::bsButton("anova_filtering_button", label="Anova Filtering", style="danger", icon=icon("exclamation-circle")),
                                                         shinyBS::bsTooltip("anova_filtering_button", "Launch a graphical window, to configure anova parameters!", placement="bottom")
                                                  )
                                                ),
                                                fluidRow(
                                                  column(12, align="center",
                                                         shinyBS::bsButton("trend_filtering_button", label="Trend Filtering", style="danger", icon=icon("exclamation-circle")),
                                                         shinyBS::bsTooltip("trend_filtering_button", "Launch a graphical window, to configure anova parameters!", placement="bottom")
                                                  )
                                                ),
                                                
                                                fluidRow(
                                                  column(12, align="center",
                                                         shinyBS::bsButton("skip_anova_filtering_button", label="Skip Filtering", style="danger", icon=icon("exclamation-circle")),
                                                         shinyBS::bsTooltip("skip_anova_filtering_button", "SkipFiltering!", placement="bottom")
                                                         
                                                  )
                                                )       
                                ),
                                bsCollapsePanel("COMPUTE BMD", style="danger",
                                                fluidRow(
                                                  column(12, align="center",
                                                         shinyBS::bsButton("bmd_button", label="Compute BMD", style="danger", icon=icon("exclamation-circle")),
                                                         shinyBS::bsTooltip("bmd_button", "Launch a graphical window to configure BMD analysis!", placement="bottom")
                                                  )
                                                )
                                ),
                                bsCollapsePanel("PATHWAY ENRICHMENT", style="danger",
                                                fluidRow(
                                                  column(12, align="center",
                                                         shinyBS::bsButton("enrich_button", label="Enrichment", style="danger", icon=icon("exclamation-circle")),
                                                         shinyBS::bsTooltip("enrich_button", "Launch a graphical window to configure enrichment analysis!", placement="bottom")
                                                  )
                                                )
                                )
                     ),
                     bsCollapse(id="bsSidebar0", open="MORE INFO",
                                bsCollapsePanel("MORE INFO", style="sample",
                                                fluidRow(
                                                  HTML(" <a style=color:blue;  target=\"_blank\"; href=\"https://github.com/Greco-Lab/BMDx\">GitHub</a>")
                                                ),
                                                fluidRow(
                                                  HTML("<a style=color:blue; target=\"_blank\";  href=\"Manual.pdf\">Manual</a>")
                                                ),
                                                fluidRow(
                                                  HTML("<a style=color:blue;  target=\"_blank\"; href=\"https://github.com/Greco-Lab/BMDx/blob/master/pheno_list_2_exp_4_TP.xlsx\">Sample pheno data</a>")
                                                ),
                                                fluidRow(
                                                  HTML("<a style=color:blue;  target=\"_blank\"; href=\"https://github.com/Greco-Lab/BMDx/blob/master/exp_mat_file_2_exp_4_tp.xlsx\">Sample expression data</a>")
                                                )
                                )
                     )
    ),
    dashboardBody(
      shinyBS::bsModal("enrichPathways", "Enrichment", "enrich_button", size="large",
                       tags$h5("1. Input gene lists"),
                       wellPanel(
                         fluidRow(
                           column(6,
                                  radioButtons("organism","1) Organisms",
                                               choices = c(human = "Human", mouse = "Mouse", rat = "Rat"),selected = "Human"),
                                  shinyBS::bsTooltip(id = "organism",title = "Note: select organism and gene ID before uploading the file",placement = "bottom")
                                  
                           ) ,
                           column(6,radioButtons("idtype","2) GeneID",
                                                 choices = c(symbols = "SYMBOL", ensemble = "ENSEMBL",entrez = "ENTREZID"),
                                                 selected = "SYMBOL"),
                                  shinyBS::bsTooltip(id = "idtype",title = "Note: select organism and gene ID before uploading the file",placement = "bottom")
                                  
                                  
                           )
                         )
                       ),
                       tags$h5("2. Functional annotation parameters"),
                       wellPanel(
                         fluidRow(
                           column(6,radioButtons("EnrichType","Select Functional Annotation",
                                                 choices = c(KEGG = "KEGG", REACTOME="REACTOME",GO = "GO"),
                                                 selected = "KEGG")
                           ),
                           column(6,radioButtons("GOType","Select GO",
                                                 choices = c(BP = "BP", CC="CC",MF = "MF"),
                                                 selected = "BP")
                           )),
                         fluidRow(
                           column(6,radioButtons("pcorrection","Correction Method",
                                                 choices = c(gSCS = "analytical", fdr = "fdr", bonferroni = "bonferroni", Nominal = "none"),
                                                 selected = "analytical"),
                                  shinyBS::bsTooltip(id = "pcorrection",
                                                     title = "Default is g:SCS. Check g:Profiler web page for more info",
                                                     placement = "top")
                           ),
                           column(6,selectInput(inputId = "pvalueTh", label = "P-value threshold:",choices = list(0.001,0.005,0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09),selected = 0.05)
                                  
                                  
                           )),
                         fluidRow(
                           column(6,checkboxInput("only_annotated", "Annotated genes only", value = TRUE)),
                           column(6, sliderInput("min_intersection", "Minumum number of genes in the intersection:",
                                                 min = 0, max = 100, value = 0))
                         )
                       ),
                       tags$h5("3. Display parameters"),
                       wellPanel(
                         fluidRow(
                           
                           column(6,radioButtons("aggregation","Aggregation Function",
                                                 choices = c(min = "min", max = "max",mean = "mean",
                                                             median = "median"),
                                                 selected = "mean")  
                           ),
                           
                           column(6,
                                  radioButtons("MapValueType","Choose Values Type",
                                               choices = c(Pvalue = "PVAL", GenesModifications="FC",GenesModifications_PValue  = "FCPV"),
                                               selected = "FC")
                           )
                         ),
                         
                         fluidRow(column(6,radioButtons("continuous","Plot modification",
                                                        choices = c(value = "continuous",
                                                                    sign = "discrete"),
                                                        selected = "continuous"))
                         )
                         
                       ),
                       tags$h5("4. Filter BMD Values"),
                       wellPanel(
                         fluidRow(
                           #column(6, sliderInput("filterBMD", "BMD values Filter:", min = 0, max = 100,value = c(0,100))),
                           column(6,  shinyBS::bsButton("enrichment_analysis", label="Run Enrichment", style="info", icon=icon("hand-o-right")))
                           
                         )
                       ),
                       tags$h5("5. Download gene lists as excel file"),
                       downloadButton("downloadFunmapponeInput", "Download")
      ),
      shinyBS::bsModal("computeBMD", "Compute BMD Value", "bmd_button", size="large",
                       fluidRow(
                         column(6, textInput("RespLev", label = "Response Level", value =1.349)),
                         column(6, checkboxInput("constantVar", "Assumption of Constant Variance", value = TRUE))
                       ),
                       fluidRow(
                         column(4,selectInput("LOOF", "Lack-of-fit PValue Th:", choices=c(0.3,0.2,0.1,0.05),selected=0.1)),
                         column(4,selectInput("min_dose_perc", "Lowest dose filter:", choices=c(0,0.1,0.2,0.3),selected=0)),
                         column(4,selectInput("max_dose_perc", "Highest dose filter:", choices=c(0,0.1,0.2,0.3),selected=0))
                       ),
                       fluidRow(
                         column(3, selectInput("BMDNCores", "Number of cores:", choices=c(1:25), selected = 3))
                       ),
                       fluidRow(
                         column(12,selectInput("BMDSettings", "Select the BMD analysis setting", choices=c("All","Regulatory","Degree of Freedom", "Custom"),selected="Custom"))
                       ),
                       fluidRow(
                         column(12,wellPanel(uiOutput("bmd_checkbox")))
                       ),
                       fluidRow(
                         column(12, align="right",shinyBS::bsButton("bmd_analysis", label="Run BMD", style="info", icon=icon("hand-o-right")))
                       ),
                       fluidRow(
                         column(12,
                                h5("Models Description"),                     
                                bsCollapse(id="BMDModHelp",
                                           bsCollapsePanel("Linear Model", style="primary",
                                                           withMathJax(),
                                                           helpText("The formula for the linear model is
                                                                    $$ f(dose) = \\beta_0 + \\beta_1 dose $$
                                                                    The linear model is a special case of the polynomial model, with $$n=1$$")
                                                           ),
                                           bsCollapsePanel("Polynomial Model (Quadratic/Cubic)", style="primary",
                                                           withMathJax(),
                                                           helpText("The formula for the polynomial model is
                                                                    $$ f(dose) = \\beta_0 + \\beta_1 dose + \\beta_2 dose^2 + \\ldots + \\beta_n dose^n $$
                                                                    Here n is the degree of the polynomial. The user can choose between $$n = 2, 3$$")
                                                           ),
                                           bsCollapsePanel("Power Model", style="primary",
                                                           withMathJax(),
                                                           helpText("The formula for the power model is
                                                                    $$ f(dose) = \\beta_0 + (dose)^\\delta $$
                                                                    The user can choose between $$\\delta = 2, 3, 4$$")                                           
                                                           ),
                                           bsCollapsePanel("Exponential Model", style="primary",
                                                           withMathJax(),
                                                           helpText("The formula for the exponential model is
                                                                    $$ f(dose) = \\beta_0 + exp(dose) $$")                                           
                                                           ),
                                           bsCollapsePanel("Hill Model", style="primary",
                                                           withMathJax(),
                                                           helpText("The formula for the hill model is 
                                                                    $$ f(dose) = \\beta_0 + \\dfrac{dose^n}{Kd + dose^n} $$
                                                                    The user can choose between $$n = 0.5,1,2,3,4,5$$ while Kd is fixed to 10.")
                                                           
                                                           ),
                                           bsCollapsePanel("Asymptotic Regression", style = "primary",
                                                           withMathJax(),
                                                           helpText("The formula for the asymptotic regression model is the following:
                                                                    $$f(dose) = c + (d-c) \\times (1-exp(-dose/e)) $$
                                                                    The parameter c is the lower limit (at x=0), the parameter d is the upper limit and 
                                                                    the parameter e>0 is determining the steepness of the increase of dose.
                                                                    The AR.3 model is the one depending from c, d and e parameters. The AR.2 model depends only on d and e parameters, while c is set to zero")
                                                           ),
                                           bsCollapsePanel("Michaelis-Menten Model", style = "primary",
                                                           withMathJax(),
                                                           helpText("The model is defined by the three-parameter model (MM.3) function
                                                                    $$f(dose, (c, d, e)) = c + \\dfrac{d-c}{1+(e/dose)}$$
                                                                    It is increasing as a function of the dose, attaining the lower limit c at dose 0 (x=0) and the upper limit d for infinitely large doses. 
                                                                    The parameter e corresponds to the dose yielding a response halfway between c and d. 
                                                                    The common two-parameter Michaelis-Menten model (MM.2) is obtained by setting c equal to 0.")
                                                           )
                                                           )
                                                           )
                                                           )
                                           ),
      shinyBS::bsModal("computeAnova", "Filter Genes by Anova", "anova_filtering_button", size="large",
                       fluidRow(
                         column(3, uiOutput("timePoint")),
                         column(3, radioButtons("adjBool", "Pvalue:", c("FDR" = TRUE,"Nominal" = FALSE),selected = FALSE)),
                         column(3, selectInput("anovaPvalTh", "Anova PValue Th:", choices=c(0.05,0.04,0.03,0.02,0.01))),
                         column(3, selectInput("anovaNCores", "Number of cores:", choices=c(1:25), selected = 3))
                       ),
                       fluidRow(
                         textOutput("anova_update")
                       ),
                       fluidRow(
                         column(12, align="right",shinyBS::bsButton("anova_analysis", label="Run Anova", style="info", icon=icon("hand-o-right")))
                       )                       
      ),
      shinyBS::bsModal("computeTrend", "Filter Genes by Trend Test", "trend_filtering_button", size="large",
                       fluidRow(
                         column(3, uiOutput("trendTimePoint")),
                         column(3, radioButtons("TrendAdjBool", "Pvalue:", c("FDR" = TRUE,"Nominal" = FALSE),selected = FALSE)),
                         column(3, selectInput("trendPvalTh", "Trend PValue Th:", choices=c(0.05,0.04,0.03,0.02,0.01))),
                         column(3, selectInput("trendNCores", "Number of cores:", choices=c(1:25), selected = 3))
                       ),
                       # fluidRow(
                       #   textOutput("anova_update")
                       # ),
                       fluidRow(
                         column(12, align="right",shinyBS::bsButton("trend_analysis", label="Run Trend Test", style="info", icon=icon("hand-o-right")))
                       )                       
      ),
      shinyBS::bsModal("importGxModal", "Import Gene Expression Table", "import_expr_submit", size="large",
                       fluidRow(
                         column(3,fileInput("gx", label="File"))
                         #column(3,uiOutput("selGxSep")),
                         #column(3,textInput("gxSepT", "Other Seperator", value=":")),
                         #column(3,uiOutput("selGxQuote"))
                       ),fluidRow(
                         column(12, align="right",shinyBS::bsButton("upload_gx_submit", label="Import", style="info", icon=icon("hand-o-right")))
                       )
      ),
      shinyBS::bsModal("importPhenoModal", "Import Phenotype Data", "import_pheno_submit", size="large",
                       fluidRow(column(3, fileInput("fPheno", label="Phenotype File"))
                                #column(3, uiOutput("selSep")),
                                #column(3, textInput("sepT", "Other Seperator", value=":")),
                                #column(3, uiOutput("selQuote"))
                       ),fluidRow(
                         column(3,actionButton("load_pheno_submit", "Preview")),
                         column(2, align="left", textOutput("phRowsText"),textOutput("phColsText"))
                       ),hr(),
                       fluidRow(
                         column(12,
                                hidden(div(id="phenoPreviewDiv",
                                           fluidRow(
                                             column(12,
                                                    rhandsontable::rHandsontableOutput("phenoTypesRH")
                                             )
                                           ),hr(),
                                           fluidRow(
                                             column(4, uiOutput("selSampleIDCol")),
                                             column(4, uiOutput("selDoseCol")),
                                             column(4, uiOutput("selTPCol"))
                                           ),fluidRow(column(12, align="right",shinyBS::bsButton("upload_pheno_submit", label="Import", style="info", icon=icon("hand-o-right"))))
                                ))
                         )
                       )
      ),
      fluidRow(column(12,
                      tabBox(id="display", title="", width=12,
                             tabPanel(value="pdTab", title="Phenotype Data",fluidRow(column(12,DT::dataTableOutput("filtered")))),
                             tabPanel(value="gExpTab", title="Gene Expression Matrix",
                                      fluidRow(column(12,DT::dataTableOutput("gExpMat")))
                             ),
                             tabPanel(value="AnovaTab", title="Filtering",
                                      fluidRow(
                                        column(4,uiOutput("ExptimeSelAnova")),
                                        column(4,uiOutput("timePointSel")),
                                        column(4,downloadButton("downloadAnovaData", "Download"))
                                      ),
                                      fluidRow(
                                        column(6,DT::dataTableOutput("Anova_table")),
                                        column(6,fluidRow(shinycssloaders::withSpinner(plotOutput("anovaPlot"), type = 6)))
                                      )
                             ),
                             tabPanel(value="BMDTab", title="BMD",
                                      fluidRow(column(12,
                                                      tabBox(id = "bmdgenes", title = "", width = 12,
                                                             #tabsetPanel(
                                                             tabPanel("Gene Level", 
                                                                      fluidRow(
                                                                        column(4,uiOutput("ExptimeSelBMD")),
                                                                        column(4,uiOutput("timePointSel2")),
                                                                        column(4,downloadButton("downloadBMDData", "Download"))
                                                                      ),
                                                                      
                                                                      fluidRow(
                                                                        column(12, DT::dataTableOutput("BMD_table"))
                                                                      ),
                                                                      fluidRow(
                                                                        column(8,fluidRow(shinycssloaders::withSpinner(plotOutput("bmd_fitting"), type = 6))),
                                                                        column(4,fluidRow(plotOutput("bmd_fitting_legend")))
                                                                      )   
                                                                      #torem
                                                                      # fluidRow(
                                                                      #   column(9, DT::dataTableOutput('x3')),
                                                                      #   column(3, verbatimTextOutput('x4'))
                                                                      # )
                                                             ),
                                                             tabPanel("Compare TP",
                                                                      tabPanel("Compare Time points",
                                                                               fluidRow(column(12,
                                                                                               bsCollapse(id="TPSidebar", open="BMD",
                                                                                                          bsCollapsePanel("BMD Values", style="primary",
                                                                                                                          #"Here i will print the BMD value distribution",
                                                                                                                          shinycssloaders::withSpinner(plotlyOutput("BMD_dist_TP"), type = 6)
                                                                                                          ),
                                                                                                          bsCollapsePanel("Lack of fit Pvalues", style="primary",
                                                                                                                          #"Here i will print the BMD value distribution",
                                                                                                                          shinycssloaders::withSpinner(plotlyOutput("BMD_pval_fitting"), type = 6)
                                                                                                                          
                                                                                                          ),
                                                                                                          bsCollapsePanel("BMD/BMDL", style="primary",
                                                                                                                          #"Here i will print the BMD value distribution",
                                                                                                                          shinycssloaders::withSpinner(plotlyOutput("BMD_BMDL"), type = 6)
                                                                                                                          
                                                                                                          ),
                                                                                                          bsCollapsePanel("Fitted models", style="primary",
                                                                                                                          #"Here i will print the BMD value distribution",
                                                                                                                          fluidRow(column(12,shinycssloaders::withSpinner(plotOutput("BMD_dist_models"), type= 6))),
                                                                                                                          fluidRow(column(12, shinycssloaders::withSpinner(plotlyOutput("BMD_BMDL_BMDU_by_model"), type = 6)))
                                                                                                          ),
                                                                                                          bsCollapsePanel("Gene by Time Point", style="primary",
                                                                                                                          #"Here i will plot the BMD for every gene at different TP",
                                                                                                                          shinycssloaders::withSpinner(plotlyOutput("NGTime"), type = 6)
                                                                                                          ),
                                                                                                          bsCollapsePanel("Venn diagram responsive genes", style="primary",
                                                                                                                          #"Here i will plot the BMD for every gene at different TP",
                                                                                                                          fluidRow(
                                                                                                                            column(6,shinycssloaders::withSpinner(plotOutput("NGVenn"),type = 6)),
                                                                                                                            column(6, DT::dataTableOutput("VennDF"))
                                                                                                                          ),
                                                                                                                          hr(),
                                                                                                                          # "Genes that show dose-response behaviour at all timepoints",
                                                                                                                          fluidRow(column(6,uiOutput("intersectionNameUI"))),
                                                                                                                          fluidRow(
                                                                                                                            column(6,uiOutput("nclustvenn")),
                                                                                                                            column(6,selectInput(inputId = "geneClustMeth", label = "Agglomerative Method",
                                                                                                                                                 choices = list("ward","average","complete","single"), selected = "ward"))
                                                                                                                          ),
                                                                                                                          fluidRow(column(12, plotlyOutput("gene_bmd_plot"))
                                                                                                                          )
                                                                                                          )
                                                                                               ) 
                                                                               )
                                                                               
                                                                               ))     
                                                                      
                                                             ),
                                                             tabPanel("Compare Experiments",
                                                                      fluidRow(
                                                                        column(3,sliderInput("range", "Range:", min = 0, max = 100, value = c(0,10))),
                                                                        column(3,sliderInput("maxInt", "Maximum intersection:", min = 1, max = 100, value = 10)),
                                                                        column(3,selectInput("groupby", "GroupBy:", c("degree"="degree", "sets" = "sets"))),
                                                                        column(3,selectInput("orderby", "OrderBy:", c("frequency"="freq", "degree" = "degree")))
                                                                      ),
                                                                      fluidRow(column(12,plotOutput("upsetplot")) )       
                                                             )
                                                      )
                                      )
                                      )), #end tabpanel
                             tabPanel(value = "enrTab",title = "Enrichment",
                                      sidebarLayout(
                                        sidebarPanel(
                                          wellPanel(
                                            tags$h5("1. Data Selection"),
                                            fluidRow(selectInput(inputId = "level", label = "Browse hierarchy: choose a level",choices = list(1,2,3))),
                                            fluidRow(
                                              column(4,uiOutput("chose_lev1")), shinyBS::bsTooltip(id = "chose_lev1",title = "Note: remove ALL from the list for specific selection.",placement = "top"),
                                              column(4,uiOutput("chose_lev2")), shinyBS::bsTooltip(id = "chose_lev2",title = "Note: remove ALL from the list for specific selection.",placement = "top"),
                                              column(3,uiOutput("chose_lev3")), shinyBS::bsTooltip(id = "chose_lev3",title = "Note: remove ALL from the list for specific selection.",placement = "top")
                                            ),
                                            fluidRow(
                                              uiOutput("selectColumn"), shinyBS::bsTooltip(id = "selectColumn",title = "Note: remove ALL from the list for specific selection.",placement = "top")
                                            )
                                          ),
                                          wellPanel(
                                            tags$h5("2. Plot section"),
                                            
                                            fluidRow(
                                              column(4,checkboxInput("doGrouping", "Show categories", value = TRUE)),
                                              column(4,checkboxInput("aspectRatio", "Keep aspect ratio", value = TRUE)),
                                              column(4,actionButton("do", "Plot Map")), shinyBS::bsTooltip(id = "do",title ="NOTE: press the Plot Mat button every time you update the map!",placement = "bottom")
                                            )
                                          ),
                                          wellPanel(
                                            tags$h5("3. Download Selection"),
                                            fluidRow(
                                              column(4,textInput(inputId ="img_width", value = 15,label = "Width")), #width
                                              column(4,textInput(inputId ="img_height", value = 30,label = "Height")),
                                              column(4,downloadButton('downloadData')), shinyBS::bsTooltip(id = "downloadData",title ="NOTE: when downloading, specify image size in inches ",placement = "bottom")
                                            )
                                          )
                                        ),
                                        mainPanel(
                                          tags$h5("Use scrollbars to navigate and see the whole map"),
                                          
                                          tabsetPanel(
                                            
                                            tabPanel("Heatmap",fluidRow(column(12,align="left",shinycssloaders::withSpinner(plotOutput(outputId="heatmap"), type=6)))),
                                            tabPanel("Cluster Bubble Plot", plotlyOutput("pathway_bubble")),
                                            tabPanel("Mean BMD for Time Point", plotlyOutput("meanBMD_timepoint")),
                                            tabPanel("Gene BMD in pathway",
                                                     fluidRow(
                                                       column(6, uiOutput("timePointSelPat"))
                                                     ),
                                                     fluidRow(
                                                       column(12,
                                                              DT::dataTableOutput("PAT_table")
                                                       )),
                                                     fluidRow(
                                                       column(12,
                                                              plotlyOutput("path_bmd_dist"))
                                                     )
                                            ),
                                            tabPanel("Pathways Table", 
                                                     #"Here I will display the enriched pathway in tabular form",
                                                     fluidRow(uiOutput("timePointSelTab")),
                                                     fluidRow(downloadButton("downloadEnrichedPathwayTables", "Download")),
                                                     fluidRow(DT::dataTableOutput("PatTable"))
                                            ),
                                            #tabPanel("Genes Heatmap", "Here I will display the gene heatmap")
                                            tabPanel("Heatmap Genes",
                                                     fluidRow(column(4,
                                                                     selectInput(inputId="levelGene", label="Choose a hierarchy level", choices=list(1,2,3))
                                                     ),column(4,
                                                              uiOutput("choosePath")
                                                     ),column(4,
                                                              selectInput(inputId="selScoreType", label="Show Values", choices=list("BMD"="lfc", "P-Value"="pval", "Combined"="comb"), selected="lfc")
                                                     )),fluidRow(column(4,
                                                                        #shinyBS::bsButton("doGeneHeatMap", label="Plot", style="danger", icon=icon("exclamation-circle"))
                                                                        actionButton("doGeneHeatMap", "Plot")
                                                     )),fluidRow(column(12,align="center",
                                                                        shinycssloaders::withSpinner(plotOutput(outputId="heatmapGenes"), type=6)
                                                     ))
                                            )
                                          )
                                        )
                                      )
                                      
                                      
                             )
                      )
      ))
                                                           )
                                                           )
                         )

