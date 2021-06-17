library(shiny)

#shinyUI(
  fluidPage(
    titlePanel("The Natural Mortality Tool: Empirical Estimators of Natural Mortality (M)"),
    h5(p(em("This tool employs various empirical estimators of natural mortality."))),
    h5(p(em("As the user enters values for the below input parameters, estimates will be displayed in the main panel."))), 
    h5(p(em("Uncertainty can also be added to the estimates by chosing a coefficient of variation (CV) value and error distribution type."))),
    h5(p(em("Downloaded R objects (.DMP) can be loaded directly into R or, if using the R terminal, drag and drop in."))),
    br(),
    h4(p("Do you have any suggested methods to add? Please submit an issue with the recommendation" ,tags$a(href="https://github.com/shcaba/Natural-Mortality-Tool/issues", "here"))),
    h4(p("References for each included method can be found",tags$a(href="javascript:window.open('References_M.html', '_blank','width=600,height=400')", "here"))),
    
    sidebarLayout(
    sidebarPanel(
    conditionalPanel(
      condition="input.conditionedPanels==1",
        h3("Parameters inputs"),
        h5(p(em("Provide a value for the CV and choose an error type to add additional uncertainty to the point estimates. The CV will be applied to all methods."))),
        h5(p(em("CV = 0 means only point estimates will be reported"))), 
        fluidRow(column(width=6,numericInput("M_CV", "CV in M", value=0,min=0, max=10, step=0.1)),
            column(width=6,selectInput("M_CV_type","Error type",c("lognormal","normal")))),    
        h5(p(em("Provide inputs below. You do not need values for all inputs; natural mortality will only be estimated for a given methods when all input requirements are met."))),
        h5(p(em("Scientific name calls the FishLife M estimator"))),
        h5(p(em("Input requirements for other methods can be found",tags$a(href="javascript:window.open('Method_inputs.html', '_blank','width=600,height=400')", "here")))),
        textInput("Genspp","Scientific name",value=""),
        fluidRow(column(width=6,numericInput("Amax", "Longevity (yrs):", value=NA,min=0.1, max=300, step=0.1)),
            column(width=6,numericInput("Linf","VBGF Linf (in cm):", value=NA,min=1, max=10000, step=0.01))),    
        fluidRow(column(6,numericInput("k_vbgf", "VBGF k:", value=NA,min = 0.001, max = 1,step=0.01)),
            column(width=6,numericInput("t0", "VBGF t0", value=NA,min = -15, max = 15,step=0.01))),    
        fluidRow(column(width=6,numericInput("Age_in","Age (yr) specified for Chen-Wat:", value=NA,min=0.1, max=300, step=0.01)),
            column(width=6,numericInput("Lt_in", "Length (cm) specified for Gislason", value=NA,min = 0, max = 10000,step=0.01))),    
        fluidRow(column(6,numericInput("Amat","Age at maturity (yrs)", value=NA,min = 0.01, max = 100,step=0.01)),
            column(width=6,numericInput("Temp","Water temp. (in C):" , value=NA,min = 0.001, max = 60,step=0.01))),
        fluidRow(column(6,numericInput("Winf","VBGF Winf (in g):", value=NA,min = 0, max = 100000,step=0.1)),
            column(width=6,numericInput("kw","VBGF kw: ", value=NA,min = 0.001, max = 5,step=0.01))),
        fluidRow(column(6,numericInput("Wdry","Total dry weight (in g):" ,value=NA,min = 0.01, max = 1000000,step=0.01)),
            column(width=6,numericInput("Wwet","Total wet weight (in g):" ,value=NA,min = 0.01, max = 1000000,step=0.01))),
        fluidRow(column(width=6,numericInput("GSI","Gonadosomatic index (GSI):",value=NA,min = 0, max = 1,step=0.001)),
             column(width=6,textInput("User_M","User M input:",value=""))), 
       ),
       
       
       conditionalPanel(
       condition="input.conditionedPanels==2",
       h3("Composite M: method weighting"),
       h5(p(em("Allows for weighting of the contribution of each method in the composite M distribution"))),
       h5("Values range from 0 to 1. A value of 0 removes the contribution; a value of 1 is full weighting."),
       h5("Default values are based on redundancies of methods using similar information. For instance,the four longevity-based methods are given a weight of 0.25, so all weighted together equal 1."),
       h5("The prior sample number generates a prior based on the number of specified samples. This value also defines the binwidth in the density plot, thus lower sample numbers will give a more diffuse prior."),
       wellPanel(
          fluidRow(
            column(6,numericInput("FishLife","FishLife",value=1,min = 0, max = 1,step=0.001))
          ),
#          h5("Uses longevity")
          fluidRow(
            column(6,numericInput("Then_nls","Then_nls",value=0.25,min = 0, max = 1,step=0.001)),
            column(6,numericInput("Then_lm","Then_lm",value=0.25,min = 0, max = 1,step=0.001))
          ),
          fluidRow(
            column(6,numericInput("Hamel_Amax","Hamel_Amax",value=0.5,min = 0, max = 1,step=0.001)),
            column(6,numericInput("Chen_Wat","Chen-Wat",value=0.333,min = 0, max = 1,step=0.001))            
          ),
          fluidRow(
            column(6,numericInput("ZM_CA_pel","ZM_CA_pel",value=0,min = 0, max = 1,step=0.001)),
            column(6,numericInput("ZM_CA_dem","ZM_CA_dem",value=0,min = 0, max = 1,step=0.001))
          ),
          fluidRow(
            column(6,numericInput("Then_VBGF","Then_VBGF",value=0.25,min = 0, max = 1,step=0.001)),
            column(6,numericInput("Hamel_VBGF","Hamel_k",value=0.5,min = 0, max = 1,step=0.001))
          ),
          fluidRow(
            column(6,numericInput("Jensen_VBGF_1","Jensen_k 1",value=0.25,min = 0, max = 1,step=0.001)),
            column(6,numericInput("Jensen_VBGF_2","Jensen_k 2",value=0.0,min = 0, max = 1,step=0.001))
          ),
          fluidRow(
            column(6,numericInput("Gislason","Gislason",value=0.333,min = 0, max = 1,step=0.001)),
            column(6,numericInput("Charnov","Charnov",value=0.333,min = 0, max = 1,step=0.001))
          ),
          fluidRow(
            column(6,numericInput("Pauly_lt","Pauly_lt",value=0.5,min = 0, max = 1,step=0.001)),
            column(6,numericInput("Roff","Roff",value=0.5,min = 0, max = 1,step=0.001))
          ),
          fluidRow(
            column(6,numericInput("Jensen_Amat","Jensen_Amat",value=0.5,min = 0, max = 1,step=0.001)),
            column(6,numericInput("Ri_Ef_Amat","Ri_Ef_Amat",value=0.5,min = 0, max = 1,step=0.001))
        ),
          fluidRow(
            column(6,numericInput("Pauly_wt","Pauly_wt",value=0.5,min = 0, max = 1,step=0.001)),
            column(6,numericInput("McGl","McC&Gil",value=0.5,min = 0, max = 1,step=0.001))
          ),

          fluidRow(
            column(6,numericInput("PnW","PnW",value=0.5,min = 0, max = 1,step=0.001)),
            column(6,numericInput("Lorenzen","Lorenzen",value=1,min = 0, max = 1,step=0.001))
        ),
         fluidRow(
            column(6,numericInput("Gonosoma","GSI",value=1,min = 0, max = 1,step=0.001)),
           column(6,numericInput("UserM_wt","User M",value=1,min = 0, max = 1,step=0.001))
        ),
        h5(p(em("M prior control parameters"))),
        fluidRow(
           column(6,numericInput("samp.num","Prior sample #",value=1000000,min = 0, max = 10000000,step=1)),
           column(6,numericInput("ad.bw","Bandwidth multiplier",value=1,min = 0.0001, max = 100,step=0.01))),
        h4(p("Composite M prior downloads")),
        h5(p(em("Downloads are DMP files that can be loaded directly into R"))),
        h5(p(em("Higher sample sizes create a composite distribution truer to the component distributions"))),
        h5(p(em("Adjusting bandwidth >1 creates a more diffuse composite distribution"))),
        fluidRow(
           column(6,downloadButton('downloadMcompositedist', 'Based on sample-sized')),
           column(6,downloadButton('downloadMcompositedistupdated', 'Bandwidth-adjusted')))
        )
       )
       ), #end sidebar
         mainPanel(
          tabsetPanel(
          tabPanel("M by method",
          h4("Natural mortality (M) estimates by method"),
          h5("Legend color indicate inputs used by each method."),
          h5("Downloadable R object contains paramter inputs, estimted M point estimates, and age-specific M values."),
          plotOutput("Mplot"),
          h4("Natural mortality (M) point estimates"),
          fluidRow(
            column(4,tableOutput("Mtable")),
            column(4,tableOutput("Mtable2")),
            column(4,tableOutput("MtableUser")),
            ),
         downloadButton('downloadMplots', 'Download M values plot'),
         downloadButton('downloadMs', 'Download M table csv file'),
         downloadButton('downloadMandPs', 'Download parameter inputs and M values R object'),
         br(),
         br(),
         br(),
         h4("Natural mortality (M) by age"),
         plotOutput("Mplot_ages"),
         downloadButton('downloadMplot_ages', 'Download M at age plot'),
         downloadButton('downloadCW_M_a', 'Download age-specific M values csv file'),
            value=1
          ),
          tabPanel("Composite M",
            h4("Method density and weights"),
            h5("When no uncertainty is expressed, point estimates with weights are shown"),
            h5("When uncertainty is expressed, the distribution of each method is given, with its frequency (determined by the weighting) plotted on the y-axis"),
            plotOutput("Mdistplots"),
            downloadButton('downloadMdensityplots', 'Download M density plot'),
            downloadButton('downloadMdistvals', 'Download M density inputs as R object'),
            br(),
            br(),
            br(),
            h4("Composite natural mortality"),
            h5(p(em("Blue vertical line indicates median value"))),
            h5(p(em("First composite M object is based on the number of specified samples"))),
            h5(p(em("Second composite M object is based on th adjusted bandwidth"))),
            plotOutput("Mcomposite"),
            downloadButton('downloadMcompositedensityplot', 'Download composite M density plot'),
            # downloadButton('downloadMcompositedist', 'Download composite M as R object'),
            # downloadButton('downloadMcompositedistupdated', 'Download bandwidth-adjusted composite M as R object'),
            value=2
          ), id="conditionedPanels"
          )
      )
    ) 
)
#)