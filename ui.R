library(shiny)

shinyUI(
  fluidPage(
    titlePanel("Estimating Natural Mortality (M)"),
    h5(p(em("This tool employs various empirical estimators of natural mortality."))),
    h5(p(em("As the user enters values for the below input parameters,"))), 
    h5(p(em("estimates will be displayed in the main panel."))),
    br(),
    h4(p("References for each method can be found",tags$a(href="javascript:window.open('References_M.html', '_blank','width=600,height=400')", "here"))),
    
    headerPanel("Input parameters"),
    sidebarLayout(
        sidebarPanel
       (
        numericInput("Amax", "Maximum age (years):", value=NA,min=1, max=300, step=0.1),    
        numericInput("Linf","Linf (in cm):", value=NA,min=1, max=1000, step=0.01),
        numericInput("k", "VBGF Growth coeff. k:", value=NA,min = 0.001, max = 1,step=0.01),
        numericInput("t0", "VBGF age at size 0 (t_0)", value=NA,min = -15, max = 15,step=0.01),
        numericInput("Amat","Age at maturity (years)", value=NA,min = 0.01, max = 100,step=0.01),
        numericInput("Winf","Asym. weight (Winf, in g):", value=NA,min = 0, max = 100000,step=0.1),
        numericInput("kw","VBGF Growth coeff. wt. (kw, in g): ", value=NA,min = 0.001, max = 5,step=0.01),
        numericInput("Temp","Water temperature (in C):" , value=NA,min = 0.001, max = 60,step=0.01),
        numericInput("Wdry","Total dry weight (in g):" ,value=NA,min = 0.01, max = 1000000,step=0.01),
        numericInput("Wwet","Total wet weight (in g):" ,value=NA,min = 0.01, max = 1000000,step=0.01),
        numericInput("Bl","Body length (cm):",value=NA,min = 0.01, max = 10000,step=0.01),
        numericInput("GSI","Gonadosomatic index:",value=NA,min = 0, max = 1,step=0.001),
        numericInput("User_M","User M input:",value=NA,min = 0, max = 10,step=0.001),
         
       br(),
       br(),
       
       h3("Composite M: method weighting"),
       h5(p(em("Allows for weighting of the contribution of each method in the composite M distribution"))),
       h5("Values range from 0 to 1. A value of 0 removes the contribution; a value of 1 is full weighting."),
       wellPanel(
          fluidRow(
            column(4,numericInput("Then_Amax_1","Then_Amax 1",value=1,min = 0, max = 1,step=0.001)),
            column(4,numericInput("Then_Amax_2","Then_Amax 2",value=1,min = 0, max = 1,step=0.001)),
            column(4,numericInput("Then_Amax_3","Then_Amax 3",value=1,min = 0, max = 1,step=0.001))
          ),
          fluidRow(
            column(4,numericInput("Hamel_Amax","Hamel_Amax",value=1,min = 0, max = 1,step=0.001)),
            column(4,numericInput("AnC","AnC",value=1,min = 0, max = 1,step=0.001)),
            column(4,numericInput("Then_VBGF","Then_VBGF",value=1,min = 0, max = 1,step=0.001))
          ),
          fluidRow(
            column(4,numericInput("Jensen_VBGF_1","Jensen_VBGF 1",value=1,min = 0, max = 1,step=0.001)),
            column(4,numericInput("Jensen_VBGF_2","Jensen_VBGF 2",value=1,min = 0, max = 1,step=0.001)),
            column(4,numericInput("Pauly_lt","Pauly_lt",value=1,min = 0, max = 1,step=0.001))
          ),
          fluidRow(
            column(4,numericInput("Gislason","Gislason",value=1,min = 0, max = 1,step=0.001)),
            column(4,numericInput("Chen_Wat","Chen-Wat",value=1,min = 0, max = 1,step=0.001)),
            column(4,numericInput("Roff","Roff",value=1,min = 0, max = 1,step=0.001))
          ),
          fluidRow(
            column(4,numericInput("Jensen_Amat","Jensen_Amat",value=1,min = 0, max = 1,step=0.001)),
            column(4,numericInput("Pauly_wt","Pauly_wt",value=1,min = 0, max = 1,step=0.001)),
            column(4,numericInput("PnW","PnW",value=1,min = 0, max = 1,step=0.001))
          ),
          fluidRow(
            column(4,numericInput("Lorenzen","Lorenzen",value=1,min = 0, max = 1,step=0.001)),
            column(4,numericInput("Gonosoma","GSI",value=1,min = 0, max = 1,step=0.001)),
            column(4,numericInput("UserM","User M",value=1,min = 0, max = 1,step=0.001)))
       )
       ),
         mainPanel(
          h4("Natural mortality (M) estimates by method"),
          plotOutput("Mplot"),
          h4("Natural mortality (M) values"),
          fluidRow(
            column(6,tableOutput("Mtable")),
            column(6,tableOutput("Mtable2")),
            downloadButton('downloadMs', 'Download M values'),
            downloadButton('downloadCW_M_a', 'Download Chen-Wat. age-specific M values'),
            br(),
            br(),
            br(),
            h4("Composite natural mortality"),
            h5(p(em("Blue vertical line indicates median value"))),
            plotOutput("Mcomposite"),
            downloadButton('downloadMcompositedensityplot', 'Download composite M denisty plot'),
            downloadButton('downloadMcompositedist', 'Download composite M for resampling')
          )
        )
    ) 
)
)