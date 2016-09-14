 library("shiny")
library(networkD3)


shinyUI(fluidPage(
  #Load D3.js
  tags$head(
    tags$script(src = 'http://d3js.org/d3.v3.min.js')
  ),

  titlePanel("SNF clustering"),
  sidebarLayout(

    #Use sidebar in order to select the files used for SNF
    sidebarPanel(
      fluidRow(
        wellPanel(
          tags$h4("GBM data set Wang et al., Nature 201"),
 
          checkboxGroupInput(inputId="data", label="Select at least two data types",c("mRNA expresion"="glio_mrna","methylation"="glio_meth","miRNA expresion"="glio_mirna"),selected = c("glio_mrna","glio_meth","glio_mirna")),

          tags$hr(),
          tags$h4("Preprocess data"),

          checkboxInput("rm", label = "Remove outliers", value = FALSE),

           conditionalPanel(
            condition="input.rm==true"
            ,
            numericInput("rmv", "remove sample and features with more than x% missing:", 20)
          ),

          checkboxInput("im", label = "Impute", value = FALSE),

          conditionalPanel(
              condition="input.im==true"
              ,
              numericInput("imv", "Select number of neighbors:", 20)
          ),

          checkboxInput("nr", label = "Normalize", value = TRUE),

          tags$hr(),

          tags$h4("Select Parameters"),

          numericInput(inputId="K",label="Number of neighbors",value=20),
          numericInput(inputId="alpha",label="hyperparameter, usually (0.3~0.8)",value=0.5),
          numericInput(inputId="T",label="Number of Iterations",value=10),
          numericInput(inputId="clusters",label="Number of Clusters",value=3,min=1,max=10,step=1),

          sliderInput(inputId="topLinks",label="Top % of Interactions",min=0,max=1,value=0.05,step=0.01),

          tags$hr()
               ),

        wellPanel(
          downloadButton('downloadWf', 'Download similarity matrix'),
          br(""),
          downloadButton('downloadGroups', 'Download group assignments'),
          tags$br(""),
          downloadButton('downloadNw', 'Download  fused network')
          #actionButton("features",label="Select Feautures")
          )

      )
      ),
    mainPanel(

      tabsetPanel(
        tabPanel("Summary", includeHTML("Desc.html")),

        tabPanel("Run SNF",

                 h2("Display Heatmap"),
                  plotOutput("heatmap"),
                  tableOutput("contents"),
                  h2("Plot Network"),
                   forceNetworkOutput('networkPlot')
                # tableOutput("features")
                  )
       )

    )

    )

))