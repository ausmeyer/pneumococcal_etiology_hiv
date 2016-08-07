library(shiny)
library(plotly)
library(shinythemes)

shinyUI(
  fluidPage(
    
    # Application title
    titlePanel("Posttest probability of pneumococcal infection in HIV patients"),
    
    sidebarLayout(
      sidebarPanel(
        radioButtons("radioBiomarker",
                     label = "Which biomarker?",
                     choices = list("C-reactive protein (CRP)" = 1, "Procalcitonin (PCT)" = 2, "Nasopharyngeal lytA" = 3), 
                     selected = 1),
#         radioButtons("radioType",
#                      label = "Continuous or dichotomous?",
#                      choices = list("Continous" = 1, "Dichotomous cutoff" = 2), 
#                      selected = 1),
        numericInput("marker.value", "Biomarker Value",  
                     min = 0, max = 100, value = 0, width = '100%'),
        sliderInput("pretest", "Prestest",  
                    min = 0, max = 1, step = 0.001, value = 0.35, width = '100%')#,
#         sliderInput("cutoff", "Cutoff",  
#                     min = 0, max = 1, step = 0.001, value = 0.5, width = '100%')
      ),
      
      mainPanel(
        
        fluidRow(column(6, plotlyOutput('continousPlot')), column(6, plotlyOutput('combinedPlot'))),
        
        br()
      )
    ),
    p("This app provides per biomarker posttest probability of S. pneumoniae infection in HIV patients with confirmed lobar pneumonia. You can input your test result to compute the posttest probability. Hover over the figure at various points to get exact values. I smooth the upper end to help with edge effects. You can also select which test you ordered and change the pretest probability. However, the continuous likelihood ratios use the pretest probability as well as logistic regression parameters from the actual data. Thus, the posttest probability will not change with pretest changes. The combined posttest probabilities will change, but they should be intrepretted cautiously. As the pretest probability is user-manipulated away from the empirical value, the greater the constraints of the dataset may be violated. I plan to add cutoffs in the near future.")
  )
)
