library(shiny)
library(plotly)

shinyUI(
  fluidPage(
    
    # Application title
    title = "Posttest probability of pneumococcal infection in HIV patients",
    
    fluidRow(column(6, plotlyOutput('continousPlot')), column(6, plotlyOutput('combinedPlot'))),
    
    br(),
    
    p("This widget provides the posttest biomarker probability of S. pneumoniae infection in HIV patients with confirmed lobar pneumonia. You can input your test result to compute the posttest probability. I smooth the upper end to help with edge effects. You can also select which test you ordered and change the pretest probability. However, the continuous likelihood ratios use the pretest probability as well as logistic regression parameters from the actual data. Thus, the posttest probability will not change with pretest changes. The combined posttest probabilities will change, but they should be intrepretted cautiously. As the pretest probability is user-manipulated away from the empirical value, the greater the constraints of the dataset may be violated."),
    
    br(),
    
    fluidRow(
      column(6,
             wellPanel(
               radioButtons("radioBiomarker",
                            label = "Which biomarker?",
                            choices = list("C-reactive protein (CRP)" = 1, "Procalcitonin (PCT)" = 2, "Nasopharyngeal lytA" = 3), 
                            selected = 1)
             )
      ),
      column(6,
             wellPanel(
               radioButtons("radioType",
                            label = "Continuous or dichotomous?",
                            choices = list("Continous" = 1, "Dichotomous cutoff" = 2), 
                            selected = 1)
             )
      )
    ),
    
    br(),
    
    fluidRow(
      column(4,
             numericInput("marker.value", "Biomarker Value",  
                         min = 0, max = 100, value = 0, width = '100%')
      ),
      column(4,
             sliderInput("pretest", "Prestest",  
                          min = 0, max = 1, step = 0.001, value = 0.35, width = '100%')
      ),
      column(4,
             sliderInput("cutoff", "Cutoff",  
                         min = 0, max = 1, step = 0.001, value = 0.5, width = '100%')
      )
    )
  )
)
