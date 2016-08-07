library(shiny)
pdf(NULL)

# Define server logic required to generate and plot a random distribution
shinyServer(function(input, output, session) {
  
  dataInput <- reactive({
    data.frame(marker.value = input$marker.value,
               pretest = input$pretest#,
               #cutoff = isolate({input$cutoff})
    )
  })
  
  dataInputCutoff <- reactive({
    data.frame(marker.value = input$marker.value,
               pretest = input$pretest,
               cutoff = input$cutoff
    )
  })
  
  ## Library Imports
  libraries.call <- c("ggplot2", "readr", 'plotly')
  lapply(libraries.call, require, character.only = TRUE)
  
  ## Utility Functions
  calc.pretest <- function(df, variable) {
    x <- sum(df[[variable]] == '1') / nrow(df)
  }
  
  calc.prob.difference <- function(df, pretest) {
    pretest.odds <- pretest / (1 - pretest)
    lr.pos <- (df$sens / (1 - df$spec))
    lr.neg <- ((1 - df$sens) / df$spec)
    pos.odds <- pretest.odds * lr.pos
    neg.odds <- pretest.odds * lr.neg
    pos.prob <- pos.odds / (1 + pos.odds)
    neg.prob <- neg.odds / (1 + neg.odds)
    
    return(data.frame(lr.pos = lr.pos, lr.neg = lr.neg, pos.prob = pos.prob, neg.prob = neg.prob))
  }
  
  generate.matrix <- function(df, cutoff, variable, response) {
    t.positive <- sum(df[[variable]] > cutoff & as.numeric(as.character(df[[response]]) == 1))
    f.positive <- sum(df[[variable]] > cutoff & as.numeric(as.character(df[[response]]) == 0))
    f.negative <- sum(df[[variable]] < cutoff & as.numeric(as.character(df[[response]]) == 1))
    t.negative <- sum(df[[variable]] < cutoff & as.numeric(as.character(df[[response]]) == 0))
    return(data.frame(t.p = t.positive, f.p = f.positive, f.n = f.negative, t.n = t.negative))
  }
  
  calculate.test.stats <- function(df) {
    sensitivity <- df$t.p / (df$t.p + df$f.n)
    specificity <- df$t.n / (df$t.n + df$f.p)
    sigma2.pos <- (1 / df$t.p) - (1 / (df$t.p + df$f.n)) + (1 / df$f.p) - (1 / (df$f.p + df$t.n))
    sigma2.neg <- (1 / df$f.n) - (1 / (df$t.p + df$f.n)) + (1 / df$t.n) - (1 / (df$f.p + df$t.n))
    
    return(data.frame(sens = sensitivity, spec = specificity, sigma2.pos = sigma2.pos, sigma2.neg = sigma2.neg))
  }
  
  calculate.empiric.lr <- function(df) {
    return(df$t.p * df$t.n / (df$t.p^2 + 2*df$f.n*df$t.p + df$f.n^2) / df$f.p / df$f.n)
  }
  
  make.test.stats.df <- function(df, which.variable) {
    df.series <- df$variable[order(df$variable)]
    df.stats <- data.frame()
    invisible(sapply(df.series, function(x) df.stats <<- rbind(df.stats, (calculate.test.stats(generate.matrix(df, cutoff = x, variable = "variable", response = which.variable))))))
    df.stats <- cbind(x.value = df.series, df.stats)
    
    return(df.stats)
  }
  
  make.empiric.lr.df <- function(df, which.variable) {
    df.series <- df$variable[order(df$variable)]
    df.stats <- data.frame()
    invisible(sapply(df.series, function(x) df.stats <<- rbind(df.stats, (calculate.empiric.lr(generate.matrix(df, cutoff = x, variable = "variable", response = which.variable))))))
    df.stats <- cbind(x.value = df.series, df.stats)
    colnames(df.stats)[2] <- 'y.value'
    
    return(df.stats)
  }
  
  make.fit <- function(df) {
    bacter.fit <- glm(factor(bacter) ~ variable, data = df, family = "binomial")
    predicted.bacteremia <- predict(bacter.fit, type="response", se = TRUE)
    
    pneumo.fit <- glm(factor(pneumo) ~ variable, data = df, family = "binomial")
    predicted.pneumo <- predict(pneumo.fit, type="response", se = TRUE)
    
    return(list(predicted = data.frame(bacter = predicted.bacteremia, pneumo = predicted.pneumo), bacter.fit = bacter.fit, pneumo.fit = pneumo.fit))
  }
  
  import.data <- function(variable) {
    df <- read_csv('/srv/shiny-server/pneumococcal_etiology_hiv/data.csv', 
                   col_types = cols(
                     study_ID = col_character(),
                     Gender = col_character(),
                     age = col_integer(),
                     Bartlett_Score = col_integer(),
                     CD4_count = col_integer(),
                     Bactrim = col_integer(),
                     HAART = col_integer(),
                     CURB65 = col_integer(),
                     MR_proANP = col_double(),
                     MR_proADM = col_double(),
                     MR_proADM = col_double(),
                     PCT = col_double(),
                     CRP = col_double(),
                     lytA = col_double(),
                     bacteremia = col_double(),
                     Pneumococcal_diagnosis = col_integer(),
                     Pneumococcal_diagnosis_expanded = col_integer()
                   ))
    
    df <- data.frame(age = df$age, variable = df[[variable]], bacter = df$bacteremia, pneumo = df$Pneumococcal_diagnosis)
    df <- na.omit(df)
    
    if(variable == 'CRP') {
      df[['variable']] <- df[['variable']] / 10
    }
    
    return(df)
  }
  
  renderCombinedPlot <- function(biomarker, xaxis) {
    
    df.raw <- import.data(variable = biomarker)
    df <- make.test.stats.df(df = df.raw, which.variable = "pneumo")
    pretest <- dataInput()$pretest
    marker.value <- dataInput()$marker.value

    prob.diff <- calc.prob.difference(df, pretest)
    df <- df[!is.na(prob.diff$pos.prob - prob.diff$neg.prob) & is.finite(prob.diff$pos.prob - prob.diff$neg.prob), ]
    prob.diff <- calc.prob.difference(df, pretest)
    df.raw <- df.raw[df.raw$variable <= max(df$x.value), ]
    fit <- make.fit(df.raw)
    
    combined.probability <- (pretest / (1 - pretest)) * prob.diff$lr.pos * prob.diff$lr.neg / 
      ((pretest / (1 - pretest)) * prob.diff$lr.pos * prob.diff$lr.neg + 1)
    
    df <- cbind(df, combined.probability)
    x.value.point <- min(which(abs(df$x.value - marker.value) == min(abs(df$x.value - marker.value))))
    
    if(x.value.point > which(df$x.value == max(df$x.value))) {
      y.value.point <- df$combined.probability[x.value.point]
      x.value.point <- max(df$x.value)
    }
    else {
      y.value.point <- df$combined.probability[x.value.point]
      x.value.point <- df$x.value[x.value.point]
    }
    
    if(which(x.value.point == df$x.value) >= length(df$x.value) - 5) {
      x.value.point <- mean(df$x.value[(which(x.value.point == df$x.value) - 5):length(df$x.value)])
      y.value.point <- df$combined.probability[min(which(abs(df$x.value - x.value.point) == min(abs(df$x.value - x.value.point))))]
    }
    
    p1 <- ggplot() + 
      geom_segment(aes(x = 0, y = pretest, xend = max(df$x.value), yend = pretest, color = 'a'), linetype = 2) +
      geom_point(data = df, aes(x = x.value, y = combined.probability, color = "b"), size = 0.5) +
      geom_point(data = data.frame(x = 0, y = pretest), aes(x = x, y = y, color = 'a')) +
      geom_point(data = data.frame(x = max(df$x.value), y = pretest), aes(x = x, y = y, color = 'a')) +
      geom_segment(aes(x = x.value.point, y = 0, xend = x.value.point, yend = y.value.point, color = 'c'), linetype = 2) +
      geom_point(data = data.frame(x = x.value.point, y = y.value.point), aes(x = x, y = y, color = 'c')) +
      
      ylim(0, 1) + 
      xlab(xaxis) +
      ylab("probability pneumococcal") +
      scale_color_discrete(name = "", labels = c(a = 'pretest probability', b = 'full posttest', c = 'posttest probability')) +
      theme_bw() +
      theme(legend.position = c(0.29, 0.86), 
            legend.title=element_blank(),
            legend.background = element_rect(fill = alpha('white', 0.0)),
            plot.margin = unit(c(0, 0, 0.5, 0.5), "lines"))

    df.raw.tmp <- df.raw
    
    pretest.odds <- pretest / (1 - pretest)
    x.1 <- (log(pretest.odds) - coef(fit[['pneumo.fit']])[[1]]) / coef(fit[['pneumo.fit']])[[2]]
    tmp.lr.continuous <- exp(coef(fit[['pneumo.fit']])[[2]] * (df.raw.tmp$variable - x.1))
    posttest.odds <- pretest.odds * tmp.lr.continuous
    posttest.prob <- posttest.odds / (1 + posttest.odds)
    df.raw.tmp <- cbind(df.raw.tmp, tmp = posttest.prob)
    colnames(df.raw.tmp)[colnames(df.raw.tmp) == 'tmp'] <- paste('p_', pretest, sep = '')
    
    df.raw.tmp <- df.raw.tmp[, c(-1, -3, -4)]
    colnames(df.raw.tmp) <- c('x.value', 'y.value')
    df.raw.tmp <- df.raw.tmp[order(df.raw.tmp$x.value), ]
    
    x.value.point.log <- min(which(abs(df.raw.tmp$x.value - marker.value) == min(abs(df.raw.tmp$x.value - marker.value))))
    if(x.value.point.log > which(df.raw.tmp$x.value == max(df.raw.tmp$x.value))) {
      y.value.point.log <- df.raw.tmp$y.value[x.value.point.log]
      x.value.point.log <- max(df.raw.tmp$x.value)
    }
    else {
      y.value.point.log <- df.raw.tmp$y.value[x.value.point.log]
      x.value.point.log <- df.raw.tmp$x.value[x.value.point.log]
    }
    
    p2 <- ggplot() + 
      geom_segment(aes(x = 0, y = pretest, xend = max(df.raw.tmp$x.value), yend = pretest, color = 'a'), linetype = 2) +
      geom_point(data = df.raw.tmp, aes(x = x.value, y = y.value, color = "b"), size = 0.5) +
      geom_point(data = data.frame(x = 0, y = pretest), aes(x = x, y = y, color = 'a')) +
      geom_point(data = data.frame(x = max(df.raw.tmp$x.value), y = pretest), aes(x = x, y = y, color = 'a')) +
      geom_segment(aes(x = x.value.point.log, y = 0, xend = x.value.point.log, yend = y.value.point.log, color = 'c'), linetype = 2) +
      geom_point(data = data.frame(x = x.value.point.log, y = y.value.point.log), aes(x = x, y = y, color = 'c')) +
      ylim(0, 1) + 
      xlab(xaxis) +
      ylab("probability pneumococcal") +
      scale_color_discrete(name = "", labels = c(a = 'pretest probability', b = 'continuous posttest', c = 'posttest probability')) +
      theme_bw() +
      theme(legend.position = c(0.29, 0.86), 
            legend.title=element_blank(),
            legend.background = element_rect(fill = alpha('white', 0.0)),
            plot.margin = unit(c(0, 0, 0.5, 0.5), "lines"))
    
    p1 <- plotly_build(p1)
    p1$data[[2]]$text <- gsub("x.value", paste(biomarker, ' Concentration', sep = ''), p1$data[[2]]$text)
    p1$data[[2]]$text <- gsub("combined.probability", "Posttest Probability", p1$data[[2]]$text)
    p1$data[[2]]$text <- gsub("b: b", "", p1$data[[2]]$text)
    p1$data[[3]]$text <- paste("Pretest Probability ", pretest, sep = '')
    p1$data[[4]]$text <- paste("Pretest Probability ", pretest, sep = '')
    p1$data[[5]]$text <- ""
    p1$data[[6]]$text <- gsub("x", paste(biomarker, ' Concentration', sep = ''), p1$data[[6]]$text)
    p1$data[[6]]$text <- gsub("y", "Posttest Probability", p1$data[[6]]$text)
    p1$data[[6]]$text <- gsub("c: c", "", p1$data[[6]]$text)
    p1 <- plot_ly(p1)
    
    p2 <- plotly_build(p2)
    p2$data[[2]]$text <- gsub("x.value", paste(biomarker, ' Concentration', sep = ''), p2$data[[2]]$text)
    p2$data[[2]]$text <- gsub("y.value", "Posttest Probability", p2$data[[2]]$text)
    p2$data[[2]]$text <- gsub("b: b", "", p2$data[[2]]$text)
    p2$data[[3]]$text <- paste("Pretest Probability ", pretest, sep = '')
    p2$data[[4]]$text <- paste("Pretest Probability ", pretest, sep = '')
    p2$data[[5]]$text <- ""
    p2$data[[6]]$text <- gsub("x", paste(biomarker, ' Concentration', sep = ''), p2$data[[6]]$text)
    p2$data[[6]]$text <- gsub("y", "Posttest Probability", p2$data[[6]]$text)
    p2$data[[6]]$text <- gsub("c: c", "", p2$data[[6]]$text)
    p2 <- plot_ly(p2)
    
    return(list(plot1 = p1, plot2 = p2))
  }
  
  runBiomarker.continuous <- function() {
    if(input$radioBiomarker == 1) {
      output$continousPlot <- renderPlotly(renderCombinedPlot(biomarker = 'CRP', 
                                                              xaxis = 'CRP (mg/dL)')[['plot1']])
      output$combinedPlot <- renderPlotly(renderCombinedPlot(biomarker = 'CRP', 
                                                             xaxis = 'CRP (mg/dL)')[['plot2']])
      
    }
    if(input$radioBiomarker == 2) {
      output$continousPlot <- renderPlotly(renderCombinedPlot(biomarker = 'PCT', 
                                                              xaxis = 'PCT (ng/mL)')[['plot1']])
      output$combinedPlot <- renderPlotly(renderCombinedPlot(biomarker = 'PCT', 
                                                             xaxis = 'PCT (ng/mL)')[['plot2']])
    }
    if(input$radioBiomarker == 3) {
      output$continousPlot <- renderPlotly(renderCombinedPlot(biomarker = 'lytA', 
                                                              xaxis = 'lytA (log10 copies/mL)')[['plot1']])
      output$combinedPlot <- renderPlotly(renderCombinedPlot(biomarker = 'lytA', 
                                                             xaxis = 'lytA (log10 copies/mL)')[['plot2']])
    }
  } 
  
  runBiomarker.cutoff <- function() {
    if(input$radioBiomarker == 1) {
      output$youdenPlot <- renderPlotly(renderCombinedPlot(biomarker = 'CRP', 
                                                              xaxis = 'CRP (mg/dL)')[['plot1']])
      output$bayesianPlot <- renderPlotly(renderCombinedPlot(biomarker = 'CRP', 
                                                             xaxis = 'CRP (mg/dL)')[['plot2']])
      
    }
    if(input$radioBiomarker == 2) {
      output$youdenPlot <- renderPlotly(renderCombinedPlot(biomarker = 'PCT', 
                                                              xaxis = 'PCT (ng/mL)')[['plot1']])
      output$bayesianPlot <- renderPlotly(renderCombinedPlot(biomarker = 'PCT', 
                                                             xaxis = 'PCT (ng/mL)')[['plot2']])
    }
    if(input$radioBiomarker == 3) {
      output$youdenPlot <- renderPlotly(renderCombinedPlot(biomarker = 'lytA', 
                                                              xaxis = 'lytA (log10 copies/mL)')[['plot1']])
      output$bayesianPlot <- renderPlotly(renderCombinedPlot(biomarker = 'lytA', 
                                                             xaxis = 'lytA (log10 copies/mL)')[['plot2']])
    }
  } 
  
  runLRtype <- function() {
    if(input$radioType == 1) {
      runBiomarker.continuous()
    }
    if(input$radioType == 2) {
      runBiomarker.cutoff()
    }
  }
  
  observeEvent(input$radioBiomarker, {
    runBiomarker.continuous()
  })
  
  observeEvent(input$radioType, {
    runLRtype()
  })
})
