## Clear workspace and set working directory
rm(list = ls())

## Library Imports
libraries.call <- c("dplyr", 'tidyr', "ggplot2", "readxl", "reshape", 'cowplot')
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
  #df.series <- seq(0, max(df$variable), length = 1000)
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

## Function to left censor the data for multiple plots
## @agm Requires ordering the data after censoring to obtain desired results
left.censor <- function(df) {
  df <- data.frame(x = df$x, y = df$y)
  new.df <- data.frame()
  lapply(1:(nrow(df)-1), function(x) new.df <<- rbind(new.df, data.frame(x = df$x[x+1], y = df$y[x])))
  df <- rbind(df, new.df)
}

import.data <- function(variable, variable2) {
  df <- read_excel('Database_BMJ_Open_25Jul2014.xlsx', col_types = c('text', 'text', 'numeric', 
                                                                     'numeric', 'numeric', 'text', 
                                                                     'text', 'numeric', 'numeric', 
                                                                     'numeric', 'numeric', 'numeric', 
                                                                     'numeric', 'numeric', 'text', 
                                                                     'text', 'text', 'text', 
                                                                     'text'))
  df[, 1:17]
  
  if(missing(variable2)) {
    df <- data.frame(age = df$age, variable = df[[variable]], bacter = df$bacteremia, pneumo = df$Pneumococcal_diagnosis)
    df <- na.omit(df)
  }
  else {
    tmp.df <- data.frame(age = df$age, variable1 = df[[variable]], variable2 = df[[variable2]], bacter = df$bacteremia, pneumo = df$Pneumococcal_diagnosis)
    tmp.df <- na.omit(tmp.df)
    df <- data.frame(age = tmp.df$age, variable = tmp.df$variable1/tmp.df$variable2, bacter = tmp.df$bacter, pneumo = tmp.df$pneumo)
  }
  return(df)
}

make.fit <- function(df) {
  bacter.fit <- glm(factor(bacter) ~ variable, data = df, family = "binomial")
  predicted.bacteremia <- predict(bacter.fit, type="response", se = TRUE)
  
  pneumo.fit <- glm(factor(pneumo) ~ variable, data = df, family = "binomial")
  predicted.pneumo <- predict(pneumo.fit, type="response", se = TRUE)
  
  return(list(predicted = data.frame(bacter = predicted.bacteremia, pneumo = predicted.pneumo), bacter.fit = bacter.fit, pneumo.fit = pneumo.fit))
}

## Plotting Functions
make.continuous.likelihoodratio.plot <- function(df, df.raw, file.name, xaxis, pretest) {
  prob.diff <- calc.prob.difference(df, pretest)
  df <- df[!is.na(prob.diff$pos.prob - prob.diff$neg.prob) & is.finite(prob.diff$pos.prob - prob.diff$neg.prob), ]
  prob.diff <- calc.prob.difference(df, pretest)
  
  df.raw <- df.raw[df.raw$variable <= max(df$x.value), ]
  fit <- make.fit(df.raw)
  
  generate.continuous.lr <- function(tmp.df, tmp.fit, tmp.pretest) {
    x.1 <- (log(tmp.pretest / (1 - tmp.pretest)) - (tmp.fit[['pneumo.fit']])[["coefficients"]][[1]]) / (tmp.fit[['pneumo.fit']])[["coefficients"]][[2]]
    tmp.lr.continuous <- exp((tmp.fit[['pneumo.fit']])[["coefficients"]][[2]] * (tmp.df$variable - x.1))
    tmp.df <- cbind(tmp.df, tmp = tmp.lr.continuous)
    colnames(tmp.df)[colnames(tmp.df) == 'tmp'] <- paste('p_', tmp.pretest, sep = '')
    return(tmp.df)
  }
  
  df.raw.tmp <- df.raw
  df.raw.tmp <- generate.continuous.lr(df.raw.tmp, fit, pretest)
  df.raw.tmp <- df.raw.tmp[, c(-1, -3, -4)]
  colnames(df.raw.tmp) <- c('x.value', 'y.value')
  
  p1 <- ggplot() + 
    geom_segment(data = df.raw.tmp, aes(x = 0, y = 1, xend = max(x.value), yend = 1), linetype = 2) +
    geom_point(data = df.raw.tmp, aes(x = x.value, y = y.value), size = 0.5) +
    ylim(0, 10) +
    xlab(xaxis) +
    ylab("likelihood ratio") +
    scale_color_discrete(name = "") +
    theme_bw() +
    theme(legend.position = "none")
  
  df.raw.tmp <- df.raw
  lapply(seq(0.05, 0.95, length = 18), function(x) df.raw.tmp <<- generate.continuous.lr(df.raw.tmp, fit, x))
  
  df.raw.tmp <- df.raw.tmp[, c(-1, -3, -4)]
  colnames(df.raw.tmp)[1] <- 'x.value'
  df.raw.tmp <- melt(df.raw.tmp, id = c('x.value'))
  
  p2 <- ggplot() + 
    geom_segment(data = df.raw.tmp, aes(x = 0, y = 1, xend = max(x.value), yend = 1), linetype = 2) +
    geom_point(data = df.raw.tmp, aes(x = x.value, y = value, color = variable), size = 0.5) +
    ylim(0, 10) +
    xlab(xaxis) +
    ylab("likelihood ratio") +
    scale_color_discrete(name = "") +
    theme_bw() +
    theme(legend.position = "none")
  
  return(list(plot1 = p1, plot2 = p2))
}

make.combined.probability.plots <- function(df.raw, df, file.name, xaxis, pretest) {
  prob.diff <- calc.prob.difference(df, pretest)
  df <- df[!is.na(prob.diff$pos.prob - prob.diff$neg.prob) & is.finite(prob.diff$pos.prob - prob.diff$neg.prob), ]
  prob.diff <- calc.prob.difference(df, pretest)
  df.raw <- df.raw[df.raw$variable <= max(df$x.value), ]
  fit <- make.fit(df.raw)
  
  combined.probability <- (pretest / (1 - pretest)) * prob.diff$lr.pos * prob.diff$lr.neg / 
    ((pretest / (1 - pretest)) * prob.diff$lr.pos * prob.diff$lr.neg + 1)
  
  df <- cbind(df, combined.probability)
  
  p1 <- ggplot() + 
    geom_segment(aes(x = 0, y = pretest, xend = max(df$x.value), yend = pretest, color = 'a'), linetype = 2) +
    geom_point(data = df, aes(x = x.value, y = combined.probability, color = "b"), size = 0.5) +
    geom_point(data = data.frame(x = 0, y = pretest), aes(x = x, y = y, color = 'a')) +
    geom_point(data = data.frame(x = max(df$x.value), y = pretest), aes(x = x, y = y, color = 'a')) +
    ylim(0, 1) + 
    xlab(xaxis) +
    ylab("probability pneumococcal") +
    scale_color_discrete(name = "", labels = c(a = 'pretest probability', b = 'combined posttest')) +
    theme_bw() +
    theme(legend.position = c(0.29, 0.86), 
          legend.title=element_blank(),
          legend.background = element_rect(fill = alpha('white', 0.0)))
  
  generate.continuous.lr <- function(tmp.df, tmp.fit, tmp.pretest) {
    x.1 <- (log(tmp.pretest / (1 - tmp.pretest)) - (tmp.fit[['pneumo.fit']])[["coefficients"]][[1]]) / (tmp.fit[['pneumo.fit']])[["coefficients"]][[2]]
    tmp.lr.continuous <- exp((tmp.fit[['pneumo.fit']])[["coefficients"]][[2]] * (tmp.df$variable - x.1))
    tmp.posttest <- (tmp.pretest / (1 - tmp.pretest)) * tmp.lr.continuous / (1 + (tmp.pretest / (1 - tmp.pretest)) * tmp.lr.continuous)
    tmp.df <- cbind(tmp.df, tmp = tmp.posttest)
    colnames(tmp.df)[colnames(tmp.df) == 'tmp'] <- paste('p_', tmp.pretest, sep = '')
    return(tmp.df)
  }
  
  df.raw.tmp <- df.raw
  df.raw.tmp <- generate.continuous.lr(df.raw.tmp, fit, pretest)
  df.raw.tmp <- df.raw.tmp[, c(-1, -3, -4)]
  colnames(df.raw.tmp) <- c('x.value', 'y.value')
  
  p2 <- ggplot() + 
    geom_segment(aes(x = 0, y = pretest, xend = max(df.raw.tmp$x.value), yend = pretest, color = 'a'), linetype = 2) +
    geom_point(data = df.raw.tmp, aes(x = x.value, y = y.value, color = "b"), size = 0.5) +
    geom_point(data = data.frame(x = 0, y = pretest), aes(x = x, y = y, color = 'a')) +
    geom_point(data = data.frame(x = max(df.raw.tmp$x.value), y = pretest), aes(x = x, y = y, color = 'a')) +
    ylim(0, 1) + 
    xlab(xaxis) +
    ylab("probability pneumococcal") +
    scale_color_discrete(name = "", labels = c(a = 'pretest probability', b = 'continuous posttest')) +
    theme_bw() +
    theme(legend.position = c(0.29, 0.86), 
          legend.title=element_blank(),
          legend.background = element_rect(fill = alpha('white', 0.0)))
  
  p3 <- ggplot() + 
    geom_segment(aes(x = 0, y = pretest, xend = max(df$x.value), yend = pretest, color = 'a'), linetype = 2) +
    geom_point(data = df, aes(x = x.value, y = combined.probability, color = "b"), size = 0.5) +
    geom_point(data = df.raw.tmp, aes(x = x.value, y = y.value, color = "c"), size = 0.5) +
    geom_point(data = data.frame(x = 0, y = pretest), aes(x = x, y = y, color = 'a')) +
    geom_point(data = data.frame(x = max(df$x.value), y = pretest), aes(x = x, y = y, color = 'a')) +
    ylim(0, 1) + 
    xlab(xaxis) +
    ylab("probability pneumococcal") +
    scale_color_discrete(name = "", labels = c(a = 'pretest', b = 'combined', c = 'continuous')) +
    theme_bw() +
    theme(legend.position = c(0.2, 0.82), 
          legend.title=element_blank(),
          legend.background = element_rect(fill = alpha('white', 0.0)))
  
  return(list(plot1 = p1, plot2 = p2, plot3 = p3))
}

## Avoid global variable by running a main function
run.main.analysis <- function() {
  raw.pneumo <- import.data(variable = "Pneumococcal_diagnosis")
  pretest.probability <- calc.pretest(raw.pneumo, 'pneumo')
  
  ## lytA setup
  lytA.filename <- 'lytA.pdf'
  lytA.xaxis <- 'lytA density (log10 copies/mL)'
  raw.lytA <- import.data(variable = 'lytA_NP')
  
  ## Plot lytA
  roc.data.lytA <- make.test.stats.df(df = raw.lytA, which.variable = "pneumo")

  p.continuous.lr.lytA <- make.continuous.likelihoodratio.plot(df = roc.data.lytA, df.raw = raw.lytA, file.name = lytA.filename, xaxis = lytA.xaxis, pretest = pretest.probability)
  p.combined.prob.lytA <- make.combined.probability.plots(df.raw = raw.lytA, df = roc.data.lytA, file.name = lytA.filename, xaxis = lytA.xaxis, pretest = pretest.probability)

}

run.main.analysis()