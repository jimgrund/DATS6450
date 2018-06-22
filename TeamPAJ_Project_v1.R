#' ---
#' title: 6450 BayesianMethods Project
#' author: AJassal
#' date: 25June2018
#' output:
#'    html_document:
#'      toc: true
#'      highlight: haddock
#' ---

# 
#
#

#'#######################################
#'
#' ### Environment Preparation
# Remove any objects in the Environment
rm(list = ls())

# Set working directory
setwd("C:/Users/akash/Desktop/GWU/6450_Bayesian_YHuang/project/project")
#source("DBDA2E-utilities.R")

# Knitr global options
library(knitr)
opts_chunk$set(eval = TRUE, echo = TRUE, warning = FALSE,
               tidy = TRUE, results = "hold", cache = TRUE)

# Set the overall seed for reproducibility
set.seed(6450)

# library
library(caret)

#'#######################################

#' ### Load Data
# 
loadData = function(testResultsFile) {
  testResultsFile <- data.frame(read.csv(testResultsFile)) # Read and imports the file 
  return(testResultsFile)
}

#' ## Driver
# Test results file
testResultsFile = "data.csv"
# Takes as input the test results file and outputs the results as a dataframe
myData <- loadData(testResultsFile)
str(myData)

#' ### Transform Data
# 
transformData1 = function(trans1) {
  trans1 <- subset(myData, select = c(2:32)) 
  return(trans1)
}

#' ## Driver
myData <- transformData(trans1)
myData <- myData[complete.cases(myData),] # cgeck if trans1 worked correctly


#'#######################################

#' ### Feature Selection
# 

#' Fit a logistic regression model
sum(myData$diagnosis == 'M') # check if dummization works
myData$diagnosis <- ifelse(myData$diagnosis == 'M',1,0)
sum(myData$diagnosis == 1) # check if dummization works
fit_glm <- glm(myData$diagnosis ~ ., family = 'binomial', data = myData)

library(caret)
glm_fi <- varImp(fit_glm)

# Fit a random forest model
library(randomForest)
fit_rf <- randomForest(myData$diagnosis ~ . , data = myData)
importance(fit_rf)

#===============================================================================









