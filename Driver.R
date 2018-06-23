#' ---
#' title: 6450 BayesianMethods Project (Driver File)
#' author: TeamPAJ
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

#' ## Environment Preparation

# Remove any objects in the Environment
rm(list = ls())

# Closes all of R's graphics windows
graphics.off()

# Knitr global options
library(knitr)
opts_chunk$set(eval = TRUE, echo = TRUE, warning = FALSE,
               tidy = TRUE, results = "hold", cache = TRUE)

# Load Necessary Libraries
suppressPackageStartupMessages(library(dplyr))

# Set the overall seed for reproducibility
set.seed(6450)

# initialize the working directory for the final project
jim_dir = "/Users/jimgrund/Documents/GWU/Bayesian_Methods/FinalProject/"
akash_dir = "C:/Users/akash/Desktop/GWU/6450_Bayesian_YHuang/project/project"
patrick_dir = "/Users/pjordan/Documents/GWU/6450/FinalProject"

for (directory in c(akash_dir, jim_dir, patrick_dir)) {
  if ( dir.exists(directory) ) {
    setwd(directory)
    break
  }
}

# Load the relevant model into R's working memory:
source("Model.R")

#' ### Graph Options

# Optional: Specify filename root and graphical format for saving output.
# Otherwise specify as NULL or leave saveName and saveType arguments 
# out of function calls.
fileNameRoot = "graph_" 
graphFileType = "png" 

#'#######################################

#' ## Data Load & Tidy

# Load the Breast Cancer dataset
myData = read.csv("data.csv")

#' ### Radius_mean
# create three bins using 0, 12.25, 14.75, 30 as the breakpoints.  
myData$radius_bin <- cut(myData$radius_mean, breaks = c(0, 12.25, 14.75, 30), labels = 1:3)
myData$tests <- 1
myData <- mutate(myData, diagnosis_code = ifelse(diagnosis == 'B',0,1))

# plot a histogram just to see the distribution of those bins
hist(as.numeric(myData$radius_bin), col = 'blue4')


#' ### Area_mean
# create three bins on area_mean using 0, 463, 680, 2600 as the breakpoints.  
myData$area_bin <- cut(myData$area_mean, breaks = c(0, 463, 680, 2600), labels = 1:3)

# plot a histogram just to see the distribution of those bins
hist(as.numeric(myData$area_bin), main="area_bin")


#' ### Compactness_mean
# create three bins on compactness_mean using 0, 0.075, 0.117, 0.5 as the breakpoints.  
myData$compactness_bin <- cut(myData$compactness_mean, breaks = c(0, 0.075, 0.117, 0.5), labels = 1:3)

# plot a histogram just to see the distribution of those bins
hist(as.numeric(myData$compactness_bin), main="compactness_bin")


#' ### Smoothness_mean
# create three bins on smoothness_mean using 0, 0.0894, 0.102, 0.17 as the breakpoints.  
myData$smoothness_bin <- cut(myData$smoothness_mean, breaks = c(0, 0.0894, 0.102, 0.17), labels = 1:3)

# plot a histogram just to see the distribution of those bins
hist(as.numeric(myData$smoothness_bin), main="smoothness_bin")


#' ### Concavity_mean
# create three bins on concavity_mean using 0, 0.039, 0.106, 0.43 as the breakpoints.  
myData$concavity_bin <- cut(myData$concavity_mean, breaks = c(0, 0.039, 0.106, 0.43), labels = 1:3)

# plot a histogram just to see the distribution of those bins
hist(as.numeric(myData$concavity_bin), main="concavity_bin")
#'#######################################

#' ## MCMC Chain

#' ### Generate

# Generate the MCMC chain:
startTime = proc.time()
mcmcCoda = genMCMC( data=myData , 
                    zName="diagnosis_code", NName="tests", sName="id", cName="radius_bin",
                    numSavedSteps=11000 , saveName=fileNameRoot ,
                    thinSteps=20 )
stopTime = proc.time()
elapsedTime = stopTime - startTime
show(elapsedTime)


#' ### Diagnostics 

# Display diagnostics of chain, for specified parameters:
parameterNames = varnames(mcmcCoda) # get all parameter names for reference
for ( parName in c("omega[1]","omegaO","kappa[1]","kappaO","theta[1]") ) { 
  diagMCMC( codaObject=mcmcCoda , parName=parName , 
            saveName=fileNameRoot , saveType=graphFileType )
}

#' ### Summary Statistics

# Get summary statistics of chain:
# summaryInfo = smryMCMC( mcmcCoda , compVal=NULL ,
#                         diffSVec=c(75,156, 159) ,
#                         diffCVec=c(1,2,3) ,
#                         compValDiff=0.0 , saveName=fileNameRoot )

#' ### Graph

# Display posterior information:
plotMCMC( mcmcCoda , data=myData ,
          zName="diagnosis_code", NName="tests", sName="id", cName="radius_bin",
          compVal=NULL ,
          diffCList=list( c(1,2), c(2,3) ) ,
          diffSList=list( c("842302","842517"), c("84501001","84458202") ) ,
          compValDiff=0.0, #ropeDiff = c(-0.05,0.05) ,
          saveName=fileNameRoot , saveType=graphFileType )
#'#######################################
