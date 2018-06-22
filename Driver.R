# Example for Jags-Ybinom-XnomSsubjCcat-MbinomBetaOmegaKappa.R 
#------------------------------------------------------------------------------- 
# Optional generic preliminaries:
graphics.off() # This closes all of R's graphics windows.
rm(list=ls())  # Careful! This clears all of R's memory!
#------------------------------------------------------------------------------- 

library(dplyr)

# initialize the working directory for the final project
jim_dir = "/Users/jimgrund/Documents/GWU/Bayesian_Methods/FinalProject/"
akash_dir = "C:/Users/akash/Desktop/GWU/6450_Bayesian/FinalProject"
patrick_dir = "/Users/pjordan/Documents/GWU/6450/FinalProject"

for (directory in c("akash_dir", "jim_dir", "patrick_dir")) {
   if ( dir.exists(directory) ) {
      setwd(directory)
      break
   }
}

# Load the Breast Cancer dataset
myData = read.csv("data.csv")

# create three bins using 0, 12.25, 14.75, 30 as the breakpoints.  
myData$bin <- cut(myData$radius_mean, breaks = c(0, 12.25, 14.75, 30), labels = 1:3)
myData$tests <- 1
myData <- mutate(myData, diagnosis_code = ifelse(diagnosis == 'B',0,1))

# plot a histogram just to see the distribution of those bins
hist(as.numeric(myData$bin))


#------------------------------------------------------------------------------- 
# Load the relevant model into R's working memory:
source("Model.R")
#------------------------------------------------------------------------------- 
# Optional: Specify filename root and graphical format for saving output.
# Otherwise specify as NULL or leave saveName and saveType arguments 
# out of function calls.
fileNameRoot = "BreastCancer-POST-" 
graphFileType = "eps" 
#------------------------------------------------------------------------------- 
# Generate the MCMC chain:
startTime = proc.time()
mcmcCoda = genMCMC( data=myData , 
                    zName="diagnosis_code", NName="tests", sName="id", cName="diagnosis",
                    numSavedSteps=11000 , saveName=fileNameRoot ,
                    thinSteps=20 )
stopTime = proc.time()
elapsedTime = stopTime - startTime
show(elapsedTime)
#------------------------------------------------------------------------------- 
# Display diagnostics of chain, for specified parameters:
parameterNames = varnames(mcmcCoda) # get all parameter names for reference
for ( parName in c("omega[1]","omegaO","kappa[1]","kappaO","theta[1]") ) { 
  diagMCMC( codaObject=mcmcCoda , parName=parName , 
                saveName=fileNameRoot , saveType=graphFileType )
}
#------------------------------------------------------------------------------- 
# Get summary statistics of chain:
summaryInfo = smryMCMC( mcmcCoda , compVal=NULL , 
                        diffSVec=c(75,156, 159,844) , 
                        diffCVec=c(1,2,3) , 
                        compValDiff=0.0 , saveName=fileNameRoot )
# Display posterior information:
plotMCMC( mcmcCoda , data=myData , 
          zName="diagnosis_code", NName="tests", sName="id", cName="diagnosis",
          compVal=NULL ,
          diffCList=list( c("M","B"), c("B","M") ) , 
          diffSList=list( c("842302","842517"), c("84501001","84458202") ) , 
          compValDiff=0.0, #ropeDiff = c(-0.05,0.05) ,
          saveName=fileNameRoot , saveType=graphFileType )
#------------------------------------------------------------------------------- 
