# mirM_functions.R
# -------------------------------------------------------------------
# Copyright 2011 Christiaan Klijn <c.klijn@gmail.com>
# Project: microRNA analyses Martine
# Description: Functions for project
# -------------------------------------------------------------------

findMinDelta <- function(samObj, fdrThres=.1) {
  
  require(siggenes)

  # Retrieve delta associated with the first FDR value
  # below a given threshold
    
  samObjFDR <- samObj@mat.fdr[samObj@mat.fdr[,'FDR'] < fdrThres,]
  minDelta <- min(samObjFDR[,'Delta'])

  return(minDelta)

}

doMirSamAnalysis <- function(array1, array2, 
  measureData, probeInfo, fdrThres=.1) {
  
  require(siggenes)

  # Wrapper for a sam analysis
  # Input is two vectors of array IDs

  # To Do:
  # Merge samProbes and samSummary@mat.sig

  # Run SAM, array1 is 1, array2 is 0 in the class vector
  samResult <- sam(measureData, colnames(measureData) %in% array1)
  sigDelta <- findMinDelta(samResult)
  samSummary <- summary(samResult, delta=sigDelta)
  samProbes <- probeInfo[samSummary@row.sig.genes,]

  # Add SAM info to the probes
  tempMat <- samSummary@mat.sig
  rownames(tempMat) <- rownames(probeInfo)[tempMat$Row]
  samProbes <- merge(samProbes, tempMat, by='row.names')
  # The merge moves the probe names to a column.
  # Restore the row names and remove the additional column
  rownames(samProbes) <- samProbes$Row.names
  samProbes <- samProbes[,!grepl('Row.names',colnames(samProbes))] 

  returnResult <- list(delta=sigDelta,
    samResult=samResult,
    samSummary=samSummary,
    samProbes=samProbes)

  return(returnResult)
  
}

plotSamFigures <- function(samList, measureData,
  dirName='Figures/', fileName='test') {
  
  # Plots SAM plot and a heatmap into a pdf file
  # input is a samResult list

  plotFile <- paste(dirName, fileName, '.pdf', )

  pdf(file=plotFile, width=5, height=5)

  plot(samList$samResult, samList$delta)

    
}

