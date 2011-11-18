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
  classVect <- colnames(measureData) %in% array1
  samResult <- sam(measureData, classVect)
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
    samProbes=samProbes,
    classVect=classVect)

  return(returnResult)
  
}

plotSamFigures <- function(samList, measureData,
  dirName='Figures/', fileName='test',
  legendText=c('class 1', 'class 2')) {
  
  # Plots SAM plot and a heatmap into a pdf file
  # input is a samResult list

  plotFile <- paste(dirName, fileName, '.pdf', sep='')

  pdf(file=plotFile, width=10, height=10)

  plot(samList$samResult, samList$delta)

  heatCols <- colorpanel(256, low='blue', high='yellow')
  classCol <- colors()[c(134, 29)]

  heatmap.2(measureData[rownames(samList$samProbes),], 
    scale='none', trace='none',
    col=heatCols, margin=c(7,10), cexRow=.8, cexCol=.8,
    ColSideColors=classCol[samList$classVect + 1])
  legend('topright', legend=legendText, fill=classCol)
  
  dev.off()
}

