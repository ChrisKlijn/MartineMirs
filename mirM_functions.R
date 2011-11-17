# mirM_functions.R
# -------------------------------------------------------------------
# Copyright 2011 Christiaan Klijn <c.klijn@gmail.com>
# Project: microRNA analyses Martine
# Description: Functions for project
# -------------------------------------------------------------------

findMinDelta <- function(samObj, fdrThres=.1) {
  
  # Retrieve delta associated with the first FDR value
  # below a given threshold
    
  samObjFDR <- samObj@mat.fdr[samObj@mat.fdr[,'FDR'] < fdrThres,]
  minDelta <- min(samObjFDR[,'Delta'])

  return(minDelta)

}


