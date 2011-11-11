# mirM_dataLoad.R
# -------------------------------------------------------------------
# Copyright 2011 Christiaan Klijn <c.klijn@gmail.com>
# Project: microRNA analyses Martine
# Description: Load data, prepocess and store in Rda
# -------------------------------------------------------------------

# Change these to the local directories
# Directories shown are mediod dirs
dataDir <- '~/data/smallproj/martineMir/'
codeDir <- '~/gitCodeChris/MartineMirs/'

setwd(dataDir)

# Read measure data, note: data contains comma's as decimal indicators
measureData <- read.delim("Data sample vs reference.txt", 
  header = TRUE,  sep = "\t",  stringsAsFactors = FALSE, dec=',')
probeInfo <- read.delim("probe data.txt", 
  header = TRUE,  sep = "\t",  stringsAsFactors = FALSE)
sampleInfo <- read.delim("sample overview.txt", 
  header = TRUE,  sep = "\t",  stringsAsFactors = FALSE)
hybInfo <- read.delim("Hybridization.txt", 
  header = TRUE,  sep = "\t",  stringsAsFactors = FALSE)

save(file='mirM_rawData.Rda', 
  list=c('measureData', 'probeInfo', 'sampleInfo', 'hybInfo'))
