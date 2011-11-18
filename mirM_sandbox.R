# mirM_sandbox.R
# -------------------------------------------------------------------
# Copyright 2011 Christiaan Klijn <c.klijn@gmail.com>
# Project: microRNA analyses Martine
# Description: Sandbox code
# -------------------------------------------------------------------

# To Do

# define comparison groups
# eliminate non-variant mirs
# possibly remove the human/rat-only mirs
# see if replicate probes are present in the list

# Change these to the local directories
# Directories shown are mediod dirs
dataDir <- '~/data/smallproj/martineMir/'
codeDir <- '~/gitCodeChris/MartineMirs/'
setwd(dataDir)

# Source functions
source(paste(codeDir, 'mirM_functions.R', sep=''))
# Use the siggenes library for the sam analysis
library(siggenes) # SAM
library(gplots) # for heatmap.2

# Load data
load('mirM_rawData.Rda')

#------------------------------------
# Sample SAM analysis

# Check all lung and LN mets vs. their donor tumors
# This analysis compared the arrays in which a specific metastasis 
# is hybridized to its primary tumor to all other arrays.

# Lung analysis
arrayLung <- subset(hybInfo, sampleType == 'Lung met')$Array
arrayNonLung <- subset(hybInfo, sampleType != 'Lung met')$Array
lungResults <- doMirSamAnalysis(arrayLung, arrayNonLung, 
  measureData, probeInfo)
plotSamFigures(lungResults, measureData, fileName='lung',
  legendText=c('Lung Met vs. Primary', 'Other Hyb'))
exportSigProbeTable(lungResults, fileName='lung')

# LN analysis
arrayLN<- subset(hybInfo, grepl('LN', sampleType))$Array
arrayNonLN <- subset(hybInfo, !grepl('LN', sampleType))$Array
LNResults <- doMirSamAnalysis(arrayLN, arrayNonLN, 
  measureData, probeInfo)
plotSamFigures(LNResults, measureData, fileName='LN',
  legendText=c('LN Met vs. Primary', 'Other Hyb'))
exportSigProbeTable(LNResults, fileName='LN')
