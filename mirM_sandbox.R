# mirM_sandbox.R
# -------------------------------------------------------------------
# Copyright 2011 Christiaan Klijn <c.klijn@gmail.com>
# Project: microRNA analyses Martine
# Description: Sandbox code
# -------------------------------------------------------------------

# To Do
# define comparison groups
# eliminate non-variant mirs
# code to plot SAM plot and a heatmap of the significant probes
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

# do a sample sam analysis
# Check all lung mets vs. their donor tumors

arrayLung <- subset(hybInfo, sampleType == 'Lung met')$Array
arrayNonLung <- subset(hybInfo, sampleType != 'Lung met')$Array

lungResults <- doMirSamAnalysis(arrayLung, arrayNonLung, 
  measureData, probeInfo)
plotSamFigures(lungResults, measureData)

