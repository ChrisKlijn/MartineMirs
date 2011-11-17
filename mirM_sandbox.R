# mirM_sandbox.R
# -------------------------------------------------------------------
# Copyright 2011 Christiaan Klijn <c.klijn@gmail.com>
# Project: microRNA analyses Martine
# Description: Sandbox code
# -------------------------------------------------------------------

# To Do
# Use SAM to find significant miRs
# Find optimal delta
# define comparison groups
# eliminate non-variant mirs

# Change these to the local directories
# Directories shown are mediod dirs
dataDir <- '~/data/smallproj/martineMir/'
codeDir <- '~/gitCodeChris/MartineMirs/'
setwd(dataDir)

# Source functions
source(paste(codeDir, 'mirM_functions.R', sep=''))
# Use the siggenes library for the sam analysis
library(siggenes)

# Functions


# Load data
load('mirM_rawData.Rda')

# do a sample sam analysis
# Check all lung mets vs. their donor tumors

arrayLung <- subset(hybInfo, sampleType == 'Lung met')$Array
arrayNonLung <- subset(hybInfo, sampleType != 'Lung met')$Array

lungSam <- sam(measureData, colnames(measureData) %in% arrayLung)
lungDelta <- findMinDelta(lungSam)
lungSamSummary <- summary(lungSam, delta=lungDelta)
lungSamProbes <- probeInfo[lungSamSummary@row.sig.genes,]

