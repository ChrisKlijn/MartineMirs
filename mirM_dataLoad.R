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

# use probeID rownames and remove the MiR column 
# from the measurement data.
# Check if the probes are in the same order before doing so

all.equal(measureData$MiR, probeInfo$description)

measureData <- as.matrix(measureData[,-1])
row.names(measureData) <- probeInfo$reporterID

# Remove _S01 from data column names
colnames(measureData) <- gsub('_S01', '', colnames(measureData))

# Split the Sample.name..user field into two separate fields
sampleInfo$type <- gsub('[0-9]{6}[ ]{1}' ,'', 
  sampleInfo$Sample.name..user.)
sampleInfo$animalID <- 
  as.numeric(strtrim(sampleInfo$Sample.name..user., 6))

# add sample info to the hybridization info
hybInfo$sampleAnimal <- sampleInfo$animalID[hybInfo$Sample]
hybInfo$sampleType <- sampleInfo$type[hybInfo$Sample]
hybInfo$refAnimal <- sampleInfo$animalID[hybInfo$Reference.sample]
hybInfo$refType <- sampleInfo$type[hybInfo$Reference.sample]

save(file='mirM_rawData.Rda', 
  list=c('measureData', 'probeInfo', 'sampleInfo', 'hybInfo'))
