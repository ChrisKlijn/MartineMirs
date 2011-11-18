#R code to analyse microRNA data
#Project together with Martine van Miltenburg
_Christiaan Klijn (c.klijn@gmail.com)_

##Overview of files
The function of the scripts:

`mirM_functions`: contains all the functions used in the analyses. This script is sources by all other scripts.
`mirM_dataLoad`: loads the data from the raw text files, performs some cleanup and saves the data in Rda format.
`mirM_sandbox`: contains work-in-progress code and short example analyses as well as a To Do list in the comments.

##Running the scripts
Required to run the scripts is a working R installation with the packages `gplots` and `siggenes` installed.

To start, put the data and the code in separate folders. In the data folder, create two subfolders: Figures and Tables. Adjust the `dataDir` and `dataDir` variables in all scripts to these locations on your own machine.

