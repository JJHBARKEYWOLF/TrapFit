#!/usr/bin/Rscript
options(warn=2) #turns warnings into errors

helpPrint <- function(){cat("
Fits a trapezoid function to block stimulus timeseries data. 
A slight variation on a method first described by Dumas A. et al., ANN NEUROL 2012 (DOI: 10.1002/ana.23566)
Input data should be like fslmeants output: one .txt file per subject,
 each containing one timeseries where each row corresponds to a timepoint.

Arguments:
Required:
    --path           - Absolute or relative path containing the timeseries .txt files (e.g. fslmeants output, 1 file per subject). Text file names will be used as the subject names.
    --TR             - TR(s)
    --onduration     - stimulus duration
    --blockduration  - duration of the block (stimulus + rest)
    --outputname     - base name of the output files (e.g. name of the experiment)

Optional:
    --CVthresh       - A coefficient of variation (CV, aka RSD) threshold above which blocks will be classified as outliers. These blocks will be removed from the analysis (defaults to no threshold).
    --restfirst      - Whether the blocks start with a rest period.'TRUE' or 'FALSE', default is FALSE. If set to TRUE, the script discards the first period of the input data.
    
    --help           - print this text

Example:
    Rscript ./trapfit.R --path='./timeseries/' --TR=3 --onduration=20 --blockduration=48 --CVthresh=3 --restfirst=FALSE --outputname='test'\n\n"
)
  
  q(save="no")
}

## Collect arguments
args <- commandArgs(TRUE)

## Help section
if("--help" %in% args) {
  helpPrint()
}

# Parse arguments (we expect the form --arg=value)
parseArgs <- function(x) strsplit(sub("^--", "", x), "=")
argsDF <- as.data.frame(do.call("rbind", parseArgs(args)))
argsL <- as.list(as.character(argsDF$V2))
names(argsL) <- argsDF$V1
argsL$TR <- as.numeric(argsL$TR)
argsL$onduration <- as.numeric(argsL$onduration)
argsL$blockduration <- as.numeric(argsL$blockduration)
if(!is.null(argsL$CVthresh)) {
  argsL$CVthresh <- as.numeric(argsL$CVthresh)
}


# check essential arguments
if(is.null(argsL$path) | is.null(argsL$TR) | is.null(argsL$onduration) | is.null(argsL$outputname) | is.null(argsL$blockduration)) {
  print('missing required arguments..')
  helpPrint()
}

# set defaults
if(is.null(argsL$restfirst)) {
  argsL$restfirst <- FALSE
}else{
  argsL$restfirst <- as.logical(argsL$restfirst)
}


#print settings
print('running with the following settings:')
print(paste0('data directory: ',argsL$path))
print(paste0('TR(s): ',argsL$TR))
print(paste0("'on'-duration(s): ", argsL$onduration))
print(paste0('block duration(s): ', argsL$blockduration))
print(paste0('blocks start with rest period: ', argsL$restfirst))
print(paste0('output names prefixed with: ', argsL$outputname))
if(!is.null(argsL$CVthresh)) {
  print(paste0('blocks with a Coefficient of Variation > ', argsL$CVthresh, '% will be removed'))
}

#import functions
source("functions.R")


###################################################### Import data
file.names <- dir(argsL$path, pattern =".txt")

print(paste0("found ",length(file.names), " input files:"))
print('importing....')
dataset <- list()
for(i in 1:length(file.names)) {
  
  rawTSmatrix <- read.table(paste(argsL$path,file.names[i],sep = ''),quote="\"", comment.char="", stringsAsFactors=FALSE)
  subjname <- unlist(strsplit(file.names[i],"[.]"))[1]
  
  #remove first three rows: these are the voxel coordinates
  dataset[[subjname]] <- rawTSmatrix
}
print('data import done')


###################################################### Process data
#process data
perblock <- list()
blockavgs <- list()
plots_blockspersub <- list()
fits <- list()
rawTSplots <- list()
outlierPlots <- list()
CV_all <- list()

for (i in 1:length(dataset)){
  subjname <- names(dataset)[i]
  print(paste0("--------------------------processing ", subjname, "....."))

  ROIavg <- unlist(dataset[i])
  
  offduration <- (argsL$blockduration - argsL$onduration)
  #when a sessions starts with rest period:
  if (argsL$restfirst == TRUE){
    restdurationInTR <- ceiling((offduration/argsL$TR))
    ROIavg <- ROIavg[-c(1:restdurationInTR)]
  }

  #plot the entire ROI averaged timeseries
  rawTSplots[[subjname]] <- plotRawTS(timeseries=ROIavg, TR = argsL$TR, subjname = subjname, onduration = argsL$onduration, blockduration = argsL$blockduration)
  
  #split blocks
  blocks <- splitBlocks(ROIavg,blockduration = argsL$blockduration, TR = argsL$TR)
  
  #assemble splitted blocks into a single df
  splittedBlocks <- assemble.df(blocks,TR = argsL$TR)
  
  if(!is.null(argsL$CVthresh)) {
    #remove outlier blocks using arbitrary criterium (cv(TS+min(TS) > CVcutoff%)
    splittedBlocks <- removeOutlierBlocks(splittedBlocks,CVcutoff = argsL$CVthresh, subjname = subjname)
    if (length(levels(as.factor(splittedBlocks$block))) <= length(blocks)/2){
      print("###WARNING####")
      print("(MORE THAN) HALF OF BLOCKS CLASSIFIED AS OUTLIER: CANCELLING FIT")
      print("###WARNING####")
      next
    }
  }

  #calculate %BOLD change using the mean of all blocks together
  splittedBlocks$value <- percBOLDchange(splittedBlocks$value)
  
  #zero mean complicates modelling; add constant to all datapoints based on minimum value over all blocks
  splittedBlocks$value<- splittedBlocks$value - min(splittedBlocks$value)
  
  #calculate the average %BOLD change per timepoint by averaging data of all blocks
  df_avgTS <- avgOverBlocks(splittedBlocks)
  df_avgTS$subject <- rep(subjname,length(df_avgTS$value))
  
  #estimate initial guess for starting parameters from the data (based on earlier experiments)
  params <- c(a = 0.01, b = 5, c = argsL$onduration, d = (argsL$onduration + 7), baseline = min(df_avgTS$value), ceiling = max(df_avgTS$value))

  print('initial guesstimate parameters:')
  print(params)

  #run curve fitting
  print(paste0("fitting curve for: ", subjname))
  #fitresult <- trapezoidFit(df_avgTS,params) # fit on avg ts
  fitresult <- trapezoidFit(splittedBlocks,params)
  annot <- paste0("RSS:",round(fitresult$value,1))

  #check for convergence
  if (fitresult$convergence != 0){
    print(paste0('convergence did not occur for',subjname))
    annot <- paste0("No convergence, double check fit! (RSS ",round(fitresult$value,1),' )')
  }
  print(paste0('final RSS: ',fitresult$value))

  #calculate best fit values  for use in plots
  bestFitValues <- generateValuesForFit(fitresult$par,xbound = max(splittedBlocks$timepoint),TR = argsL$TR, xresolution = 0.2)
  
  #plot best fit values over the timeseries of all blocks + their average
  plots_blockspersub[[subjname]] <- plotResponsePerBlock(splittedBlocks,avgTS = df_avgTS,onduration=argsL$onduration,subjname=subjname,
                                                          bestfit=bestFitValues, annotation = annot)
  
  #save timeseries and fit for reuse
  fits[[subjname]] <- fitresult
  
  perblock[[subjname]] <- splittedBlocks
  blockavgs[[subjname]] <- df_avgTS
  
  print(paste0("DONE for ", subjname))
}


################################################################### OUTPUTS
print('writing outputs to files....')
#make table with results
options(scipen=999) #effectively disable scientific notation
df <- data.frame(do.call(rbind, fits))
outputtable <- data.frame(do.call(rbind,lapply(df$par,unlist)))
outputtable <- cbind(outputtable,RSS = do.call(rbind,lapply(df$value,unlist)))



#calculate physiologically relevant parameters
outputtable$TTB <- outputtable$d - outputtable$c
outputtable$TTP <- outputtable$b - outputtable$a
outputtable$amplitude <- outputtable$ceiling - outputtable$baseline

#combine all average timeseries in a dataframe for future use in plotting (only useful when running interactively)
allavgs <- do.call(rbind,blockavgs)

#output table and settings
write.csv(outputtable, file = paste0(argsL$outputname,"_fitparameters.csv"))
write.csv(data.frame(argsL), file = paste0(argsL$outputname,"_trapfit_settings.csv"))

#collect CV values
if(!is.null(argsL$CVthresh)) {
  CV_df <- do.call(rbind, CV_all)
  rownames(CV_df) <- seq_along(CV_df$value)
  write.csv(CV_df, file = paste0(argsL$outputname,"CV_per_block.csv"))
}

#render all plots and output to pdf
pdf(file = paste0(argsL$outputname,"_allfits.pdf"),paper = 'a4r', width = 11.5, height = 8)
plotBatches <- split(plots_blockspersub, ceiling(seq_along(plots_blockspersub)/20))
for (i in 1:length(plotBatches)){
  do.call(grid.arrange,c(plotBatches[i][[1]],ncol = 5, nrow=4))
}
dev.off()

pdf(file = paste0(argsL$outputname,"_raw_inputs.pdf"),paper = 'a4', width = 11.5, height = 8)
plotBatches <- split(rawTSplots, ceiling(seq_along(rawTSplots)/5))
for (i in 1:length(plotBatches)){
  do.call(grid.arrange,c(plotBatches[i][[1]],ncol=1,nrow=5))
}
dev.off()

if (length(outlierPlots) > 0){
  print(paste0(length(outlierPlots),' outlier blocks detected, see *_outlierblocks.pdf. '))
  pdf(file = paste0(argsL$outputname,"_outlierblocks.pdf"),paper = 'a4',width = 11.5, height = 8)
  plotBatches <- split(outlierPlots, ceiling(seq_along(outlierPlots)/3))
  for (i in 1:length(plotBatches)){
    do.call(grid.arrange,c(plotBatches[i][[1]],ncol=1,nrow=3))
  }
  dev.off()
}else{
  print('no outlier blocks detected at this CV threshold or no threshold set.')
}

pdf(file = paste0(argsL$outputname,"_blockaverages.pdf"),paper = 'a4r', width = 11.5, height = 8)
print(plotAllAvgResponses(blockavgs,20))


print('done!')
