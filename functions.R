require(ggplot2)
require(gridExtra)

######################################### PREPROCESSING
avgOverROI <- function(TSmatrix){
  averageTS <- rowMeans(TSmatrix)
  return(averageTS)
}

percBOLDchange <- function(TS_vector){
  #use timecourse average as baseline
  baseline <- mean(TS_vector)
  normalised <- ((TS_vector - baseline)/baseline)*100
  return(normalised)
}

avgOverBlocks <- function(splittedblocks,fun="mean"){
  blockavg_df <- aggregate(list(value=splittedblocks$value), FUN = fun, by = list(timepoint = splittedblocks$timepoint))
  return(blockavg_df)
}

splitBlocks <- function(TS_vector,blockduration,TR){
  #blocksize is (pause duration + stimulus duration)/TR,assuming each dynamic belongs to a block.
  blocksize <- blockduration / TR #size of block in volumes
  
  blocks <- split(TS_vector, ceiling(seq_along(TS_vector)/blocksize))
  return(blocks)
}

removeOutlierBlocks <- function(splittedBlocks, CVcutoff, subjname = "", plotOutliers=TRUE,onduration = argsL$onduration){
  ####does not work on data with mean = 0!!!
  
  #calculate coefficient of variation
  CVs <- aggregate(list(value = splittedBlocks$value), by=list(block=splittedBlocks$block), FUN = function(x){(sd(x)/mean(x))*100})
  CV_all[[subjname]] <<- CVs
  
  print("CV per block:")
  print(CVs)

  CVs$sub <- rep(subjname,length(CVs$block))
  
  #cut and paste
  illegalblocks <- CVs[CVs$value > CVcutoff,]$block
  
  illegaldata <- splittedBlocks[splittedBlocks$block %in% illegalblocks,]
  legaldata <- splittedBlocks[!splittedBlocks$block %in% illegalblocks,]
  
  
  if (plotOutliers == TRUE & length(illegalblocks) > 0){
    print(paste0("The following blocks are classified as outliers using a CV cutoff of ", CVcutoff, ":"))
    print(illegalblocks)
    
    if(is.null(outlierPlots)) {outlierPlots <<- list()} #make global
    plt <- plotResponsePerBlock(illegaldata,subjname = paste0("Outlier blocks (CV > ",CVcutoff,"%) of subj ", subjname),
                                legend=TRUE,blockfacet = TRUE, ytitle = "BOLD signal", onduration = onduration, zeroline=FALSE)
    outlierPlots[[subjname]] <<- plt
  }
  
  return(legaldata)
}

assemble.df <- function(blocks,TR){
  #add timepoint and blocknr and assemble dataframe
  useful <- list()
  for (i in 1:length(blocks)){
    blockTS <- as.numeric(unname(blocks[i])[[1]])
    blockdf <- data.frame(cbind(value = blockTS,volume = 0:(length(blockTS)-1)))
    blockdf$block <- rep(as.character(i),length(blockTS))
    blockdf$timepoint <- blockdf$volume * TR
    useful[[i]] <- blockdf
  }
  df <- do.call(rbind,useful)
  return(df)
}

Denoise <- function(TS_vector,freq=3){
  #double check encoding of period. Essential for sizing of moving averages/loess windows. Governs smoothness.
  TSobj <- ts(TS_vector, deltat = 1/freq)
  
  
  #two methods: pick one
  #decompose by LOESS
  #stled <- stl(TSobj, "periodic")
  #denoised <- as.numeric(stled$time.series[,2])
  
  #decompose using the Moving Averages model (http://www.itl.nist.gov/div898/handbook/pmc/section4/pmc422.htm)
  dcmpsd <- decompose(TSobj,"additive")
  denoised <- dcmpsd$trend
  
  return(denoised)
}

######################################### PLOTTING
plotRawTS <- function(timeseries,TR,subjname,onduration,blockduration){
  maxVol <- (length(timeseries)-1)
  maxTime <- maxVol*TR
  shadestart <- seq(from = 0, to = maxTime, by = blockduration)
  shade_end <- sapply(shadestart,function(x){min(x+onduration,maxTime)})
  shading <- data.frame(xstart = shadestart, xend = shade_end)
  df <- data.frame(cbind(value = timeseries,volume = 0:maxVol))
  df$timepoint <- df$volume * TR
  
  tsplot <- ggplot() + geom_line(data=df,aes(timepoint,value),size = 1) +
    #    geom_hline(yintercept = 0, color = "black") + 
    labs(title = subjname, y = 'BOLD signal', x = "Time(s)") +
    geom_rect(data = shading, aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf), alpha = 0.4) +
    theme(legend.position="none")
  return(tsplot)
}

plotResponsePerBlock <- function(splittedblocks,avgTS=NULL,onduration=NULL, zeroline = TRUE,
                                 subjname="",ytitle = "% BOLD Change", bestfit=NULL, legend=FALSE, blockfacet = FALSE,annotation = ''){
  #splittedblocks,bestfit, and avgTS must be data frames with at least timepoint(in s) and the boldperc value as variables.
  
  ####plot all blocks for one subject
  df <- splittedblocks
  
  tsplot <- ggplot() +
    geom_line(data=df,aes(timepoint,value, color = block),size = 0.6) +
    labs(title = subjname, y = ytitle, x = "Time(s)")
  
  if (zeroline == TRUE){
    tsplot <- tsplot + geom_hline(yintercept = 0, color = "black") 
  }
  
  if (blockfacet == TRUE){
    tsplot <- tsplot + facet_wrap(~ block, scales = "free_y")
  }
  
  if (!is.null(bestfit)){
    tsplot <- tsplot + geom_line(data=avgTS,aes(timepoint,value, color = "average"), colour="#000099", size = 2)
  }
  
  if (!is.null(onduration)){
    offduration <- (max(df$timepoint) - onduration)
    shading <- data.frame(xstart = 0, xend = onduration)
    tsplot <- tsplot + geom_rect(data = shading, aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf), alpha = 0.4)
  }
  
  if (!is.null(bestfit)){
    tsplot <- tsplot + geom_line(data=bestfit,aes(timepoint,value), colour="#f90505", size = 1)
  }
  
  if (legend == FALSE){
    tsplot <- tsplot + theme(legend.position="none")
  }
  
  #annotate with params
  xMax <- max(df$timepoint)
  yMax <- max(df$value)
  tsplot <- tsplot + annotate("text", x = c(xMax), y=c(yMax), label = c(annotation), hjust=1,size = 1.8)
  
  return(tsplot)
}

plotAllAvgResponses <- function(blockAverages,onduration){
  #plot the average timeseries for all subjects
  tmp <- list()
  for (i in 1:length(blockAverages)){
    subjname <- names(blockAverages[i])
    df <- blockAverages[i][[1]]
    df$subj <- rep(as.character(subjname),length(df$value))
    tmp[[subjname]] <- df
  }
  
  folded <- do.call(rbind,tmp)
  folded$subj <- as.factor(folded$subj)
  
  offduration <- (max(df$timepoint) - onduration)
  shading <- data.frame(xstart = 0, xend = onduration)

  avgPlot <- ggplot() +
    geom_line(data=folded,aes(timepoint,value, color = subj),size = 1) +
    #geom_line(data=df_avgTS,aes(timepoint,x, color = "average"), colour="#000099", size = 2) +
    geom_hline(yintercept = 0, color = "black") + 
    geom_rect(data = shading, aes(xmin = xstart, xmax = xend, ymin = -Inf, ymax = Inf), alpha = 0.4) +
    labs(title = "", y = "% BOLD Change", x = "Time(s)")
  
  return(avgPlot)
}

######################################## Curve Fitting
trapFunc <- function(x,a,b,c,d,ceiling,baseline){
  if (0 <= x & x < a){
    y <- baseline
    return(y)
  }
  if(a <= x & x < b){
    y <- baseline + (1-((b-x)/(b-a)))*(ceiling - baseline)
    return(y)
  }
  if (b <= x & x < c){
    y <- ceiling
    return(y)
  }
  if (c <= x & x <= d){
    y <- baseline + ((d-x)/(d-c))*(ceiling-baseline)
    return(y)
  }
  if (d < x){
    y <- baseline
    return(y)
  }
}


evalFunc <- function(par,xvals,yvals){
  #par is supplied by optimisation routine
  #xvals, yvals have to be provided to calculate RSS
  a <- par['a']
  b <- par['b']
  c <- par['c']
  d <- par['d']
  ceiling <- par['ceiling']
  baseline <- par['baseline']
  
  #fit
  fittedvals <- sapply(xvals,trapFunc,a,b,c,d,ceiling,baseline)
  
  #calculate RSS (note different lengths of yvals and fittedvals)
  RSS <- sum((yvals-fittedvals)^2)
 
  return(RSS)
}

grad <- function(f, x0, step = .Machine$double.eps^(1/3), ...) {
  #Computes the gradient numerically using the central difference formula.
  if (!is.numeric(x0))
    stop("Argument 'x0' must be a numeric value.")
  
  fun <- match.fun(f)
  f <- function(x) fun(x, ...)
  
  if (length(f(x0)) != 1)
    stop("Function 'f' must be a univariate function of >2 variables.")
  n <- length(x0)
  
  hh <- rep(0, n)
  gr <- numeric(n)
  for (i in 1:n) {
    hh[i] <- step
    gr[i] <- (f(x0 + hh) - f(x0 - hh)) / (2*step)
    hh[i] <- 0
  }
  return(gr)
}

gradientFunc <- function(par,xvals,yvals){
  # xvals, yvals are ignored but present for coding purposes (required by the constrOptim API)
  return(grad(evalFunc,par,xvals = xvals,yvals = yvals))
  
}

trapezoidFit <- function(data_df, startingparams){
  #fit the trapezoid function
  #volumes is the total number of volumes in a block to constrain (a+b+c+d < Tmax)
  Tmax <- max(data_df$timepoint)
  
  #formulate constraints so that constrMatrix %*% params - constrVector >= 0
  #completely depends on the order of the params and constraint vectors!
  constrMatrix <- rbind(c(1,0,0,0,0,0),c(-1,1,0,0,0,0),c(0,-1,1,0,0,0),c(0,0,-1,1,0,0),c(0,0,0,-1,0,0))
  constrVector <- c(0,0.1,0.1,0.1,-Tmax) #0.1 should be zeroes, but forcing distance between breakpoints works better
  
  # print(Tmax)
  # print(constrMatrix %*% startingparams - constrVector)


  ##no gradient (Nelder-Mead)
  # result <- constrOptim(theta = startingparams,f = evalFunc, grad = NULL,outer.iterations = 500, ui = constrMatrix,ci = constrVector,
  #                       xvals = data_df$timepoint, yvals = data_df$value)
  
  ##with gradient (BFGS)
  # result <- constrOptim(theta = startingparams, f = evalFunc, grad = gradientFunc,outer.iterations = 500, ui = constrMatrix,ci = constrVector,
  #                      xvals = data_df$timepoint, yvals = data_df$value)
  
  #with gradient (conjugate-gradient method)
   result <- constrOptim(theta = startingparams, f = evalFunc, grad = gradientFunc,outer.iterations = 500, method = "CG", ui = constrMatrix,ci = constrVector,
                         xvals = data_df$timepoint, yvals = data_df$value)

  return(result)
}


generateValuesForFit <- function(par,xbound,TR,xresolution){
  #generate values from the function and a set of parameters for use in plotting
  a <- par['a']
  b <- par['b']
  c <- par['c']
  d <- par['d']
  ceiling <- par['ceiling']
  baseline <- par['baseline']
  xvals = seq(from = 0, to = xbound, by = xresolution)
  fittedVals <- sapply(xvals,trapFunc,a,b,c,d,ceiling,baseline)
  bestfit_df <- data.frame(value = fittedVals, volume = xvals/TR ,timepoint = xvals)
  return(bestfit_df)
}
