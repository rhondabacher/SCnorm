#' quickReg
#'
#' Perform the single gene regressions using quantile regression.
#' 
#' @details Perform the single gene regressions using quantile regression.
#' 
#' @param x gene to perform the regression on.
#' @param InputData list of data needed for the regression.
#'   
#' @return gene slope.


quickReg <- function(x,InputData) {
  
    LogData = InputData[[1]]
    SeqDepth = InputData[[2]]
    X = InputData[[3]][x]
    Tau = InputData[[4]]
    ditherFlag = InputData[[5]]

    if(ditherFlag == TRUE) {
      slope <- try(quantreg::rq(quantreg::dither(LogData[X, ], type="symmetric", value=.01) ~ 
                  log(SeqDepth), tau = Tau, method="br")$coef[2], silent=TRUE)
      if(methods::is(slope, "try-error")){
        slope <- try(quantreg::rq(quantreg::dither(LogData[X, ], type="symmetric", value=.01) ~ 
                    log(SeqDepth), tau = Tau, method="fn")$coef[2], silent=TRUE)
        if(methods::is(slope, "try-error")){
          slope <- NA
        }    
      } 
    } else {
      slope <- try(quantreg::rq(LogData[X, ] ~ log(SeqDepth), tau = Tau, 
              method="br")$coef[2], silent=TRUE)
      if(methods::is(slope, "try-error")){
        slope <- try(quantreg::rq(LogData[X, ] ~ log(SeqDepth), tau = Tau, 
                method="fn")$coef[2], silent=TRUE)
        if(methods::is(slope, "try-error")){
          slope <- NA
        }    
      }
    }
    names(slope) <- X

    return(slope)
}

#' redoBox
#' 
#' @details Function to log data and turn zeros to NA to mask/ignore in 
#'  functions.
#' 
#' @param DATA data set to.
#' @param smallc, what value to ignore, typically is zero.
#'   
#' @return the dataset has been logged with values below smallc masked. 


redoBox <- function(DATA, smallc) {

    DATA[DATA <= smallc] <- NA #truncate some values first, usually just zero
    y <- log(DATA)

    return(y)
}

#' splitGroups
#' 
#' @details helper function to get split a vector into a specified number of groups
#' 
#' @param DATA vector to be splot.
#' @param NumGroups number of groups
#'   
#' @return list, length is equal to NumGroups

splitGroups <- function(DATA, NumGroups=10) {
  splitby <- sort(DATA) 
  grps <- length(splitby) / NumGroups
  sreg <- split(splitby, ceiling(seq_along(splitby) / grps))
  return(sreg)
}

#' getDens
#' 
#' @details get density of slopes in different expression groups
#' 
#' @param ExprGroups expression groups already split.
#' @param byGroup factor (usually slopes) to get density based on ExprGroups.
#' @param RETURN whether to return Mode or Height of density.
#' 
#' @return list, length is equal to NumGroups

getDens <- function(ExprGroups, byGroup, RETURN=c("Mode", "Height")) {
 NumExpressionGroups <- length(ExprGroups)
 Mode <- rep(NA, NumExpressionGroups)
 Height <- rep(NA, NumExpressionGroups)
 for (i in seq_len(NumExpressionGroups)) {
   useg <- names(ExprGroups[[i]])
   rqdens <- density(na.omit(byGroup[useg]))
   peak <- which.max(rqdens$y)
   Mode[i] <- rqdens$x[peak]
   Height[i] <- rqdens$y[peak]
  }
  RETURN <- eval(parse(text = match.arg(RETURN)))
  return(RETURN)
}

#' @title getCounts
#' 
#' @usage getCounts(DATA)
#'
#' @param DATA An object of class \code{SummarizedExperiment} that contains 
#' single-cell expression and metadata
#' 
#' @description Convenient helper function to extract the normalized expression
#' matrix from the SummarizedExperiment
#' 
#' @return A \code{matrix} which contains the count data
#'  where genes are in rows and cells are in columns
#'
#' @export
#'
#' @importFrom SummarizedExperiment assays
#' @examples 
#' data(ExampleSimSCData)
#' ExampleData <- SummarizedExperiment::SummarizedExperiment(assays=list("Counts"=ExampleSimSCData))
#' myData <- getCounts(ExampleData)

getCounts <- function(DATA){
  return(SummarizedExperiment::assays(DATA)[["Counts"]])
}



#' @title results
#'   
#' @param DATA An object of class \code{SummarizedExperiment} that contains 
#' normalized single-cell expression and other metadata, and the output of the
#' \code{SCnorm} function.
#' 
#' @param type A character variable specifying which output is desired, 
#'  with possible values "NormalizedData", "ScaleFactors", and "GenesFilteredOut". 
#'  By default results() will
#'  return type="NormalizedData", which is the matrix of normalized counts from SCnorm.
#'  By specifiying type="ScaleFactors" a matrix of scale factors (only returned if 
#'  reportSF=TRUE when running \code{SCnorm()}) can be obtained.
#'  type="GenesFilteredOut" returns a list of genes that were not normalized using 
#'  SCnorm, these are genes that
#'  did not pass the filter critiera.
#'    
#' @description Convenient helper function to extract the results (
#' normalized data, list of genes filtered out, or scale factors). Results 
#' data.frames/matrices are stored in the 
#' \code{metadata} slot and can also be accessed without the help of this 
#' convenience function by calling \code{metadata(DataNorm)}.
#'
#' @return A \code{data.frame} containing output as detailed in the
#'  description of the \code{type} input parameter
#'
#' @export
#'
#' @importFrom S4Vectors metadata
#' @examples 
#' data(ExampleSimSCData)
#' Conditions = rep(c(1), each= 90)
#' #NormData <- SCnorm(Data=ExampleSimSCData, Conditions=Conditions)
#' #normDataMatrix <- results(NormData)

results <- function(DATA, type=c("NormalizedData", "ScaleFactors", "GenesFilteredOut")){
  type <- match.arg(type)
  return(S4Vectors::metadata(DATA)[[type]])
}




#' correctWithin
#'
#' Perform the correction within each sample (See loess normalization in 
#' original publication Risso et al., 2011 (BMC Bioinformatics)).
#' Similar to function in EDAseq v2.8.0.
#'  
#' @details Performs within sample normalization.
#' 
#' @param y gene to perform the regression on.
#' @param correctFactor list of data needed for the regression.
#'   
#' @return within-cell normalized expression estimates

correctWithin <- function(y, correctFactor) {
  
  #don't use zeros or outliers
  useg <- which(y > 0 & y <= quantile(y, probs=0.995)) 
  X <- correctFactor[useg]
  Y <- log(y[useg])

  calcL <- loess(Y ~ X)
  counts.fit <- predict(calcL, newdata = correctFactor)
  names(counts.fit) <- names(y)
  counts.fit[is.na(counts.fit)] <- 0

  scaleC <- y / exp(counts.fit - median(Y)) #correct
  return(scaleC)
} ##from EDAseq v2.8.0