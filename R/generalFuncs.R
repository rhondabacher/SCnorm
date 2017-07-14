#' quickreg
#'
#' Perform the single gene regressions using quantile regression.
#' 
#' @details Perform the single gene regressions using quantile regression.
#' 
#' @param x gene to perform the regression on.
#' @param InputData list of data needed for the regression.
#'   
#' @return gene slope.


quickreg<-function(x,InputData)
{
    LogData = InputData[[1]]
    SeqDepth = InputData[[2]]
    X = InputData[[3]][x]
    Tau = InputData[[4]]
    ditherFlag = InputData[[5]]

    if(ditherFlag == TRUE) {
      slope <- try(rq(dither(LogData[X, ], type="symmetric", value=.01) ~ 
                  log(SeqDepth), tau = Tau, method="br")$coef[2], silent=TRUE)
      if(class(slope) == "try-error"){
        slope <- try(rq(dither(LogData[X, ], type="symmetric", value=.01) ~ 
                    log(SeqDepth), tau = Tau, method="fn")$coef[2], silent=TRUE)
        if(class(slope) == "try-error"){
          slope <- NA
        }    
      } 
    } else {
      slope <- try(rq(LogData[X, ] ~ log(SeqDepth), tau = Tau, 
              method="br")$coef[2], silent=TRUE)
      if(class(slope) == "try-error"){
        slope <- try(rq(LogData[X, ] ~ log(SeqDepth), tau = Tau, 
                method="fn")$coef[2], silent=TRUE)
        if(class(slope) == "try-error"){
          slope <- NA
        }    
      }
    }
    names(slope) <- X

    return(slope)
}

#' redobox
#' 
#' @details Function to log data and turn zeros to NA to mask/ignore in 
#'  functions.
#' 
#' @param DATA data set to.
#' @param smallc, what value to ignore, typically is zero.
#'   
#' @return the dataset has been logged with values below smallc masked. 


redobox <- function(DATA, smallc) {

    DATA[DATA <= smallc] <- NA #truncate some values first, usually just zero
    y <- log(DATA)

    return(y)
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
#' data(ExampleData)
#' ExampleData <- SummarizedExperiment::SummarizedExperiment(assays=list("Counts"=ExampleData),
#'                                      colData=data.frame(colnames(ExampleData)))
#' myData <- getCounts(ExampleData)

getCounts <- function(DATA){
  return(SummarizedExperiment::assays(DATA)[["Counts"]])
}



#' @title getResults
#'
#' @usage getResults(DATA, type=c("NormalizedData", "ScaleFactors", "GenesFilteredOut"))
#'
#'   
#' @param DATA An object of class \code{SummarizedExperiment} that contains 
#' normalized single-cell expression and other metadata, and the output of the
#' \code{SCnorm} function.
#' 
#' @param type A character variable specifying which output is desired, 
#'  with possible values "NormalizedData", "ScaleFactors", and "GenesFilteredOut". 
#'  By default getResults() will
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
#' data(ExampleData)
#' Conditions = rep(c(1), each= 90)
#' #runNorm <- SCnorm(Data=ExampleData, Conditions=Conditions)
#' #normData <- getresults(runNorm)

getResults <- function(DATA, type=c("NormalizedData", "ScaleFactors", "GenesFilteredOut")){
  type <- match.arg(type)
  return(S4Vectors::metadata(DATA)[[type]])
}