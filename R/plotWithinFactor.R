#' @title Evaluate gene-specific factors in the the data.

#' @param Data can be a matrix of single-cell expression with cells  
#'   where rows are genes and columns are samples. Gene names should
#'   not be a column in this matrix, but should be assigned to rownames(Data).
#'   Data can also be an object of class \code{SummarizedExperiment} that contains 
#'   the single-cell expression matrix and other metadata. The \code{assays} 
#'   slot contains the expression matrix and is named \code{"Counts"}.  
#'   This matrix should have one row for each gene and one sample for each column.  
#'   The \code{colData} slot should contain a data.frame with one row per 
#'   sample and columns that contain metadata for each sample.  This data.frame
#'   should contain a variable that represents biological condition 
#'   in the same order as the columns of \code{NormCounts}). 
#'   Additional information about the experiment can be contained in the
#'   \code{metadata} slot as a list.
#' @param withinSample a vector of gene-specific features. 
#' @param Conditions vector of condition labels, this should correspond to 
#'  the columns of the un-normalized expression matrix. If provided the cells
#' will be colored by Condition instead of individually.
#' @param FilterExpression exclude genes having median of non-zero expression
#'  below this threshold.
#' @param NumExpressionGroups the number of groups to split the within 
#'   sample factor into, e.g genes will be split into equally sized groups
#'   based on their GC content/Gene length/etc.

#' @description This function can be used to evaluate the extent of 
#'   gene-specific biases in the data. If a bias exists, the plots provided 
#'   here will identify whether it affects cells equally or not. Correction 
#'   for such features may be considered especially if the bias is different 
#'   between conditions (see SCnorm vignette for details).

#' @return produces a plot and returns the data the plot is based on.
#' @export

#' @author Rhonda Bacher
#' @import graphics
#' @import ggplot2
#' @import grDevices
#' @importFrom methods is as
#' @importFrom SummarizedExperiment SummarizedExperiment assayNames assays colData
#' @importFrom S4Vectors metadata
#' @importFrom data.table melt
#'  
#' @examples 
#'  
#' data(ExampleSimSCData)
#' Conditions = rep(c(1,2), each= 90) 
#' exampleFactor = runif(dim(ExampleSimSCData)[1], 0, 1)
#' names(exampleFactor) = rownames(ExampleSimSCData)
#' #plotWithinFactor(Data = ExampleSimSCData,  
#'   #withinSample=exampleFactor, Conditions = Conditions)

plotWithinFactor <- function(Data, withinSample=NULL, Conditions = NULL, 
                              FilterExpression = 0, 
                              NumExpressionGroups = 4) {
      

  #Checks
  if (methods::is(Data, "SummarizedExperiment") | methods::is(Data, "SingleCellExperiment")) {
      Data <- methods::as(Data, "SingleCellExperiment")
    if (is.null(SummarizedExperiment::assayNames(Data)) || SummarizedExperiment::assayNames(Data)[1] != "counts") {
      message("Renaming the first element in assays(Data) to 'counts'")
        SummarizedExperiment::assayNames(Data)[1] <- "counts"

    if (is.null(colnames(SingleCellExperiment::counts(Data)))) {stop("Must supply sample/cell names!")}


    }
    Data <- as.matrix(SingleCellExperiment::counts(Data))
  }
      
    Data <- data.matrix(Data)
    if(anyNA(Data)) {stop("Data contains at least one value of NA. 
                           Unsure how to proceed.")}
    if(is.null(withinSample)) {stop("Please provide a vector of gene-specific features!")}
    if(is.null(names(withinSample))) {stop("Make sure feature names match the 
                                            row names of input data matrix")}
      
    if (is.null(rownames(Data))) {stop("Must supply gene/row names!")}
    if (is.null(colnames(Data))) {stop("Must supply sample/cell names!")}
   
    if(is.null(Conditions)) {Conditions <- rep("1", dim(Data)[2])}

    if(ncol(Data) != length(Conditions)) {stop("Number of columns in 
         expression matrix must match length of conditions vector!")}

    
    SeqDepthList <- colSums(Data)
    MedExpr <- apply(Data, 1, function(x) median(x[x != 0]))
    GeneFilter <- names(which(MedExpr >= FilterExpression))
    
    if (length(GeneFilter) < 0) {stop(" Less than 100 genes pass the filter specified! 
         Try lowering thresholds or perform more QC on your data.")}
    
    withinSplit <- splitGroups(withinSample, NumExpressionGroups)
    
    SplitMedExprList <- lapply(seq_len(ncol(Data)), function(x) {
        vapply(seq_len(NumExpressionGroups), function(y) {
            useg <- intersect(names(withinSplit[[y]]), GeneFilter)
            Y <- Data[useg,x]
            median(Y[Y!=0])
        }, FUN.VALUE=numeric(1))
      })
        


    withinCells <- data.frame(Sample = colnames(Data), Condition = factor(Conditions),
                             do.call(rbind, SplitMedExprList))
    rownames(withinCells) <- colnames(Data)
    colnames(withinCells) <- c("Sample", "Condition", 
                             paste0("Group",1:NumExpressionGroups))
    longdata <- data.table::melt(withinCells, id=c("Sample", "Condition"))

    if (length(unique(Conditions))==1) { 
  
       QQ <- ggplot2::ggplot(data=longdata, aes_(x=forcats::fct_inorder(factor(longdata$variable)),
                                                 y=longdata$value,
                                                 colour=longdata$Sample, 
                                                 group=longdata$Sample))+
                    geom_line(size=.2,show.legend = FALSE)
    	
     } else {
      QQ <- ggplot2::ggplot(data=longdata, aes_(x=forcats::fct_inorder(factor(longdata$variable)),
                                                y=longdata$value,
                                                colour=longdata$Condition, 
                                                group=longdata$Sample))+
                    geom_line(size=.2,show.legend = TRUE)
                     
    
      }
   plot(QQ  +  labs(x="GC content group", y="Log non-zero median expression") + 
               guides(colour=guide_legend(title="Conditions"))+
               geom_point(show.legend = FALSE) + 
               theme_bw())
  
    
  return(withinCells)
}
