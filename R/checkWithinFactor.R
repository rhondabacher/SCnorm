#' @title Evaluation within sample factors in the the data.

#' @usage checkWithinFactor(Data, withinSample=NULL, Conditions = NULL, 
#'     FilterExpression = 0, 
#'     NumExpressionGroups = 4) 

#' @param Data matrix of un-normalized expression counts. Rows are genes
#'  and columns are samples.
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
#'   for such feature may be considered especially if the bias if different 
#'   between conditions (see SCnorm vignette for details).

#' @return produces a plot and returns the data the plot is based on.
#' @export

#' @author Rhonda Bacher
#' @import graphics
#' @import ggplot2
#' @import grDevices
#' @importFrom SummarizedExperiment SummarizedExperiment assayNames assays colData
#' @importFrom S4Vectors metadata
#' @importFrom data.table melt
#'  
#' @examples 
#'  
#' data(ExampleData)
#' Conditions = rep(c(1,2), each= 90) 
#' exampleFactor = runif(dim(ExampleData)[1], 0, 1)
#' #checkWithinFactor(Data = ExampleData,  
#'   #withinSample=exampleFactor, Conditions = Conditions)

checkWithinFactor <- function(Data, withinSample=NULL, Conditions = NULL, 
    FilterExpression = 0, 
    NumExpressionGroups = 4) {
      

  if("SummarizedExperiment" %in% class(Data)) {
    if(is.null(SummarizedExperiment::assayNames(Data)) || SummarizedExperiment::assayNames(Data)[1] != "Counts") {
      message("renaming the first element in assays(Data) to 'Counts'")
      SummarizedExperiment::assayNames(Data)[1] <- "Counts"
    }

    Data = SCnorm::getCounts(Data)
  }
      
    Data <- data.matrix(Data)
    if(anyNA(Data)) {stop("Data contains at least one value of NA. 
      Unsure how to proceed.")}
    if(is.null(withinSample)) {stop("Please provide a vector of gene-specific features!")}
    if(is.null(names(withinSample))) {stop("Make sure feature names match the rownames of input data matrix")}
      
    if(is.null(rownames(Data))) {rownames(Data) <- as.vector(sapply("X_", 
        paste0, 1:dim(Data)[1]))}
    if(is.null(colnames(Data))) {stop("Must supply sample/cell names!")}
    if(is.null(Conditions)) {Conditions <- rep("1", dim(Data)[2])}

    if(dim(Data)[2] != length(Conditions)) {stop("Number of columns in 
         expression matrix must match length of conditions vector!")}

    
    SeqDepthList <- colSums(Data)
    MedExpr <- apply(Data, 1, function(c) median(c[c != 0]))
    GeneFilter <- names(which(MedExpr >= FilterExpression))
    
    if(length(GeneFilter) == 0) {stop("No genes pass the filter specified! 
         Try lowering thresholds or perform more QC on your data.")}

    splitby <- sort(withinSample) 
    grps <- length(splitby) / NumExpressionGroups
    withinSplit <- split(splitby, ceiling(seq_along(splitby) / grps))
  
    SplitMedExprList <- lapply(1:dim(Data)[2], function(x) {
        sapply(1:NumExpressionGroups, function(y) {
            useg <- intersect(names(withinSplit[[y]]), GeneFilter)
            median(exp(redobox(Data[useg,x],0)), na.rm=TRUE)
        })
      })
        


    withinCells <- data.frame(Sample = colnames(Data), Condition = factor(Conditions),
     do.call(rbind, SplitMedExprList))
    rownames(withinCells) <- colnames(Data)
    colnames(withinCells) <- c("Sample", "Condition", 
        paste0("Group",1:NumExpressionGroups))
    longdata <- data.table::melt(withinCells, id=c("Sample", "Condition"))

if(length(unique(Conditions))==1) { 
  
    QQ <- ggplot(data=longdata, aes_(x=forcats::fct_inorder(factor(longdata$variable)),y=longdata$value,
        colour=longdata$Sample, group=longdata$Sample))+
         geom_point(show.legend = FALSE) + geom_line(size=.2,show.legend = FALSE)
    	
      } else {
    QQ <- ggplot(data=longdata, aes_(x=forcats::fct_inorder(factor(longdata$variable)),y=longdata$value,
        colour=longdata$Condition, group=longdata$Sample))+
         geom_point(show.legend = FALSE) + geom_line(size=.2,show.legend = TRUE) 
    
      }
    plot(QQ +  labs(x="GC content group", y="Log non-zero median expression") + 
    guides(colour=guide_legend(title="Conditions"))+
    	 theme_bw())
    
    
    return(withinCells)
}
