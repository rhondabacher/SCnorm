#' @title Evaluate the count-depth relationship before (or after) normalizing
#' the data.
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
#' @param NormalizedData matrix of normalized expression counts. Rows are
#'  genesand columns are samples. Only input this if evaluating already
#'  normalized data.
#' @param Conditions vector of condition labels, this should correspond to 
#'  the columns of the un-normalized expression matrix. If not provided data
#'  is assumed to come from same condition/batch.
#' @param Tau value of quantile for the quantile regression used to
#'  estimate gene-specific slopes (default is Tau = .5 (median)). 
#' @param FilterCellProportion the proportion of non-zero expression estimates
#'  required to include the genes into the evaluation. Default is .10, and 
#'  will not go below a proportion which uses less than 10 total cells/samples.
#' @param FilterExpression exclude genes having median of non-zero expression
#'  below this threshold from count-depth plots (default = 0).
#' @param NumExpressionGroups the number of groups to split the data into,
#'  genes are split into equally sized groups based on their non-zero median
#'  expression. 
#' @param NCores number of cores to use, default is detectCores() - 1. 
#' This will be used to set up a parallel environment using either MulticoreParam (Linux, Mac) 
#' or SnowParam (Windows) with NCores using the package BiocParallel. 
#' @param ditherCounts whether to dither/jitter the counts, may be used for
#'  data with many ties, default is FALSE. 

#' @description Quantile regression is used to estimate the dependence of read
#'  counts on sequencing depth for every gene. If multiple conditions are
#'  provided, a separate plot is provided for each and the filters are 
#'  applied within each condition separately. The plot can be used to evaluate
#'  the extent of the count-depth relationship in the dataset or can be be
#'  used to evaluate data normalized by alternative methods.

#' @return returns a data.frame containing each gene's slope (count-depth relationship)
#' and it's associated expression group. A plot will be output if showPlot=TRUE.
#' @export

#' @author Rhonda Bacher
#' @importFrom parallel detectCores
#' @import graphics
#' @import grDevices
#' @importFrom methods is
#' @importFrom BiocParallel bplapply  
#' @importFrom BiocParallel register
#' @importFrom BiocParallel MulticoreParam
#' @importFrom BiocParallel bpparam
#' @importFrom parallel detectCores
#' @importFrom S4Vectors metadata
#' @importFrom SummarizedExperiment SummarizedExperiment assayNames assays colData
#' @examples 
#'  
#' data(ExampleSimSCData)
#' Conditions = rep(c(1,2), each= 90) 
#' #plotCountDepth(Data = ExampleSimSCData, Conditions = Conditions, 
#'   #FilterCellProportion = .1)

plotCountDepth <- function(Data, NormalizedData= NULL, Conditions = NULL, 
                           Tau = .5, FilterCellProportion = .10, 
                           FilterExpression = 0, NumExpressionGroups = 10, NCores=NULL, 
                           ditherCounts = FALSE) {
    
    # Checks  
    if (methods::is(Data, "SummarizedExperiment")) {
      if (is.null( SummarizedExperiment::assayNames(Data)) ||  SummarizedExperiment::assayNames(Data)[1] != "Counts") {
        message("renaming the first element in assays(Data) to 'Counts'")
        SummarizedExperiment::assayNames(Data)[1] <- "Counts"
    }

    Data = SCnorm::getCounts(Data)
    }
    
    Data <- data.matrix(Data)
    if (anyNA(Data)) {stop("Data contains at least one value of NA. 
                            Unsure how to proceed.")}

    if (is.null(NCores)) {NCores <- max(1, parallel::detectCores() - 1)}
    if (.Platform$OS.type == "windows") {
      prll=BiocParallel::SnowParam(workers=NCores)
      BiocParallel::register(BPPARAM = prll, default=TRUE)
    } else {   
      prll=BiocParallel::MulticoreParam(workers=NCores)
      BiocParallel::register(BPPARAM = prll, default=TRUE)
    }
    
    if (is.null(rownames(Data))) {stop("Must supply gene/row names!")}
    if (is.null(colnames(Data))) {stop("Must supply sample/cell names!")}

    if (is.null(Conditions)) {Conditions <- rep("1", ncol(Data))}
    if (ncol(Data) != length(Conditions)) {stop("Number of columns in 
         expression matrix must match length of conditions vector!")} 
    if (ditherCounts == TRUE) {RNGkind("L'Ecuyer-CMRG");set.seed(1);
         message("Jittering values introduces some randomness, for 
           reproducibility set.seed(1) has been set")}
  
    Levels <- unique(Conditions)
    DataList <- lapply(seq_along(Levels), function(x) {
                      Data[,which(Conditions == Levels[x])]}) # split conditions

    # Restrict this to not go below less than 10 cells/samples
    FilterCellProportion <-  lapply(seq_along(Levels), function(x) {
        max(FilterCellProportion, 10 / dim(DataList[[x]])[2])})

    SeqDepthList <- lapply(seq_along(Levels), function(x) {
        colSums(Data[,which(Conditions == Levels[x])])})

    PropZerosList <- lapply(seq_along(Levels), function(x) {
         rowSums(DataList[[x]] != 0) / ncol(DataList[[x]])})
     
    MedExprList <- lapply(seq_along(Levels), function(x) {
        apply(DataList[[x]], 1, function(c) median(c[c != 0])) })

    BeforeNorm <- TRUE
    # Switch to the normalized data:
    if(!is.null(NormalizedData)) {  
       DataList <- lapply(seq_along(Levels), function(x) {
          NormalizedData[,which(Conditions == Levels[x])]})
       BeforeNorm <- FALSE
    }

    GeneFilterList <- lapply(seq_along(Levels), function(x) {
        names(which(PropZerosList[[x]] >= FilterCellProportion[[x]] & 
          MedExprList[[x]] >= FilterExpression))})
    NM <- unlist(lapply(seq_along(Levels), function(x) {
            length(GeneFilterList[[x]] )}))
   if(any(NM < 100)) {stop("Less than 100 genes pass the filter specified! 
               Try lowering thresholds or perform more QC on your data.")}

    # Get median quantile regr. slopes.
    SlopesList <- lapply(seq_along(Levels), function(x) {
          getSlopes(Data = DataList[[x]][GeneFilterList[[x]],], 
                    SeqDepth = SeqDepthList[[x]], 
                    Tau = Tau, FilterCellNum = 10,
                    ditherCounts = ditherCounts)})

    plotSlots <- lapply(seq_along(Levels), function(x) {
            generateEvalPlot(MedExpr = MedExprList[[x]][GeneFilterList[[x]]], 
                             SeqDepth = SeqDepthList[[x]], 
                             Slopes = SlopesList[[x]], 
                             Name = Levels[[x]], 
                             NumExpressionGroups = NumExpressionGroups, 
                             BeforeNorm = BeforeNorm)})
    names(plotSlots) <- Levels
    
    return(plotSlots)

}
