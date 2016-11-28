#' @title Evaluation the count-depth relationship before (or after) normalizing the data.

#' @usage checkCountDepth(Data, Conditions, OutputName, PLOT = T, Tau = .5, FilterCellProportion = .10, NumExpressionGroups = 10, NCores)

#' @param Data matrix of un-normalized expression counts. Rows are genes and columns are samples.
#' @param NormalizedData matrix of normalized expression counts. Rows are genes and columns are samples. Only input this if evaluating already normalized data.

#' @param Conditions vector of condition labels, this should correspond to the columns of the un-normalized expression matrix. If not provided data is assumed to come from same condition/batch.
#' @param OutputName specify the path and/or name of output files.
#' @param PLOT whether to save the evaluation plots in addition to printing to screen.
#' @param Tau value of quantile for the quantile regression used to estimate gene-specific slopes (default is median, Tau = .5 ). 
#' @param FilterCellProportion the proportion of non-zero expression estimates required to include the genes into the evaluation. Default is .10. 
#' @param NumExpressionGroups the number of groups to split the data into, groups are split into equally sized groups based on non-zero median expression. 
#' @param NCores number of cores to use, default is detectCores() - 1.

#' @description Quantile regression is used to estimate the dependence of read counts on sequencing depth for every gene. If multiple conditions are provided, a separate plot is provided for each. Can be used to evaluate the extent of the count-depth relationship in the dataset or can be be used to evaluate data normalized by alternative methods.

#' @return outputs a plot.
#' @export

#' @author Rhonda Bacher


checkCountDepth <- function(Data, NormalizedData= NULL, Conditions = NULL, OutputName, PLOT = T, Tau = .5, FilterCellProportion = .10, 
	FilterExpression = 0, NumExpressionGroups = 10, NCores=NULL) {
  
  Data <- data.matrix(Data)
  ## checks
  if (is.null(rownames(Data))) {rownames(Data) <- as.vector(sapply("X_", paste0, 1:dim(Data)[1]))}
  if (is.null(colnames(Data))) {stop("Must supply sample/cell names!")}
  if(is.null(Conditions)) {Conditions <- rep("1", dim(Data)[2])}
  if (dim(Data)[2] != length(Conditions)) {stop("Number of columns in expression matrix must match length of conditions vector!")}
  if (is.null(NCores)) {NCores <- max(1, detectCores() - 1)}
	  

  Levels <- levels(as.factor(Conditions)) # Number of conditions
  
  DataList <- lapply(1:length(Levels), function(x) Data[,which(Conditions == Levels[x])]) # split conditions

  FilterCellProportion <-  lapply(1:length(Levels), function(x) max(FilterCellProportion, 10 / dim(DataList[[x]])[2])) 
  
  SeqDepthList <- lapply(1:length(Levels), function(x) colSums(Data[,which(Conditions == Levels[x])]))

  PropZerosList <- lapply(1:length(Levels), function(x) { apply(DataList[[x]], 1, function(c) sum(c != 0)) / length(SeqDepthList[[x]]) })
  
  MedExprAll <- apply(Data, 1, function(c) median(c[c != 0]))
	
  MedExprList <- lapply(1:length(Levels), function(x) { apply(DataList[[x]], 1, function(c) median(c[c != 0])) })
  
  BeforeNorm <- TRUE
  #switch to the normalized data:
  if(!is.null(NormalizedData)) {  
	  DataList <- lapply(1:length(Levels), function(x) NormalizedData[,which(Conditions == Levels[x])])
	  BeforeNorm <- FALSE
  }
 
  GeneFilterList <- lapply(1:length(Levels), function(x) names(which(PropZerosList[[x]] >= FilterCellProportion[[x]] & MedExprAll >= FilterExpression)))
 
  
  # Get median quantile regr. slopes.
  SlopesList <- lapply(1:length(Levels), function(x) GetSlopes(DataList[[x]][GeneFilterList[[x]],], SeqDepthList[[x]], Tau, NCores))
  
  
  # Data, SeqDepth, Slopes, CondNum, PLOT = TRUE, PropToUse, outlierCheck, Tau
  ROWS <- round(length(Levels) / 2)
  par(mfrow=c(ROWS,2))
  lapply(1:length(Levels), function(x) {
 	initialEvalPlot(MedExpr = MedExprList[[x]][GeneFilterList[[x]]], SeqDepth = SeqDepthList[[x]], 
                                                        Slopes = SlopesList[[x]], Name = Levels[[x]], NumExpressionGroups, BeforeNorm = BeforeNorm)
	})
  
   if(PLOT==TRUE) {  
	 dev.copy(pdf, file=paste0(OutputName, "_count-depth_evaluation.pdf"), height=5*ROWS, width=10)
     dev.off()
     }

  
}
