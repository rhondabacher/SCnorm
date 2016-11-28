#' @title SCnorm

#' @usage SCnorm(Data, Conditions, K, OutputName = "SCnorm", PLOT = TRUE, PropToUse = .25, outlierCheck = TRUE, Tau = .5, reportSF = FALSE)

#' @param Data matrix of un-normalized expression counts. Rows are genes and columns are samples.
#' @param Conditions vector of condition labels, this should correspond to the columns of the un-normalized expression matrix.
#' @param OutputName specify the path and/or name of output files.
#' @param PLOT whether to save all output the evaluation plots while determining the optimal K.
#' @param PropToUse proportion of genes closest to the slope mode used for the group fitting, default is set at .25. This number mainly affects speed. 
#' @param outlierCheck cells/samples with relatively small sequencing depths may recieve very small scaling factors, to ensure these 
#' small scaling factors do not create outliers, genes with potential outliers are first flagged and then if necessary, their normalized expression 
#' values are corrected. First, gene expression counts in cells/samples having the three smallest sequencing depths are flagged if they contain
#' expression counts larger than 3 times the predicted counts (from the group regression). For each gene, if the flagged values
#' are the largest normalized expression value, then the count is corrected to be the median of the non-zero normalized values (default = 
#' TRUE)
#' @param Tau value of quantile for the quantile regression used to estimate gene-specific slopes (default is median, Tau = .5 ). 
#' @param reportSF whether to provide a matrix of scaling counts in the output (default = FALSE).
#' @param FilterCellNumber the number of non-zero expression estimate required to include the genes into the SCnorm fitting (default = 10). THe initial 
#' grouping fits a quantile regression to each gene, making this value too low gives unstable fits.
#' @param K the number of groups for normalizing. If left unspecified, an evaluation procedure will determine the optimal value of K
#' (recommended). If you're sure about specifiyng K, then a vector equal to the number of conditions may be used.
#' @param NCores number of cores to use, default is detectCores() - 1.

#' @description Quantile regression is used to estimate the dependence of read counts on sequencing depth for every gene. Genes
#' with similar dependence are then grouped, and a second quantile regression is used to estimate scale factors within each 
#' group. 
#' Within-group adjustment for sequencing depth is then performed using the estimated scale factors to provide normalized 
#' estimates of expression. If multiple conditions are provided, normalization is performed within condition and then
#' normalized estimates are scaled between conditions.

#' @return List containing matrix of normalized expression (and optionally a matrix of size factors if reportSF = TRUE ).
#' @export

#' @author Rhonda Bacher


SCnorm <- function(Data, Conditions, OutputName, PLOT = T, PropToUse = .25, outlierCheck= T, Tau = .5, 
                   reportSF = F, FilterCellNum = 10, K = NULL, NCores = NULL, FilterExpression = 0) {
  
  Data <- data.matrix(Data)
  if(anyNA(Data)) {stop("Data contains at least one value of NA. Unsure how to proceed.")}
  ## checks
  if (is.null(rownames(Data))) {rownames(Data) <- as.vector(sapply("X_", paste0, 1:dim(Data)[1]))}
  if (is.null(colnames(Data))) {stop("Must supply sample/cell names!")}
  if (dim(Data)[2] != length(Conditions)) {stop("Number of columns in expression matrix must match length of conditions vector!")}
  if (!is.null(K)) {warning(paste0("SCnorm will normalize assuming ", K, " is the optimal number of groups. It is not advised to set this."))}
  if (is.null(NCores)) {NCores <- max(1, detectCores() - 1)}
  
  Levels <- levels(as.factor(Conditions)) # Number of conditions
  
  DataList <- lapply(1:length(Levels), function(x) Data[,which(Conditions == Levels[x])]) # split conditions
  Genes <- rownames(Data) 
  
  #options: correctGC=FALSE, methodGC="loess", 
  # if(correctGC==TRUE) {
  # 	  library(EDASeq)
  # 	  withinLaneNormalization(Data, "gc", which=methodGC)
  #
  
  
  
  SeqDepthList <- lapply(1:length(Levels), function(x) colSums(Data[,which(Conditions == Levels[x])]))
  
  # Get median quantile regr. slopes.
  SlopesList <- lapply(1:length(Levels), function(x) GetSlopes(DataList[[x]], SeqDepthList[[x]], Tau, NCores))
  
  
  NumZerosList <- lapply(1:length(Levels), function(x) { apply(DataList[[x]], 1, function(c) sum(c != 0)) })
  
  
  GeneFilterList <- lapply(1:length(Levels), function(x) names(which(NumZerosList[[x]] >= FilterCellNum)))
  
  GeneFilterOUT <- lapply(1:length(Levels), function(x) names(which(NumZerosList[[x]] < FilterCellNum)))
  
  print("Gene filter is applied within each condition.")
  
  lapply(1:length(Levels), function(x) print(paste0(length(GeneFilterOUT[[x]]), 
         " genes were not included in the normalization due to having less than ", FilterCellNum, 
         " non-zero values.")))
  
  print("A list of these genes can be accessed in output, see vignette for example.") 
  
  # Data, SeqDepth, Slopes, CondNum, PLOT = TRUE, PropToUse, outlierCheck, Tau
  
  if(PLOT==TRUE) {  pdf(paste0(OutputName, "_k_evaluation.pdf"), height=10, width=10)
    par(mfrow=c(2,2))}
  #   if k is NOT provided
  if (is.null(K)) {
    NormList <- lapply(1:length(Levels), function(x) {
      Normalize(Data = DataList[[x]], 
                SeqDepth = SeqDepthList[[x]], Slopes = SlopesList[[x]],
                CondNum = Levels[x], OutputName= OutputName,
                PLOT = PLOT,
                PropToUse = PropToUse,
                outlierCheck = outlierCheck,
                Tau = Tau, NCores= NCores)
    }) 
    
  }
  
  
  # if specific k then do:
  # if length of k is less than number of conditions.
  if (!is.null(K) ) {
    if (length(K) == length(Levels)) {
      NormList <- lapply(1:length(Levels), function(x) {
        SCnorm_fit(Data = DataList[[x]], 
                   SeqDepth = SeqDepthList[[x]], Slopes = SlopesList[[x]],
                   K = K[x], NCores = NCores)
      })
    } else if (length(K) == 1) {
      K <- rep(K, length(Levels))
      NormList <- lapply(1:length(Levels), function(x) {
        SCnorm_fit(Data = DataList[[x]], 
                   SeqDepth = SeqDepthList[[x]], Slopes = SlopesList[[x]],
                   K = K[x], NCores = NCores)
      }) 
    } else (stop("Check that the specification of K is correct!"))
  }	
  
  if (PLOT==TRUE) {  dev.off() }
  
  FilterCellProportion = lapply(1:length(Levels), function(x) FilterCellNum / dim(DataList[[x]])[2])
  
  NORMDATA <- do.call(cbind, lapply(1:length(Levels), function(x) {NormList[[x]]$NormData}))
  
  ## plot the normalized data to screen
  
  checkCountDepth(Data = Data, NormalizedData = NORMDATA,
                  Conditions = Conditions, OutputName = "SCnorm_NormalizedData_FinalK", PLOT = PLOT,
                  FilterCellProportion = FilterCellProportion, FilterExpression = FilterExpression)
  
  
  
  
  if (length(Levels) > 1) {
    
    # Scaling
    # Genes = Reduce(intersect, GeneFilterList)
    
    ScaledNormData <- scaleNormMultCont(NormList, Data, Genes)
    names(ScaledNormData) <- c("NormalizedData", "ScaleFactors")
    ScaledNormData <- c(ScaledNormData, GenesFilteredOut = GeneFilterOUT)
    if(reportSF == T) {
      return(ScaledNormData) 
    } else {
      ScaledNormData$ScaleFactors <- NULL
      return(ScaledNormData) 
    }
  } else {
    NormDataFull <- NormList[[x]]$NormData
    ScaleFactorsFull <- NormList[[x]]$ScaleFactors
    
    if(reportSF == T) {
      FinalNorm <- list(NormalizedData = NormDataFull, ScaleFactors = ScaleFactorsFull, GenesFilteredOut = GeneFilterOUT)
      return(FinalNorm) 
    } else {
      FinalNorm <-list(NormalizedData = NormDataFull, GenesFilteredOut = GeneFilterOUT)
      return(FinalNorm) 
    }
  }
  
  
  
  
  
  
  
}


