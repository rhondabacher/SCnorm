#' @title SCnorm
#'
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
#' @param Conditions vector of condition labels, this should correspond to
#'    the columns of the expression matrix.
#' @param PrintProgressPlots whether to automatically produce plot as SCnorm 
#'    determines the optimal number of groups (default is FALSE, highly 
#'    suggest using TRUE). Plots will be printed to the current device.
#' @param reportSF whether to provide a matrix of scaling counts in the
#'    output (default = FALSE).
#' @param FilterCellNum the number of non-zero expression estimate required
#'    to include the genes into the SCnorm fitting
#' (default = 10). The initial grouping fits a quantile regression to each
#'    gene, making this value too low gives unstable fits.
#' @param FilterExpression exclude genes having median of non-zero expression
#'    from the normalization.
#' @param Thresh threshold to use in evaluating the sufficiency of K, default
#'    is .1.
#' @param K the number of groups for normalizing. If left unspecified, an
#'    evaluation procedure will determine the optimal value of K
#'    (recommended). 
#' @param NCores number of cores to use, default is detectCores() - 1. 
#' This will be used to set up a parallel environment using either MulticoreParam (Linux, Mac) 
#' or SnowParam (Windows) with NCores using the package BiocParallel. 
#' @param ditherCounts whether to dither/jitter the counts, may be used for
#'    data with many ties, default is FALSE.
#' @param PropToUse proportion of genes closest to the slope mode used for
#'    the group fitting, default is set at .25. This number #' mainly affects
#'    speed.
#' @param Tau value of quantile for the quantile regression used to estimate
#'    gene-specific slopes (default is median, Tau = .5 ).
#' @param withinSample a vector of gene-specific features to correct counts
#'    within a sample prior to SCnorm. If NULL(default) then no correction will
#'    be performed. Examples of gene-specific features are GC content or gene
#'    length.
#' @param useSpikes whether to use spike-ins to perform across condition
#'    scaling (default=FALSE). SCnorm will assume spike-in names start with "ERCC-".
#' @param useZerosToScale whether to use zeros when scaling across conditions (default=FALSE).
#'
#' @description Quantile regression is used to estimate the dependence of
#'    read counts on sequencing depth for every gene. Genes with similar
#'     dependence are then grouped, and a second quantile regression is used to
#'    estimate scale factors within each group. Within-group adjustment for
#'    sequencing depth is then performed using the estimated scale factors to
#'    provide normalized estimates of expression. If multiple conditions are
#'    provided, normalization is performed within condition and then
#'    normalized estimates are scaled between conditions. If withinSample=TRUE
#'    then the method from Risso et al. 2011 will be implemented.


#' @return List containing matrix of normalized expression (and optionally a
#'    matrix of size factors if reportSF = TRUE ).
#' @export


#' @importFrom parallel detectCores
#' @import graphics
#' @import grDevices
#' @import stats
#' @importFrom BiocParallel bplapply  
#' @importFrom BiocParallel register
#' @importFrom BiocParallel MulticoreParam
#' @importFrom BiocParallel bpparam
#' @importFrom parallel detectCores
#' @importFrom S4Vectors metadata
#' @importFrom SummarizedExperiment SummarizedExperiment assayNames assays colData
#' @author Rhonda Bacher
#' @examples 
#'  
#'  data(ExampleData)
#'    Conditions = rep(c(1,2), each= 90)
#'    #DataNorm <- SCnorm(ExampleData, Conditions, 
#'    #FilterCellNum = 10)
#'    #str(DataNorm)

SCnorm <- function(Data=NULL, Conditions=NULL,
                    PrintProgressPlots=FALSE, reportSF=FALSE,
                    FilterCellNum=10, FilterExpression=0, Thresh=.1, 
                    K=NULL, NCores=NULL, ditherCounts=FALSE, 
                    PropToUse=.25, Tau=.5, 
                    withinSample=NULL, useSpikes=FALSE, useZerosToScale=FALSE) {
  
    if (is.null(Conditions)) {stop("Must supply conditions.")}
    
    if (is(Data, "SummarizedExperiment")) {
      if (is.null(  SummarizedExperiment::assayNames(Data)) || SummarizedExperiment::assayNames(Data)[1] != "Counts") {
        message("Renaming the first element in assays(Data) to 'Counts'")
          SummarizedExperiment::assayNames(Data)[1] <- "Counts"
  
      if (is.null(colnames(SCnorm::getCounts(Data)))) {stop("Must supply sample/cell names!")}

      }
    }
    
      
    if (!(is(Data, "SummarizedExperiment"))) {
      Data <- data.matrix(Data)
      Data <- SummarizedExperiment(assays=list("Counts"=Data))
     }
       
    ## Checks
    
    if(any(colSums(SCnorm::getCounts(Data)) == 0)) {stop("Data contains at least one 
            column will all zeros. Please remove these columns before 
            calling SCnorm(). Performing quality control on your data is highly recommended prior
            to running SCnorm!")}
      
    if(anyNA(SCnorm::getCounts(Data))) {stop("Data contains at least one value of NA. SCnorm is unsure how to proceed.")}
     
    if (is.null(NCores)) {NCores <- max(1, parallel::detectCores() - 1)}
    if (.Platform$OS.type == "windows") {NCores = 1}
    options(mc.cores = NCores)
    message(paste0("Setting up parallel computation using ", 
                      NCores, " cores" ))

    
    if (is.null(rownames(SCnorm::getCounts(Data)))) {stop("Must supply gene/row names!")}
    if (is.null(colnames(SCnorm::getCounts(Data)))) {stop("Must supply sample/cell names!")}

    if (ncol(SCnorm::getCounts(Data)) != length(Conditions)) {stop("Number of columns in 
      expression matrix must match length of conditions vector!")}
    
    if (!is.null(K)) {message(paste0("SCnorm will normalize assuming ",
      K, " is the optimal number of groups. It is not advised to set this."))}
    

    if (ditherCounts == TRUE) {RNGkind("L'Ecuyer-CMRG");
      set.seed(1);message("Jittering values introduces some randomness, 
        for reproducibility set.seed(1) has been set.")}
      
    Levels <- unique(Conditions) # Number of conditions
  
  
    # Option to normalize within samples:
    if(!is.null(withinSample)) {
        if(length(withinSample) == nrow(SCnorm::getCounts(Data))) {
          message("Using loess method described in ''GC-Content Normalization 
          for RNA-Seq Data'', Risso et al. to perform within-sample 
          normalization. For other options see the original publication and 
          package EDASeq." )
        
        S4Vectors::metadata(Data)[["OriginalData"]] <- Data
        SummarizedExperiment::assays(Data)[["Counts"]] = apply(SCnorm::getCounts(Data), 2, SCnorm::correctWithin, correctFactor = withinSample)
        
        } else{
          message("Length of withinSample should match the number of 
            genes in Data!")
        }
    }

    DataList <- lapply(seq_len(Levels), function(x) {
        SCnorm::getCounts(Data)[,which(Conditions == Levels[x])]}) # split conditions
    Genes <- rownames(SCnorm::getCounts(Data)) 
    
    SeqDepthList <- lapply(seq_len(Levels), function(x) {
        colSums(SCnorm::getCounts(Data)[,which(Conditions == Levels[x])])})

    if(any(do.call(c, SeqDepthList) <= 10000)) {
       warning("At least one cell/sample has less than 10,000 counts total. 
       Check the quality of your data or filtering criteria. 
       SCnorm may not be appropriate for your data (see vignette for details).")
     }
  
     NumZerosCellList <- lapply(seq_len(Levels), function(x) {
         colSums(DataList[[x]]!= 0) })
     if(any(do.call(c, NumZerosCellList) <= 100)) {
        warning("At least one cell/sample has less than 100 genes detected (non-zero). 
        Check the quality of your data or filtering criteria. 
        SCnorm may not be appropriate for your data (see vignette for details).")
      }
      
    message("Gene filter is applied within each condition.")
    GeneZerosList <- lapply(seq_len(Levels), function(x) {
         rowSums(DataList[[x]]!= 0) })
    MedExprList <- lapply(seq_len(Levels), function(x) {
        apply(DataList[[x]], 1, function(c) median(c[c != 0])) })
    GeneFilterList <- lapply(seq_len(Levels), function(x) {
        names(which(GeneZerosList[[x]] >= FilterCellNum & MedExprList[[x]] >= FilterExpression))})
  
    checkGeneFilter <- vapply(seq_len(Levels), function(x) {
             length(GeneFilterList[[x]])}, FUN.VALUE=numeric(1))
    if(any(checkGeneFilter < 100)) {
       stop("At least one condition has less then 100 genes that pass the specified filter. Check the quality of your data or filtering criteria. 
       SCnorm may not be appropriate for your data (see vignette for details).")
     }
       

    GeneFilterOUT <- lapply(seq_len(Levels), function(x) {
        names(which(GeneZerosList[[x]] < FilterCellNum | MedExprList[[x]] < FilterExpression))})
    names(GeneFilterOUT) <- paste0("GenesFilteredOutGroup", unique(Conditions))
  
 
    NM <- lapply(seq_len(Levels), function(x) {
        message(paste0(length(GeneFilterOUT[[x]]), 
           " genes in condition ", Levels[x]," will not be included in the normalization due to 
             the specified filter criteria."))})
  
    message("A list of these genes can be accessed in output, 
    see vignette for example.") 
    
    # Get median quantile regr. slopes.
    SlopesList <- lapply(seq_len(Levels), function(x) {
            getSlopes(Data = DataList[[x]][GeneFilterList[[x]],], 
                      SeqDepth = SeqDepthList[[x]], 
                      Tau=Tau, 
                      FilterCellNum=FilterCellNum, 
                      ditherCounts=ditherCounts)})
  
 
    # If k is NOT provided
    if (is.null(K)) {
        NormList <- lapply(seq_len(Levels), function(x) {
          normWrapper(Data = DataList[[x]], 
                      SeqDepth = SeqDepthList[[x]], 
                      Slopes = SlopesList[[x]],
                      CondNum = Levels[x], 
                      PrintProgressPlots = PrintProgressPlots,
                      PropToUse = PropToUse,
                      Tau = Tau, 
                      Thresh = Thresh, 
                      ditherCounts=ditherCounts)
      }) 
    }
    # If specific K then do:
    # If length of K is less than number of conditions, assume the same K.
    if (!is.null(K) ) {
      if (length(K) == length(Levels)) {
        NormList <- lapply(seq_len(Levels), function(x) {
          SCnormFit(Data = DataList[[x]], 
                    SeqDepth = SeqDepthList[[x]], Slopes = SlopesList[[x]],
                    K = K[x], PropToUse = PropToUse, 
                    ditherCounts=ditherCounts)
        })
      } else if (length(K) == 1) {
        K <- rep(K, length(Levels))
        NormList <- lapply(seq_len(Levels), function(x) {
          SCnormFit(Data = DataList[[x]], 
                    SeqDepth = SeqDepthList[[x]], Slopes = SlopesList[[x]],
                    K = K[x], PropToUse = PropToUse,
                    ditherCounts=ditherCounts)
        }) 
      } else (stop("Check that the specification of K is correct!"))
    }    
  
   FilterCellProportion = lapply(seq_len(Levels), function(x) {
        FilterCellNum / ncol(DataList[[x]])})
  
    NORMDATA <- do.call(cbind, lapply(seq_len(Levels), function(x) {
        NormList[[x]]$NormData}))
  
   
    if (length(Levels) > 1) {
    
      # Scaling
      ### Genes = Reduce(intersect, GeneFilterList)
      message("Scaling data between conditions...")
      ScaledNormData <- scaleNormMultCont(NormData = NormList, 
                        OrigData = SCnorm::getCounts(Data), 
                        Genes = Genes, useSpikes = useSpikes, 
                        useZerosToScale = useZerosToScale)
      names(ScaledNormData) <- c("NormalizedData", "ScaleFactors")
      
      NormDataFull <- ScaledNormData$NormalizedData
      ScaleFactorsFull <- ScaledNormData$ScaleFactors

      } else {
        NormDataFull <- NormList[[1]]$NormData
        ScaleFactorsFull <- NormList[[1]]$ScaleFactors
        }
    if(reportSF == FALSE) {
        ScaleFactorsFull <- NULL
        }
    
    
    # Return 
    S4Vectors::metadata(Data)[["NormalizedData"]] <- NormDataFull
    S4Vectors::metadata(Data)[["ScaleFactors"]] <- ScaleFactorsFull
    S4Vectors::metadata(Data)[["GenesFilteredOut"]] <- GeneFilterOUT

  
    message("Done!")
  
return(Data)  
  
}


