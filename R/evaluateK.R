#' @title Evaluate normalization using K slope groups

#' @param Data matrix of normalized expression counts. Rows are genes and
#'    columns are samples.
#' @param SeqDepth vector of sequencing depths estimated as columns sums of
#'    un-normalized expression matrix.
#' @param OrigData matrix of un-normalized expression counts. Rows are genes
#'    and columns are samples.
#' @param Slopes vector of slopes estimated in the GetSlopes() function. Only
#'    used here to obtain the names of genes considered in the
#' normalization.
#' @param Name plot title
#' @param Tau value of quantile for the quantile regression used to estimate
#'    gene-specific slopes (default is median, Tau = .5 ).
#' @param PrintProgressPlots whether to automatically produce plot as SCnorm 
#'    determines the optimal number of groups (default is FALSE, highly 
#'    suggest using TRUE). Plots will be printed to the current device.
#' @param ditherCounts whether to dither/jitter the counts, may be used for
#'    data with many ties, default is FALSE.

#' @description Median quantile regression is fit for each gene using the
#'    normalized gene expression values. A slope near zero indicate the
#'    sequencing depth effect has been successfully removed. 
#'  Genes are divided into ten equally sized groups based on their non-zero
#'    median expression. Slope densities are plot for each group and estimated
#'    modes are calculated. If any of the ten group modes is larger than .1, the
#'    K is not sufficient to normalize the data.
#' @return value of largest mode and a plot of the ten normalized slope
#'    densities. 
#' @author Rhonda Bacher



evaluateK <- function(Data, SeqDepth, OrigData, Slopes, Name, Tau, PrintProgressPlots, ditherCounts) {
    

    Genes <- names(Slopes) #Genes for normalizing
    NormSlopes <- getSlopes(Data[Genes,], SeqDepth, Tau=.5, FilterCellNum = 0, ditherCounts)
    
    colors <- colorRampPalette(c("#00C3FF", "blue","black", "#FF0700"), 
            bias=2)(n = 10)
        
    MedExpr <- apply(OrigData, 1, function(x) median(x[x != 0])) 

    ExprGroups <- splitGroups(MedExpr[intersect(names(MedExpr), names(NormSlopes))], 10)
    Modes <- getDens(ExprGroups, NormSlopes, "Mode") 

    if (PrintProgressPlots == TRUE) {
      XX <- generateEvalPlot(MedExpr = MedExpr, 
         SeqDepth = SeqDepth, Slopes = NormSlopes, 
         Name = Name, NumExpressionGroups = 10, BeforeNorm = FALSE)
      }
      
    return(Modes)
}
    