#' @title Estimate gene specific count-depth relationships
#' @usage GetSlopes(Data, SeqDepth, Tau, FilterCellNum, NCores, ditherCounts)
#' @param Data matrix of un-normalized expression counts. Rows are genes and
#'    columns are samples.
#' @param SeqDepth vector of sequencing depths estimated as columns sums of
#'    un-normalized expression matrix.
#' @param FilterCellNum the number of non-zero expression estimate required to
#'    include the genes into the SCnorm fitting (default = 10). The initial 
#' @param Tau value of quantile for the quantile regression used to estimate
#'    gene-specific slopes (default is median, Tau = .5 ). 
#' @param NCores number of cores to use, default is detectCores() - 1.
#' @param ditherCounts whether to dither/jitter the counts, may be used for data
#'    with many ties, default is FALSE.

#' @description This is the gene-specific fitting function, where a median 
#'    (Tau = .5) quantile regression is fit for each gene. Only genes having at
#'    least 10 non-zero expression values are considered.
#' @return vector of estimated slopes.
#' @author Rhonda Bacher
#' @export
#' @importFrom quantreg rq
#' @importFrom BiocParallel bplapply
#' @examples 
#'  data(ExampleData)
#'  myslopes <- GetSlopes(ExampleData)

GetSlopes <- function(Data, SeqDepth = 0, Tau = .5, FilterCellNum = 10, 
    NCores = 1, ditherCounts=FALSE) {

    if(any(SeqDepth==0))
        {SeqDepth = colSums(Data)}
    
    NumNonZeros <- apply(Data, 1, function(x) sum(x != 0))

    Genes <- names(which(NumNonZeros >= FilterCellNum)) ##filter for now    
    LogData <- redobox(Data, 0) #log data
    
    AllReg <- unlist(bplapply(X = 1:length(Genes), FUN = quickreg, 
                InputData = list(LogData, SeqDepth, Genes, Tau, ditherCounts)))
    
    return(AllReg)
}



