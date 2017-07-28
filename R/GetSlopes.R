#' @title Estimate gene specific count-depth relationships
#' @param Data matrix of un-normalized expression counts. Rows are genes and
#'    columns are samples.
#' @param SeqDepth vector of sequencing depths estimated as columns sums of
#'    un-normalized expression matrix.
#' @param FilterCellNum the number of non-zero expression estimate required to
#'    include the genes into the SCnorm fitting (default = 10). The initial 
#' @param Tau value of quantile for the quantile regression used to estimate
#'    gene-specific slopes (default is median, Tau = .5 ). 
#' @param ditherCounts whether to dither/jitter the counts, may be used for data
#'    with many ties, default is FALSE.

#' @description This is the gene-specific fitting function, where a median 
#'    (Tau = .5) quantile regression is fit for each gene. Only genes having at
#'    least 10 non-zero expression values are considered.
#' @return vector of estimated slopes.
#' @author Rhonda Bacher
#' @export
#' @importFrom quantreg rq dither
#' @importFrom BiocParallel bplapply
#' @examples 
#'  data(ExampleSimSCData)
#'  myslopes <- getSlopes(ExampleSimSCData)

getSlopes <- function(Data, SeqDepth = 0, Tau = .5, FilterCellNum = 10, ditherCounts=FALSE) {

    if(any(SeqDepth==0))
        {SeqDepth = colSums(Data)}
    
    NumNonZeros <- rowSums(Data!=0)
    Genes <- names(which(NumNonZeros >= FilterCellNum)) ##filter for now    
    LogData <- redoBox(Data, 0) #log data
    
    AllReg <- unlist(BiocParallel::bplapply(X = seq_len(length(Genes)), FUN = quickReg, 
                InputData = list(LogData, SeqDepth, Genes, Tau, ditherCounts)))
    
    return(AllReg)
}




