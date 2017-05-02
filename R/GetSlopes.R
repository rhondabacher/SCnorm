#' @title Estimate gene specific count-depth relationships
#' @usage GetSlopes(Data, SeqDepth, FilterCellNum = 10, Tau, NCores, ditherCounts)
#' @param Data matrix of un-normalized expression counts. Rows are genes and columns are samples.
#' @param SeqDepth vector of sequencing depths estimated as columns sums of un-normalized expression matrix.
#' @param FilterCellNumber the number of non-zero expression estimate required to include the genes into the SCnorm fitting (default = 10). The initial 
#' @param Tau value of quantile for the quantile regression used to estimate gene-specific slopes (default is median, Tau = .5 ). 
#' @param NCores number of cores to use, default is detectCores() - 1.
#' @param ditherCounts whether to dither/jitter the counts, may be used for data with many ties, default is FALSE.

#' @description This is the gene-specific fitting function, where a median (Tau = 
#' .5) quantile regression is fit for each gene. Only genes having at least 10 non-zero expression values are considered.
#' @return vector of estimated slopes.
#' @author Rhonda Bacher
#' @export


GetSlopes <- function(Data, SeqDepth = 0, Tau = .5, FilterCellNum = 10, NCores, ditherCounts) {

	if(any(SeqDepth==0))
		{SeqDepth = colSums(Data)}
	
	NumNonZeros <- apply(Data, 1, function(x) sum(x != 0))

	Genes <- names(which(NumNonZeros >= FilterCellNum)) ##filter for now	
	LogData <- redobox(Data, 0) #log data
    
	if (.Platform$OS.type == "windows") {
		NCores = 1
	}
	
	AllReg <- unlist(mclapply(1:length(Genes), function(x) {
        quickreg(InputData = list(LogData[Genes[x],], SeqDepth, Genes[x], Tau, ditherCounts))}, mc.cores = NCores))
	
	return(AllReg)
}



