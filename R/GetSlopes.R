#' @title Estimate gene specfic count-depth relationships
#' @usage GetSlopes(Data, SeqDepth, Tau)
#' @param Data matrix of un-normalized expression counts. Rows are genes and columns are samples.
#' @param SeqDepth vector of sequencing depths estimated as columns sums of un-normalized expression matrix.
#' @param Tau value of quantile for the quantile regression used to estimate gene-specific slopes (default is median, Tau = .5 ). 

#' @description This is the gene-specific fitting function, where a median (Tau = 
#' .5) quantile regression is fit for each gene. Only genes having at least 10 non-zero expression values are considered.
#' @return vector of estimated slopes.
#' @author Rhonda Bacher
#' @export


GetSlopes <- function(Data, SeqDepth = 0, Tau = .5, NCores) {

	if(any(SeqDepth==0))
		{SeqDepth = colSums(Data)}
	
	NumNonZeros <- apply(Data, 1, function(x) sum(x != 0))

	Genes <- names(which(NumNonZeros >= 10)) ##filter for now	
	LogData <- redobox(Data, 0) #log data
    
	if (.Platform$OS.type == "windows") {
		NCores = 1
	}
	AllReg <- unlist(mclapply(X = 1:length(Genes), FUN = quickreg, InputData = list(LogData, SeqDepth, Genes, Tau), mc.cores = NCores))
	
	return(AllReg)
}



