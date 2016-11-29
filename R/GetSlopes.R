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



quickreg<-function(x,InputData)
{
	LogData = InputData[[1]]
	SeqDepth = InputData[[2]]
	X = InputData[[3]][x]
	Tau = InputData[[4]]
	
	slope <- rq(LogData[X, ] ~ log(SeqDepth), tau = Tau)$coef[2] 
	names(slope) <- X
	
	return(slope)
}

redobox <- function(DATA, smallc) {

  DATA[DATA <= smallc] <- NA #truncate some values first, usually just zero
  y <- log(DATA)
  
return(y)
}
