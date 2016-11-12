#' @title Iteratively fit group regression and evaluate to choose optimal K
#' @usage Normalize(Data, SeqDepth, Slopes, CondNum, PLOT, PropToUse, outlierCheck, Tau)
#' @param Data matrix of un-normalized expression counts. Matrix rows are genes and columns are samples. 
#' @param SeqDepth vector of sequencing depths estimated as columns sums of un-normalized expression matrix.
#' @param Slopes vector of slopes estimated in the GetSlopes() function.
#' @param Condnum which condition is being normalized, used for plotting only.
#' @param PLOT whether to output the evaluation plots.

#' @param PropToUse proportion of genes closest to the slope mode used for the group fitting, default is set at .25. This number mainly affects speed. 
#' @param outlierCheck cells/samples with relatively small sequencing depths may recieve very small scaling factors, to ensure these 
#' small scaling factors do not create outliers, genes with potential outliers are first flagged and then if necessary, their normalized expression 
#' values are corrected. First, gene expression counts in cells/samples having the three smallest sequencing depths are flagged if they contain
#' expression counts larger than 3 times the predicted counts (from the group regression). For each gene, if the flagged values
#' are the largest normalized expression value, then the count is corrected to be the median of the non-zero normalized values (default = 
#' TRUE)
#' @param Tau value of quantile for the quantile regression used to estimate gene-specific slopes (default is median, Tau = .5 ). 



#' @description This function iteratively normalizes using K groups and then evaluates whether K is sufficient. If the maximum mode received 
#' from the GetK() function is larger than .1, K is increased to K + 1. 
#' @return matrix of normalized and scaled expression values for all conditions. If PLOT = TRUE then the evaluation plots are 
#' output 
#' for each attempted value of K.
#' @author Rhonda Bacher
#' @export

Normalize <- function(Data, SeqDepth, Slopes, CondNum, PLOT = TRUE, OutputName = OutputName, PropToUse, outlierCheck, Tau, NCores) {
	# Set up
	GetMax = 1
	i = 0
	
		while(GetMax > .1) {
			i = i + 1 

			NormDataList <- SCnorm_fit(Data = Data, SeqDepth = SeqDepth, Slopes = Slopes, K = i, 
								PropToUse = PropToUse, outlierCheck = outlierCheck, Tau = Tau, NCores = NCores)
			
			NAME = paste0("Condition: ", CondNum, "\n K = ", i)
		
			GetMax <- GetK(NormDataList$NormData, SeqDepth, Data, Slopes, NAME, PLOT = PLOT, Tau=Tau, NCores = NCores)
		}
	
return(NormDataList)
}

