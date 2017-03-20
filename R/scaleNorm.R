#' @title Scale multiple conditions
#' @usage scaleNormMultCont(NormData, OrigData, Genes)
#' @param NormData list of matrices of normalized expression counts and scale factors for each condition. Matrix rows are genes and columns are samples. 
#' @param OrigData list of matrices of un-normalized expression counts. Matrix rows are genes and columns are samples. Each item in list is a different condition.
#' @param Genes vector of genes that will be used to scale conditions, only want to use genes that were normalized.
#' @param useSpikes whether to use spike-ins to perform between condition scaling (default=FALSE). Assumes spike-in names start with "ERCC-".

#' @description After conditions are independtly normalized with the 
#' count-depth 
#' effect removed, conditions need to be additionally scaled prior to further analysis. Genes that were normalized in both
#' conditions are split into quartiles based on their un-normalized non-zero medians. Genes in each quartile are scaled to the 
#' median fold change of 
#' condition specific gene means and overall gene means.
#' @return matrix of normalized and scaled expression values for all conditions.
#' @author Rhonda Bacher
#' @export

scaleNormMultCont <- function(NormData, OrigData, Genes, useSpikes)
{
sreg <- list()
NumCond <- length(NormData)
AllGenes <- Genes
K <- 4
avgexp <- log(apply(OrigData[Genes,], 1, function(x) median(x[x!=0]))) #conditional median
groups <- K
splitby <- sort(avgexp)
splitS <- length(splitby)/groups
sreg <- split(splitby, ceiling(seq_along(splitby) / splitS))

# Need to put a check here later on to make sure the rownames are in the same order
OC <- do.call(cbind, lapply(1:NumCond, function(x) {cbind(NormData[[x]]$NormData)}))

ScaleMat<-list()
ScaleFacs<-list()
for(i in 1:NumCond){

	C1 <- NormData[[i]]$NormData
	SF <- NormData[[i]]$ScaleFactors
	
	for(r in 1:groups){
		 qgenes <- names(sreg[[r]])
		 scalegenes <- qgenes
		 if(useSpikes==TRUE) {
			 scalegenes <- qgenes[grep("ERCC-", qgenes)] #which are spikes
			 if(length(scalegenes) <= 5) {
			 	stop("Not enough spike-ins or spike-ins do not span range of expression to perform reasonably well.")
			 }
		 }

		 ss1 <- apply(C1[scalegenes,], 1, function(x) mean(x[x != 0]))
		 os <- apply(OC[scalegenes,], 1, function(x) mean(x[x != 0]))
	
		 rr <- median(ss1/os, na.rm=T); #print(rr)
		 C1[qgenes,] <- round((C1[qgenes,] / rr),2 )
		 SF[qgenes,] <- SF[qgenes,] * rr
	}

	ScaleMat[[i]] <- C1[AllGenes,] #ensures order remains the same here
	
	ScaleFacs[[i]] <- SF[AllGenes,] #ensures order remains the same here
}


ScaleMatAll <- do.call(cbind, ScaleMat)
ScaleFacsAll <- do.call(cbind, ScaleFacs)

ScaledData <- list(ScaledData = ScaleMatAll, ScaleFactors = ScaleFacsAll)
return(ScaledData)
}

