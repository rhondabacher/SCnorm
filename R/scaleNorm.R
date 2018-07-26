#' @title Scale multiple conditions
#' @param NormData list of matrices of normalized expression counts and scale
#'    factors for each condition. Matrix rows are genes and 
#'    columns are samples. 
#' @param OrigData list of matrices of un-normalized expression counts. Matrix
#'    rows are genes and columns are samples. Each item in list is a different
#'    condition.
#' @param Genes vector of genes that will be used to scale conditions, only
#'    want to use genes that were normalized.
#' @param useSpikes whether to use spike-ins to perform between condition
#'    scaling (default=FALSE). Assumes spike-in names start with "ERCC-".
#' @param useZerosToScale whether to use zeros when scaling across conditions (default=FALSE).
#'
#' @description After conditions are independently normalized with the 
#' count-depth 
#' effect removed, conditions need to be additionally scaled prior to further
#'    analysis. Genes that were normalized in both
#' conditions are split into quartiles based on their un-normalized non-zero
#'    medians. Genes in each quartile are scaled to the 
#' median fold change of 
#' condition specific gene means and overall gene means.
#' @return matrix of normalized and scaled expression values for all
#'    conditions.
#' @author Rhonda Bacher
#' @import SingleCellExperiment



scaleNormMultCont <- function(NormData, OrigData, Genes, useSpikes, useZerosToScale) {

    NumCond <- length(NormData)
    NumGroups <- 4
    MedExpr <- apply(counts(OrigData)[Genes,], 1, function(x) median(x[x != 0]))
    ExprGroups <- splitGroups(MedExpr, NumGroups = NumGroups)
    
    FullNormMat <- do.call(cbind, lapply(seq_len(NumCond), function(x) {cbind(NormData[[x]]$NormData)[Genes,]}))

    ScaledDataList <- vector("list", NumCond)
    ScaledFacsList <- vector("list", NumCond)
    
    for(i in seq_len(NumCond)){

        CondNormMat <- NormData[[i]]$NormData
        CondScaleMat <- NormData[[i]]$ScaleFactors
    
        for(r in seq_len(NumGroups)){
             qgenes <- names(ExprGroups[[r]])
             scalegenes <- qgenes
             if(useSpikes==TRUE) {
                 scalegenes <- qgenes[Genes[which(SingleCellExperiment::isSpike(OrigData))]] #which are spikes
                 if(length(scalegenes) <= 5) {
                     warnings("Not enough spike-ins or spike-ins do not 
                     span range of expression to perform reasonably well. 
                     Using all genes instead. SCnorm expects spike-ins to
                     contain 'ERCC-' in the gene name.")
                     scalegenes <- qgenes
                 }
             }

             if (useZerosToScale) {
               spCond <- apply(CondNormMat[scalegenes,], 1, function(x) mean(x))
               allCond <- apply(FullNormMat[scalegenes,], 1, function(x) mean(x))
              } else {
                spCond <- apply(CondNormMat[scalegenes,], 1, function(x) mean(x[x != 0]))
                allCond <- apply(FullNormMat[scalegenes,], 1, function(x) mean(x[x != 0]))
              }
             condFac <- median(spCond / allCond, na.rm=TRUE)
             CondNormMat[qgenes,] <- round((CondNormMat[qgenes,] / condFac), 2)
             CondScaleMat[qgenes,] <- CondScaleMat[qgenes,] * condFac
        }

        ScaledDataList[[i]] <- CondNormMat[Genes,] #ensures order remains the same here
    
        ScaledFacsList[[i]] <- CondScaleMat[Genes,] #ensures order remains the same here
    }


    ScaleNormMatAll <- do.call(cbind, ScaledDataList)
    ScaleNormFacsAll <- do.call(cbind, ScaledFacsList)

    ScaledDataList <- list(ScaledData = ScaleNormMatAll, ScaleFactors = ScaleNormFacsAll)
    return(ScaledDataList)
}

