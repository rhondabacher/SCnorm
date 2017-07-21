#' @title Iteratively fit group regression and evaluate to choose optimal K
#' @inheritParams SCnorm
#' @param SeqDepth sequencing depth for each cell/sample.
#' @param Slopes per gene estimates of the count-depth relationship.
#' @param CondNum name of group being normalized, just for printing messages.

#' @description This function iteratively normalizes using K groups and then
#'    evaluates whether K is sufficient. If the maximum mode received 
#' from the GetK() function is larger than .1, K is increased to K + 1. Uses
#'    params sent from SCnorm.
#' @return matrix of normalized and scaled expression values for all conditions
#'    and the evaluation plots are 
#' output for each attempted value of K.
#' @author Rhonda Bacher

normWrapper <- function(Data, SeqDepth=NULL, Slopes=NULL, CondNum=NULL, 
                        PrintProgressPlots, PropToUse, Tau,
                        Thresh, ditherCounts) {
    
    # Set up
    MaxMode = 1
    i = 0
    
    message(paste0("Finding K for Condition ", CondNum))
        while(MaxMode > Thresh) {
            i = i + 1 
            message(paste0("Trying K = ", i))
            NormDataList <- SCnormFit(Data = Data, SeqDepth = SeqDepth, 
                Slopes = Slopes, K = i, PropToUse = PropToUse, Tau = Tau,
                ditherCounts=ditherCounts)
            
            NAME = paste0("Condition: ", CondNum, "\n K = ", i)
        
            Modes <- evaluateK(Data = NormDataList$NormData, 
                               SeqDepth = SeqDepth, OrigData = Data, 
                               Slopes = Slopes, Name = NAME, 
                               Tau=Tau, PrintProgressPlots= PrintProgressPlots,
                               ditherCounts=ditherCounts)

             MaxMode <- max(abs(Modes))
             if(i > 25) {stop("SCnorm is unable to converge. 
                         Consider altering the filter criteria such as FilterExpression. 
                         See vigenette for additional details.")} 
          }        
    
    return(NormDataList)
}

