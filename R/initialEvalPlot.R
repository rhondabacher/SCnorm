#' @title Evaluate data before normalization.
#' @usage initialEvalPlot(MedExpr, SeqDepth, Slopes, Name, 
#'    NumExpressionGroups = 10, BeforeNorm = TRUE)

#' @inheritParams checkCountDepth
#' @param MedExpr non-zero median expression for all genes.
#' @param SeqDepth sequencing depth for each cell/sample.
#' @param Slopes per gene estimates of the count-depth relationship.
#' @param Name name for plot title.
#' @param BeforeNorm whether dat have already been normalized.

#' @description Genes are divided into NumExpressionGroups = 10 equally sized
#'    groups based on their non-zero median expression. Slope densities are plot
#'    for each group. 
#' @return a plot of the un-normalized slope densities. 
#' @author Rhonda Bacher
#' 
#' @import grDevices

initialEvalPlot <- function(MedExpr, SeqDepth, Slopes, Name, 
    NumExpressionGroups = 10, BeforeNorm = TRUE) {
    
    sreg <- list()
        
    colors <- colorRampPalette(c("#00C3FF", "blue","black", "#FF0700"), 
                bias=2)(n = NumExpressionGroups)
        
    splitby <- sort(MedExpr[intersect(names(MedExpr), names(Slopes))])
    grps <- length(splitby) / NumExpressionGroups
    sreg <- split(splitby, ceiling(seq_along(splitby) / grps))
        
    Mode <- c()
    DensH <- c()
    for (i in 1:NumExpressionGroups) {
      useg <- names(sreg[[i]])
      rqdens <- density(na.omit(Slopes[useg]))
      peak <- which.max(rqdens$y)
      Mode[i] <- rqdens$x[peak]
      DensH[i] <- rqdens$y[peak]
      }

    YMax <- pmin(round(max(DensH), 2) + .2, 10) #just for plotting

    plot(density(na.omit(Slopes), from = -3, to = 3), xlim = c(-3,3), 
                ylim = c(0,YMax), lwd = 3, col = "white", xlab = "Slope",
                cex.axis = 1.5, main = paste0(Name), cex.lab = 1.5, cex.main = 1.7)
    for (i in 1:length(sreg)) {
        useg <- names(sreg[[i]])
        lines(density(na.omit(Slopes[useg]), from=-3, to=3, adjust=1), 
                    lwd=3, col=colors[i])
    }
    if(BeforeNorm == TRUE) {abline(v=1, lwd=3, col="black")}
    if(BeforeNorm == FALSE) {abline(v=0, lwd=3, col="black")}

}
    