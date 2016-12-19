#' @title Evaluate data before normalization.
#' @usage initialEvalPlot(Data, SeqDepth, Slopes, Name, NumExpressionGroups = 10)
#' @param Data matrix of normalized expression counts. Rows are genes and columns are samples.
#' @param SeqDepth vector of sequencing depths estimated as columns sums of un-normalized expression matrix.
#' @param Slopes vector of slopes estimated in the GetSlopes() function. 
#' @param Name plot title
#' @param NumExpressionGroups number of (non-zero median) expression groups to split the slope densities into. 
#' @param BeforeNorm whether data are already normalized or not.


#' @description Genes are divided into NumExpressionGroups = 10 equally sized groups based on their non-zero median expression. Slope densities are plot for each group. 
#' @return a plot of the un-normalized slope densities. 
#' @author Rhonda Bacher
#' @export


initialEvalPlot <- function(MedExpr, SeqDepth, Slopes, Name, NumExpressionGroups = 10, BeforeNorm = TRUE) {
	
	sreg <- list()
		
	colors <- colorRampPalette(c("#00C3FF", "blue","black", "#FF0700"), bias=2)(n = NumExpressionGroups)
		
	splitby <- sort(MedExpr)
	grps <- length(splitby) / NumExpressionGroups
	sreg <- split(splitby, ceiling(seq_along(splitby) / grps))
		
	Mode <- c()
	DensH <- c()
	for (i in 1:NumExpressionGroups) {
	  useg <- names(sreg[[i]])
	  rqdens <- density(Slopes[useg])
	  peak <- which.max(rqdens$y)
	  Mode[i] <- rqdens$x[peak]
	  DensH[i] <- rqdens$y[peak]
  	}

	# if point mass occurs, arbitrarily 10
	# this is to avoid plotting it
	if(any(DensH > 10)) {
		XX = which(DensH > 10)
		TOKEEP <- setdiff(1:10, XX)

		DensH <- DensH[TOKEEP]	
		Mode <- Mode[TOKEEP]
		colors = colors[TOKEEP]

		sreg <- sreg[TOKEEP]
	}
	
	
	YMax <- round(max(DensH), 2) + .2 #just for plotting

	plot(density(na.omit(Slopes), from = -3, to = 3), xlim = c(-3,3), ylim = c(0,YMax), lwd = 3, col = "white", 
		xlab = "Slope",	cex.axis = 2, main = paste0(Name), cex.lab = 2, cex.main = 2)
	for (i in 1:NumExpressionGroups) {
		useg <- names(sreg[[i]])
		lines(density(na.omit(Slopes[useg]), from=-3, to=3, adjust=1), lwd=3, col=colors[i])
	}
	if(BeforeNorm == TRUE) {abline(v=1, lwd=3, col="black")}
	if(BeforeNorm == FALSE) {abline(v=0, lwd=3, col="black")}

}
	