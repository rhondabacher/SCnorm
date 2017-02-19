#' @title Evaluate normalization using K slope groups
#' @usage GetK(Data, SeqDepth, OrigData, Slopes, Name, PLOT = TRUE)
#' @param Data matrix of normalized expression counts. Rows are genes and columns are samples.
#' @param SeqDepth vector of sequencing depths estimated as columns sums of un-normalized expression matrix.
#' @param OrigData matrix of un-normalized expression counts. Rows are genes and columns are samples.
#' @param Slopes vector of slopes estimated in the GetSlopes() function. Only used here to obtain the names of genes considered in the
#' normalization.
#' @param Name plot title
#' @param PLOT whether to output the evaluation plots.

#' @description Median quantile regresion is fit for each gene using the normalized gene expression values. 
#' A slope near zero indicate the sequencing depth effect has been successfully removed. 
#' Genes are divided into ten equally sized groups based on their non-zero median expression. Slope densities are plot for each 
#' group (if PLOT = TRUE) and estimated modes are calculated. If any of the ten group modes is larger than .1, the K is not sufficient
#' to normalize the data.
#' @return value of largest mode and (if PLOT = TRUE) a plot of the ten normalized slope densities. 
#' @author Rhonda Bacher
#' @export


GetK <- function(Data, SeqDepth, OrigData, Slopes, Name, PLOT = TRUE, Tau=Tau, NCores) {
	
	sreg <- list()
	
	Genes <- names(Slopes) #Genes for normalizing
	
	LogData <- redobox(Data, 0) #LOG data
		
	NormSlopes <- unlist(mclapply(X = 1:length(Genes), FUN = quickreg, 
					InputData = list(LogData, SeqDepth, Genes, Tau), mc.cores = NCores))
	
	colors <- colorRampPalette(c("#00C3FF", "blue","black", "#FF0700"), bias=2)(n = 10)
		
	MedExp <- apply(OrigData, 1, function(x) median(x[x != 0])) 
	splitby <- sort(MedExp[Genes]) 
	grps <- length(splitby) / 10
	sreg <- split(splitby, ceiling(seq_along(splitby) / grps))
		
	Mode <- c()
	DensH <- c()
	for (i in 1:10) {
	  useg <- names(sreg[[i]])
	  rqdens <- density(NormSlopes[useg])
	  peak <- which.max(rqdens$y)
	  Mode[i] <- rqdens$x[peak]
	  DensH[i] <- rqdens$y[peak]
  	}

	# # if point mass occurs, arbitrarily 10
# 	# this is to avoid plotting it
# 	if(any(DensH > 10)) {
# 		XX = which(DensH > 10)
# 		TOKEEP <- setdiff(1:10, XX)
#
# 		DensH <- DensH[TOKEEP]
# 		Mode <- Mode[TOKEEP]
# 		colors = colors[TOKEEP]
#
# 		sreg <- sreg[TOKEEP]
# 	}
	
	
	
	YMax <- round(max(DensH), 2) + .2 #just for plotting
	if(PLOT == TRUE) {	
		plot(density(na.omit(NormSlopes), from = -3, to = 3), xlim = c(-3,3), ylim = c(0,YMax), lwd = 3, col = "white", 
			xlab = "Slope",	cex.axis = 2, main = paste0(Name), cex.lab = 2, cex.main = 2)
		for (i in 1:length(sreg)) {
			useg <- names(sreg[[i]])
			lines(density(na.omit(NormSlopes[useg]), from=-3, to=3, adjust=1), lwd=3, col=colors[i])
		}
		abline(v=0, lwd=3, col="black") }

	MAX <- max(abs(Mode))
	return(MAX)
}
	