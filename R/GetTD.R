#' @title Fit group regression for specific quantile and degree
#' @usage GetTD(x, InputData)
#' @param x specifies a column of the grid matrix of tau and d.
#' @param InputData contains the expression values, sequencing depths to fit the group regression, and the quantile used in the individual 
#' gene regression for grouping.

#' @description This is an internal fitting of the group regression. For a single combination of possible tau and d values the group regression is fist fit, then predicted values are obtained and regressed against the original sequencing depths. The estimates slope is passed back to the SCnorm_fit() function.
#' @author Rhonda Bacher
#' @export

GetTD <- function(x, InputData) {


	TauGroup <- seq(.05, .95, by =.05)
	D <- x
	
	O <- InputData[[1]]
	Y <- InputData[[2]]
	SeqDepth <- InputData[[3]]
	Tau <- InputData[[4]]
	ditherFlag <- InputData[[5]]
	
	polyX <- try(poly(O, degree = D, raw = FALSE), silent=T)
	
	if(!is.null(dim(polyX))){
		Xmat <- data.table(model.matrix( ~ polyX ))
	
		polydata <- data.frame(Y = Y, Xmat = Xmat[,-1])
	
		if(ditherFlag == TRUE) {
			rqfit <- try(rq(dither(Y, type="symmetric", value=.01) ~ ., data = polydata, na.action = na.exclude, tau = TauGroup, method="fn"), silent=T)
		} else {
			rqfit <- try(rq(Y ~ ., data = polydata, na.action = na.exclude, tau = TauGroup, method="fn"), silent=T)
		}
		if(class(rqfit) != "try-error"){
			revX <- data.frame(predict(polyX, SeqDepth))
					
			colnames(revX) <- colnames(polydata[-1])
			pdvalsrq <- predict(rqfit, newdata=data.frame(revX))
			rownames(pdvalsrq) <- colnames(SeqDepth)
	
			checkS <- apply(pdvalsrq, 2, function(y) {
						if (min(y) > 0) {
								S <- rq(y ~ SeqDepth, tau = Tau)$coef[2]
							} else {S <- -50} })
		} else {checkS <- rep(-50, length(TauGroup)) }
	} else {checkS <- rep(-50, length(TauGroup))}

return(as.numeric(checkS))
}

