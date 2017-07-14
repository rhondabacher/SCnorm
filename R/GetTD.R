#' @title Fit group regression for specific quantile and degree
#' @usage GetTD(x, InputData)
#' @param x specifies a column of the grid matrix of tau and d.
#' @param InputData contains the expression values, sequencing depths to fit 
#'  the group regression, and the quantile used in the individual 
#'  gene regression for grouping.

#' @description This is an internal fitting of the group regression. For a
#' single combination of possible tau and d values the group regression is
#' fist fit, then predicted values are obtained and regressed against the
#' original sequencing depths. The estimates slope is passed back to the
#' SCnorm_fit() function.
#' @return estimated count-depth relationship of predicted values for one
#' value of tau and degree.
#' @author Rhonda Bacher
#' @importFrom data.table data.table

GetTD <- function(x, InputData) {


    TauGroup <- InputData[[4]][x,1]
    DG <- InputData[[4]][x,2]
    
    O <- InputData[[1]]
    Y <- InputData[[2]]
    SeqDepth <- InputData[[3]]
    Tau <- InputData[[5]]
    ditherFlag <- InputData[[6]]
    
    polyX <- try(poly(O, degree = DG, raw = FALSE), silent=TRUE)
    
    if(!is.null(dim(polyX))){
        Xmat <- data.table(model.matrix( ~ polyX ))
    
        polydata <- data.frame(Y = Y, Xmat = Xmat[,-1])
    
        if(ditherFlag == TRUE) {
            rqfit <- try(rq(dither(Y, type="symmetric", value=.01) ~ ., 
            data = polydata, na.action = na.exclude, tau = TauGroup, 
            method="fn"), silent=TRUE)
        } else {
            rqfit <- try(rq(Y ~ ., data = polydata, na.action = na.exclude,
                 tau = TauGroup, method="fn"), silent=TRUE)
        }
        if(class(rqfit) != "try-error"){
            revX <- data.frame(predict(polyX, SeqDepth))
                    
            colnames(revX) <- colnames(polydata[-1])
            pdvalsrq <- predict(rqfit, newdata=data.frame(revX))

            names(pdvalsrq) <- colnames(SeqDepth)
    
            if (min(pdvalsrq) > 0) { 
                S <- rq(pdvalsrq ~ SeqDepth, tau = Tau)$coef[2]
            } else {S <- -50}
        } else {S <- -50}    
    } else {S <- -50}

    return(as.numeric(S))
}
