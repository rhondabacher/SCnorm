#' quickreg
#'
#' Perform the single gene regressions using quantile regression.
#' 
#' @details Perform the single gene regressions using quantile regression.
#' 
#' @param x gene to perform the regression on.
#' @param InputData list of data needed for the regression.
#'   
#' @return gene slope.


quickreg<-function(x,InputData)
{
    LogData = InputData[[1]]
    SeqDepth = InputData[[2]]
    X = InputData[[3]][x]
    Tau = InputData[[4]]
    ditherFlag = InputData[[5]]

    if(ditherFlag == TRUE) {
      slope <- try(rq(dither(LogData[X, ], type="symmetric", value=.01) ~ 
                  log(SeqDepth), tau = Tau, method="br")$coef[2], silent=TRUE)
      if(class(slope) == "try-error"){
        slope <- try(rq(dither(LogData[X, ], type="symmetric", value=.01) ~ 
                    log(SeqDepth), tau = Tau, method="fn")$coef[2], silent=TRUE)
        if(class(slope) == "try-error"){
          slope <- NA
        }    
      } 
    } else {
      slope <- try(rq(LogData[X, ] ~ log(SeqDepth), tau = Tau, 
              method="br")$coef[2], silent=TRUE)
      if(class(slope) == "try-error"){
        slope <- try(rq(LogData[X, ] ~ log(SeqDepth), tau = Tau, 
                method="fn")$coef[2], silent=TRUE)
        if(class(slope) == "try-error"){
          slope <- NA
        }    
      }
    }
    names(slope) <- X

    return(slope)
}

#' redobox
#' 
#' @details Function to log data and turn zeros to NA to mask/ignore in 
#'  functions.
#' 
#' @param DATA data set to.
#' @param smallc, what value to ignore, typically is zero.
#'   
#' @return the dataset has been logged with values below smallc masked. 


redobox <- function(DATA, smallc) {

    DATA[DATA <= smallc] <- NA #truncate some values first, usually just zero
    y <- log(DATA)

    return(y)
}
