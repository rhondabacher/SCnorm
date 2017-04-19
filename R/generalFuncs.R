

quickreg<-function(x,InputData)
{
  LogData = InputData[[1]]
  SeqDepth = InputData[[2]]
  X = InputData[[3]][x]
  Tau = InputData[[4]]
  ditherFlag = InputData[[5]]
  
  if(ditherFlag == TRUE) {
    slope <- try(rq(dither(LogData[X, ], type="symmetric", value=.01) ~ log(SeqDepth), tau = Tau, method="br")$coef[2], silent=T) 
    if(class(slope) == "try-error"){
      slope <- try(rq(dither(LogData[X, ], type="symmetric", value=.01) ~ log(SeqDepth), tau = Tau, method="fn")$coef[2], silent=T)  
      if(class(slope) == "try-error"){
		slope <- NA
	  }	
	} 
  } else {
    slope <- try(rq(LogData[X, ] ~ log(SeqDepth), tau = Tau, method="br")$coef[2], silent=T)
    if(class(slope) == "try-error"){
      slope <- try(rq(LogData[X, ] ~ log(SeqDepth), tau = Tau, method="fn")$coef[2], silent=T)
      if(class(slope) == "try-error"){
		slope <- NA
	  }	
    }
  }
  names(slope) <- X
  
  return(slope)
}


redobox <- function(DATA, smallc) {

  DATA[DATA <= smallc] <- NA #truncate some values first, usually just zero
  y <- log(DATA)
  
return(y)
}
