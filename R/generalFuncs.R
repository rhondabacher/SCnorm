
quickreg<-function(x,InputData)
{
	LogData = InputData[[1]]
	SeqDepth = InputData[[2]]
	X = InputData[[3]][x]
	Tau = InputData[[4]]
	
	slope <- rq(LogData[X, ] ~ log(SeqDepth), tau = Tau)$coef[2] 
	names(slope) <- X
	
	return(slope)
}

redobox <- function(DATA, smallc) {

  DATA[DATA <= smallc] <- NA #truncate some values first, usually just zero
  y <- log(DATA)
  
return(y)
}
