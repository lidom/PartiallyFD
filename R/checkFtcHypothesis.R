checkFtcHypothesis 	<- function(hypDec) {
	# hypDec <- c(TRUE, FALSE, TRUE)
	# Decides whether a sequence of hypotheses belongs to
	# F1: All FALSE
	# F2: At least H1 TRUE and another Hypothesis H_j true
	# C: Only H1 true
	if( sum(hypDec) == 0 ) { 
		return("Null") # if all are false (if independent)
	} else if( isTRUE(hypDec[1]) && !any(hypDec[-1]) ) { 
		return("(V)") # if first true and rest false
	} else if( isTRUE(hypDec[1]) ) { 
		return("Other with (V)") # if first true and some others
	} else {
		return("Other")
	} 
}

