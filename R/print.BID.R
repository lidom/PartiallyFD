print.BID <-
function(x){
	cat(paste0("Selected Basis Length: ", x$selBasis, "\n"))
	cat(paste0("Correlation of d and 1. Score: ", round(x$cor,2), "\n"))
	hypDec 	<- checkFtcHypothesis(x$romWolf$ent)
	cat(paste0("RomWolf scenario ", hypDec, "\n"))
	cat(paste0("Romano Wolf p-Value: ", x$romWolf$pval, "\n"))
}
