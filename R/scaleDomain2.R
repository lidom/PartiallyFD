scaleDomain2 <-
function(x, a=0.5, b=1){ # x <- pnorm(scores1, mean, mean = xi_1_mean, sd = sqrt(ev) ); plot(x); hist(x)
	# Cuts a [0,1] into [0.5, 1] by transforming all values bigger than 0.5 uniformly to [0.5,1] and taking all values smaller than 0.5 to 0.5 and 
	x[which(x >= a) ] <- ((b-a)*(x - min(x)))/(max(x) - min(x)) + a
	x[which(x < a) ] <- a
	return(x)
 }
