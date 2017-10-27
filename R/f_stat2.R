f_stat2 <-
function(X, d, k){
	# calculates the univariate F statistics from a multipe linear regression
	# d ~ a_1 x_1 + ... + a_b x_b
    # dim(data.frame(d, X))
    sum <- summary(lm <- lm(d ~ 1 + ., data = data.frame(d, X)))
    # coefs(sum)
    # if there are NA coefs, sum responds only non NA and therefore the length of sum = k+1 is less!
    if(dim(sum$coefficients)[1] != (k+1)){
        ret <- numeric(k)
        ret[which(is.na(coef(lm)[-1]))] <- 0
        ret[which(!is.na(coef(lm)[-1]))] <- sum$coefficients[2:(dim(sum$coefficients)[1]) , "t value"]^2
        return(ret)
    } else {
        return(sum$coefficients[2:(k+1) , "t value"]^2)
    }
}
