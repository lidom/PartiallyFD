integrate_xy <-
function (x,y) {
## Summary: Calulates the function F(z) = int_0^z f(x) dx
## Input:
# Two vectors x, y of same length  
## Output:
# Gives a vector of same length with evaluated numbers at each point
	if(length(x) != length(y)) stop("Both vectors should have same length")
    idx = 2:length(x)
	int_right 				<- cumsum(as.double((x[idx] - x[idx - 1]) * (y[idx] + y[idx - 1]))/2)
	return(c(0,int_right))
}
