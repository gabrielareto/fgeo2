# Author: Gabriel Arellano
# gabriel.arellano.torres@gmail.com

# Funtion to generate many random torus translations of a set of points in n dimensions.
# 'data' is a matrix with one row per point and one column per dimension.
# 'mins' and 'maxs' are the lower and upper limits in any dimension. If not provided,
# the function will use the range of the data as lower and upper limits.
# Returns one matrix per dimension, with points in rows and simulations in columns.
# For a plot with known bounds, supply mins and maxs explicitly.
translate_around_torus <- function(data, mins = rep(NA, ncol(data)), maxs = rep(NA, ncol(data)), nsim = 999)
{
	sapply(1:ncol(data), function(i) {
		x <- data[,i]
		minx = mins[i]
		maxx = maxs[i]
		if(is.na(minx)) minx = min(x, na.rm = TRUE)
		if(is.na(maxx)) maxx = max(x, na.rm = TRUE)
		steps <- runif(nsim, 0, maxx - minx)
		minx + (outer(x, steps, '+') - minx) %% (maxx - minx)
	}, simplify = FALSE)
}

# Examples are disabled when sourcing this file. Run either block explicitly.
if(FALSE)
{
	set.seed(123)
	# Easy example for visualization:
	x = c(1, 2, 2, 2, 3, 4)
	y = c(1, 2, 3, 3, 3, 4)
	data <- cbind(x, y)

	tt1 <- translate_around_torus(data = data, nsim = 1, mins = c(0, 0), maxs = c(5, 5))
	tt2 <- translate_around_torus(data = data, nsim = 1, mins = c(0, 0), maxs = c(5, 5))

	par(mfrow=c(3, 3), mar = c(0, 0, 0, 0))
	for(i in 1:9) 
	{
		plot(data, xlim = c(0, 5), ylim = c(0, 5), xaxs = "i", yaxs = "i", type = "b")
		points(tt1[[1]], tt1[[2]], col = "red")
		points(tt2[[1]], tt2[[2]], col = "green")
	}
}

if(FALSE)
{
	set.seed(123)
	# Example with bigger dataset:
	n = 200000 # 200000 is roughly the sample size of CTFS plots
	x = runif(n, 0, 1000)
	y = runif(n, 0, 500)
	z = runif(n, 700, 1500)
	data <- cbind(x, y, z)
	head(data)

	system.time({tt <- translate_around_torus(data = data, nsim = 999)})
	str(tt)
}
