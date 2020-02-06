# Doubly robust uniform confidence band for the conditional average treatment effect funtion
# Sokbae Lee, Ryo Okui & Yoon-Jae Whang
# Dependecies : locpol, KernSmooth, ggplot2 (not yet)

# function argument should be changed
# Y : depedent variable
# X : covariate of interest
# V : other covariates
# D : treatment 
# alpha : confidence level, default = 0.05
# ps : propensity score, default = "logit"

drcate <- function(Y, X, V, D, alpha = 0.05, ps = "logit"){

require("locpol")
require("KernSmooth")
require("ggplot2")

### initialization 

Y <- as.matrix(Y)
X <- as.matrix(X)
D <- as.matrix(D)
V <- as.matrix(V)
  
Z <- cbind(X, V)                    # Z : covariates
Y_1 <- Y[D==1]                      # Y_1 : depedent variable of treated samples
Y_0 <- Y[D==0]  		            # Y_0 : depedent variable of controlled samples
Z_1 <- Z[D==1,]                     # Z_1 : covariates of treated samples
Z_0 <- Z[D==0,]		            # Z_0 : covariates of controlled samples
k <- ncol(Z)                        # k : number of covariates
n <- nrow(Z)                        # n : number of samples

N = 101					# N : number of grid	
x_min <- sort(X)[floor(0.1 * n)]
x_max <- sort(X)[ceiling(0.9 * n)]
x_grid <- seq(x_min, x_max, length.out = N)

### propensity score
pi_hat <- fitted(glm(D~Z, family = binomial(link = ps)))  # pi_hat : fitted value of propensity score

### regression
mu_1 <- cbind(1, Z) %*% as.matrix(lm(Y_1~Z_1)$coef)   # mu_1 : fitted value of regression using treated samples
mu_0 <- cbind(1, Z) %*% as.matrix(lm(Y_0~Z_0)$coef)   # mu_2 : fitted value of regression using controlled samples

### conditional ATE given all Zs
psi_1 <- D * Y / pi_hat - (D - pi_hat) * mu_1 / pi_hat                    # psi_1 : augmented inverse probability weighting of treated
psi_0 <- (1 - D) * Y / (1 - pi_hat) + (D - pi_hat) * mu_0 / (1 - pi_hat)  # psi_0 : augmented inverse probability weighting of controlled
psi <- psi_1 - psi_0                                            

### unconditional ATE
ate <- mean(psi)

ate_func <- Vectorize(FUN = function(x){
	result <- ate
}, vectorize.args = "x")

### Gaussian Kernel
ker <- gaussK            # library("locpol")
rk <- computeRK(ker)
lambda <- 1 - (1 / (4 * sqrt(pi))) / rk

### bandwidth selection : Ruppert, Sheather and Wand (1995)
h_RSW <- dpill(x = X, y = psi)     # library("KernSmooth")
if (h_RSW == "NaN"){
	h <- 0.01
}else{
	h <- h_RSW * n^(1/5) * n^(-2/7) 
}

### conditional average treatment effect
drcate_func <- Vectorize(FUN = function(x){
	ghat <- locPolSmootherC(x = X, y = psi, xeval = x, bw = h, deg = 1, kernel = ker)$beta0
	result <- ghat
}, vectorize.args = "x")

### level alpha critical value
a2 <- 2*log(h^(-1)*(x_max-x_min))+2*log(lambda^(1/2)/(2*pi))
if (a2 >= 0){
	a <- sqrt(a2)
}else{
	a <- 0
}
critical <- sqrt(a^2 - 2*log(log((1 - alpha)^(-1/2))))

### confidence interval
drcate_ci_func <- Vectorize(FUN = function(x){
   	fX_hat <- mean(ker((X - x)/h))/h
	sigmasq_hat <- mean((psi-drcate_func(x))^2*ker((X-x)/h))/fX_hat/h
	sigmasq_hat <- sigmasq_hat * n / (n - 3 * k - 3)
	s_hat <- sqrt(rk * sigmasq_hat / fX_hat)
	sg_hat <- s_hat/sqrt(n*h)
	l_ci <- drcate_func(x) - critical * sg_hat
	u_ci <- drcate_func(x) + critical * sg_hat

	ci <- cbind(l_ci, u_ci)
	colnames(ci) <- c("95% CI lower", "95% CI upper") 
	result <- ci
}, vectorize.args = "x")

### plot
x_lim <- cbind(x_min, x_max)
ghat <- drcate_func(x_grid)
ci <- drcate_ci_func(x_grid)
drcate_ci <- t(ci)

y_max <- max(drcate_ci)
y_min <- min(drcate_ci)
y_max <- (y_max - y_min) / 2 + y_max
y_min <- y_min - (y_max - y_min) / 3 

drcate_plot <- ggplot(data = data.frame(x = x_grid), aes(x = x_grid))
drcate_plot <- drcate_plot + geom_line(aes(x = x_grid, y = ghat, colour = "CATE")) 
drcate_plot <- drcate_plot + geom_line(aes(x = x_grid, y = ate * rep(1, N), colour = "ATE"), linetype = "dashed")
drcate_plot <- drcate_plot + scale_color_manual(values = c('CATE' = 'black','ATE' = 'darkblue')) + labs(color = 'Treatment Effect')
drcate_plot <- drcate_plot + geom_ribbon(aes(x = x_grid, ymin = drcate_ci[,1], ymax = drcate_ci[,2]), alpha = 0.1)
drcate_plot <- drcate_plot + labs(x = "x", y = "") + ylim(y_min, y_max)

### results 
result <- list(ate = ate_func(0), bwidth = h, drcate_func = drcate_func, drcate_ci_func = drcate_ci_func, drcate_plot = drcate_plot)

return(result)
}