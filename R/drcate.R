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
N <- 100                            # N : number of grids
x_min <- sort(X)[floor(0.1 * n)]
x_max <- sort(X)[ceiling(0.9 * n)]
x_axis <- seq(x_min, x_max, length.out = N)      # x_axis : discretization of X axis


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
ate <- mean(psi) * rep(1, N)

### Gaussian Kernel
ker <- gaussK            # library("locpol")
rk <- computeRK(ker)
lambda <- 1 - (1 / (4 * sqrt(pi))) / rk

### bandwidth selection : Ruppert, Sheather and Wand (1995)
h_RSW <- robustdpill(x = X, y = psi)     # library("KernSmooth")
h <- h_RSW * n^(1/5) * n^(-2/7)
if (h == "NaN"){
h <- 0.01
}
ghat <- locPolSmootherC(x = X, y = psi, xeval = x_axis, bw = h, deg = 1, kernel = ker)$beta0

### standard error
fX_hat <- numeric(N)
sigmasq_hat <- numeric(N)
s_hat <- numeric(N)

for(i in 1 : N){
fX_hat[i] <- mean(ker((X - x_axis[i])/h))/h
sigmasq_hat[i] <- mean((psi-ghat[i])^2*ker((X-x_axis[i])/h))/fX_hat[i]/h
sigmasq_hat[i] <- sigmasq_hat[i] * n / (n - 3 * k - 3)
s_hat[i] <- sqrt(rk * sigmasq_hat[i] / fX_hat[i])
}

### level alpha critical value
a2 <- 2*log(h^(-1)*(x_max-x_min))+2*log(lambda^(1/2)/(2*pi))
if (a2 >= 0){
a <- sqrt(a2)
}else{
a <- 0
}
critical <- sqrt(a^2 - 2*log(log((1 - alpha)^(-1/2))))


### confidence interval
sg_hat <- s_hat/sqrt(n*h)
cblower <- numeric(N)
cbupper <- numeric(N)
cblower <- ghat - critical * sg_hat
cbupper <- ghat + critical * sg_hat

y_max <- max(cbupper)
y_min <- min(cblower)
y_max <- (y_max - y_min) / 2 + y_max
y_min <- y_min - (y_max - y_min) / 3  

png("result.jpg")

plot(x_axis, ghat, ylim = c(y_min, y_max), xlab ="X", ylab="", type="l", lwd = 2)
par(new = T)
plot(x_axis, cblower, ylim = c(y_min, y_max), xlab ="", ylab="", type="l", lty =2)
par(new = T)
plot(x_axis, cbupper, ylim = c(y_min, y_max), xlab ="",ylab="", type="l", lty=2)
par(new = T)
plot(x_axis, ate, ylim = c(y_min, y_max), xlab ="",ylab="", type="l", lty=4)

legend("topright", c("CATEF", "CB", "ATE"), lty = c(1,2,4), lwd = c(2,1,1), cex =1.2)
dev.off()

}

robustdpill<- function(x, y){
tryCatch(dpill(x, y), warning = function(w) {Nan}, error = function(e) {NaN})
}
