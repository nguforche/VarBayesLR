#' SynData
#' 
#'  Generate two types of artificial binary classification data: 
#'  two simple multivariate normal distributions with shifted centers 
#'  and two multivariate normal distributions a noisy hyperplane. 
#' 
#' @param n sample size 
#' @param p dimension   
#' @param type type of artificial data.
#' @param class.prior proportion o positive examples 
#' @param noise numeric constant in (0,1) for noise level
#' @param ... Further arguments passed to or from other methods.
#' @return  A data frame  
#' @author  Che Ngufor Ngufor.Che@@mayo.edu
#' @examples 
#' head(SynData())  
#'  
#' @export
#' @import mvtnorm
#' 
SynData <- function(n=100, p=10, type = c("two normal", "hyperplane")[1], 
                 class.prior = 0.3, noise = 0.01, ...){
u1 = rep(0, p) 
if(type == "two normal"){
sigma =  diag(x = 1.5, nrow = p, ncol = p)
u2 = rep(0.5, p) 
N1 = floor(n*class.prior)
N2 = n-N1 
X1 =  cbind.data.frame(key = rep(0, N1), rmvnorm(N1, mean = u1, sigma = sigma))
X2 =  cbind.data.frame(key = rep(1, N2), rmvnorm(N2, mean = u2, sigma = sigma))
dat = rbind.data.frame(X1, X2) 
} else if(type == "hyperplane") {
Xy = rmvnorm(n, mean = u1, sigma = sigma)
a = runif(p, -1, 1)
a0 = sum(a)*class.prior
Ha = apply(t(apply(Xy, 1, "*", a)), 1, sum) + rnorm(n, 0, noise) 
ix = which(Ha > a0)
y = rep(1, n)
y[ix] = 0
dat = data.frame(key = y, Xy)
} else stop("type must be 1 or 2") 
names(dat)[-1] <- paste0("V", 1:p)
dat = dat[sample(nrow(dat), nrow(dat)), ]
return(dat)
}

