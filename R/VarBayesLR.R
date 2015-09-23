#' Variational Bayes Logistic Regression 
#' 
#' Train Variational Bayesian Logistic Regression  
#'   
#'
#' @name VarBayesLR
#'
#' @param X data.frame/matrix of predictors.
#' @param Y binary (0/1) vector of response.
#' @param dat data.frame containing predictors and a binary (0/1) response 
#'  variable in formula.
#' @param form formula 
#' @param max.iter maximum number of iteration  
#' @param tol convergence tolerance 
#' @param type.prob type of approximation for posterior probabilities: 
#'  Mckay = Mckay approximation or logistic = sigmoid function the default
#' @param w estimated predictor weights 
#' @param V covariance matrix of weights
#' @param \dots Further arguments passed to or from other methods.
#' @return An object of class \code{VarBayesLR}; a list with items 
#' \item{w}{weights}
#' \item{V}{weight covariance matrix}
#' \item{inV}{inverse of V} 
#' \item{VarImp}{Variable Relevance: These are mean of gamma distribution: ratio of hyperparameters}
#' \item{L}{variational lower bound}
#' \item{prob}{fitted probabilities}
#' \item{thresh}{optimal classification threshold}
#' \item{type}{type of approximation for posterior probabilities}
#'
#' @references
#' MacKay, David JC. "The evidence framework applied to classification networks." 
#' Neural computation 4.5 (1992): 720-736.
#'
#' Jaakkola, T., and Michael I. Jordan. "A variational approach to Bayesian logistic 
#' regression models and their extensions." Sixth International Workshop on Artificial 
#' Intelligence and Statistics. 1997.

#' @author  Che Ngufor Ngufor.Che@@mayo.edu
#' @import matrixcalc 
#' 
NULL 
#'
#'
#' @rdname VarBayesLR  
#' @export
VarBayesLR  <- function(X, ...) UseMethod("VarBayesLR")
#'
#' @rdname VarBayesLR  
#' @export
#' @examples
#' \dontrun{
#' set.seed(12345)
#' dat <- SynData()
#' ix = sample(nrow(dat), floor(nrow(dat)*0.75))
#' dat.trn = dat[ix, ]
#' dat.tst = dat[-ix, ]
#' form <- as.formula(paste("key ~ ", paste(names(dat)[!names(dat)%in%"key"], collapse = "+")))
#' mod <- VarBayesLR(form, dat.trn)
#' pred <- predict(mod, dat.tst)
#' perf <- Performance.measures(pred[,2], dat.tst$key, mod$thresh)
#' perf
#' }
#'
VarBayesLR.default <- function(X, Y, max.iter=100, tol=1e-6, 
               type.prob = c("Mckay", "logistic")[2], ...){
               
n = nrow(X)
D = ncol(X) 
Y <- as.numeric(factor(Y))-1

## hyperprior parameters
a0 = 1e-2;
b0 = 1e-4;
an = a0 + 0.5; 
 
D_gammaln_an_an = D * (lgamma(an) + an);
t_w = matrix(0.5*colSums ((2*Y-1)*X, na.rm = TRUE), ncol = 1)  
 
## start first iteration kind of here, with xi = 0 -> lam.xi = 1/8 
lam.xi = rep(1/8, n); 
E_a = rep(a0/b0, D);
invV = diag(E_a) + 2 * t(X)%*% (lam.xi*X)
V = solve(invV);
w = V%*%t_w;
bn = as.numeric(b0 + 0.5 * (w^2 + diag(V)));

L_last = -n * log(2) + 0.5 * (t(w)%*%invV%*%w - logdet(invV)) - sum((b0 * an)/bn) 
         - sum(an * log(bn)) + D_gammaln_an_an;

for(ii in 1:max.iter){
## update xi by EM-algorithm
    xi = sqrt(rowSums(hprod(X, X %*%(V + w%*%t(w)))));
    lam.xi = mytanh(xi);
# update posterior parameters of a based on xi
    bn = as.numeric(b0 + 0.5 * (w^2 + diag(V)));
    E_a = an/bn;
# recompute posterior parameters of w
#    % invV = diag(E_a) + 2 * X' * (X .* repmat(lam.xi, 1, D))
invV = diag(E_a) + 2 * t(X)%*% (X*lam.xi);
V = solve(invV);
logdetV = - logdet(invV);
w = V%*%t_w;
# variational bound, ingnoring constant terms for now
LowB = - sum(log(1 + exp(-xi))) + sum(lam.xi* xi^2) + 0.5 * (t(w)%*%invV%*%w + logdetV - sum(xi))- 
        sum(b0 * E_a) - sum(an * log(bn)) + D_gammaln_an_an;
        
err = abs(L_last-LowB)
if(err <= tol){
#print(sprintf("Variational Bayes EM Converged at %d iteration with error  %f.3 ", ii, err))
break 
}
L_last = LowB;  
}
if(err > tol) warning("Variational Bayes EM did not converge")
LowB = LowB - D * (lgamma(a0) - a0 * log(b0));


x = X%*%w
switch(type.prob,
      logistic={prob = sigmoid(x)
          thresh <- opt.thresh(prob, Y)},
      Mckay={prob <- Mackay.Approximation(X, w, V)
           thresh <-  opt.thresh(prob, Y)}, 
      stop("Wrong probability approximation choice: choices are Mckay or logistic")
)
names(E_a) = colnames(X)
m1 <- min(E_a)
m2 <- max(E_a)
VarImp <- (((E_a - m1)*(99))/(m2-m1) + 1)
VarImp <- data.frame(VarImp)
names(VarImp)[1] <- "importance"
VarImp <- VarImp[order(VarImp[, "importance"], decreasing = TRUE), ,drop=FALSE]   

res <- list(w = w, V = V, invV = invV, VarImp = VarImp, LowB= LowB, prob = prob, thresh = thresh, 
type.prob = type.prob)
class(res) <- "VarBayesLR"
return(res)
}
#'
#' @rdname VarBayesLR
#' @export
VarBayesLR.formula <- function(form, dat, max.iter=100, tol=1e-6, 
                           type.prob = c("Mckay", "logistic")[2], ...){                                                         
Y = dat[, lhs.form(form)]
X <- model.matrix(object = form, data = dat)
X <- X[, -1, drop = FALSE]
res <- VarBayesLR(X, Y, max.iter, tol, type.prob)
res$form <- form
class(res) <- "VarBayesLR"
return(res)
}
#'
#' @rdname VarBayesLR
#' @export
#'
### Mackay approximation for continuous outputs 
Mackay.Approximation <- function(X, w, V){
mu <- X%*%w 
s2 <- rowSums(hprod(X, X %*%V))  
l <- mu/sqrt(1+ pi*s2/8)
post = sigmoid(l)
return(post)
}
