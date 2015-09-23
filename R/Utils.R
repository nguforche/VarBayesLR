#' Utility functions.
#' 
#' Useful functions  
#' Some are used internally by other functions  
#'
#' @name Utils 
#' @param x,y A numeric or character vector depending on the function 
#' @param formula formula object  
#' @param pred predicted probabilities 
#' @param obs observed outcomes 
#' @param threshold classification cut-off 
#' @return  A numeric or character vector depending on the function 
#'
#' @importFrom PresenceAbsence optimal.thresholds
#' @importFrom PresenceAbsence presence.absence.accuracy
#'
NULL 
#' @rdname  Utils 
sigmoid <- function(x) 1.0 / (1 + exp(-x))
#' @rdname  Utils  
sigminus <- function(x) 1.0 + exp(-x) 
#' @rdname  Utils  
sigplus <- function(x) 1.0 + exp(x)
#' @rdname  Utils  
tr <- function(x) matrix.trace(x)
#' @rdname  Utils  
hprod <- function(x,y){hadamard.prod(x, y)}
#'
#' @rdname Utils
### log determinant 
logdet <- function(x) as.numeric(determinant(x)$modulus)   
#'
#' @rdname Utils
## tangent in approximation of logistic function for variational Bayes logistic 
## regression 
mytanh <- function(x) {
tn <- tanh(x/2)/(4*x)
tn[is.nan(tn)]<-1/8
return(tn)
} 
#' @rdname Utils
# get left hand side of formula
lhs.form <- function(formula) {
    tt <- terms(formula)
    vars <- as.character(attr(tt, "variables"))[-1] ## [1] is the list call
    response <- attr(tt, "response") # index of response var
    vars[response] 
}
#' @rdname  Utils
# get right hand side of formula 
rhs.form <- function(formula) {
    tt <- terms(formula)
    vars <- as.character(attr(tt, "variables"))[-1] ## [1] is the list call
    response <- attr(tt, "response") # index of response var
    vars[-response] 
}
#' @rdname  Utils
#' @export
opt.thresh <- function(pred, obs){
thresh = 0.5 
if(length(unique(obs)) > 1){
obs <- as.numeric(as.factor(obs))-1 
SIMDATA = cbind.data.frame(plotID = 1:length(obs), Observed = obs, Predicted = pred)
thresh <- optimal.thresholds(SIMDATA, threshold = 101, which.model = 1, opt.methods = 9)
thresh <- ifelse(length(thresh["Predicted"]) >= 1,as.numeric(thresh["Predicted"]), 0.5)
}
return(thresh)
}
#' @rdname  Utils 
#' @export
## performance measures 
Performance.measures <- function(pred, obs, threshold=NULL){
obs <- as.numeric(factor(obs))-1 
## get best cut-off 
if(is.null(threshold))
threshold <- opt.thresh(pred, obs)
### get these performance measures
nme = c("PCC", "PCC.sd", "AUC", "AUC.sd", "sensitivity", "sensitivity.sd", 
"specificity", "specificity.sd")
xx = cbind.data.frame(plotID = 1:length(pred), Observed = obs, Predicted = pred)
accuracy <- presence.absence.accuracy(xx, threshold = threshold, st.dev = TRUE)[, nme]
accuracy$G.mean <- sqrt(as.numeric(accuracy$sensitivity)*as.numeric(accuracy$specificity))
accuracy$BER <- 1 - 0.5*(as.numeric(accuracy$sensitivity) + as.numeric(accuracy$specificity))  
return(accuracy)
}




