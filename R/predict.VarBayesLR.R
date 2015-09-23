#' predict.VarBayesLR 
#' 
#' Make predictions using a fitted VarBayesLR model
#'  
#' @param object Fitted model from \code{\link{VarBayesLR}}.
#' @param newdata A new input data frame.
#' @param type of prediction: "prop" for probabilities and "class" for 
#' class labels.
#' @param ... Further arguments passed to or from other methods.
#' @return predictions depending on type 
#'
#' @author  Che Ngufor Ngufor.Che@@mayo.edu
#' @export
predict.VarBayesLR <- function(object, newdata, type = c("class", "prob")[2], ...){
w <- object$w 
V <- object$V
thresh <- object$thresh
type.prob <- object$type.prob 

if(is.null(object$form)) {
cat("Direct call to the default method assumes data matrix is all numeric ! ") 
X <- data.matrix(newdata)    
} else {
X <- model.matrix(object = object$form, data = newdata)
X <- X[, -1, drop = FALSE]
}
x = X%*%w
switch(type.prob,
      logistic={prob = sigmoid(x)},
      Mckay={prob <- Mackay.Approximation(X, w, V)} 
)
if(type == "class") 
pred = ifelse(prob >= thresh, 1, 0)
else if(type == "prob") {
pred =  cbind(1-prob, prob)
colnames(pred) = c("0", "1")
} else stop("type unknown") 
return(pred)
}


